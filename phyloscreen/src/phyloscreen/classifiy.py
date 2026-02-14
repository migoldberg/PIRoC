"""
classifiy.py
Provides core functions and classification logic.
"""

import os 
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

from .metadata import parse_species_name
from .support import get_node_bootstrap
from .branch_length import is_on_long_branch

def compute_clade_composition(node, species_to_group):
    species_counts = Counter()
    group_counts = Counter()

    for leaf in node.get_leaves():
        if not leaf.name:
            continue
        species_id = parse_species_name(leaf.name)
        species_counts[species_id] += 1

        group = species_to_group.get(species_id, "UNKNOWN")
        group_counts[group] += 1

    total_leaves = sum(species_counts.values())
    return species_counts, group_counts, total_leaves

def analyze_sister_clades(focal_mrca, species_to_group, focal_group, outgroup_names):
    if outgroup_names is None:
        outgroup_names = {"Outgroup"}

    parent_node = focal_mrca.up

    if parent_node is None:
        # focal_mrca is root, so there is no sister clade
        return {
            "sister_groups": Counter(),
            "sister_species": Counter(),
            "has_outgroup_sister": False,
            "has_same_group_sister": False,
            "sister_is_pure_outgroup": False
        }
    
    sisters = [child for child in parent_node.children if child != focal_mrca]

    sister_groups = Counter()
    sister_species = Counter()

    for sister in sisters:
        for leaf in sister.get_leaves():
            if not leaf.name:
                continue
            species = parse_species_name(leaf.name)
            group = species_to_group.get(species, "UNKNOWN")
            sister_species[species] += 1
            sister_groups[group] += 1

    has_outgroup = any(g in outgroup_names for g in sister_groups.keys())
    has_same_group = focal_group in sister_groups

    # check if sister is purely outgroup
    sister_is_pure_outgroup = (
        len(sister_groups) > 0 and
        all(g in outgroup_names for g in sister_groups.keys())
    )

    return {
        "sister_groups": sister_groups,
        "sister_species": sister_species,
        "has_outgroup_sister": has_outgroup,
        "has_same_group_sister": has_same_group,
        "sister_is_pure_outgroup": sister_is_pure_outgroup
    }

def calculate_clade_target_group_fraction(node, species_to_group, focal_group):
    species_counts, group_counts, total_leaves = compute_clade_composition(node, species_to_group)
    if total_leaves == 0:
        return 0.0
    else:
        return group_counts.get(focal_group, 0) / total_leaves

def mrca_clade_contains_outgroup_intruders(group_counts, outgroup_names):
    return any(g in outgroup_names for g in group_counts.keys())

def other_focal_group_leaves_in_clade(species_counts, species_to_group, focal_group, focal_species):
    focal_leavs_in_clade = sum(
        count for sp, count in species_counts.items()
        if species_to_group.get(sp) == focal_group and sp != focal_species
    )
    return focal_leavs_in_clade > 0

def classify_sequence(
    sequence, 
    focal_species, 
    focal_group, 
    tree, 
    species_to_group, 
    branch_length_stats, 
    min_support, 
    min_target_purity, 
    max_contaminant_purity, 
    outgroup_names
    ):

    if outgroup_names is None:
        outgroup_names = {"Outgroup"}

    sequence_name = sequence.name
    parent_node = sequence.up

    if parent_node is None:
        return "FLAG", {
            "sequence_name": sequence_name,
            "bootstrap": "NA",
            "clade_target_group_fraction": 0.0,
            "total_leaves": 1,
            "classification_notes": ["sequence_at_root"]
        }
    
    species_counts, group_counts, total_leaves = compute_clade_composition(parent_node, species_to_group)
    bootstrap = get_node_bootstrap(parent_node)
    clade_target_group_fraction = calculate_clade_target_group_fraction(parent_node, species_to_group, focal_group)
    sequence_on_long_branch = is_on_long_branch(sequence, tree, branch_length_stats)
    sister_clade_metrics = analyze_sister_clades(parent_node, species_to_group, focal_group, outgroup_names)
    has_outgroup_in_clade = mrca_clade_contains_outgroup_intruders(group_counts, outgroup_names)
    has_other_focal_group_leaves_in_clade = other_focal_group_leaves_in_clade(species_counts, species_to_group, focal_group, focal_species)
    has_other_focal_species_leaves_in_clade = species_counts.get(focal_species, 0)

    # ------ Classify ------
    classification = "FLAG"
    classification_notes = []

    high_support = bootstrap is not None and bootstrap >= min_support

    if sister_clade_metrics["sister_is_pure_outgroup"] and high_support:
        classification = "CONTAMINANT"
        classification_notes.append("sister_clade_pure_outgroup")
    elif sequence_on_long_branch and has_outgroup_in_clade:
        classification = "CONTAMINANT"
        classification_notes.append("long_branch_with_outgroup_in_clade")
    elif high_support and clade_target_group_fraction <= max_contaminant_purity:
        classification = "CONTAMINANT"
        classification_notes.append("low_ctgf_high_support")

    elif high_support and clade_target_group_fraction >= min_target_purity:
        if not sequence_on_long_branch and not has_outgroup_in_clade:
            classification = "TARGET"
            classification_notes.append("high_ctgf_high_support")

            if has_other_focal_group_leaves_in_clade:
                classification_notes.append("other_focal_group_in_clade")
            if has_other_focal_species_leaves_in_clade:
                classification_notes.append("other_focal_species_in_clade")
        elif sequence_on_long_branch:
            classification = "FLAG"
            classification_notes.append("high_ctgf_but_long_branch")
        elif has_outgroup_in_clade:
            classification = "FLAG"
            classification_notes.append("high_ctgf_but_outgroup_in_clade")
    
    else:
        if bootstrap is None or bootstrap < min_support:
            classification_notes.append("low_or_missing_support")
        if clade_target_group_fraction > max_contaminant_purity and clade_target_group_fraction < min_target_purity:
            classification_notes.append("intermediate_target_group_fraction")
        if sequence_on_long_branch:
            classification_notes.append("long_branch_detected")
    
    metrics = {
        "sequence_name": sequence_name,
        "bootstrap": bootstrap if bootstrap is not None else "NA",
        "clade_target_group_fraction": clade_target_group_fraction,
        "total_leaves": total_leaves,
        "focal_group": focal_group,
        "focal_group_leaves": group_counts.get(focal_group, 0),
        "focal_species_in_clade": species_counts.get(focal_species, 0),
        "has_other_focal_group_leaves_in_clade": has_other_focal_group_leaves_in_clade,
        "has_other_focal_species_leaves_in_clade": has_other_focal_species_leaves_in_clade,
        "has_outgroup_in_clade": has_outgroup_in_clade,
        "is_on_long_branch": sequence_on_long_branch,
        "classification_notes": classification_notes
    }

    return classification, metrics
   