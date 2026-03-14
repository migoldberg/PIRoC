"""
classifiy.py
Provides core functions and classification logic.
"""

import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

from .metadata import parse_species_name
from .support import get_node_bootstrap
from .branch_length import is_on_long_branch

def compute_clade_composition(node, species_to_group):
    """
    Computes the composition of a clade by finding the number of differnet specise and groups in the clade.
    Returns the number of species, groups, and total leaves.
    """

    # initializes empty counters to store the number of species and groups
    species_counts = Counter()
    group_counts = Counter()

    # loops through each leaf in the clade
    for leaf in node.get_leaves():
        if not leaf.name:
            continue

        # parses the species name from the leaf 
        species_id = parse_species_name(leaf.name)

        # increments the species count for that species
        species_counts[species_id] += 1

        # gets the metadata provided group for the species
        group = species_to_group.get(species_id, "UNKNOWN")

        # increments the group count for that group
        group_counts[group] += 1

    # sums the total number of leaves in the clade
    total_leaves = sum(species_counts.values())

    # returns the species counts, group counts, and total leaves
    return species_counts, group_counts, total_leaves

def analyze_sister_clades(focal_mrca, species_to_group, focal_group, contaminant_names):
    """
    Finds and analyzes the composition of the sister clade(s) of a focal clade.
    """

    # if a custom contaminant names list is not provided, uses the default group name "Contaminant"
    # this is from the metadata file
    if contaminant_names is None:
        contaminant_names = {"Contaminant"}

    # gets the parent node of MRCA of the focal clade
    parent_node = focal_mrca.up

    # if the parent node is None, then the focal clade is the root and there are no sister clades
    if parent_node is None:
        return {
            "sister_groups": Counter(),
            "sister_species": Counter(),
            "has_contaminant_sister": False,
            "has_same_group_sister": False,
            "sister_is_pure_contaminant": False
        }
    
    # gets the sister clades of the focal clade
    sisters = [child for child in parent_node.children if child != focal_mrca]

    # initializes empty counters to store the number of sister groups and species
    sister_groups = Counter()
    sister_species = Counter()

    # loops through each sister clade
    for sister in sisters:
        # loops through each leaf in the sister clade
        for leaf in sister.get_leaves():
            if not leaf.name:
                continue

            # parses the species name from the leaf
            species = parse_species_name(leaf.name)

            # parse the group from the species name
            group = species_to_group.get(species, "UNKNOWN")

            # increments the sister species and group counts
            sister_species[species] += 1
            sister_groups[group] += 1

    #  checks if any of the sister clades contain contaminant sequences
    has_contaminant = any(g in contaminant_names for g in sister_groups.keys())

    # checks if any of hte sister clades contain the focal group
    has_same_group = focal_group in sister_groups

    # checks if the sister clade is purely composed of contaminant sequences
    sister_is_pure_contaminant = (
        len(sister_groups) > 0 and
        all(g in contaminant_names for g in sister_groups.keys())
    )

    # returns the sister clade metrics
    return {
        "sister_groups": sister_groups,
        "sister_species": sister_species,
        "has_contaminant_sister": has_contaminant,
        "has_same_group_sister": has_same_group,
        "sister_is_pure_contaminant": sister_is_pure_contaminant
    }

def calculate_clade_target_group_fraction(node, species_to_group, focal_group):
    """
    Calculates the "Clade Target Group Fraction" (CTGF) for a clade.
    The math here is the number of focal group sequences divided by the total number of sequences in the clade.
    Essentially, this function determines what percentage of the clade is composed of the focal group. 
    The group is used instead of species so that clades of closely related species are not identified as contaminants.
    """

    # gets the composition of the clade
    species_counts, group_counts, total_leaves = compute_clade_composition(node, species_to_group)

    if total_leaves == 0:
        return 0.0
    else:
        return group_counts.get(focal_group, 0) / total_leaves

def mrca_clade_contains_contaminant_intruders(group_counts, contaminant_names):
    """
    Checks if the clade of the MRCA contains any contaminant intruders.
    """
    return any(g in contaminant_names for g in group_counts.keys())

def other_focal_group_leaves_in_clade(species_counts, species_to_group, focal_group):
    """
    Checks if the clade contains any other focal group sequences other than the one being classified.
    """

    # sums the number of focal group sequences in the clade other than the one being classified
    focal_leaves_in_clade = sum(
        count for sp, count in species_counts.items()
        if species_to_group.get(sp) == focal_group
    )

    # returns True if there are any other focal group sequences in the clade
    return focal_leaves_in_clade > 0

def classify_sequence(
    sequence, 
    focal_group, 
    tree, 
    species_to_group, 
    branch_length_stats, 
    min_support, 
    min_target_purity, 
    max_contaminant_purity, 
    contaminant_names
    ):
    """
    Main classification function for PIRoC.
    This function is called for each sequence in the tree individually. 
    """

    # if a custom set of contaminant names is not provided in the metadata file as the "groups" column, uses the default group name "Contaminant"
    if contaminant_names is None:
        contaminant_names = {"Contaminant"}

    # gets the name of the sequence
    sequence_name = sequence.name

    # gets the parent node of the sequence
    parent_node = sequence.up

    # if there is no parent node, then the sequence is at the root and is not classified
    if parent_node is None:
        return "FLAG", {
            "sequence_name": sequence_name,
            "bootstrap": "NA",
            "clade_target_group_fraction": 0.0,
            "total_leaves": 1,
            "classification_notes": ["sequence_at_root"]
        }
    
    # gets the composition of the clade the sequence is in
    species_counts, group_counts, total_leaves = compute_clade_composition(parent_node, species_to_group)

    # gets the bootstrap support for the clade via the parent node
    bootstrap = get_node_bootstrap(parent_node)

    # calcualtes the clade target group fraction of the clade the sequence is in
    clade_target_group_fraction = calculate_clade_target_group_fraction(parent_node, species_to_group, focal_group)

    # checks if the sequence is on a long branch
    sequence_on_long_branch = is_on_long_branch(sequence, tree, branch_length_stats)

    # analyzes the sister clade(s) of the clade the sequence is in
    sister_clade_metrics = analyze_sister_clades(parent_node, species_to_group, focal_group, contaminant_names)

    # checks if the clade containing the sequence has any contaminant intruders
    has_contaminant_in_clade = mrca_clade_contains_contaminant_intruders(group_counts, contaminant_names)

    # checks if the clade containing the sequence has any other focal group sequences other than the one being classified
    has_other_focal_group_leaves_in_clade = other_focal_group_leaves_in_clade(species_counts, species_to_group, focal_group)

    # ------ Classification ------
    # starts by assuming the sequence should be flagged
    classification = "FLAG"

    # initializes an empty list to store the classification notes
    classification_notes = []

    # creates a boolean for if a node is highly supported based on bootstrap support and the provided minimum support threshold argument
    high_support = bootstrap is not None and bootstrap >= min_support


    if sister_clade_metrics["sister_is_pure_contaminant"] and high_support:
        # if the sister clade(s) are purely contaminants AND the node is highly supported = CONTAMINANT
        classification = "CONTAMINANT"
        classification_notes.append("sister_clade_pure_contaminant")
    elif sequence_on_long_branch and has_contaminant_in_clade:
        # if the sequence is on a long branch AND the clade contains contaminant intruders = CONTAMINANT
        classification = "CONTAMINANT"
        classification_notes.append("long_branch_with_contaminant_in_clade")
    elif high_support and clade_target_group_fraction <= max_contaminant_purity:
        # if the node is highly supported AND the clade target group fraction is less than the maximum contaminant purity argument = CONTAMINANT
        classification = "CONTAMINANT"
        classification_notes.append("low_ctgf_high_support")

    elif high_support and clade_target_group_fraction >= min_target_purity:
        if not sequence_on_long_branch and not has_contaminant_in_clade:
            # if the node is highly supported AND the clade target group fraction is greater than the minimum target purity argument = TARGET
            classification = "TARGET"
            classification_notes.append("high_ctgf_high_support")

            # adds notes supporting the TARGET classification if the clade contains other focal group or species sequences
            if has_other_focal_group_leaves_in_clade:
                classification_notes.append("other_focal_group_in_clade")
        elif sequence_on_long_branch:
            # if the node is highly supported AND the clade target group fraction is greater than the minimum target purity argument 
            # AND the sequence is on a long branch = FLAG
            classification = "FLAG"
            classification_notes.append("high_ctgf_but_long_branch")
        elif has_contaminant_in_clade:
            # if the node is highly supported AND the clade target group fraction is greater than the minimum target purity argument 
            # AND the clade contains contaminant intruders = FLAG
            classification = "FLAG"
            classification_notes.append("high_ctgf_but_contaminant_in_clade")
    
    else:
        if bootstrap is None or bootstrap < min_support:
            # if the node is not highly supported = FLAG
            classification_notes.append("low_or_missing_support")
        if clade_target_group_fraction > max_contaminant_purity and clade_target_group_fraction < min_target_purity:
            # if the clade target group fraction is intermediate between the maximum contaminant purity and minimum target purity = FLAG
            classification_notes.append("intermediate_target_group_fraction")
        if sequence_on_long_branch:
            # if the sequence is on a long branch = FLAG
            classification_notes.append("long_branch_detected")
    
    # creates a dictionary of metrics for the sequence
    metrics = {
        "sequence_name": sequence_name,
        "bootstrap": bootstrap if bootstrap is not None else "NA",
        "clade_target_group_fraction": clade_target_group_fraction,
        "total_leaves": total_leaves,
        "focal_group": focal_group,
        "focal_group_leaves": group_counts.get(focal_group, 0),
        "has_other_focal_group_leaves_in_clade": has_other_focal_group_leaves_in_clade,
        "has_contaminant_in_clade": has_contaminant_in_clade,
        "is_on_long_branch": sequence_on_long_branch,
        "classification_notes": classification_notes
    }

    # returns the classification and metrics
    return classification, metrics
   