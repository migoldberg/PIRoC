# Developed by Mark Goldberg

"""
A tool for identifiction and removal of contaminants in inverterbrate trasncriptomes using phylogenetic trees.
This script is the main entry point and workflow controller.
"""

import os 
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

from .cli import init_cli
from .metadata import parse_species_name, get_group_from_species_name
from .root import root_tree_at_outgroup
from .support import collapse_low_support_nodes
from .branch_length import compute_branch_length_stats
from .classifiy import classify_sequence
from .output import write_summary, write_sequence_classifications, write_sequence_lists
from .clean_trees import clean_trees

def main():
    params = init_cli()

    logger = params["logger"]
    sys.stdout = logger

    tree_dir = params["tree_dir"]
    tree_suffix = params["tree_suffix"]
    output_dir = params["output_dir"]
    focal_species = params["focal_species"]
    focal_group = params["focal_group"]
    species_to_group = params["species_to_group"]
    min_support = params["min_support"]
    min_target_purity = params["min_target_purity"]
    max_contaminant_purity = params["max_contaminant_purity"]
    collapse_threshold = params["collapse_threshold"]
    outgroup_names = params["outgroups"]
    remove_contaminants = params["remove_contaminants"]
    debug = params["debug"]

    sequence_classifications = {} # og_id::sequence_name -> classification
    sequence_metrics = {} # og_id::sequence_name -> metrics

    run_metrics = {
        'no_focal_species_in_tree': 0,
        'total_trees': 0,
        'total_errors': 0,
        'total_sequences_classified': 0,
        'total_nodes_collapsed': 0,
    }

    for filename in os.listdir(tree_dir):
        if not filename.endswith(tree_suffix):
            continue
        
        run_metrics['total_trees'] += 1
        tree_path = os.path.join(tree_dir, filename)
        og_id = filename.replace(tree_suffix, "")

        try:
            with open(tree_path, 'r') as f:
                tree_string = f.read().strip()

            t = Tree(tree_string, format=0)
           
            outgroup_leaves = [
                l for l in t.get_leaves()
                if get_group_from_species_name(l.name, species_to_group) in outgroup_names 
            ]

            tree_metrics = {
                'total_leaves': len(t.get_leaves()),
                'total_outgroup_leaves': len(outgroup_leaves),
                'rooted': False,
                'nodes_collapsed': 0,
                'branch_length_mean': None,
                'branch_length_std': None,
                'branch_length_long_branch_threshold': None,
                'focal_leaves': 0,
            }

            t, rooted = root_tree_at_outgroup(t, outgroup_leaves)
            
            
            nodes_collapsed = collapse_low_support_nodes(t, collapse_threshold)
            branch_length_stats = compute_branch_length_stats(t)

            focal_leaves = [
                l for l in t.get_leaves()
                if l.name and parse_species_name(l.name) == focal_species
            ]

            if not focal_leaves:
                run_metrics['no_focal_species_in_tree'] += 1
                continue

            tree_metrics['rooted'] = rooted
            tree_metrics['nodes_collapsed'] = nodes_collapsed
            tree_metrics['branch_length_mean'] = branch_length_stats['mean']
            tree_metrics['branch_length_std'] = branch_length_stats['std']
            tree_metrics['branch_length_long_branch_threshold'] = branch_length_stats['long_branch_threshold']
            tree_metrics['focal_leaves'] = len(focal_leaves)
            run_metrics['total_nodes_collapsed'] += nodes_collapsed

            # classify each sequence
            for sequence in focal_leaves:
                classification, metrics = classify_sequence(
                    sequence = sequence,
                    focal_species = focal_species,
                    focal_group = focal_group,
                    tree = t,
                    species_to_group = species_to_group,
                    branch_length_stats = branch_length_stats,
                    min_support = min_support,
                    min_target_purity = min_target_purity,
                    max_contaminant_purity = max_contaminant_purity,
                    outgroup_names = outgroup_names
                )

                sequence_name = sequence.name
                sequence_id = f"{og_id}::{sequence_name}"

                sequence_classifications[sequence_id] = classification
                sequence_metrics[sequence_id] = metrics

                run_metrics['total_sequences_classified'] += 1

                notes = ",".join(metrics["classification_notes"]) if metrics["classification_notes"] else "-"
                print(
                    f"[{og_id}] {sequence_name:40s} | {classification:12s} | "
                    f"bs={metrics['bootstrap']:>5} | "
                    f"ctgf={metrics['clade_target_group_fraction']:.2f} | "
                    f"longbr={metrics['is_on_long_branch']} | "
                    f"notes={notes}"
                )

                if debug:
                    print(tree_metrics)

        except Exception as e:
            print(f"error: {e}")
            run_metrics['total_errors'] += 1

    if remove_contaminants:
        clean_trees(tree_dir, tree_suffix, sequence_classifications, output_dir)

    # Write output summary
    summary_file = os.path.join(output_dir, "classification_summary.txt")
    write_summary(summary_file, focal_species, focal_group, min_support, min_target_purity, max_contaminant_purity, collapse_threshold, outgroup_names, run_metrics, sequence_classifications)

    sequence_classifications_file = os.path.join(output_dir, "sequence_classifications.tsv")
    write_sequence_classifications(sequence_classifications_file, sequence_classifications, sequence_metrics)

    sequence_lists_dir = os.path.join(output_dir, "sequence_classifications")
    write_sequence_lists(sequence_lists_dir, sequence_classifications)  

    sys.stdout = logger.terminal
    logger.close()

if __name__ == "__main__":
    main()
