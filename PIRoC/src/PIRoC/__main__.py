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
from .root import root_tree_at_contaminant
from .branch_length import compute_branch_length_stats
from .classify import classify_sequence
from .output import write_summary, write_sequence_classifications, write_sequence_lists
from .clean_orthogroups import clean_orthogroups
from .review_flags import start_review_flags

def main():
    params = init_cli()

    logger = params["logger"]
    sys.stdout = logger

    tree_dir = params["tree_dir"]
    tree_suffix = params["tree_suffix"]
    output_dir = params["output_dir"]
    focal_group = params["focal_group"]
    species_to_group = params["species_to_group"]
    min_support = params["min_support"]
    min_clean_purity = params["min_clean_purity"]
    max_contaminant_purity = params["max_contaminant_purity"]
    contaminant_names = params["contaminants"]
    remove_contaminants = params["remove_contaminants"]
    loud = params["loud"]
    review_flags = params["review_flags"]
    
    # create empty dictionaries to store sequence classifications and metrics
    sequence_classifications = {} # og_id::sequence_name -> classification
    sequence_metrics = {} # og_id::sequence_name -> metrics

    # metrics for the run
    run_metrics = {
        'no_focal_group_in_tree': 0,
        'total_trees': 0,
        'total_errors': 0,
        'total_sequences_classified': 0,
    }

    # count total tree files before processing
    tree_files = [f for f in os.listdir(tree_dir) if f.endswith(tree_suffix)]
    total_tree_files = len(tree_files)

    # loop through each tree in the tree directory
    for filename in tree_files:
        run_metrics['total_trees'] += 1
        tree_path = os.path.join(tree_dir, filename)
        og_id = filename.replace(tree_suffix, "")

        # try to load the tree
        try:
            with open(tree_path, 'r') as f:
                tree_string = f.read().strip()

            # loads tree using ete3
            t = Tree(tree_string, format=0)
           
            # gets the contaminant leaves
            contaminants_leaves = [
                l for l in t.get_leaves()
                if get_group_from_species_name(l.name, species_to_group) in contaminant_names 
            ]

            # creates a dictionary to store tree metrics
            tree_metrics = {
                'total_leaves': len(t.get_leaves()),
                'total_contaminants_leaves': len(contaminants_leaves),
                'rooted': False,
                'branch_length_mean': None,
                'branch_length_std': None,
                'branch_length_long_branch_threshold': None,
                'focal_leaves': 0,
            }

            # attempts to root the tree at a monophyletic contaminant clade and returns the rooted (or not) tree
            t, rooted = root_tree_at_contaminant(t, contaminants_leaves)
            
            # computes branch length information for later classification of sequences on long branches
            branch_length_stats = compute_branch_length_stats(t)

            # gets the focal leaves
            focal_leaves = [
                l for l in t.get_leaves()
                if l.name and species_to_group.get(parse_species_name(l.name), "UNKNOWN") == focal_group
            ]

            # if no focal leaves are found, increment the total number of trees with no focal species and continue
            if not focal_leaves:
                run_metrics['no_focal_group_in_tree'] += 1
                continue

            # update tree metrics
            tree_metrics['rooted'] = rooted
            tree_metrics['branch_length_mean'] = branch_length_stats['mean']
            tree_metrics['branch_length_std'] = branch_length_stats['std']
            tree_metrics['branch_length_long_branch_threshold'] = branch_length_stats['long_branch_threshold']
            tree_metrics['focal_leaves'] = len(focal_leaves)

            # loops through each focal leaf in the tree
            for sequence in focal_leaves:
                # classifies the sequence and returns the classification and metrics
                classification, metrics = classify_sequence(
                    sequence = sequence,
                    focal_group = focal_group,
                    tree = t,
                    species_to_group = species_to_group,
                    branch_length_stats = branch_length_stats,
                    min_support = min_support,
                    min_clean_purity = min_clean_purity,
                    max_contaminant_purity = max_contaminant_purity,
                    contaminant_names = contaminant_names
                )

                # stores the sequence name and id
                sequence_name = sequence.name
                sequence_id = f"{og_id}::{sequence_name}"

                # stores the classification and metrics
                sequence_classifications[sequence_id] = classification
                sequence_metrics[sequence_id] = metrics

                # increment the total number of sequences classified
                run_metrics['total_sequences_classified'] += 1

                # joins the classification notes with a comma so they can be printed together later
                notes = ",".join(metrics["classification_notes"]) if metrics["classification_notes"] else "-"

                # prints quick classification information for the sequence to the console during the run               
                if loud:
                    print(
                        f"[{og_id}] {sequence_name:40s} | {classification:12s} | "
                        f"bs={metrics['bootstrap']:>5} | "
                        f"ctgf={metrics['clade_target_group_fraction']:.2f} | "
                        f"longbr={metrics['is_on_long_branch']} | "
                        f"notes={notes}"
                    )
                else:
                    print(f"\r  Classifying Sequences: [{run_metrics['total_trees']}/{total_tree_files}]", end="", flush=True)
                
        except Exception as e:
            print(f"error: {e}")
            run_metrics['total_errors'] += 1

    if not loud:
        print()

    if review_flags:
        start_review_flags(tree_dir, tree_suffix, sequence_classifications, sequence_metrics)

    # if the remove-contaminants flag is enabled, the function to produce clean trees is called
    if remove_contaminants:
        clean_orthogroups(tree_dir, tree_suffix, sequence_classifications, contaminant_names, species_to_group, output_dir, loud)
        print()
        

    print()
    print("\033[4mWriting Output\033[0m")

    # writes the summary file
    summary_file = os.path.join(output_dir, "classification_summary.txt")
    write_summary(summary_file, focal_group, min_support, min_clean_purity, max_contaminant_purity, contaminant_names, run_metrics, sequence_classifications)

    sequence_classifications_file = os.path.join(output_dir, "sequence_classifications.tsv")
    write_sequence_classifications(sequence_classifications_file, sequence_classifications, sequence_metrics)

    sequence_lists_dir = os.path.join(output_dir, "sequence_classifications")
    write_sequence_lists(sequence_lists_dir, sequence_classifications)  

    print("\033[4mRun Complete\033[0m")
    print(f"  Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  https://github.com/migoldberg/PIRoC\n")

    sys.stdout = logger.terminal
    logger.close()

if __name__ == "__main__":
    main()
