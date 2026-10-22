"""
output.py
Provides the functions for writing the summary and sequence-level classification results to output files.
"""

import os 
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict

def write_summary(summary_file, focal_group, min_support, min_clean_purity, max_contaminant_purity, contaminant_names, run_metrics, sequence_classifications, distance_params):
    with open(summary_file, 'w') as f:
        f.write("PIRoC Classification Summary\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Parameters:\n")
        f.write(f"  Focal group:       {focal_group}\n")
        f.write(f"  Min support:         {min_support}\n")
        f.write(f"  Min clean purity:   {min_clean_purity}\n")
        f.write(f"  Max contam purity:   {max_contaminant_purity}\n")
        f.write(f"  Contaminants:        {', '.join(contaminant_names)}\n\n")
        f.write(f"  Max contam distance: {distance_params['max_contaminant_distance'] if distance_params['max_contaminant_distance'] is not None else 'off'}\n")

        if distance_params['max_contaminant_distance'] is not None:
            f.write(f"  Min background dist: {distance_params['min_average_distance']}\n")
            f.write(f"  Max shared fraction: {distance_params['max_shared_fraction']}\n")
            f.write("\n")

        f.write(f"  Resolve low support: {distance_params['resolve_low_support']}\n\n")
        
        f.write("Processing Stats:\n")
        f.write(f"  Total trees:             {run_metrics['total_trees']}\n")
        f.write(f"  Trees with no focal group:        {run_metrics['no_focal_group_in_tree']}\n")
        f.write(f"  Errors:                  {run_metrics['total_errors']}\n")
        f.write(f"  Total sequences classified:  {run_metrics['total_sequences_classified']}\n\n")

        if sequence_classifications:
            sequence_counts = Counter(sequence_classifications.values())
            f.write("Classification Summary:\n")
            f.write("-" * 40 + "\n")
            for label in ["CLEAN", "CONTAMINANT", "FLAG"]:
                count = sequence_counts.get(label, 0)
                pct = (count / len(sequence_classifications)) * 100.0
                f.write(f"  {label:12s}: {count:5d} ({pct:5.1f}%)\n")
            f.write("\n")

    print(f"  Summary written to: {summary_file}")

def format_distance_helper(value):
    if value is None or value == float('inf'):
        return "NA"
    return f"{value:.5f}"

def write_sequence_classifications(sequence_classifications_file, sequence_classifications, sequence_metrics):
    with open(sequence_classifications_file, 'w') as f:
        header = [
            "OG_ID", "Sequence_Name",
            "Classification", "Confidence_Notes",
            "Bootstrap", "Is_Long_Branch",
            "Total_Leaves", "n_Species_In_Clade",
            "Focal_Group", "n_Focal_Group_Leaves_In_Clade", "Has_Other_Focal_Group_In_Clade",
            "n_Contaminant_Leaves_In_Clade", "Has_Contaminant_In_Clade",
            "Clade_Target_Group_Fraction", "Clade_Group_Composition",
            "Sister_Total_Leaves", "Sister_Contains_Focal_Group",
            "Sister_Contains_Contaminants", "Sister_Is_Pure_Contaminant",
            "Sister_Group_Composition",
            "Pre_Distance_Classification","Distance_Rule_Status", "Nearest_Contaminant_Distance", "Nearest_Contaminant", "Average_Contaminant_Distance",
            "Nearest_Other_Species", "Nearest_Other_Species_Distance",
        ]
        f.write("\t".join(header) + "\n")
        
        for sequence_key in sorted(sequence_classifications.keys()):
            og_id, sequence_name = sequence_key.split("::", 1)
            m = sequence_metrics[sequence_key]
            row = [
                og_id,
                sequence_name,
                sequence_classifications[sequence_key],
                ";".join(m["classification_notes"]),
                str(m["bootstrap"]),
                str(m["is_on_long_branch"]),
                str(m["total_leaves"]),
                str(m["n_species_in_clade"]),
                str(m["focal_group"]),
                str(m["focal_group_leaves"]),
                str(m["has_other_focal_group_leaves_in_clade"]),
                str(m["n_contaminant_leaves_in_clade"]),
                str(m["has_contaminant_in_clade"]),
                f"{m['clade_target_group_fraction']:.4f}",
                str(m["clade_group_composition"]),
                str(m["sister_total_leaves"]),
                str(m["sister_contains_focal_group"]),
                str(m["sister_contains_contaminants"]),
                str(m["sister_is_pure_contaminant"]),
                str(m["sister_group_composition"]),
                str(m["pre_distance_classification"]),
                str(m["distance_rule_status"]),
                format_distance_helper(m["nearest_contaminant_distance"]),
                str(m["nearest_contaminant"] or "NA"),
                format_distance_helper(m["average_distance"]),
                str(m["nearest_neighbor"] or "NA"),
                format_distance_helper(m["nearest_neighbor_distance"]),
            ]
            f.write("\t".join(row) + "\n")

    print(f"  Classification results written to: {sequence_classifications_file}\n")

def write_sequence_lists(sequence_lists_dir, sequence_classifications):
    os.makedirs(sequence_lists_dir, exist_ok=True)

    for label in ["CLEAN", "CONTAMINANT", "FLAG"]:
        list_file = os.path.join(sequence_lists_dir, f"{label.lower()}_sequences.txt")
        with open(list_file, 'w') as f:
            for sequence_key in sorted(sequence_classifications.keys()):
                if sequence_classifications[sequence_key] == label:
                    og_id, sequence_name = sequence_key.split("::", 1)
                    f.write(f"{og_id}\t{sequence_name}\n")