"""
cli.py
Provides the command line interface (CLI) functions.
"""

import os 
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

from .metadata import load_species_metadata, parse_species_name

class Logger:
    """
    Writes output to a log file and prints
    """
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log_file = open(log_path, 'w')
        self.log_file.flush()

    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()

    def flush(self):
        self.terminal.flush()
        self.log_file.flush()

    def close(self):
        self.log_file.close()

def init_cli() -> None:
    """
    Returns/Loads all informatoin from the CLI for the run.
    """
    parser = argparse.ArgumentParser(
        description="PIRoC: Phylogeny-based orthogroup classification (TARGET / CONTAMINANT / FLAG)"
    )
    parser.add_argument(
        "-t", "--tree_dir", required=True,
        help="Directory containing Newick tree files"
    )
    parser.add_argument(
        "-s", "--suffix", default=".tre",
        help="Tree filename suffix (default: .tre)"
    )
    parser.add_argument(
        "-o", "--output_dir", default="PIRoC_output",
        help="Output directory (default: PIRoC_output)"
    )
    parser.add_argument(
        "-f", "--focal_species", required=True,
        help="Species ID of the focal species (e.g., obim)"
    )
    parser.add_argument(
        "-m", "--metadata", required=True,
        help="TSV file with species_id and group columns"
    )
    parser.add_argument(
        "--contaminants", type=str, default="Contaminant",
        help="Comma-separated list of known contaminant group names (default: Contaminant)"
    )
    parser.add_argument(
        "-rm", "--remove_contaminants", action="store_true",
        help="After classification remove contaminant sequences from the orthogroups."
    )
    parser.add_argument(
        "--min_support", type=float, default=70.0,
        help="Minimum bootstrap support for confident classification (default: 70)"
    )
    parser.add_argument(
        "--min_target_purity", type=float, default=0.8,
        help="Minimum clade target group fraction for TARGET classification (default: 0.8)"
    )
    parser.add_argument(
        "--max_contaminant_purity", type=float, default=0.5,
        help="Maximum clade target group fraction for CONTAMINANT classification (default: 0.5)"
    )
    parser.add_argument(
        "--collapse_threshold", type=float, default=50.0,
        help="Collapse nodes with bootstrap support below this value (default: 50)"
    )
    parser.add_argument(
        "--review-flags", action="store_true",
        help="At the end of the run, review the flagged sequences via CLI and make final decisions"
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Enable quiet mode (minimul output to the console)"
    )

    # parses the arguments
    args = parser.parse_args()

    tree_dir = args.tree_dir
    tree_suffix = args.suffix
    output_dir = args.output_dir
    focal_species = args.focal_species
    metadata_path = args.metadata
    contaminants = set(args.contaminants.split(","))
    remove_contaminants = args.remove_contaminants
    quiet = args.quiet
    review_flags = args.review_flags

    # creates the output directory
    os.makedirs(output_dir, exist_ok=True)

    # creates the directory for sequence classifications
    sequence_classifications_dir = os.path.join(output_dir, "sequence_classifications")
    os.makedirs(sequence_classifications_dir, exist_ok=True)

    # creates the log file
    log_path = os.path.join(output_dir, "PIRoC.log")
    logger = Logger(log_path)
    sys.stdout = logger

    # loads the species to group dictionary from the metadata file
    species_to_group = load_species_metadata(metadata_path)

    # if the species to group dictionary is empty or improperly formatted then an error is raised
    if not species_to_group:
        raise ValueError("Metadata file is empty or improperly formatted.")

    # if the focal species provided in the arguments is not found in the species to group dictionary then an error is raised
    if focal_species not in species_to_group:
        raise ValueError(f"Focal species '{focal_species}' not found in metadata.")

    # gets the focal group for the focal species
    focal_group = species_to_group[focal_species]

    # prints quick stats before the run begins
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print()
    print("\033[4mPIRoC - Phylogeny-Informed Removal of Contaminants\033[0m")
    print(f"Version 1.0.0")
    print(f"Run Date: {timestamp}\n")

    print("\033[4mRun Parameters\033[0m")
    print("  [~] required   [ ] default   [X] user set\n")
    params_display = [
        ("tree_dir",             tree_dir,                    True,  None),
        ("suffix",               tree_suffix,                 False, ".tre"),
        ("output_dir",           output_dir,                  False, "PIRoC_output"),
        ("focal_species",        focal_species,               True,  None),
        ("focal_group",          focal_group,                 True,  None),
        ("metadata",             metadata_path,               True,  None),
        ("contaminants",         args.contaminants,           False, "Contaminant"),
        ("remove_contaminants",  remove_contaminants,         False, False),
        ("min_support",          args.min_support,            False, 70.0),
        ("min_target_purity",    args.min_target_purity,      False, 0.8),
        ("max_contaminant_purity", args.max_contaminant_purity, False, 0.5),
        ("collapse_threshold",   args.collapse_threshold,     False, 50.0),
        ("quiet",                quiet,                       False, False),
        ("review_flags",         review_flags,                False, False),
    ]
    max_name_len = max(len(name) for name, *_ in params_display)
    for name, value, required, default in params_display:
        if required:
            marker = "~"
        elif value == default:
            marker = " "
        else:
            marker = "X"
        print(f"  [{marker}] {name:<{max_name_len}}  {value}")
        
    print()
    print("\033[4mRun Progress\033[0m")

    # returns a dictionary of the arguments and information for the run
    return {
        "tree_dir": tree_dir,
        "tree_suffix": tree_suffix,
        "output_dir": output_dir,
        "sequence_classifications_dir": sequence_classifications_dir,
        "focal_species": focal_species,
        "focal_group": focal_group,
        "species_to_group": species_to_group,
        "min_support": args.min_support,
        "min_target_purity": args.min_target_purity,
        "max_contaminant_purity": args.max_contaminant_purity,
        "collapse_threshold": args.collapse_threshold,
        "contaminants": contaminants,
        "logger": logger,
        "remove_contaminants": remove_contaminants,
        "quiet": quiet,
        "review_flags": review_flags,
    }