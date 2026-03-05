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
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log_file.write(f"PIRoC Log - {timestamp}\n")
        self.log_file.write("=" * 50 + "\n\n")
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
        help="Comma-separated list of outgroup names (default: Outgroup)"
    )
    parser.add_argument(
        "-rm", "--remove_contaminants", action="store_true",
        help="Produce new tree files with contaminants removed"
    )
    parser.add_argument(
        "--min_support", type=float, default=70.0,
        help="Minimum bootstrap support for confident classification (default: 70)"
    )
    parser.add_argument(
        "--min_target_purity", type=float, default=0.8,
        help="Minimum focal-group purity for TARGET (default: 0.8)"
    )
    parser.add_argument(
        "--max_contaminant_purity", type=float, default=0.5,
        help="Maximum focal-group purity for CONTAMINANT (default: 0.5)"
    )
    parser.add_argument(
        "--collapse_threshold", type=float, default=50.0,
        help="Collapse nodes with bootstrap below this value (default: 50)"
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable debug mode"
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
    debug = args.debug

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
    print(f"Focal species: {focal_species} (Group: {focal_group})")
    print(f"Contaminant group(s) names: {contaminants}")
    print(f"Collapse Threshold: {args.collapse_threshold}")
    print()

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
        "debug": debug,
    }