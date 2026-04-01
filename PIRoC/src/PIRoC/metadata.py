"""
metadata.py
Provides the functions and helper functions for loading and using the species metadata file.
"""

import os 
import re
import sys
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

def load_species_metadata(metadata_path):
    """
    Loads the species metadata file and returns a dictionary of species to group.
    """

    # initializes an empty dictionary to store the species to group mapping
    species_to_group = {}

    # opens the metadata file
    with open(metadata_path, 'r') as f:
        # reads the header line and splits it into a list of column names
        header = f.readline().strip().split('\t')

        # creates an index for the species_id and group columns
        try:
            sp_idx = header.index("species_id")
            grp_idx = header.index("group")
        except ValueError:
            raise ValueError("Metadata file must have 'species_id' and 'group' columns")

        # loops through each line in the metadata file
        for line in f:
            # skips empty lines
            if not line.strip():
                continue

            # splits the line into columns
            parts = line.strip().split('\t')

            # skips lines that don't have enough columns
            if len(parts) <= max(sp_idx, grp_idx):
                continue

            # gets the species id and group from the columns
            species_id = parts[sp_idx]
            group = parts[grp_idx]

            # adds the species id and group to the dictionary
            if species_id:
                species_to_group[species_id] = group

    # returns the species to group dictionary
    return species_to_group


def parse_species_name(leaf_name):
    """
    Extract species ID from a leaf name of form 'species|gene'.
    """
    return leaf_name.split("|")[0]


def get_group_from_species_name(leaf_name, species_to_group):
    """
    Gets the group from the species name using the species to group dictionary.
    """
    species_name = parse_species_name(leaf_name)
    return species_to_group.get(species_name, "UNKNOWN")

def verify_metadata(species_to_group, tree_dir, tree_suffix, focal_group):
    # if the species to group dictionary is empty or improperly formatted then an error is raised
    if not species_to_group:
        raise ValueError("Metadata file is empty or improperly formatted.")

    # if the focal group provided in the arguments is not found in the species to group dictionary then an error is raised
    if focal_group not in species_to_group.values():
        raise ValueError(f"Focal group '{focal_group}' not found in metadata.")

    missing_species = []

    for tree_file in os.listdir(tree_dir):
        if tree_file.endswith(tree_suffix):
            tree = Tree(os.path.join(tree_dir, tree_file))
            for leaf in tree.get_leaves():
                if parse_species_name(leaf.name) not in species_to_group:
                    missing_species.append(parse_species_name(leaf.name))

    if missing_species:
        unique_missing_species = sorted(set(missing_species))
        print(f"WARNING: {len(unique_missing_species)} species not found in metadata:")
        for species in unique_missing_species:
            print(f"  - {species}")