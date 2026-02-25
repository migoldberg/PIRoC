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