"""
metadata.py
Provides the functions and helper functions for loading and using the species metadata file.
"""

import os 
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree


def load_species_metadata(metadata_path):
    species_to_group = {}
    with open(metadata_path, 'r') as f:
        header = f.readline().strip().split('\t')
        try:
            sp_idx = header.index("species_id")
            grp_idx = header.index("group")
        except ValueError:
            raise ValueError("Metadata file must have 'species_id' and 'group' columns")

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) <= max(sp_idx, grp_idx):
                continue
            species_id = parts[sp_idx]
            group = parts[grp_idx]
            if species_id:
                species_to_group[species_id] = group

    return species_to_group


def parse_species_name(leaf_name):
    """
    Extract species ID from a leaf name of form 'species|gene'.
    """
    return leaf_name.split("|")[0]


def get_group_from_species_name(leaf_name, species_to_group):
    species_name = parse_species_name(leaf_name)
    return species_to_group.get(species_name, "UNKNOWN")