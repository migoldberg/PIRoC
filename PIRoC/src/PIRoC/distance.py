"""
distance.py

"""

from collections import defaultdict
import statistics
import math
from .metadata import parse_species_name

def get_nearest_contaminants(tree, contaminant_leaves):
    contaminants = set(contaminant_leaves)
    no_contaminant = (float('inf'), None)

    below_node = {}
    for node in tree.traverse("postorder"):
        if node in contaminants:
            below_node[node] = (0.0, node)
        elif node.is_leaf():
            below_node[node] = no_contaminant
        else:
            below_node[node] = min(
                ((below_node[child][0] + child.dist, below_node[child][1]) for child in node.children),
                key=lambda x: x[0]
            )

    nearest_node = {}
    for node in tree.traverse("preorder"):
        if node.up is None:
            nearest_node[node] = below_node[node]
        else:
           via_parent = (nearest_node[node.up][0] + node.dist, nearest_node[node.up][1])
           nearest_node[node] = min(below_node[node], via_parent, key=lambda x: x[0])

    return nearest_node


def get_nearby_nodes(tree):
    nodes = defaultdict(list)
    for node in tree.traverse():
        if node.up is not None:
            nodes[node].append((node.up, node.dist))
            nodes[node.up].append((node, node.dist))
    return nodes

def get_distance_from_leaf(leaf, nearby_nodes):
    distances = {leaf: 0.0}
    order = [leaf]
    while order:
        node = order.pop()
        for nearby_node, branch_length in nearby_nodes[node]:
            if nearby_node not in distances:
                distances[nearby_node] = distances[node] + branch_length
                order.append(nearby_node)
    return distances

def get_nearest_other_species_leaves(tree, focal_leaves):
    nearby_nodes = get_nearby_nodes(tree)
    leaves = tree.get_leaves()
    leaf_species = {leaf: parse_species_name(leaf.name) for leaf in leaves}

    nearest = {}
    for focal_leaf in focal_leaves:
        distances = get_distance_from_leaf(focal_leaf, nearby_nodes)
        nearest[focal_leaf] = min(
            ((distances[leaf], leaf) for leaf in leaves if leaf_species[leaf] != leaf_species[focal_leaf]),
            key=lambda x: x[0],
            default=(float('inf'), None)
        )
    return nearest

def get_average_distance(tree, focal_leaves, contaminant_leaves):
    sources, targets = (
        (contaminant_leaves, focal_leaves) if len(contaminant_leaves) <= len(focal_leaves) else (focal_leaves, contaminant_leaves)
    )

    nearby_nodes = get_nearby_nodes(tree)

    pairwise_distances = []
    for source in sources:
        distances = get_distance_from_leaf(source, nearby_nodes)
        pairwise_distances.extend(distances[target] for target in targets)

    return statistics.median(pairwise_distances)

def compute_distance_rule(tree, focal_leaves, contaminant_leaves, max_contaminant_distance, min_average_distance, max_shared_fraction, n_focal_species):
    if max_contaminant_distance is None:
        return {leaf: {"distance_rule": "skipped"} for leaf in focal_leaves}

    if not contaminant_leaves:
        return {leaf: {"distance_rule": "no_contaminants_in_tree"} for leaf in focal_leaves}

    nearest_contaminants = get_nearest_contaminants(tree, contaminant_leaves)
    average_distance = get_average_distance(tree, focal_leaves, contaminant_leaves)

    nearby_focal_species = {
        parse_species_name(leaf.name) for leaf in focal_leaves
        if nearest_contaminants[leaf][0] <= max_contaminant_distance
    }

    max_shared_species = max(3, math.ceil(max_shared_fraction * n_focal_species))

    if average_distance < min_average_distance:
        distance_rule = "skipped_low_average"
    elif len(nearby_focal_species) >= max_shared_species:
        distance_rule = "skipped_conserved_gene"
    else:
        distance_rule = "applied"

    return {
        leaf: {
            "distance_rule": distance_rule,
            "max_contaminant_distance": max_contaminant_distance,
            "nearest_contaminant_distance": nearest_contaminants[leaf][0],
            "nearest_contaminant": nearest_contaminants[leaf][1].name,
            "average_distance": average_distance,
            "n_nearby_focal_species": len(nearby_focal_species)
        }
        for leaf in focal_leaves
    }

