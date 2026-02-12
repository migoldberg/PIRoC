# PhyloScreen 2026-01-03
# Mark Goldberg
# Kocot Lab, The University of Alabama

import os 
import re
import sys
import argparse
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

# ----------------------------
# Logging utility
# ---------------------------- 

class Logger:
    """
    Dual-output logger that writes to both stdout and a log file
    """
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log_file = open(log_path, 'w')
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log_file.write(f"PhyloScreen Log - {timestamp}\n")
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

# ----------------------------
# Metadata loading and parsing
# ---------------------------- 

def load_species_metadata(metadata_path):
    """
    Loads species to group mapping from metadata file (tested with .TSV)
    """
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

# ----------------------------
# Tree pre-processing
# ---------------------------- 

def collapse_low_support_nodes(tree, threshold=50.0):
    """
    Collapses internal nodes with bootstrap support below the provided threshold
    """
    nodes_collapsed = 0

    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue

        support = getattr(node, 'support', None)

        # normalize support value if a 0-1 range
        if support is not None and 0.0 <= support <= 1.0:
            support = support * 100.0
        
        if support is not None and support < threshold:
            node.delete()
            nodes_collapsed += 1
        
    return nodes_collapsed

# ----------------------------
# Long branch detection
# ---------------------------- 

def compute_branch_length_stats(tree):
    """
    Compute branch length statistics for a tree
    """
    lengths = []
    for node in tree.traverse():
        if node.dist is not None and node.dist > 0 and not node.is_root():
            lengths.append(node.dist)

    if not lengths:
        return None

    lengths_sorted = sorted(lengths)
    n = len(lengths_sorted)
    median = lengths_sorted[n // 2]
    mean = sum(lengths) / n

    # calculate standard deviation
    variance = sum((x - mean) ** 2 for x in lengths) / n
    std = variance ** 0.5

    # long branch threshold is median + 3*std OR 5x the median, whichever is smaller
    # this is conservative, may be adjusted in the future. 0.05 is set to prevent branches from being classified as long if they are very short.
    threshold = max(min(median + 3 * std, median * 5), 0.05)

    return {
        "median": median,
        "mean": mean,
        "std": std,
        "long_branch_threshold": threshold
    }

def check_long_branches(focal_leaves, branch_stats):
    """
    Checks if any focal species leaves are on a long branch
    """
    if branch_stats is None:
        return {
            "has_long_branch": False,
            "max_relative_length": None,
            "long_branch_leaves": []
        }
    
    threshold = branch_stats["long_branch_threshold"]
    median = branch_stats["median"]

    long_branch_leaves = []
    max_relative = 0.0

    for leaf in focal_leaves:
        if leaf.dist is not None and leaf.dist > 0:
            relative = leaf.dist / median if median > 0 else 0
            max_relative = max(max_relative, relative)

            if leaf.dist > threshold:
                long_branch_leaves.append(leaf.name)

    return {
        "has_long_branch": len(long_branch_leaves) > 0,
        "max_relative_length": max_relative,
        "long_branch_leaves": long_branch_leaves
    }

# ----------------------------
# Monophyly Check
# ---------------------------- 

def check_focal_monophyly(tree, focal_leaves, species_to_group, focal_group):
    """
    Checks if a focal species leaves form a monophyletic group
    """

    if len(focal_leaves) <= 1:
        # single leaf is trivially monophyletic
        return {
            "is_monophyletic": True,
            "is_group_monophyletic": True,
            "num_intruders": 0,
            "intruder_groups": Counter(),
            "intruder_species": Counter()
        }

    mrca = tree.get_common_ancestor(focal_leaves)
    mrca_leaf_names = set(mrca.get_leaf_names())
    focal_names = set(l.name for l in focal_leaves)

    intruder_names = mrca_leaf_names - focal_names

    intruder_groups = Counter()
    intruder_species = Counter()

    for name in intruder_names:
        species = parse_species_name(name)
        group = species_to_group.get(species, "UNKNOWN")
        intruder_species[species] += 1
        intruder_groups[group] += 1

    is_group_monophyletic = all(
        g == focal_group for g in intruder_groups.keys()
    ) if intruder_groups else True

    return {
        "is_monophyletic": len(intruder_names) == 0,
        "is_group_monophyletic": is_group_monophyletic,
        "num_intruders": len(intruder_names),
        "intruder_groups": intruder_groups,
        "intruder_species": intruder_species
    }

# ----------------------------
# Sister clade analysis
# ---------------------------- 

def analyze_sister_clades(focal_mrca, species_to_group, focal_group, outgroup_names=None):
    """
    Analyze the sister clade of the focal MRCA to find phylogenetic context
    If the sister group contains outgroup taxa, this is concerning
    If the sister group contains same-group taxa, this is good
    """

    if outgroup_names is None:
        outgroup_names = {"Outgroup"}

    parent = focal_mrca.up

    if parent is None:
        # focal_mrca is root, so there is no sister clade
        return {
            "sister_groups": Counter(),
            "sister_species": Counter(),
            "has_outgroup_sister": False,
            "has_same_group_sister": False,
            "sister_is_pure_outgroup": False
        }
    
    # if not root, then find the sister nodes
    sisters = [child for child in parent.children if child != focal_mrca]

    sister_groups = Counter()
    sister_species = Counter()

    for sister in sisters:
        for leaf in sister.get_leaves():
            if not leaf.name:
                continue
            species = parse_species_name(leaf.name)
            group = species_to_group.get(species, "UNKNOWN")
            sister_species[species] += 1
            sister_groups[group] += 1

    has_outgroup = any(g in outgroup_names for g in sister_groups.keys())
    has_same_group = focal_group in sister_groups

    # check if sister is purely outgroup
    sister_is_pure_outgroup = (
        len(sister_groups) > 0 and
        all(g in outgroup_names for g in sister_groups.keys())
    )

    return {
        "sister_groups": sister_groups,
        "sister_species": sister_species,
        "has_outgroup_sister": has_outgroup,
        "has_same_group_sister": has_same_group,
        "sister_is_pure_outgroup": sister_is_pure_outgroup
    }

# ----------------------------
# Bootstrap Support
# ---------------------------- 

def get_node_bootstrap(node, tree_string=None):
    """
    Get bootstrap support for a node
    """

    # Previously used tree_string as a fallback but extensive testing showed node.support is reliable enough

    if hasattr(node, "support") and node.support is not None:
        bs = float(node.support)
        # normalize
        if 0.0 <= bs <= 1.0:
            return bs * 100.0
        return bs
    else:
        return 0

# ----------------------------
# Clade composition
# ---------------------------- 

def compute_clade_composition(node, species_to_group):
    """
    Given a node (clade), finds the compositions of species/groups
    """
    species_counts = Counter()
    group_counts = Counter()

    for leaf in node.get_leaves():
        if not leaf.name:
            continue
        species_id = parse_species_name(leaf.name)
        species_counts[species_id] += 1

        group = species_to_group.get(species_id, "UNKNOWN")
        group_counts[group] += 1

    total_leaves = sum(species_counts.values())
    return species_counts, group_counts, total_leaves

# ----------------------------
# Classification logic
# ---------------------------- 

def classify_orthogroup(
    focal_species,
    focal_group,
    focal_leaves,
    mrca_node,
    tree,
    species_to_group,
    tree_string,
    branch_stats,
    min_support=70.0,
    min_target_purity=0.8,
    max_contaminant_purity=0.5,
    outgroup_names=None
):
    """
    Classifies an orthogroup using multiple lines of evidence:
        1. Bootstrap support of MRCA
        2. Purity of MRCA clade
        3. Monophyly of focal species
        4. Long branch detection
        5. Sister clade composition

    Classification logic:
        - CONTAMINANT: Strong evidence focal species doesn't belong
            * Sister clade is pure outgroup, OR
            * Low purity + high support, OR
            * Long branch + outgroup intruders
        - TARGET: Strong evidence focal species belongs
            * High purity + high support
            * No long branches
            * No outgroup in immediate sister
        - FLAG: Ambiguous or insufficient evidence
    """
    if outgroup_names is None:
        outgroup_names = {"Outgroup"}
    
    # finds metrics about the clade's composition
    species_counts, group_counts, total_leaves = compute_clade_composition(mrca_node, species_to_group)
    bootstrap = get_node_bootstrap(mrca_node, tree_string=tree_string)

    # calculates "purity" value
    if total_leaves == 0:
        purity = 0.0
    else:
        focal_group_count = group_counts.get(focal_group, 0)
        purity = focal_group_count / float(total_leaves)

    # check monophyly
    monophyly = check_focal_monophyly(tree, focal_leaves, species_to_group, focal_group)

    # check for long branch
    long_branch = check_long_branches(focal_leaves, branch_stats)

    # sister clade analysis
    sister = analyze_sister_clades(mrca_node, species_to_group, focal_group, outgroup_names)

    # count focal species leaves
    num_focal_species_leaves = sum(
        count for sp, count in species_counts.items()
        if species_to_group.get(sp, None) == focal_group
    )

    # check for outgroup intruders in the MRCA
    has_outgroup_intruders = any(
        g in outgroup_names for g in monophyly["intruder_groups"].keys()
    )

    # ------ Classificatoin ------

    classification = "FLAG" # default
    confidence_notes = []

    high_support = bootstrap is not None and bootstrap >= min_support

    # Contaminant signals ordered by priority 
    if sister["sister_is_pure_outgroup"] and high_support:
        classification = "CONTAMINANT"
        confidence_notes.append("sister_clade_pure_outgroup")
    elif has_outgroup_intruders and long_branch["has_long_branch"]:
        classification = "CONTAMINANT"
        confidence_notes.append("long_branch_with_outgroup_intruders")
    elif high_support and purity <= max_contaminant_purity:
        classification = "CONTAMINANT"
        confidence_notes.append("low_purity_high_support")

    # Target signals
    elif high_support and purity >= min_target_purity:
        if not long_branch["has_long_branch"] and not has_outgroup_intruders:
            classification = "TARGET"
            confidence_notes.append("high_purity_high_support")

            if monophyly["is_monophyletic"]:
                confidence_notes.append("focal_monophyletic")
            elif monophyly["is_group_monophyletic"]:
                confidence_notes.append("group_monophyletic")
        elif long_branch["has_long_branch"]:
            classification = "FLAG"
            confidence_notes.append("high_purity_but_long_branch")
        elif has_outgroup_intruders:
            classification = "FLAG"
            confidence_notes.append("high_purity_but_outgroup_intruders")

    # Ambigious cases for FLAG signal
    else:
        if bootstrap is None or bootstrap < min_support:
            confidence_notes.append("low_or_missing_support")
        if purity > max_contaminant_purity and purity < min_target_purity:
            confidence_notes.append("intermediate_purity")
        if long_branch["has_long_branch"]:
            confidence_notes.append("long_branch_detected")

    metrics = {
        # Basic metrics
        "bootstrap": bootstrap if bootstrap is not None else "NA",
        "total_leaves": total_leaves,
        "focal_group": focal_group,
        "focal_group_leaves": group_counts.get(focal_group, 0),
        "purity": purity,
        "num_focal_species_leaves": num_focal_species_leaves,
        "num_species": len(species_counts),
        "num_groups": len(group_counts),
        
        # Monophyly metrics
        "is_monophyletic": monophyly["is_monophyletic"],
        "is_group_monophyletic": monophyly["is_group_monophyletic"],
        "num_intruders": monophyly["num_intruders"],
        "intruder_groups": dict(monophyly["intruder_groups"]),
        
        # Long branch metrics
        "has_long_branch": long_branch["has_long_branch"],
        "max_relative_branch_length": long_branch["max_relative_length"],
        "long_branch_leaves": long_branch["long_branch_leaves"],
        
        # Sister clade metrics
        "sister_groups": dict(sister["sister_groups"]),
        "has_outgroup_sister": sister["has_outgroup_sister"],
        "has_same_group_sister": sister["has_same_group_sister"],
        "sister_is_pure_outgroup": sister["sister_is_pure_outgroup"],
        
        # Classification rationale
        "confidence_notes": confidence_notes
    }

    return classification, metrics

# ----------------------------
# Main
# ---------------------------- 


def main():
    parser = argparse.ArgumentParser(
        description="PhyloScreen: Phylogeny-based orthogroup classification (TARGET / CONTAMINANT / FLAG)"
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
        "-o", "--output_dir", default="phyloscreen_output",
        help="Output directory (default: phyloscreen_output)"
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
        "--outgroups", type=str, default="Outgroup",
        help="Comma-separated list of outgroup names (default: Outgroup)"
    )

    args = parser.parse_args()

    tree_dir = args.tree_dir
    tree_suffix = args.suffix
    output_dir = args.output_dir
    focal_species = args.focal_species
    metadata_path = args.metadata
    outgroup_names = set(args.outgroups.split(","))

    os.makedirs(output_dir, exist_ok=True)

    # Create OG_Classifications subfolder
    og_class_dir = os.path.join(output_dir, "OG_Classifications")
    os.makedirs(og_class_dir, exist_ok=True)

    # Set up dual logging (terminal + log file)
    log_path = os.path.join(output_dir, "log.txt")
    logger = Logger(log_path)
    sys.stdout = logger

    species_to_group = load_species_metadata(metadata_path)

    if not species_to_group:
        raise ValueError(f"No species metadata found in {metadata_path}")

    if focal_species not in species_to_group:
        raise ValueError(f"Focal species '{focal_species}' not found in metadata file")

    focal_group = species_to_group[focal_species]
    print(f"Focal species: {focal_species} (group: {focal_group})")
    print(f"Outgroup names: {outgroup_names}")
    print(f"Collapse threshold: {args.collapse_threshold}")
    print()

    classifications = {}
    metrics_by_og = {}

    no_focal_count = 0
    error_count = 0
    total_trees = 0
    total_nodes_collapsed = 0

    for filename in os.listdir(tree_dir):
        if not filename.endswith(tree_suffix):
            continue

        total_trees += 1
        tree_path = os.path.join(tree_dir, filename)
        og_id = filename.replace(tree_suffix, "")

        try:
            with open(tree_path, 'r') as f:
                tree_string = f.read().strip()

            t = Tree(tree_string, format=0)

            # find outgroup leaves
            outgroup_leaves = [
                l for l in t.get_leaves() 
                if species_to_group.get(parse_species_name(l.name)) in outgroup_names
            ]
            
            # try to root the tree with the outgroup leaves
            if outgroup_leaves:
                try:
                    # root at a outgroup leaf
                    t.set_outgroup(outgroup_leaves[0])
                except Exception as e:
                    # if failed, root at the midpoint
                    print(f"[{og_id}] Could not root on outgroup, using midpoint: {e}")
                    t.midpoint_root()
            else:
                # if no outgroup leaves, roots at the midpoint
                t.midpoint_root()

            # collapse low-support nodes
            nodes_collapsed = collapse_low_support_nodes(t, threshold=args.collapse_threshold)
            total_nodes_collapsed += nodes_collapsed

            # compute branch length statistics for this tree
            branch_stats = compute_branch_length_stats(t)

            # find leaves from focal species
            focal_leaves = [
                leaf for leaf in t.get_leaves()
                if leaf.name and parse_species_name(leaf.name) == focal_species
            ]

            if not focal_leaves:
                no_focal_count += 1
                continue

            # get MRCA of focal leaves
            if len(focal_leaves) == 1:
                if focal_leaves[0].up is None:
                    error_count += 1
                    continue
                mrca = focal_leaves[0].up
            else:
                mrca = t.get_common_ancestor(focal_leaves)

            classification, metrics = classify_orthogroup(
                focal_species=focal_species,
                focal_group=focal_group,
                focal_leaves=focal_leaves,
                mrca_node=mrca,
                tree=t,
                species_to_group=species_to_group,
                tree_string=tree_string,
                branch_stats=branch_stats,
                min_support=args.min_support,
                min_target_purity=args.min_target_purity,
                max_contaminant_purity=args.max_contaminant_purity,
                outgroup_names=outgroup_names
            )

            classifications[og_id] = classification
            metrics_by_og[og_id] = metrics

            # compact output
            notes = ",".join(metrics["confidence_notes"]) if metrics["confidence_notes"] else "-"
            print(
                f"[{og_id}] {classification:12s} | "
                f"bs={metrics['bootstrap']:>5} | "
                f"purity={metrics['purity']:.2f} | "
                f"mono={metrics['is_monophyletic']} | "
                f"longbr={metrics['has_long_branch']} | "
                f"notes={notes}"
            )

        except Exception as e:
            print(f"[{og_id}] Error: {e}")
            error_count += 1

    # Write summary
    summary_file = os.path.join(output_dir, "classification_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("PhyloScreen Classification Summary\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Parameters:\n")
        f.write(f"  Focal species:       {focal_species} ({focal_group})\n")
        f.write(f"  Min support:         {args.min_support}\n")
        f.write(f"  Min target purity:   {args.min_target_purity}\n")
        f.write(f"  Max contam purity:   {args.max_contaminant_purity}\n")
        f.write(f"  Collapse threshold:  {args.collapse_threshold}\n")
        f.write(f"  Outgroups:           {', '.join(outgroup_names)}\n\n")
        
        f.write("Processing Stats:\n")
        f.write(f"  Total trees:           {total_trees}\n")
        f.write(f"  Successfully classified: {len(classifications)}\n")
        f.write(f"  No focal species:      {no_focal_count}\n")
        f.write(f"  Errors:                {error_count}\n")
        f.write(f"  Total nodes collapsed: {total_nodes_collapsed}\n\n")

        if classifications:
            counts = Counter(classifications.values())
            f.write("Classification Breakdown:\n")
            f.write("-" * 30 + "\n")
            for label in ["TARGET", "CONTAMINANT", "FLAG"]:
                count = counts.get(label, 0)
                pct = (count / len(classifications)) * 100.0
                f.write(f"  {label:12s}: {count:5d} ({pct:5.1f}%)\n")

    print(f"\nSummary written to {summary_file}")

    # Write detailed results
    results_file = os.path.join(output_dir, "classification_results.tsv")
    with open(results_file, 'w') as f:
        header = [
            "OG_ID", "Classification", "Bootstrap", "Purity", "Total_Leaves",
            "Focal_Group", "Focal_Group_Leaves",
            "Is_Monophyletic", "Is_Group_Monophyletic", "Num_Intruders",
            "Has_Long_Branch", "Max_Relative_Branch",
            "Has_Outgroup_Sister", "Sister_Pure_Outgroup",
            "Confidence_Notes"
        ]
        f.write("\t".join(header) + "\n")
        
        for og_id in sorted(classifications.keys()):
            m = metrics_by_og[og_id]
            row = [
                og_id,
                classifications[og_id],
                str(m["bootstrap"]),
                f"{m['purity']:.4f}",
                str(m["total_leaves"]),
                m["focal_group"],
                str(m["focal_group_leaves"]),
                str(m["is_monophyletic"]),
                str(m["is_group_monophyletic"]),
                str(m["num_intruders"]),
                str(m["has_long_branch"]),
                f"{m['max_relative_branch_length']:.2f}" if m["max_relative_branch_length"] else "NA",
                str(m["has_outgroup_sister"]),
                str(m["sister_is_pure_outgroup"]),
                ";".join(m["confidence_notes"])
            ]
            f.write("\t".join(row) + "\n")

    print(f"Detailed results written to {results_file}")

    # Write per-classification lists for easy downstream filtering
    for label in ["TARGET", "CONTAMINANT", "FLAG"]:
        list_file = os.path.join(og_class_dir, f"{label.lower()}_orthogroups.txt")
        with open(list_file, 'w') as f:
            for og_id in sorted(classifications.keys()):
                if classifications[og_id] == label:
                    f.write(og_id + "\n")
        print(f"{label} list written to {list_file}")

    # Restore stdout and close logger
    sys.stdout = logger.terminal
    logger.close()
    print(f"Log written to {log_path}")


if __name__ == "__main__":
    main()
