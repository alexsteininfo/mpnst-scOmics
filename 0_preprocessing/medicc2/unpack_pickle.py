#!/usr/bin/env python3
"""
Unpack MEDICC2 bootstrap tree pickles into multi-tree Newick files readable by ape.

Each pickle holds a DataFrame with unique tree topologies and their counts across
bootstrap replicates. This script expands by count so that the output Newick file
contains one tree per replicate (e.g. 100 lines for 100 bootstrap replicates),
which is the format expected by ape::read.tree() / prop.clades().

Usage (with micromamba python env):
    micromamba run -n python python unpack_pickle.py
"""

import io
import os
import pickle
import glob

from Bio import Phylo

INPUT_DIR  = "/srv/home/aste0033/projects/CNA-PhyloAnalysis/MEDICC2/10X_DLP/per_cluster_bootstrap/minseg_1.5"
OUTPUT_DIR = "/srv/home/aste0033/projects/CNA-PhyloAnalysis/MEDICC2/10X_DLP/per_cluster_bootstrap/unpacked_minseg_1.5"

os.makedirs(OUTPUT_DIR, exist_ok=True)

for sample_dir in sorted(os.listdir(INPUT_DIR)):
    sample_path = os.path.join(INPUT_DIR, sample_dir)
    if not os.path.isdir(sample_path):
        continue

    pickles = glob.glob(os.path.join(sample_path, "*_bootstrap_trees_df.pickle"))
    if not pickles:
        continue

    with open(pickles[0], "rb") as f:
        df = pickle.load(f)

    # Expand unique topologies by their replicate count
    newick_lines = []
    for _, row in df.iterrows():
        buf = io.StringIO()
        Phylo.write(row["tree"], buf, "newick")
        newick = buf.getvalue().strip()
        newick_lines.extend([newick] * int(row["count"]))

    out_path = os.path.join(OUTPUT_DIR, f"{sample_dir}_bootstrap_trees.new")
    with open(out_path, "w") as f:
        f.write("\n".join(newick_lines) + "\n")

    print(f"{sample_dir}: {len(newick_lines)} trees written -> {out_path}")
