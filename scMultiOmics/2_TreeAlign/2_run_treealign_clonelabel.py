#!/usr/bin/env python3
"""Run TreeAlign's clone-label mode: assign scRNA cells directly to the 18
existing scDNA k-means clusters (no phylogenetic tree needed).

Input (from 0_prepare_scdna_cn.R / 1_prepare_scrna_expr.R):
    <IN_DIR>/gene_by_cell_cnv.csv.gz    gene x cell scDNA copy number
    <IN_DIR>/cell_clone_labels.csv      scDNA cell_id -> clone_id
    <IN_DIR>/expr_counts.mtx.gz + expr_genes.tsv.gz + expr_cells.tsv.gz
                                        gene x cell raw scRNA counts

Output (under OUT_DIR, see config below):
    clone_assign_df.csv       scRNA cell_id -> assigned clone_id (+ probability)
    gene_type_score_df.csv    per-gene CN-dosage score (0-1)
    allele_assign_prob_df.csv empty (allele-specific submodel not used, see README)

Requires: the `treealign` micromamba env (see scMultiOmics/2_TreeAlign/README.md
for exact creation commands — torch==1.11.0, pyro-ppl==1.8.1, TreeAlign itself
cloned from github.com/AlexHelloWorld/TreeAlign).
Usage:
    micromamba run -n treealign python scMultiOmics/2_TreeAlign/2_run_treealign_clonelabel.py
"""

import gzip
import os
import time

import pandas as pd
import pyro
import scipy.io
from treealign import CloneAlignClone

# ── Configuration ─────────────────────────────────────────────────────────────

MODE = os.environ.get("TREEALIGN_MODE", "test")  # "test" or "full"
assert MODE in ("test", "full")

SEED = 0
REPEAT = 1 if MODE == "test" else 10  # library default is 10; 1 for a fast smoke test

IN_DIR = f"/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/inputs/{MODE}"
OUT_DIR = f"/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/outputs/{MODE}/clonelabel"
RESULTS_DIR = "scMultiOmics/2_TreeAlign/results"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_sparse_expr(mtx_gz_path, genes_gz_path, cells_gz_path):
    """gene x cell raw counts, as a dense pandas.DataFrame (what TreeAlign expects)."""
    with gzip.open(mtx_gz_path, "rt") as f:
        mat = scipy.io.mmread(f).tocsr()
    genes = pd.read_csv(genes_gz_path, header=None)[0].tolist()
    cells = pd.read_csv(cells_gz_path, header=None)[0].tolist()
    return pd.DataFrame(mat.toarray(), index=genes, columns=cells)


print(f"[2_run_treealign_clonelabel.py] MODE={MODE}, repeat={REPEAT}, input={IN_DIR}")

pyro.set_rng_seed(SEED)

print("Loading inputs...")
cnv = pd.read_csv(os.path.join(IN_DIR, "gene_by_cell_cnv.csv.gz"), index_col=0)
clone = pd.read_csv(os.path.join(IN_DIR, "cell_clone_labels.csv"))
expr = load_sparse_expr(
    os.path.join(IN_DIR, "expr_counts.mtx.gz"),
    os.path.join(IN_DIR, "expr_genes.tsv.gz"),
    os.path.join(IN_DIR, "expr_cells.tsv.gz"),
)
print(f"cnv: {cnv.shape[0]} genes x {cnv.shape[1]} scDNA cells")
print(f"clone: {clone.shape[0]} scDNA cells across {clone['clone_id'].nunique()} clones")
print(f"expr: {expr.shape[0]} genes x {expr.shape[1]} scRNA cells")

print("Running CloneAlignClone...")
start = time.time()
obj = CloneAlignClone(clone=clone, expr=expr, cnv=cnv, repeat=REPEAT)
obj.assign_cells_to_clones()
clone_assign_df, gene_type_score_df, allele_assign_prob_df = obj.generate_output()
print(f"Done in {time.time() - start:.1f}s")

clone_assign_df.to_csv(os.path.join(OUT_DIR, "clone_assign_df.csv"), index=False)
gene_type_score_df.to_csv(os.path.join(OUT_DIR, "gene_type_score_df.csv"), index=False)
allele_assign_prob_df.to_csv(os.path.join(OUT_DIR, "allele_assign_prob_df.csv"), index=False)

qc = pd.DataFrame([{
    "mode": MODE,
    "n_scrna_cells": expr.shape[1],
    "n_scdna_cells": cnv.shape[1],
    "n_clones": clone["clone_id"].nunique(),
    "n_genes": cnv.shape[0],
    "n_assigned": clone_assign_df["clone_id"].notna().sum(),
    "n_unassigned": clone_assign_df["clone_id"].isna().sum(),
    "runtime_sec": round(time.time() - start, 1),
}])
qc.to_csv(os.path.join(RESULTS_DIR, f"2_run_treealign_clonelabel_qc_{MODE}.tsv"), sep="\t", index=False)
print(qc.to_string(index=False))
