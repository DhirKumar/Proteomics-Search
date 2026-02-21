#!/usr/bin/env python3
"""
benchmark/compare_results.py
─────────────────────────────
Compare pipeline PSM output against a known answer key to evaluate:
  • Sensitivity   — % of true PSMs recovered
  • Precision     — % of our PSMs that are correct
  • F1 score      — harmonic mean of sensitivity and precision
  • FDR accuracy  — actual FDR vs nominal FDR threshold
  • Speed         — spectra per second (from pipeline log if available)

Usage
─────
  python benchmark/compare_results.py \\
      --our-psms   pipeline_output/identified_psms.csv \\
      --answer-key benchmark/answer_key.csv \\
      --key-col    Peptide \\
      --our-col    Peptide

Answer key format
─────────────────
  A CSV file with at least one column containing correct peptide sequences.
  The iPRG2012 answer key can be downloaded from:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3764048/

Example with iPRG2012
─────────────────────
  python benchmark/compare_results.py \\
      --our-psms   pipeline_output/identified_psms.csv \\
      --answer-key benchmark/iprg2012_correct_psms.csv \\
      --key-col    peptide \\
      --our-col    Peptide
"""

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
    import numpy as np
except ImportError:
    sys.exit("Install pandas and numpy:  pip install pandas numpy")


# ─────────────────────────────────────────────────────────────────────────────
def load_csv(path, label):
    p = Path(path)
    if not p.exists():
        sys.exit(f"[ERROR] {label} file not found: {path}")
    df = pd.read_csv(p)
    print(f"  Loaded {label}: {len(df):,} rows  ({p.name})")
    return df


def normalise_peptide(seq):
    """Strip modifications and convert to uppercase for matching."""
    import re
    return re.sub(r"[^A-Z]", "", str(seq).upper())


# ─────────────────────────────────────────────────────────────────────────────
def run_benchmark(our_psms_path, answer_key_path,
                  our_col, key_col, score_col, fdr_nominal):

    print("\n" + "═" * 60)
    print("  MS/MS PIPELINE BENCHMARK")
    print("═" * 60)

    df_our = load_csv(our_psms_path,   "Our PSMs")
    df_key = load_csv(answer_key_path, "Answer Key")

    # validate columns
    for col, df, label in [(our_col, df_our, "our PSMs"),
                           (key_col, df_key, "answer key")]:
        if col not in df.columns:
            available = list(df.columns)
            sys.exit(f"[ERROR] Column '{col}' not in {label}. "
                     f"Available: {available}")

    # normalise sequences
    df_our["_pep_norm"] = df_our[our_col].apply(normalise_peptide)
    df_key["_pep_norm"] = df_key[key_col].apply(normalise_peptide)

    our_set = set(df_our["_pep_norm"])
    key_set = set(df_key["_pep_norm"])

    # ── core metrics ─────────────────────────────────────────────────────────
    true_positives  = our_set & key_set
    false_positives = our_set - key_set
    false_negatives = key_set - our_set

    n_tp = len(true_positives)
    n_fp = len(false_positives)
    n_fn = len(false_negatives)
    n_our = len(our_set)
    n_key = len(key_set)

    sensitivity = n_tp / n_key       if n_key  else 0.0
    precision   = n_tp / n_our       if n_our  else 0.0
    f1          = (2 * sensitivity * precision
                   / (sensitivity + precision)
                   if (sensitivity + precision) else 0.0)
    actual_fdr  = n_fp / n_our       if n_our  else 0.0

    # ── score distribution of TP vs FP ───────────────────────────────────────
    score_stats = {}
    if score_col and score_col in df_our.columns:
        tp_mask = df_our["_pep_norm"].isin(true_positives)
        fp_mask = df_our["_pep_norm"].isin(false_positives)
        tp_scores = df_our.loc[tp_mask, score_col]
        fp_scores = df_our.loc[fp_mask, score_col]
        if len(tp_scores) and len(fp_scores):
            score_stats = dict(
                tp_median=tp_scores.median(),
                fp_median=fp_scores.median(),
                tp_mean=tp_scores.mean(),
                fp_mean=fp_scores.mean(),
            )

    # ── print report ─────────────────────────────────────────────────────────
    bar = "─" * 60
    print(f"\n{bar}")
    print(f"  {'Metric':<30} {'Value':>12}")
    print(bar)
    print(f"  {'Our unique peptides':<30} {n_our:>12,}")
    print(f"  {'Answer key peptides':<30} {n_key:>12,}")
    print(bar)
    print(f"  {'True Positives (TP)':<30} {n_tp:>12,}")
    print(f"  {'False Positives (FP)':<30} {n_fp:>12,}")
    print(f"  {'False Negatives (FN)':<30} {n_fn:>12,}")
    print(bar)
    print(f"  {'Sensitivity (Recall)':<30} {sensitivity:>11.1%}")
    print(f"  {'Precision':<30} {precision:>11.1%}")
    print(f"  {'F1 Score':<30} {f1:>11.3f}")
    print(bar)
    print(f"  {'Nominal FDR':<30} {fdr_nominal:>11.1%}")
    print(f"  {'Actual FDR (FP / total)':<30} {actual_fdr:>11.1%}")
    fdr_ok = "✔ within 2×" if actual_fdr <= fdr_nominal * 2 else "⚠ inflated"
    print(f"  {'FDR status':<30} {fdr_ok:>12}")

    if score_stats:
        print(bar)
        print(f"  {'TP median score':<30} {score_stats['tp_median']:>12.2f}")
        print(f"  {'FP median score':<30} {score_stats['fp_median']:>12.2f}")
        sep = score_stats['tp_median'] - score_stats['fp_median']
        print(f"  {'Score separation':<30} {sep:>12.2f}")

    print(bar)

    # ── save detailed results ─────────────────────────────────────────────────
    out_dir = Path("pipeline_output")
    out_dir.mkdir(exist_ok=True)

    # TP peptides
    tp_df = df_our[df_our["_pep_norm"].isin(true_positives)].copy()
    tp_df.drop(columns=["_pep_norm"], inplace=True)
    tp_df.to_csv(out_dir / "benchmark_true_positives.csv", index=False)

    # FP peptides
    fp_df = df_our[df_our["_pep_norm"].isin(false_positives)].copy()
    fp_df.drop(columns=["_pep_norm"], inplace=True)
    fp_df.to_csv(out_dir / "benchmark_false_positives.csv", index=False)

    # FN peptides (from answer key only)
    fn_df = df_key[df_key["_pep_norm"].isin(false_negatives)].copy()
    fn_df.drop(columns=["_pep_norm"], inplace=True)
    fn_df.to_csv(out_dir / "benchmark_false_negatives.csv", index=False)

    # summary CSV
    summary = pd.DataFrame([{
        "our_peptides":   n_our,
        "answer_peptides": n_key,
        "true_positives":  n_tp,
        "false_positives": n_fp,
        "false_negatives": n_fn,
        "sensitivity":     round(sensitivity, 4),
        "precision":       round(precision, 4),
        "f1_score":        round(f1, 4),
        "nominal_fdr":     fdr_nominal,
        "actual_fdr":      round(actual_fdr, 4),
    }])
    summary.to_csv(out_dir / "benchmark_summary.csv", index=False)

    print(f"\n  Results saved to pipeline_output/")
    print(f"    benchmark_summary.csv")
    print(f"    benchmark_true_positives.csv   ({n_tp:,} rows)")
    print(f"    benchmark_false_positives.csv  ({n_fp:,} rows)")
    print(f"    benchmark_false_negatives.csv  ({n_fn:,} rows)")
    print()

    return dict(sensitivity=sensitivity, precision=precision,
                f1=f1, actual_fdr=actual_fdr)


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Benchmark pipeline PSMs against a known answer key.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    parser.add_argument("--our-psms",   required=True,
                        help="Path to identified_psms.csv from the pipeline")
    parser.add_argument("--answer-key", required=True,
                        help="Path to answer key CSV with correct PSMs")
    parser.add_argument("--our-col",    default="Peptide",
                        help="Column name in our PSMs file (default: Peptide)")
    parser.add_argument("--key-col",    default="Peptide",
                        help="Column name in answer key file (default: Peptide)")
    parser.add_argument("--score-col",  default="Score",
                        help="Score column for TP/FP distribution (default: Score)")
    parser.add_argument("--fdr",        type=float, default=0.01,
                        help="Nominal FDR used in search (default: 0.01)")

    args = parser.parse_args()
    run_benchmark(
        our_psms_path=args.our_psms,
        answer_key_path=args.answer_key,
        our_col=args.our_col,
        key_col=args.key_col,
        score_col=args.score_col,
        fdr_nominal=args.fdr,
    )


if __name__ == "__main__":
    main()
