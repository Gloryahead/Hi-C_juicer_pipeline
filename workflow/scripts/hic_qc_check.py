#!/usr/bin/env python3
"""
hic_qc_check.py — parse Juicer's own inter_30.txt alignment/contact
statistics and flag values outside the tutorial's stated expected ranges
(inter-chromosomal %, valid read pairs). These ranges are the tutorial
author's own rule-of-thumb, not a peer-reviewed statistical benchmark — see
docs/KEY_CONCEPTS.md before treating a WARN as a hard failure.

This script never exits non-zero: it is a checkable QC report wired into the
DAG (workflow/rules/02_juicer_core.smk's hic_qc_summary rule), not a gate
that blocks the pipeline — matching how every other QC checkpoint in this
pipeline family (FastQC, MultiQC, Juicer's own validate step) works.

Usage:
    python3 hic_qc_check.py --stats inter_30.txt --sample hic1 \
        --inter-min 5 --inter-max 40 --pairs-min 10000000 --pairs-high-res 100000000 \
        --out qc/hic1_hic_qc_summary.txt
"""
import argparse
import re
import sys


def parse_stats(path):
    """Extract the handful of fields Juicer's inter_30.txt always reports.
    Field names are matched loosely (case-insensitive, tolerant of the
    ':'/'%' formatting Juicer uses) since the exact wording has drifted
    slightly across Juicer versions.
    """
    text = open(path).read()
    fields = {}

    def find_int(label):
        m = re.search(rf"{label}\s*:\s*([\d,]+)", text, re.IGNORECASE)
        return int(m.group(1).replace(",", "")) if m else None

    def find_pct(label):
        m = re.search(rf"{label}\s*:[^%]*\(\s*([\d.]+)\s*%\s*\)", text, re.IGNORECASE)
        return float(m.group(1)) if m else None

    fields["sequenced_pairs"] = find_int("Sequenced Read Pairs") or find_int("Total Read Pairs")
    fields["alignable_pairs"] = find_int(r"Alignable \(Normal\+Chimeric Paired\)") or find_int("Alignable")
    fields["duplicates_pct"]  = find_pct("Duplicates") or find_pct("PCR Duplicates")
    fields["inter_pct"]       = find_pct("Inter-chromosomal")
    fields["intra_pct"]       = find_pct("Intra-chromosomal")
    return fields


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--inter-min", type=float, required=True)
    ap.add_argument("--inter-max", type=float, required=True)
    ap.add_argument("--pairs-min", type=int, required=True)
    ap.add_argument("--pairs-high-res", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    f = parse_stats(args.stats)
    lines = [f"Hi-C QC summary — sample: {args.sample}", "=" * 50]
    warnings = []

    valid_pairs = f["alignable_pairs"]
    if valid_pairs is None:
        lines.append("Valid (alignable) read pairs: UNKNOWN — could not parse inter_30.txt")
        warnings.append("Could not parse alignable-pairs count from inter_30.txt")
    else:
        lines.append(f"Valid (alignable) read pairs: {valid_pairs:,}")
        if valid_pairs < args.pairs_min:
            warnings.append(
                f"Valid pairs ({valid_pairs:,}) below the basic-resolution floor "
                f"({args.pairs_min:,}) — TADs/compartments may be unreliable."
            )
        elif valid_pairs < args.pairs_high_res:
            lines.append(
                f"  (basic-resolution range; below the {args.pairs_high_res:,} floor "
                f"typically needed for confident loop calling — see docs/KEY_CONCEPTS.md)"
            )

    inter = f["inter_pct"]
    if inter is None:
        lines.append("Inter-chromosomal contact %: UNKNOWN — could not parse inter_30.txt")
        warnings.append("Could not parse inter-chromosomal %% from inter_30.txt")
    else:
        lines.append(f"Inter-chromosomal contact %: {inter:.2f}%")
        if not (args.inter_min <= inter <= args.inter_max):
            warnings.append(
                f"Inter-chromosomal %% ({inter:.2f}%%) outside the expected "
                f"{args.inter_min}-{args.inter_max}%% range — check for a failed "
                f"in-situ ligation step or excessive noise."
            )

    if f["duplicates_pct"] is not None:
        lines.append(f"PCR duplicate %: {f['duplicates_pct']:.2f}%")

    lines.append("")
    if warnings:
        lines.append(f"STATUS: WARN ({len(warnings)} item(s) outside expected range)")
        for w in warnings:
            lines.append(f"  - {w}")
    else:
        lines.append("STATUS: PASS")

    with open(args.out, "w") as out:
        out.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    sys.exit(0)  # advisory QC, never fails the DAG — see module docstring


if __name__ == "__main__":
    main()
