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

    # Real Juicer 2.0 / juicer_tools 2.17.00 lines look like:
    #   Sequenced Read Pairs:  193035634
    #   Total Duplicates: 18692703 (10.30%, 9.68%)
    #   Hi-C Contacts: 129,150,544 (66.91% / 79.36%)
    #   Inter-chromosomal: 20,646,363 (10.70% / 12.69%)
    # The two percentages are of all sequenced pairs and of unique pairs,
    # NOT of Hi-C contacts, so the inter-chromosomal share is computed from
    # the counts below instead of read off a percentage.
    def find_int(label):
        m = re.search(rf"^\s*{label}\s*:\s*([\d,]+)", text, re.IGNORECASE | re.MULTILINE)
        return int(m.group(1).replace(",", "")) if m else None

    def find_first_pct(label):
        m = re.search(rf"^\s*{label}\s*:\s*[\d,]+\s*\(\s*([\d.]+)\s*%", text, re.IGNORECASE | re.MULTILINE)
        return float(m.group(1)) if m else None

    fields["sequenced_pairs"] = find_int("Sequenced Read Pairs")
    fields["valid_pairs"]     = find_int("Hi-C Contacts")
    fields["inter_count"]     = find_int("Inter-chromosomal")
    fields["duplicates_pct"]  = find_first_pct("Total Duplicates")
    if fields["valid_pairs"] and fields["inter_count"] is not None:
        fields["inter_pct"] = 100.0 * fields["inter_count"] / fields["valid_pairs"]
    else:
        fields["inter_pct"] = None
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

    if f["sequenced_pairs"] is not None:
        lines.append(f"Sequenced read pairs: {f['sequenced_pairs']:,}")

    valid_pairs = f["valid_pairs"]
    if valid_pairs is None:
        lines.append("Valid read pairs (Hi-C contacts, MAPQ>=30, de-duplicated): UNKNOWN — could not parse inter_30.txt")
        warnings.append("Could not parse the Hi-C Contacts count from inter_30.txt")
    else:
        lines.append(f"Valid read pairs (Hi-C contacts, MAPQ>=30, de-duplicated): {valid_pairs:,}")
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
        warnings.append("Could not parse the inter-chromosomal count from inter_30.txt")
    else:
        lines.append(f"Inter-chromosomal contact % (of Hi-C contacts): {inter:.2f}%")
        if not (args.inter_min <= inter <= args.inter_max):
            warnings.append(
                f"Inter-chromosomal % ({inter:.2f}%) outside the expected "
                f"{args.inter_min}-{args.inter_max}% range — check for a failed "
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
