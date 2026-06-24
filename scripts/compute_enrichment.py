#!/usr/bin/env python3
"""Quantify enrichment of LSR variants from Nanopore alignments, normalised to the BxB1 dead control.

Barcode-to-construct assignments are read from an external CSV (--barcodes) that is kept local and
not committed, so no proprietary sequences appear in the published code.

CSV format (header required, one row per construct):
    construct_id,barcode,role
    construct_28,<barcode_sequence>,standard       # count reads carrying this barcode
    construct_357,,gblock                           # control: count all aligned reads
    construct_359,<barcode_sequence>,dead_control   # BxB1 dead control, used for normalization

Usage:
    python compute_enrichment.py --input-bam input.sorted.bam --output-bam output.sorted.bam \
        --barcodes barcodes.csv --output-csv enrichment_results.csv
"""

import argparse

import numpy as np
import pandas as pd
import pysam


def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))


# Load the barcode map from an external CSV (keeps proprietary sequences out of the code)
def load_constructs(path):
    """Read construct_id, barcode and role for every construct from the barcode CSV"""
    df = pd.read_csv(path, dtype=str).fillna("")
    constructs = []
    for _, row in df.iterrows():
        constructs.append({
            "construct_id": row["construct_id"].strip(),
            "barcode": row["barcode"].strip().upper(),
            "role": row["role"].strip().lower(),
        })
    return constructs


def count_reads(samfile, construct, barcode, role):
    """Count reads for a construct: all reads for gBlocks, otherwise reads containing the given barcode"""
    count = 0
    # Construct may not be present in this BAM's reference list
    try:
        reads = samfile.fetch(construct)
    except ValueError:
        return 0
    for read in reads:
        seq = read.query_sequence
        if not seq:
            continue
        # gBlocks: count all reads
        if role == "gblock":
            count += 1
        # Everything else: only count reads that carry the (original or inverted) barcode
        elif barcode and barcode in seq:
            count += 1
    return count


def main():
    # Set up argument parser for command-line options
    parser = argparse.ArgumentParser(description="Compute log2 enrichment of LSR variants, normalised to the BxB1 dead control.")
    parser.add_argument("--input-bam", required=True, help="Sorted, indexed input (pre-selection) BAM file")
    parser.add_argument("--output-bam", required=True, help="Sorted, indexed output (post-selection) BAM file")
    parser.add_argument("--barcodes", required=True, help="CSV of construct_id,barcode,role (kept local, not committed)")
    parser.add_argument("--output-csv", default="enrichment_results.csv", help="Path for the results CSV")
    args = parser.parse_args()

    # Load construct definitions (all 367 constructs) from the barcode file
    constructs = load_constructs(args.barcodes)

    # Use the BxB1 dead control as normalization reference (exactly one row must have this role)
    dead_controls = [c["construct_id"] for c in constructs if c["role"] == "dead_control"]
    if len(dead_controls) != 1:
        raise ValueError(f"Expected exactly one row with role 'dead_control', found {len(dead_controls)}.")
    dead_control_id = dead_controls[0]

    # Load sam files
    print("Loading SAM files...")
    samfile_inp = pysam.AlignmentFile(args.input_bam, "rb")
    samfile_out = pysam.AlignmentFile(args.output_bam, "rb")

    # Count reads
    print(f"\nCounting reads for {len(constructs)} constructs...")
    input_dict = {}
    output_dict = {}
    for i, c in enumerate(constructs):
        if i % 50 == 0:
            print(f"  Progress: {i}/{len(constructs)}")
        cid, barcode, role = c["construct_id"], c["barcode"], c["role"]
        # INPUT reads use the original barcode
        input_dict[cid] = count_reads(samfile_inp, cid, barcode, role)
        # OUTPUT reads use the INVERTED barcode (only present if the recombinase worked)
        inverted_barcode = reverse_complement(barcode) if barcode else ""
        output_dict[cid] = count_reads(samfile_out, cid, inverted_barcode, role)

    samfile_inp.close()
    samfile_out.close()

    # BxB1 dead control counts used for normalization
    ref_input = input_dict.get(dead_control_id, 0)
    ref_output = output_dict.get(dead_control_id, 0)
    print(f"\nReference (BxB1 dead) - Input: {ref_input}, Output: {ref_output}")

    # Calculate enrichment
    results = []
    for c in constructs:
        cid = c["construct_id"]
        inp = input_dict[cid]
        out = output_dict[cid]
        # Normalize to BxB1 dead reference
        norm_inp = inp / ref_input if ref_input > 0 else 0
        norm_out = out / ref_output if ref_output > 0 else 0
        enrichment = norm_out / norm_inp if norm_inp > 0 else 0
        log2_enrichment = np.log2(enrichment) if enrichment > 0 else -np.inf
        results.append({
            'construct': cid,
            'role': c["role"],
            'input_counts': inp,
            'output_counts': out,
            'normalized_input': norm_inp,
            'normalized_output': norm_out,
            'enrichment': enrichment,
            'log2_enrichment': log2_enrichment
        })

    # Create DataFrame and sort by enrichment
    df = pd.DataFrame(results)
    df_sorted = df.sort_values('enrichment', ascending=False)

    # Save to CSV
    df_sorted.to_csv(args.output_csv, index=False)
    print(f"\n✓ Saved to {args.output_csv}")

    # Print summary statistics
    print("\n=== SUMMARY ===")
    active = df_sorted[df_sorted['enrichment'] > 2]
    print(f"Constructs with enrichment > 2: {len(active)}")


if __name__ == "__main__":
    main()
