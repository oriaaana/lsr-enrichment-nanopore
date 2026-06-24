#!/usr/bin/env python3
"""Filter Nanopore reads by length and average quality across a directory of BAM or FASTQ files."""

import os
import gzip
import argparse

import pysam
from Bio import SeqIO


# Helper function for calculating mean quality from a list of qualities
def mean_quality(qual):
    return sum(qual) / len(qual) if qual else 0


# Filter a BAM file by read length and average quality
def filter_bam_file(input_path, output_path, min_len, max_len, min_qual):
    print(f"Filtering BAM file: {input_path} -> {output_path}")
    bam_in = pysam.AlignmentFile(input_path, "rb")
    bam_out = pysam.AlignmentFile(output_path, "wb", template=bam_in)
    for read in bam_in:
        # Skip reads without sequence or qualities
        if read.query_sequence is None or read.query_qualities is None:
            continue
        seq_len = len(read.query_sequence)
        avg_qual = sum(read.query_qualities) / seq_len if seq_len > 0 else 0
        # Write read if it passes filters
        if (min_len <= seq_len <= max_len) and (avg_qual >= min_qual):
            bam_out.write(read)
    bam_in.close()
    bam_out.close()


# Filter a FASTQ or FASTQ.GZ file by read length and average quality
def filter_fastq_file(input_path, output_path, min_len, max_len, min_qual):
    print(f"Filtering FASTQ file: {input_path} -> {output_path}")
    # Open input file (gzipped or plain text)
    handle_in = gzip.open(input_path, "rt") if input_path.endswith('.gz') else open(input_path, "rt")
    # Open output file (gzipped or plain text)
    handle_out = gzip.open(output_path, "wt") if output_path.endswith('.gz') else open(output_path, "wt")
    for record in SeqIO.parse(handle_in, "fastq"):
        seq_len = len(record.seq)
        avg_qual = mean_quality(record.letter_annotations["phred_quality"])
        # Write record if it passes filters
        if (min_len <= seq_len <= max_len) and (avg_qual >= min_qual):
            SeqIO.write(record, handle_out, "fastq")
    handle_in.close()
    handle_out.close()


def main():
    # Set up argument parser for command-line options
    parser = argparse.ArgumentParser(description="Filter BAM/FASTQ reads by length and average quality in a directory.")
    parser.add_argument("-d", "--directory", required=True, help="Input directory containing BAM or FASTQ(.gz) files")
    parser.add_argument("-t", "--type", required=True, choices=["bam", "fastq"], help="Input file type: bam or fastq")
    parser.add_argument("-O", "--output_dir", required=True, help="Output directory for filtered files")
    parser.add_argument("--min_len", type=int, default=1350, help="Minimum read length (default: 1350)")
    parser.add_argument("--max_len", type=int, default=1700, help="Maximum read length (default: 1700)")
    parser.add_argument("--min_qual", type=int, default=20, help="Minimum average quality (default: 20)")
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Process all files in the input directory
    for fname in os.listdir(args.directory):
        fpath = os.path.join(args.directory, fname)
        out_path = os.path.join(args.output_dir, f"filtered_{fname}")
        # Filter BAM files
        if args.type == "bam" and fname.endswith(".bam"):
            filter_bam_file(fpath, out_path, args.min_len, args.max_len, args.min_qual)
            print(f"Filtered BAM: {fname} -> {out_path}")
        # Filter FASTQ or FASTQ.GZ files
        elif args.type == "fastq" and (fname.endswith(".fastq") or fname.endswith(".fastq.gz")):
            filter_fastq_file(fpath, out_path, args.min_len, args.max_len, args.min_qual)
            print(f"Filtered FASTQ: {fname} -> {out_path}")


if __name__ == "__main__":
    main()
