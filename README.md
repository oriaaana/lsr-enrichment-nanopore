# LSR Enrichment — Nanopore Analysis

Code and pipeline used to quantify the enrichment of computationally mined
large serine recombinase (LSR) variants following selection. Pre-selection
("input") and post-selection ("output") libraries were sequenced on a MinION
Mk1D, and enrichment was calculated per candidate variant from the resulting
read counts.

## Enrichment metric

For each screened LSR *j*, enrichment is calculated relative to the BxB1 dead
control:

```
enrichment_j = log2( (count_j,output / count_dead,output)
                     / (count_j,input  / count_dead,input) )
```

## Pipeline overview

| Step | Tool | Custom code? |
|------|------|--------------|
| 1. Basecalling | Dorado v1.4.0 (r10.4.1 model) | external |
| 2. Demultiplexing | Dorado demux (SQK-NBD114-24) | external |
| 3. Read QC | NanoPlot / FastQC | external |
| 4. Length & quality filtering | `scripts/filter_reads.py` | **yes** |
| 5. Alignment | minimap2 v2.26 (`-ax map-ont`) | external |
| 6. SAM → sorted/indexed BAM | samtools v1.17 | external |
| 7. Read counting & enrichment | `scripts/compute_enrichment.py` | **yes** |

## Repository structure

```
lsr-enrichment-nanopore/
├── README.md
├── environment.yml
├── data/
│   ├── reference/
│   │   └── candidate_lsr_reference.fasta   # 367 candidate LSRs + BxB1 dead control
│   ├── barcodes_template.csv               # format only; no real sequences
│   └── README.md
├── scripts/
│   ├── filter_reads.py
│   └── compute_enrichment.py
└── results/
    └── enrichment_results.csv              # example output
```

## Setup

Create the analysis environment:

```bash
mamba env create -f environment.yml
mamba activate nanopore
```

Dorado is installed separately from Oxford Nanopore. Use the build that
matches your operating system (the example below is for Apple Silicon macOS;
see https://github.com/nanoporetech/dorado for other builds):

```bash
curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-osx-arm64.zip" -o dorado.zip
unzip dorado.zip
export PATH="$PWD/dorado-1.4.0-osx-arm64/bin:$PATH"
dorado --version

# One-time basecalling model download
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0
```

## Running the pipeline

Replace `run1` with your own run identifier throughout.

### 1. Basecalling

```bash
mkdir -p ~/nanopore_data/run1
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.2.0 path/to/pod5/ \
    --kit-name SQK-NBD114-24 \
    > ~/nanopore_data/run1/calls.bam
```

### 2. Demultiplexing

Split the basecalled reads into one FASTQ per barcode (one per sample):

```bash
dorado demux --kit-name SQK-NBD114-24 --emit-fastq \
    --output-dir ~/nanopore_data/run1/demux/ \
    ~/nanopore_data/run1/calls.bam
```

### 3. Read QC

```bash
NanoPlot --fastq ~/nanopore_data/run1/demux/SAMPLE.fastq \
    --outdir ~/nanopore_data/run1/nanoplot
```

### 4. Length & quality filtering

```bash
python scripts/filter_reads.py \
    -d ~/nanopore_data/run1/demux \
    -t fastq \
    -O ~/nanopore_data/run1/filtered \
    --min_len 1600 --max_len 2200 --min_qual 20
```

Optionally re-run FastQC on the filtered output:

```bash
fastqc ~/nanopore_data/run1/filtered/filtered_SAMPLE.fastq.gz
```

### 5. Alignment + sorted/indexed BAM

Align each library (input and output) to the candidate reference, then sort
and index. Shown here for the input library; repeat for the output library.

```bash
minimap2 -ax map-ont data/reference/candidate_lsr_reference.fasta \
    ~/nanopore_data/run1/filtered/filtered_input.fastq.gz > input.sam
samtools sort input.sam -o input.sorted.bam
samtools index input.sorted.bam
```

### 6. Read counting & enrichment

```bash
python scripts/compute_enrichment.py \
    --input-bam input.sorted.bam \
    --output-bam output.sorted.bam \
    --barcodes barcodes.csv \
    --output-csv results/enrichment_results.csv
```

## Barcode file

`compute_enrichment.py` reads barcode-to-construct assignments from a CSV.
`data/barcodes_template.csv` shows the required format (with placeholder IDs
and sequences). To run the analysis, copy it to `barcodes.csv`, fill in the
real construct IDs, barcode sequences, and roles (`standard`, `gblock`, or
`dead_control`), and keep that file local — it is git-ignored so that no
proprietary sequences are published.

## Data availability

Raw POD5 files and basecalled/aligned reads are not included in this
repository due to IP reasons
