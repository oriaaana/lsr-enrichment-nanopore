# Nanopore Analysis of Mined LSR 

code used to quantify the enrichment of computationally mined
large serine recombinase (LSR) variants following selection. Pre-selection
("input") and post-selection ("output") libraries were sequenced on a MinION
Mk1D. Enrichment was calculated per variant.

## Enrichment 

for each screened LSR *j*, enrichment is calculated relative to the BxB1 dead
control:

```
enrichment_j = log2( (count_j,output / count_dead,output)
                     / (count_j,input  / count_dead,input) )
```

## Setup

Create the analysis environment:

```bash
mamba env create -f environment.yml
mamba activate nanopore
```

Install Dorado:

```bash
curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-osx-arm64.zip" -o dorado.zip
unzip dorado.zip
export PATH="$PWD/dorado-1.4.0-osx-arm64/bin:$PATH"
dorado --version

# download basecalling model
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0
```

## Pipeline

replace `run1` with run identifier.

### 1. basecall

```bash
mkdir -p ~/nanopore_data/run1
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.2.0 path/to/pod5/ \
    --kit-name SQK-NBD114-24 \
    > ~/nanopore_data/run1/calls.bam
```

### 2. demultiplex

split the basecalled reads into one FASTQ per barcode:

```bash
dorado demux --kit-name SQK-NBD114-24 --emit-fastq \
    --output-dir ~/nanopore_data/run1/demux/ \
    ~/nanopore_data/run1/calls.bam
```

### 3. read QC

```bash
NanoPlot --fastq ~/nanopore_data/run1/demux/SAMPLE.fastq \
    --outdir ~/nanopore_data/run1/nanoplot
```

### 4. filter Reads 

```bash
python scripts/filter_reads.py \
    -d ~/nanopore_data/run1/demux \
    -t fastq \
    -O ~/nanopore_data/run1/filtered \
    --min_len 1600 --max_len 2200 --min_qual 20
```

### 5. align

align each library (input and output) to the candidate reference, then sort
and index. Here shown for input ->> repeat for output.

```bash
minimap2 -ax map-ont data/reference/candidate_lsr_reference.fasta \
    ~/nanopore_data/run1/filtered/filtered_input.fastq.gz > input.sam
samtools sort input.sam -o input.sorted.bam
samtools index input.sorted.bam
```

### 6. count reads and compute enrichment

```bash
python scripts/compute_enrichment.py \
    --input-bam input.sorted.bam \
    --output-bam output.sorted.bam \
    --barcodes barcodes.csv \
    --output-csv enrichment_results.csv
```

## barcode file

`compute_enrichment.py` reads barcode-to-construct assignments from a CSV.
`data/barcodes_template.csv` shows the required format (with placeholder IDs
and sequences).

## Data availability

raw POD5 files, basecalled/aligned reads, candidate lsr reference amino acid fastas are not included in this
repository due to IP reasons.
