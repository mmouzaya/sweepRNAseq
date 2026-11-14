# sweepRNAseq

sweepRNAseq is a Python pipeline for cleaning up bulk RNA-seq data. It removes rRNA reads with SortMeRNA, checks quality with FastQC before and after cleaning, and aligns the cleaned reads to a reference genome with STAR to check for leftover genomic contamination.

We built this for BME 160 (Research Programming in the Life Sciences) at UC Santa Cruz in Fall 2025.

## How it works

For each FASTQ file, the pipeline:

1. Runs FastQC on the raw reads
2. Removes rRNA reads with SortMeRNA
3. Runs FastQC on the cleaned reads
4. Aligns the cleaned reads to the reference genome with STAR
5. Runs FastQC on the aligned reads to confirm the cleanup worked

Each tool is wrapped in its own Python class (`FastQC`, `SortMeRNA`, `STAR`), and a `Pipeline` class runs them in order.

## Requirements

- Python 3.10+
- FastQC
- SortMeRNA (we used v4.3.7)
- STAR

We installed everything in a conda environment:

```bash
conda create -n sweepRNAseq-env -c bioconda -c conda-forge python=3.10 fastqc sortmerna=4.3.7 star
conda activate sweepRNAseq-env
```

## Usage

```bash
python sweepRNAseq.py sample1.fastq sample2.fastq annotation.gtf reference-genome.fa rrna-database.fa
```

The last two arguments always have to be the reference genome and the rRNA database. The GTF annotation file is optional, but if you include it, it goes right before those two. Everything before that is treated as a FASTQ file.

## Output

Everything is saved in an `output/` folder:

- `fastqc/` – FastQC reports for the raw, cleaned, and aligned reads
- `sortmerna/` – the rRNA reads that were removed
- `cleaned/` – the cleaned FASTQ files
- `star/` – the STAR genome index and sorted BAM files

## Testing

We tested the pipeline on bulk RNA-seq data from six brain organoid samples from our lab. The data isn't included here.

## Authors

Matilda Mouzaya, Sofia Naderi, and Sophia Oriaku
