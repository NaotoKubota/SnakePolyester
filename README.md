# SnakePolyester (v0.1.0)

[![GitHub Release](https://img.shields.io/github/v/release/NaotoKubota/SnakePolyester?style=flat)](https://github.com/NaotoKubota/SnakePolyester/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/NaotoKubota/SnakePolyester)](https://github.com/NaotoKubota/SnakePolyester/releases)
[![Create Release](https://github.com/NaotoKubota/SnakePolyester/actions/workflows/release.yaml/badge.svg)](https://github.com/NaotoKubota/SnakePolyester/actions/workflows/release.yaml)

A snakemake workflow for simulating RNA-seq reads by [Polyester](https://github.com/alyssafrazee/polyester). You can generate RNA-seq reads for multiple scenarios with different coverage, replicates, and size factors.

This workflow generates RNA-seq read files in fastq format and maps them to the reference genome using [STAR](https://github.com/alexdobin/STAR).

## Install

- Make sure you have [Singularity](https://docs.sylabs.io/guides/3.7/user-guide/index.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed.
- Clone this repository.

```bash
git clone https://github.com/NaotoKubota/SnakePolyester.git
```

## Usage

```bash
snakemake -s /path/to/SnakePolyester/SnakePolyester.smk --configfile config.yaml --cores 32 --use-singularity --singularity-args "--bind $HOME:$HOME" --rerun-incomplete
```

## Input files

- `transcript_table.tsv`: A table containing transcript IDs and their fold-changes. Columns are `transcript_id`, `fold-change in condition1`, `fold-change in condition2`.

```text
ENSMUST00000000003	1.0	1.4
ENSMUST00000000028	1.0	0.6
ENSMUST00000000080	1.0	1.0
ENSMUST00000000090	1.0	1.4
ENSMUST00000000153	1.0	1.0
ENSMUST00000000161	1.0	0.6
ENSMUST00000000284	1.0	1.4
ENSMUST00000000299	1.0	0.6
ENSMUST00000000304	1.0	1.0
ENSMUST00000000335	1.0	0.6
```

- `config.yaml`: Configuration file.

```yaml
workdir:
    /path/to/outdir
transcript_table:
    /path/to/transcript_table.tsv
cdna_fasta:
    /path/to/Mus_musculus.GRCm38.cdna.all.fa # cDNA fasta file
star_index:
    /path/to/star_index # STAR index directory

scenarios:
    scenario1: # Scenario name
        coverage: # Coverage of the simulated reads for each transcript
            30
        replicates: # Number of replicates for each condition
            1
        size_factor: # controls the per-transcript mean/variance relationship. Increase the value of size factor to introduce more variance into your simulations.
            3 # low-variance situation
    scenario2:
        coverage:
            30
        replicates:
            2
        size_factor:
            3
    scenario3:
        coverage:
            30
        replicates:
            3
        size_factor:
            3
```

## Output files

- `benchmark`: [Benchmarking](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html) results.
- `data`: Simulated RNA-seq read files.
- `logs`: Log files.

```bash
outdir/
├── benchmark
├── data
└── logs
```

### An example of `data` directory structure

```bash
outdir/data/
├── scenario1
│   ├── bam_table.tsv
│   ├── fastq_table.tsv
│   ├── flag.txt
│   ├── polyester
│   │   ├── sample_01_1.fastq.gz
│   │   ├── sample_01_2.fastq.gz
│   │   ├── sample_02_1.fastq.gz
│   │   ├── sample_02_2.fastq.gz
│   │   ├── sim_counts_matrix.rda
│   │   ├── sim_rep_info.txt
│   │   └── sim_tx_info.txt
│   └── star
│       ├── sample_01
│       │   ├── sample_01_Aligned.out.bam
│       │   ├── sample_01_Aligned.out.bam.bai
│       │   ├── sample_01_Log.final.out
│       │   ├── sample_01_Log.out
│       │   ├── sample_01_Log.progress.out
│       │   └── sample_01_SJ.out.tab
│       └── sample_02
│           ├── sample_02_Aligned.out.bam
│           ├── sample_02_Aligned.out.bam.bai
│           ├── sample_02_Log.final.out
│           ├── sample_02_Log.out
│           ├── sample_02_Log.progress.out
│           └── sample_02_SJ.out.tab
├── scenario2
│   ├── bam_table.tsv
│   ├── fastq_table.tsv
│   ├── flag.txt
│   ├── polyester
│   │   ├── sample_01_1.fastq.gz
│   │   ├── sample_01_2.fastq.gz
│   │   ├── sample_02_1.fastq.gz
│   │   ├── sample_02_2.fastq.gz
│   │   ├── sample_03_1.fastq.gz
│   │   ├── sample_03_2.fastq.gz
│   │   ├── sample_04_1.fasta
│   │   ├── sample_04_2.fasta
│   │   ├── sim_counts_matrix.rda
│   │   ├── sim_rep_info.txt
│   │   └── sim_tx_info.txt
│   └── star
│       ├── sample_01
│       │   ├── sample_01_Aligned.out.bam
│       │   ├── sample_01_Aligned.out.bam.bai
│       │   ├── sample_01_Log.final.out
│       │   ├── sample_01_Log.out
│       │   ├── sample_01_Log.progress.out
│       │   └── sample_01_SJ.out.tab
│       ├── sample_02
│       │   ├── sample_02_Aligned.out.bam
│       │   ├── sample_02_Aligned.out.bam.bai
│       │   ├── sample_02_Log.final.out
│       │   ├── sample_02_Log.out
│       │   ├── sample_02_Log.progress.out
│       │   └── sample_02_SJ.out.tab
│       └── sample_03
│           ├── sample_03_Aligned.out.bam
│           ├── sample_03_Aligned.out.bam.bai
│           ├── sample_03_Log.final.out
│           ├── sample_03_Log.out
│           ├── sample_03_Log.progress.out
│           └── sample_03_SJ.out.tab
├── scenario3
│   ├── bam_table.tsv
│   ├── fastq_table.tsv
│   ├── flag.txt
│   ├── polyester
│   │   ├── sample_01_1.fastq.gz
│   │   ├── sample_01_2.fastq.gz
│   │   ├── sample_02_1.fastq.gz
│   │   ├── sample_02_2.fastq.gz
│   │   ├── sample_03_1.fastq.gz
│   │   ├── sample_03_2.fastq.gz
│   │   ├── sample_04_1.fastq.gz
│   │   ├── sample_04_2.fastq.gz
│   │   ├── sample_05_1.fasta
│   │   ├── sample_05_2.fasta
│   │   ├── sample_06_1.fasta
│   │   ├── sample_06_2.fasta
│   │   ├── sim_counts_matrix.rda
│   │   ├── sim_rep_info.txt
│   │   └── sim_tx_info.txt
│   └── star
│       ├── sample_01
│       │   ├── sample_01_Aligned.out.bam
│       │   ├── sample_01_Aligned.out.bam.bai
│       │   ├── sample_01_Log.final.out
│       │   ├── sample_01_Log.out
│       │   ├── sample_01_Log.progress.out
│       │   └── sample_01_SJ.out.tab
│       ├── sample_02
│       │   ├── sample_02_Aligned.out.bam
│       │   ├── sample_02_Aligned.out.bam.bai
│       │   ├── sample_02_Log.final.out
│       │   ├── sample_02_Log.out
│       │   ├── sample_02_Log.progress.out
│       │   └── sample_02_SJ.out.tab
│       ├── sample_03
│       │   ├── sample_03_Aligned.out.bam
│       │   ├── sample_03_Aligned.out.bam.bai
│       │   ├── sample_03_Log.final.out
│       │   ├── sample_03_Log.out
│       │   ├── sample_03_Log.progress.out
│       │   └── sample_03_SJ.out.tab
│       └── sample_04
│           ├── sample_04_Aligned.out.bam
│           ├── sample_04_Aligned.out.bam.bai
│           ├── sample_04_Log.final.out
│           ├── sample_04_Log.out
│           ├── sample_04_Log.progress.out
│           └── sample_04_SJ.out.tab
└── transcript_cdna.fa
```
