# 🌾 Barley Cold-Stress RNA-seq Analysis

### Transcriptomic analysis of *Hordeum vulgare* under cold stress using a Linux-based RNA-seq workflow.

---

## 🔬 Project Overview

This project investigates the transcriptional response of barley
(*Hordeum vulgare*) to cold stress using publicly available paired-end
RNA-seq data.

The analysis is being performed in a Linux/WSL environment using
command-line bioinformatics tools.

### Objectives

- Assess raw sequencing quality
- Perform adapter and quality trimming
- Align reads to the barley reference genome
- Process alignment files
- Quantify gene expression
- Identify differentially expressed genes
- Perform functional/pathway analysis

---

## 🧬 Dataset

| Feature | Information |
|---|---|
| Organism | *Hordeum vulgare* |
| Condition | Cold stress |
| Data type | Paired-end RNA-seq |
| Reference genome | Barley MorexV3 |
| Data source | ENA / SRA |

Raw FASTQ files are not included in this repository because of their
large size.

---

## 🧪 RNA-seq Workflow

```text
Raw FASTQ
    │
    ▼
FastQC
    │
    ▼
MultiQC
    │
    ▼
fastp
    │
    ▼
Trimmed FASTQ
    │
    ▼
FastQC + MultiQC
    │
    ▼
HISAT2
    │
    ▼
SAM/BAM
    │
    ▼
SAMtools
    │
    ▼
featureCounts
    │
    ▼
Gene Count Matrix
    │
    ▼
Differential Expression
    │
    ▼
GO / Pathway Analysis
