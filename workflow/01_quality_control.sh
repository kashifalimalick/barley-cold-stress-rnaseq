#!/bin/bash

# ==========================================
# Step 01: Raw RNA-seq Quality Control
# ==========================================

mkdir -p results/qc/raw

fastqc raw_data/*.fastq.gz \
    --outdir results/qc/raw
