#!/bin/bash

source ~/.bashrc
conda activate mass-spec-analysis
module load singularity
snakemake
