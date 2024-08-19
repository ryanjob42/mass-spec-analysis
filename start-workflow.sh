#!/bin/bash

conda activate mass-spec-analysis
module load singularity
snakemake
