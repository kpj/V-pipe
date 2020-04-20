#!/usr/bin/env bash

bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -q normal.24h \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  -s data_survey.smk \
  --use-conda \
  "$@"
