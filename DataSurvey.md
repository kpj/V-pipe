# Data Survey

## Usage

```bash
$ snakemake -s data_survey.smk -pr --cores 2 --use-conda
$ python prepare_vpipe_input.py
$ snakemake -s vpipe.snake -pr --cores 2 --use-conda
$ python gather_vcf_files.py
```
