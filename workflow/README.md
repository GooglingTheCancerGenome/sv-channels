run from the root repo directory (not this workflow directory)
```
snakemake --profile workflow/ -s workflow/Snakefile --configfile workflow/sm.config.yml --jobs 8 train
```
