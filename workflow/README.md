
run from the root repo directory (not this workflow directory)

## Train

```
snakemake --profile workflow/ -s workflow/Snakefile --configfile workflow/sm.config.yml --jobs 8 --use-conda train
```


## Predict

You will need to modify `workflow/sm.config.yml`
```
snakemake --profile workflow/ -s workflow/Snakefile --configfile workflow/sm.config.yml --jobs 8 predict
```
