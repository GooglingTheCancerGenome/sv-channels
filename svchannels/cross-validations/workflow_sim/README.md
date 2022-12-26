workflow_sim: workflow to process simulated data

Run the snakemake workflow with:
```
snakemake -s Snakefile --use-conda -p --profile profile
```

Extract chromosomes chr10 and chr12 from file "exclude.N-and-encodebl.bed"
and rename them 10 and 12. Coordinates are the same because both hs37d5
```commandline
grep -e 'chr10' -e 'chr12' exclude.N-and-encodebl.bed | \
 sed 's/chr//' > exclude.N-and-encodebl.10_12_hs37d5.bed
```
