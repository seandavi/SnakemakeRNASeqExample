To start the pipeline:

In your ~/.bashrc script:

```{sh}
export MODULEPATH=/data/ngs/modules:$MODULEPATH
```

To be able to use snakemake:

```{sh}
module load snakemake
```

To submit a snakemake job, change into the directory with the Snakefile and execute:

```{sh}
qsub -l nodes=1:c16 submit.sh
```
