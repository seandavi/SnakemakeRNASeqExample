To start the pipeline:

In your ~/.bashrc script:

```{sh}
export MODULEPATH=/data/ngs/modules:$MODULEPATH
```

Then, to submit the pipeline:

```{sh}
qsub -l nodes=1:c16 submit.sh
```
