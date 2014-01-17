#!/bin/bash
#
cd ${PBS_O_WORKDIR}
module load snakemake
snakemake --js /data/CCRBioinfo/jobscript.sh -k --stats snakemake.stats --debug -T --rerun-incomplete -j 400 --cluster 'qsub {params.batch}' >&  snakemake.log
