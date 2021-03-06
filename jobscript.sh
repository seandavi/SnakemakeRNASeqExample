#!/bin/sh
#rule: {job}
#input: {job.input}
#output: {job.output}
cd $PBS_O_WORKDIR
module load snakemake
snakemake --snakefile {self.workflow.snakefile} \
--force -j{self.cores} \
--directory {workdir} --nocolor --notemp --quiet --nolock {job.output} \
> /dev/null && touch "{jobfinished}" || touch "{jobfailed}"
exit 0
