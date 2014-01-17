import os
import config

from snakemake.utils import report
from snakemake.utils import R

COUNTFILES = expand('bam/{source}/{source}.md.counts',source=set([res['source'] for res in config.studyresults]))

include: 'bampost.snk'

localrules: final, bamindex

rule final:
    input: COUNTFILES, "reports/index.html"
    params: batch='-l nodes=1:gpfs'

rule fai:
    input: '{base}.fa'
    output: '{base}.fa.fai'
    params: batch='-l nodes=1:gpfs'
    threads: 2
    shell:  "samtools faidx {input}"
    
rule sampleslist:
    input: "samplesheet.txt"
    output: "samples.list"
    params: batch='-l nodes=1:gpfs'
    threads: 2
    run:
        with open(output,'w') as f:
            f.write('Sample ID\tBam File\tNotes\n')
            for res in config.studyresults:
                f.write('%s\tbam/%s/%s.md.bam\tNone\n'.format(res['source'],res['source'],res['source']))

rule rnaseqqc:
    input: samplelist="samples.list",
           gtf = config.GTF, 
           fasta = config.REFERENCE,
           fai = config.REFERENCE + '.fai'
    output: "reports/index.html"
    params: batch='-l nodes=1:gpfs',gtf = config.GTF
    threads: 2
    shell:  "java -jar RNA-SeQC_v1.1.7.jar -o reports -r {input.fasta} -s {input.samplelist} -t {params.gtf} -ttype 2 "
        

rule rnapostalignment:
    input: bam="bam/{sample,\w+}/aligned.bam",chimeric="bam/{sample}/Chimeric.out.sam"
    output: bam=temp("bam/{sample}/{sample}.bam"),chimeric=temp("bam/{sample}/{sample}.chimeric.bam")
    params: batch='-l nodes=1:gpfs'
    threads: 2
    shell: """
clearscratch
samtools view -bS {input.chimeric} \
    | samtools sort - bam/{wildcards.sample}/{wildcards.sample}.chimeric
java -jar /usr/local/picard/MergeSamFiles.jar \
    I={output.chimeric} \
    I={input.bam} \
    O=/scratch/tmp.bam \
    AS=true
java -jar /usr/local/picard/AddOrReplaceReadGroups.jar \
    RGID={wildcards.sample} \
    RGLB={wildcards.sample} \
    RGPL=illumina \
    RGPU={wildcards.sample} \
    RGSM={wildcards.sample} \
    INPUT=/scratch/tmp.bam \
    OUTPUT={output.bam}
clearscratch
"""


rule rnacount:
    input: bam='{base}.bam',bai='{base}.bam.bai',gtf=config.GTF
    output: "{base}.counts"
    params: batch="-l nodes=1:gpfs"
    threads: 16
    shell: """
module load subread
featureCounts -a {input.gtf} -i {input.bam} -o {output} -b -T {threads}
"""

rule starindex:
    input: fasta = config.REFERENCE,gtf = config.GTF
    output: [os.path.join(config.STAR_INDEX_DIRECTORY,fname) for fname in ["chrLength.txt","chrNameLength.txt", "chrName.txt", "chrStart.txt", "Genome", "genomeParameters.txt", "SA", "SAindex", "sjdbInfo.txt"]]
    params: batch     = "-l nodes=1 -q ccr",
            overhang  = str(config.READLENGTH-1),
            genomeDir = config.STAR_INDEX_DIRECTORY
    threads: 32
    shell: """
mkdir -p {params.genomeDir}
/data/CCRBioinfo/biowulf/local/STAR_2.3.1n/STAR --runMode genomeGenerate \
    --genomeDir {params.genomeDir} \
    --sjdbGTFfile {input.gtf} \
    --genomeFastaFiles {input.fasta} \
    --sjdbOverhang {params.overhang} \
    --runThreadN {threads}
"""

rule rnaalignment:
    input:  config.source2fastq,indexfiles = [os.path.join(config.STAR_INDEX_DIRECTORY,fname) for fname in ["chrLength.txt","chrNameLength.txt", "chrName.txt", "chrStart.txt", "Genome", "genomeParameters.txt", "SA", "SAindex", "sjdbInfo.txt"]]
    output: temp("bam/{source}/aligned.bam"), temp("bam/{source}/Chimeric.out.sam")
    params:  batch="-q ccr -l nodes=1:gpfs",
             outdir="bam/{source}",
             genomeDir = config.STAR_INDEX_DIRECTORY
    priority: 10
    threads: 32
    version: "2.3.1n"
    log: "bam/{source}/star.log}"
    run:
        readfiles = input[0:(len(input)-9)]
        read1 = readfiles[0::2]
        read2 = readfiles[1::2]
        shell("""mkdir -p {params.outdir}
/data/CCRBioinfo/biowulf/local/STAR_2.3.1n/STAR \
   --genomeDir {params.genomeDir} \
   --outFilterType BySJout \
   --readFilesIn <( zcat {read1} ) <( zcat {read2} ) \
   --runThreadN {threads} \
   --outSAMattributes Standard \
   --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
   --alignIntronMax 100000 \
   --outSAMstrandField intronMotif \
   --outFileNamePrefix {params.outdir}/ \
   --outSAMunmapped Within \
   --chimSegmentMin 25 \
   --chimJunctionOverhangMin 25 \
   --outStd SAM \
      | samtools view -bS - \
      | samtools sort -m 15000000000 - {params.outdir}/aligned
""")
