import os
import csv
from itertools import chain

GTF                  = 'Mus_musculus.GRCm38.74.gtf'
STAR_INDEX_DIRECTORY = 'genomeDir'
REFERENCE            = 'Mus_musculus.GRCm38.74.dna.toplevel.fa'
READLENGTH           = 100
PICARD_BASE          = '/usr/local/apps/picard/1.94'
studyresults = []

with open('samplesheet.txt') as ss:
    reader = csv.DictReader(ss,delimiter="\t")
    for row in reader:
        studyresults.append(row)
    
def source2fastq(wc):
    # for rna-seq alignment with STAR, need all fastqfiles for the source.  
    # Turns out these are all paired-end, so form read1 and read2 lists for return
    sourceres = filter(lambda res: res['source']==wc.source,studyresults)
    zippedlist = list(chain(*list([[res['read1'],res['read2']] for res in sourceres])))
    # zippedlist is a list of size 2 (assuming all paired-end reads)
    # with the first element being the first read and the second
    # being the second read
    return(zippedlist)
















