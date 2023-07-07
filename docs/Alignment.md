# Cancer genomics course, EMBL-EBI, 2023

## Alignment workflow

Let's define some variables before starting to run things

```
# Tools to be used
SEQPURGE=SeqPurge
BWA=bwa
SAMTOOLS=samtools
GATK=gatk
FASTQC=fastqc

# Main folder
FOLDER=/home/training/Documents/DNAseq_short_reads/

# Variant calling (practical)
REFGENOME=${FOLDER}/source_data/chr7.fa

# dbSNP vcf with known germline sites
DBSNP=${FOLDER}/source_data/Homo_sapiens_assembly38.dbsnp138.vcf.gz

# Define and create out directory
OUT_FOLDER=${FOLDER}/practical_results/alignment
mkdir -p $OUT_FOLDER

# Fastq reads tumour
read1_T=${FOLDER}/fastq_files/COLO829T.R1.fastq.gz
read2_T=${FOLDER}/fastq_files/COLO829T.R2.fastq.gz

# Fastq reads normal
read1_N=${FOLDER}/fastq_files/COLO829BL.R1.fastq.gz
read2_N=${FOLDER}/fastq_files/COLO829BL.R2.fastq.gz
```
