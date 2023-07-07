# Cancer genomics course, EMBL-EBI, 2023

# Alignment workflow

## Set the environment for running this workflow

Let's define some variables before starting to run things.

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

## Initial Quality Control (QC) step

```
OUT_FASTQC=$OUT_FOLDER/fastqc_results
mkdir -p $OUT_FASTQC

$FASTQC --noextract --nogroup -o $OUT_FASTQC ${FOLDER}/fastq_files/*.fastq.gz
```
## Trimming adapters

```
OUT_TRIM=$OUT_FOLDER/trimmed_reads
mkdir -p $OUT_TRIM

## For tumour
trimmed1_T=$OUT_TRIM/COLO829T.R1.trimmed.fastq.gz
trimmed2_T=$OUT_TRIM/COLO829T.R2.trimmed.fastq.gz
${SEQPURGE} -in1 $read1_T -in2 $read2_T -out1 $trimmed1_T -out2 $trimmed2_T -qcut 0 -ncut 0 -threads 4

## For normal
trimmed1_N=$OUT_TRIM/COLO829BL.R1.trimmed.fastq.gz
trimmed2_N=$OUT_TRIM/COLO829BL.R2.trimmed.fastq.gz
${SEQPURGE} -in1 $read1_N -in2 $read2_N -out1 $trimmed1_N -out2 $trimmed2_N -qcut 0 -ncut 0 -threads 4
```

## Alignment
```
OUT_BAMS=$OUT_FOLDER/bams
mkdir -p $OUT_BAMS

## For tumour
$BWA mem -t 4 ${REFGENOME} $trimmed1_T $trimmed2_T | $SAMTOOLS view -Shb - > $OUT_BAMS/COLO829T.bam
$SAMTOOLS sort -o $OUT_BAMS/COLO829T.sorted.bam $OUT_BAMS/COLO829T.bam
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.bam

## For normal
$BWA mem -t 4 ${REFGENOME} $trimmed1_N $trimmed2_N | $SAMTOOLS view -Shb - > $OUT_BAMS/COLO829BL.bam
$SAMTOOLS sort -o $OUT_BAMS/COLO829BL.sorted.bam $OUT_BAMS/COLO829BL.bam
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.bam
```

## Preparing bam files for variant calling
```
## For tumour
# 3.1 Mark duplicates
$GATK MarkDuplicates -I $OUT_BAMS/COLO829T.sorted.bam -O $OUT_BAMS/COLO829T.sorted.MD.bam -M $OUT_BAMS/COLO829T.marked_dup_metrics.txt -REMOVE_DUPLICATES true --TMP_DIR $OUT_BAMS
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.MD.bam

# 3.2 Add read groups
$GATK AddOrReplaceReadGroups -I $OUT_BAMS/COLO829T.sorted.MD.bam -O $OUT_BAMS/COLO829T.sorted.MD.RG.bam -LB COLO829T -PL illumina -PU COLO829T -SM COLO829T
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.MD.RG.bam

# 3.3 Base quality score recalibration
$GATK BaseRecalibrator -I $OUT_BAMS/COLO829T.sorted.MD.RG.bam  -R $REFGENOME --known-sites $DBSNP -O $OUT_BAMS/COLO829T.recal_data.table
$GATK ApplyBQSR  -R $REFGENOME -I $OUT_BAMS/COLO829T.sorted.MD.RG.bam --bqsr-recal-file $OUT_BAMS/COLO829T.recal_data.table -O $OUT_BAMS/COLO829T.sorted.MD.RG.BQSR.bam
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.MD.RG.BQSR.bam

## For normal

# 3.1 Mark duplicates
$GATK MarkDuplicates -I $OUT_BAMS/COLO829BL.sorted.bam -O $OUT_BAMS/COLO829BL.sorted.MD.bam -M $OUT_BAMS/COLO829BL.marked_dup_metrics.txt -REMOVE_DUPLICATES true --TMP_DIR $OUT_BAMS
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.MD.bam

# 3.2 Add read groups
$GATK AddOrReplaceReadGroups -I $OUT_BAMS/COLO829BL.sorted.MD.bam -O $OUT_BAMS/COLO829BL.sorted.MD.RG.bam -LB COLO829BL -PL illumina -PU COLO829BL -SM COLO829BL
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.MD.RG.bam

# 3.3 Base quality score recalibration
$GATK BaseRecalibrator -I $OUT_BAMS/COLO829BL.sorted.MD.RG.bam  -R $REFGENOME --known-sites $DBSNP -O $OUT_BAMS/COLO829BL.recal_data.table
$GATK ApplyBQSR  -R $REFGENOME -I $OUT_BAMS/COLO829BL.sorted.MD.RG.bam --bqsr-recal-file $OUT_BAMS/COLO829BL.recal_data.table -O $OUT_BAMS/COLO829BL.sorted.MD.RG.BQSR.bam
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.MD.RG.BQSR.bam
```

