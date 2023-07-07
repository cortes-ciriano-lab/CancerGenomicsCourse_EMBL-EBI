# Variant calling
Somatic (short) variant calling is the process of identifying genetic variations specific to tumour cells by comparing the genomic data from the tumour and matched normal samples. It involves detecting somatic mutations, such as single-nucleotide variants (SNVs) and/or small insertions/deletions (indels), that are acquired during tumorigenesis and distinguishing them from germline variants present in both tumour and normal cells.

In this section, we will explain how to run different variant calling algorithms ([Strelka2](https://www.nature.com/articles/s41592-018-0051-x) and [MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-)), as well as how to filter the raw variants to get a high-quality set of somatic mutations. 

This workflow continues with the results obtained from the [alignment exercise](https://github.com/cortes-ciriano-lab/CancerGenomicsCourse_EMBL-EBI/blob/main/docs/Alignment.md#alignment-workflow---cancer-genomic-course-2023---embl-ebi), where we generated the COLO829T (tumour) and COLO829BL (matched-normal) bam files. 

## Strelka2 - Somatic workflow

```
STRELKA2=/home/training/strelka-2.9.10.centos6_x86_64/
REFGENOME=${FOLDER}/source_data/chr7.fa

# Main folder
FOLDER=/home/training/Documents/DNAseq_short_reads/

# COLO289 bam files
BAM_N=${FOLDER}/practical_results/alignment/bams/COLO829BL.sorted.MD.RG.BQSR.bam
BAM_T=${FOLDER}/practical_results/alignment/bams/COLO829T.sorted.MD.RG.bam

# Out folder
OUT_FOLDER1=${FOLDER}/practical_results/variant_calling/strelka2_results
mkdir -p $OUT_FOLDER1
rm -r $OUT_FOLDER1/*
```

```
${STRELKA2}/bin/configureStrelkaSomaticWorkflow.py \
--normalBam $BAM_N \
--tumorBam $BAM_T \
--referenceFasta $REFGENOME \
--runDir $OUT_FOLDER1
```
```
${OUT_FOLDER1}/runWorkflow.py  -m local -j 4 -g 15
```
```
# Keep only pass mutations
# For SNVs
VCF_SNV_STRELKA2=$OUT_FOLDER1/results/variants/somatic.snvs.vcf.gz
VCF_SNV_STRELKA2_PASS=$(echo $VCF_SNV_STRELKA2 | sed 's/.vcf.gz/.pass.vcf/g')
zgrep '#\|PASS' ${VCF_SNV_STRELKA2} > ${VCF_SNV_STRELKA2_PASS}
```
```
# For indels
VCF_INDEL_STRELKA2=$OUT_FOLDER1/results/variants/somatic.indels.vcf.gz
VCF_INDEL_STRELKA2_PASS=$(echo $VCF_INDEL_STRELKA2 | sed 's/.vcf.gz/.pass.vcf/g')
zgrep '#\|PASS' ${VCF_INDEL_STRELKA2} > ${VCF_INDEL_STRELKA2_PASS}
```


## GATK4 - MuTect2

```
GATK=gatk

REFGENOME=${FOLDER}/source_data/chr7.fa

# Main folder
FOLDER=/home/training/Documents/DNAseq_short_reads/

SNP1000G=${FOLDER}/source_data/1000G_phase1.snps.high_confidence.hg38.vcf.gz
BED=${FOLDER}/source_data/Targeted_regions.chr7.bed

OUT_FOLDER2=${FOLDER}/practical_results/variant_calling/mutect2_results
mkdir -p $OUT_FOLDER2
```

```
#---
# 1. Estimate Normal/tumour contamination
#---

## Generate GATK pileup tables

# Pileup table for the tumor sample
${GATK} GetPileupSummaries \
   -I ${BAM_T} \
   -V ${SNP1000G} \
   -L chr7 \
   -O $OUT_FOLDER2/COLO829T.pileups.table

# Pileup table for the normal sample 
${GATK} GetPileupSummaries \
   -I ${BAM_N} \
   -V ${SNP1000G} \
   -L chr7 \
   -O $OUT_FOLDER2/COLO829BL.pileups.table

# Esitmate contamination
${GATK} CalculateContamination \
   -I $OUT_FOLDER2/COLO829T.pileups.table \
   -matched $OUT_FOLDER2/COLO829BL.pileups.table \
   -O $OUT_FOLDER2/contamination.table


#---
# 2. Get coverage per sample
#---

## Get depth of coverage
${GATK} DepthOfCoverage \
   -R ${REFGENOME} \
   -O ${OUT_FOLDER2}/Coverage \
   -I ${BAM_T} \
   -I ${BAM_N} \
   -L ${BED} \
   --omit-interval-statistics \
   --omit-depth-output-at-each-base

#---
# 3. Run MuTect2
#---

## Run MuTect2 
# Variants MuTecT2
${GATK} Mutect2 \
   -R ${REFGENOME} \
   -I ${BAM_T} \
   -I ${BAM_N} \
   -normal COLO829BL \
   -tumor COLO829T \
   --germline-resource ${SNP1000G} \
   --panel-of-normals ${FOLDER}/source_data/1000g_pon.hg38.vcf.gz \
   -O ${OUT_FOLDER2}/mutect2.vcf \
   -L chr7

# Filtering
${GATK} FilterMutectCalls \
   -R ${REFGENOME} \
   -V ${OUT_FOLDER2}/mutect2.vcf \
   --contamination-table $OUT_FOLDER2/contamination.table \
   --min-allele-fraction 0.01 \
   -O ${OUT_FOLDER2}/mutect2.filtered.vcf

# Get only PASS calls
grep '#\|PASS' ${OUT_FOLDER2}/mutect2.filtered.vcf > ${OUT_FOLDER2}/mutect2.filtered.pass.vcf
```
