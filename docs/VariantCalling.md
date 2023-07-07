# Variant calling
Somatic (short) variant calling is the process of identifying genetic variations specific to tumour cells by comparing the genomic data from the tumour and matched normal samples. It involves detecting somatic mutations, such as single-nucleotide variants (SNVs) and/or small insertions/deletions (indels), that are acquired during tumorigenesis and distinguishing them from germline variants present in both tumour and normal cells.

In this section, we will explain how to run different variant calling algorithms ([Strelka2](https://www.nature.com/articles/s41592-018-0051-x) and [MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-)), as well as how to filter the raw variants to get a high-quality set of somatic mutations. 

This workflow continues with the results obtained from the [alignment exercise](https://github.com/cortes-ciriano-lab/CancerGenomicsCourse_EMBL-EBI/blob/main/docs/Alignment.md#alignment-workflow---cancer-genomic-course-2023---embl-ebi), where we generated the COLO829T (tumour) and COLO829BL (matched-normal) bam files. 

## Strelka2 - Somatic workflow
[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized to analyse germline variation in small cohorts and somatic variation in tumour/normal sample pairs. The somatic calling model improves on the original Strelka liquid and late-stage tumour analysis method by accounting for possible tumour cell contamination in the normal sample. A final empirical variant re-scoring step using random forest models trained on various call quality features has been added to both callers to improve precision further.

### Set the environment for running this workflow

- Path to Strelka2
```
STRELKA2=/home/training/strelka-2.9.10.centos6_x86_64/
```

- Main working directory and out folder for the results of this algorithm
```
# Main folder
FOLDER=/home/training/Documents/DNAseq_short_reads/
# Out folder
OUT_FOLDER1=${FOLDER}/practical_results/variant_calling/strelka2_results
mkdir -p $OUT_FOLDER1
rm -r $OUT_FOLDER1/*
```
- Human reference genome (intersection from Hg38)
```
# Reference genome
REFGENOME=${FOLDER}/source_data/chr7.fa
```

- Tumour and matched-normal bam files
```
# Matched-normal sample
BAM_N=${FOLDER}/practical_results/alignment/bams/COLO829BL.sorted.MD.RG.BQSR.bam
# Tumour sample
BAM_T=${FOLDER}/practical_results/alignment/bams/COLO829T.sorted.MD.RG.bam
```

### Running Strelka2 somatic workflow
Strelka2 is run in two steps: (1) configuration (specifying input data and options) and (2) workflow execution (specifying parameters on how Strelka2 is executed). 

#### Configuration
Configuration (specifying input data and options):
```
${STRELKA2}/bin/configureStrelkaSomaticWorkflow.py \
--normalBam $BAM_N \
--tumorBam $BAM_T \
--referenceFasta $REFGENOME \
--runDir $OUT_FOLDER1
```

#### Workflow execution
Here we specify parameters on how Strelka2 is executed. This step can also be interrupted and restarted without changing the final result of the workflow.
```
${OUT_FOLDER1}/runWorkflow.py  -m local -j 4 -g 15
```
- **QUESTION**: Check the results (vcf files), trying to understand the different fields of each column:
```
zless -S $OUT_FOLDER1/results/variants/somatic.snvs.vcf.gz
```

#### Keep PASS mutations
The final output of Strelka2 has all the potential somatic variants detected by the tool. However, many of them have been filtered by the tool and represent low-quality (or false) variants that we want to remove. You can check this in the column FILTER of the vcf file. For these reasons, we only want to keep the *PASS* (high-quality) mutations. To do so, we need to run: 

- For SNVs:
```
# For SNVs
VCF_SNV_STRELKA2=$OUT_FOLDER1/results/variants/somatic.snvs.vcf.gz
VCF_SNV_STRELKA2_PASS=$(echo $VCF_SNV_STRELKA2 | sed 's/.vcf.gz/.pass.vcf/g')
zgrep '#\|PASS' ${VCF_SNV_STRELKA2} > ${VCF_SNV_STRELKA2_PASS}
```

- For indels
```
# For indels
VCF_INDEL_STRELKA2=$OUT_FOLDER1/results/variants/somatic.indels.vcf.gz
VCF_INDEL_STRELKA2_PASS=$(echo $VCF_INDEL_STRELKA2 | sed 's/.vcf.gz/.pass.vcf/g')
zgrep '#\|PASS' ${VCF_INDEL_STRELKA2} > ${VCF_INDEL_STRELKA2_PASS}
```
- **QUESTION**: Compare the number of variants before and after keeping the PASS mutations using `grep`, `zgrep` and `wc -l` commands.  

## GATK4 - MuTect2
[MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) call short somatic mutations via local assembly of haplotypes. The caller uses a Bayesian somatic genotyping model and uses the assembly-based machinery of [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) (GATK4 germline variant caller). 

MuTect2 is split into different steps to understand better the background noise in the samples and obtain high-quality somatic mutations. Computationally, this tool is much slowlier than Strelka2, especially with WGS data. However, it is one of the most used tools in the field.

### Set the environment for running this workflow

- GATK4 tool
```
GATK=gatk
```

- Main working directory and out folder for the results of this algorithm
```
# Main folder
FOLDER=/home/training/Documents/DNAseq_short_reads/
OUT_FOLDER2=${FOLDER}/practical_results/variant_calling/mutect2_results
mkdir -p $OUT_FOLDER2
```

- Human reference genome (intersection from Hg38)
```
REFGENOME=${FOLDER}/source_data/chr7.fa
```
- Bed file (genomic coordinates) of the regions sampled for this exercise
```
BED=${FOLDER}/source_data/Targeted_regions.chr7.bed
```
- High confidence set of 1000G variants (vcf file)
```
SNP1000G=${FOLDER}/source_data/1000G_phase1.snps.high_confidence.hg38.vcf.gz
```
- Panel Of Normals (PoN) already generated by the GATK community (vcf file)
```
PON=${FOLDER}/source_data/1000g_pon.hg38.vcf.gz
```

### Estimating normal/tumour contamination

- Generating GATK pileup tables in known variant sites
```
# Pileup table for the tumour sample
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
```
- Estimating contamination
```
${GATK} CalculateContamination \
   -I $OUT_FOLDER2/COLO829T.pileups.table \
   -matched $OUT_FOLDER2/COLO829BL.pileups.table \
   -O $OUT_FOLDER2/contamination.table
```

### Getting coverage per sample
```
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
   --panel-of-normals ${PON} \
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
