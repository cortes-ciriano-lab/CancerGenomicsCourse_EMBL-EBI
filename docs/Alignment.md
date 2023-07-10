## Alignment workflow - Cancer Genomic Course 2023 - EMBL-EBI
This section will outline the various stages of the alignment workflow before the variant calling step. This includes quality control (QC) and alignment steps to produce bam files suitable for the variant calling. 

The example data utilized is derived from [Valle-Inclan et al 2022](https://www.sciencedirect.com/science/article/pii/S2666979X22000726?via%3Dihub) and represents a subset of reads from the WGS data generated for the [COLO829](https://depmap.org/portal/cell_line/COLO829_SKIN?tab=overview) cancer cell line (small regions from chromosome 7). This example has a tumour and matched normal sample, so it is a good example of a standard tumour-matched-normal (paired) somatic variant calling approach.

<div align="center">
<img src="/docs/VC_workflow.png" width="30%">
</div>

### Set the environment for running this workflow

Let's start defining some variables that will be important for the workflow:

- Tools to be used
```
SEQPURGE=SeqPurge
BWA=bwa
SAMTOOLS=samtools
GATK=gatk
FASTQC=fastqc
```

- Main working directory
```
FOLDER=/home/training/Documents/DNAseq_short_reads/
```

- Reference genome used (intersected from Hg38 human reference genome)
```
REFGENOME=${FOLDER}/source_data/chr7.fa
```

- DbSNP vcf file with high-frequency germline variants
```
DBSNP=${FOLDER}/source_data/Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

- Out folder where we will generate and save all our results:
```
OUT_FOLDER=${FOLDER}/practical_results/alignment
mkdir -p $OUT_FOLDER
```

- Paths to the example fastq files
```
# Fastq reads tumour
read1_T=${FOLDER}/fastq_files/COLO829T.R1.fastq.gz
read2_T=${FOLDER}/fastq_files/COLO829T.R2.fastq.gz

# Fastq reads normal
read1_N=${FOLDER}/fastq_files/COLO829BL.R1.fastq.gz
read2_N=${FOLDER}/fastq_files/COLO829BL.R2.fastq.gz
```

### Initial Quality Control (QC) step
The first step of this workflow is to check the quality of the data we will analyse. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a tool that provides a quick assessment of the quality of high-throughput sequencing data, highlighting potential issues such as adapter contamination, sequencing errors, or overrepresented sequences.

```
OUT_FASTQC=$OUT_FOLDER/fastqc_results
mkdir -p $OUT_FASTQC

$FASTQC --noextract --nogroup -o $OUT_FASTQC ${FOLDER}/fastq_files/*.fastq.gz
```
- **QUESTION**: Could you quickly check the quality of your data (open the *html files in the browser)? Compare the generated results to the high- and low-quality examples shown [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Trimming adapters
The adapter trimming step is essential in sequencing workflows to remove adapter sequences introduced during library preparation. Adapters are short DNA sequences used to ligate sequencing primers to the target DNA fragments. Still, if left untrimmed, they can lead to false-positive alignments, decreased mapping efficiency, and affect downstream analyses. By accurately trimming adapters, researchers can improve the accuracy of sequence alignment, enhance mapping rates, and reduce the potential for bias or artefacts in downstream analysis pipelines. Here we use [SeqPurge](https://pubmed.ncbi.nlm.nih.gov/27161244/).

We do it first for the tumour reads:
```
OUT_TRIM=$OUT_FOLDER/trimmed_reads
mkdir -p $OUT_TRIM

## For tumour
trimmed1_T=$OUT_TRIM/COLO829T.R1.trimmed.fastq.gz
trimmed2_T=$OUT_TRIM/COLO829T.R2.trimmed.fastq.gz
${SEQPURGE} -in1 $read1_T -in2 $read2_T -out1 $trimmed1_T -out2 $trimmed2_T -qcut 0 -ncut 0 -threads 4
```

And now for the matched normal reads
```
## For normal
trimmed1_N=$OUT_TRIM/COLO829BL.R1.trimmed.fastq.gz
trimmed2_N=$OUT_TRIM/COLO829BL.R2.trimmed.fastq.gz
${SEQPURGE} -in1 $read1_N -in2 $read2_N -out1 $trimmed1_N -out2 $trimmed2_N -qcut 0 -ncut 0 -threads 4
```

### Alignment
Read alignment is the process of aligning short DNA sequences, called reads, to a reference genome or transcriptome. It involves finding the best match or mapping position on the reference for each read, enabling researchers to determine the origin and location of the reads within the genome and facilitating further analysis such as variant calling or gene expression quantification.

As we already prepared our raw fastq files, we will align all of them against a small fraction of the human reference genome (GRCh38 - hg38). For that, we use the tool [BWA-mem](https://github.com/lh3/bwa).

- First, the tumour reads:
```
OUT_BAMS=$OUT_FOLDER/bams
mkdir -p $OUT_BAMS

## For tumour
$BWA mem -t 4 ${REFGENOME} $trimmed1_T $trimmed2_T | $SAMTOOLS view -Shb - > $OUT_BAMS/COLO829T.bam
$SAMTOOLS sort -o $OUT_BAMS/COLO829T.sorted.bam $OUT_BAMS/COLO829T.bam
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.bam
```
- And now, the matched normal reads:
```
## For normal
$BWA mem -t 4 ${REFGENOME} $trimmed1_N $trimmed2_N | $SAMTOOLS view -Shb - > $OUT_BAMS/COLO829BL.bam
$SAMTOOLS sort -o $OUT_BAMS/COLO829BL.sorted.bam $OUT_BAMS/COLO829BL.bam
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.bam
```


- **QUESTION**: Take a look at the bam files and get familiar with this [format](https://samtools.github.io/hts-specs/SAMv1.pdf):
```
$SAMTOOLS view -h $OUT_BAMS/COLO829T.sorted.bam | less -S
```
- **QUESTION**: Open the BAM file and check some read [flag](https://samtools.github.io/hts-specs/SAMv1.pdf) values. You can use this [link](https://broadinstitute.github.io/picard/explain-flags.html) to understand the meaning of each value.   

### Preparing bam files for variant calling

#### Marking duplicates
Marking duplicates in a BAM file refers to the process of identifying and flagging PCR or optical duplicates, which are multiple reads derived from the same original DNA fragment. This information is useful in downstream analysis, such as variant calling, as it allows for more accurate estimation of library complexity and reduces potential biases introduced by PCR amplification or sequencing artefacts. This step is run using [GATK4](https://gatk.broadinstitute.org/hc/en-us)

- For tumour:
```
## For tumour
$GATK MarkDuplicates -I $OUT_BAMS/COLO829T.sorted.bam -O $OUT_BAMS/COLO829T.sorted.MD.bam -M $OUT_BAMS/COLO829T.marked_dup_metrics.txt -REMOVE_DUPLICATES true --TMP_DIR $OUT_BAMS
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.MD.bam
```

- For matched normal:
```
## For normal
$GATK MarkDuplicates -I $OUT_BAMS/COLO829BL.sorted.bam -O $OUT_BAMS/COLO829BL.sorted.MD.bam -M $OUT_BAMS/COLO829BL.marked_dup_metrics.txt -REMOVE_DUPLICATES true --TMP_DIR $OUT_BAMS
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.MD.bam
```

#### Adding read groups
This step assigns all the reads in a file to a single new read-group. It is usually required by the [GATK4](https://gatk.broadinstitute.org/hc/en-us) tools that we will use in the downstream analysis. 

- For tumour:
```
## For tumour
# Add read groups
$GATK AddOrReplaceReadGroups -I $OUT_BAMS/COLO829T.sorted.MD.bam -O $OUT_BAMS/COLO829T.sorted.MD.RG.bam -LB COLO829T -PL illumina -PU COLO829T -SM COLO829T
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.MD.RG.bam
```

- For matched-normal:
```
## For matched-normal
$GATK AddOrReplaceReadGroups -I $OUT_BAMS/COLO829BL.sorted.MD.bam -O $OUT_BAMS/COLO829BL.sorted.MD.RG.bam -LB COLO829BL -PL illumina -PU COLO829BL -SM COLO829BL
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.MD.RG.bam
```

#### Base quality score recalibration (BQSR)
Base quality score recalibration (BQSR, [GATK4](https://gatk.broadinstitute.org/hc/en-us)) is a process in which we apply machine learning to model potential errors empirically and adjust the quality scores accordingly. It is not mandatory but highly recommended for GATK variant calling pipelines. It requires a set of known variant sites, and in this case, we use the variants found in the [dbSNP](https://en.wikipedia.org/wiki/DbSNP) database.

- For tumour:
```
## For tumour
$GATK BaseRecalibrator -I $OUT_BAMS/COLO829T.sorted.MD.RG.bam  -R $REFGENOME --known-sites $DBSNP -O $OUT_BAMS/COLO829T.recal_data.table
$GATK ApplyBQSR  -R $REFGENOME -I $OUT_BAMS/COLO829T.sorted.MD.RG.bam --bqsr-recal-file $OUT_BAMS/COLO829T.recal_data.table -O $OUT_BAMS/COLO829T.sorted.MD.RG.BQSR.bam
$SAMTOOLS index $OUT_BAMS/COLO829T.sorted.MD.RG.BQSR.bam
```

- For matched-normal
```
## For matched-normal
$GATK BaseRecalibrator -I $OUT_BAMS/COLO829BL.sorted.MD.RG.bam  -R $REFGENOME --known-sites $DBSNP -O $OUT_BAMS/COLO829BL.recal_data.table
$GATK ApplyBQSR  -R $REFGENOME -I $OUT_BAMS/COLO829BL.sorted.MD.RG.bam --bqsr-recal-file $OUT_BAMS/COLO829BL.recal_data.table -O $OUT_BAMS/COLO829BL.sorted.MD.RG.BQSR.bam
$SAMTOOLS index $OUT_BAMS/COLO829BL.sorted.MD.RG.BQSR.bam
```

- **QUESTION**: Look at all the bam files generated in each workflow step. The size is very different, while the number of reads and sequences does not change. Play with samtools to check for differences and run *samtools flagstat* to observe the characteristics of the reads in the bam file:
```
$SAMTOOLS flagstat $OUT_BAMS/COLO829T.sorted.MD.RG.BQSR.bam
$SAMTOOLS flagstat $OUT_BAMS/COLO829BL.sorted.MD.RG.BQSR.bam
```

