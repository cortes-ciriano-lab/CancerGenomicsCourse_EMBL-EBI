# Variant annotation
Variant annotation refers to the process of providing additional contextual information and functional interpretation for identified genetic variants. It involves assigning biological and clinical annotations to variants, such as their location in the genome, impact on gene function, population frequency and predicted pathogenicity.

In this section, we will learn to annotate somatic variant calls using [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to find pathogenic variants in cancer samples. 

This workflow continues with the results obtained from the [variant calling exercise](https://github.com/cortes-ciriano-lab/CancerGenomicsCourse_EMBL-EBI/blob/main/docs/VariantCalling.md#variant-calling), where we detected and filtered somatic variants in the COLO829 sample with two algorithms ([Strelka2](https://www.nature.com/articles/s41592-018-0051-x) and [MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-))

### Set the environment for running this workflow
- Annovar tool and databases to be used
```
ANNOVAR=/home/training/annovar
hummandb=/home/training/annovar/humandb
```

- Main working directory and out folder for the annotation files
```
# Main folder
FOLDER=/home/training/Documents/DNAseq_short_reads/
# Out folder
OUT_FOLDER3=${FOLDER}/practical_results/variant_calling/variant_annotation
mkdir -p $OUT_FOLDER3
```

## Strelka2 - somatic calls

### Strelka2 - SNV calls

- Path to vcf file with the PASS SNVs
```
VCF1=${FOLDER}/practical_results/variant_calling/strelka2_results/results/variants/somatic.snvs.pass.vcf
```
- Prepare vcf format to ANNOVAR format
```
perl $ANNOVAR/convert2annovar.pl -format vcf4old $VCF1 > ${OUT_FOLDER3}/Strelka2.snvs.pass.avinput
```
- Annotate variants
```
perl $ANNOVAR/table_annovar.pl ${OUT_FOLDER3}/Strelka2.snvs.pass.avinput \
   $hummandb -buildver hg38 \
   -out ${OUT_FOLDER3}/Strelka2.snvs.annovar \
   -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
   -nastring . -csvout -polish --otherinfo
```

### Strelka2 - Indel calls
- Path to vcf file with the PASS indels
```
VCF2=${FOLDER}/practical_results/variant_calling/strelka2_results/results/variants/somatic.indels.pass.vcf
```
- Prepare vcf format to ANNOVAR format
```
perl $ANNOVAR/convert2annovar.pl -format vcf4old $VCF2 > ${OUT_FOLDER3}/Strelka2.indels.pass.avinput
```
- Annotate variants
```
perl $ANNOVAR/table_annovar.pl ${OUT_FOLDER3}/Strelka2.indels.pass.avinput \
   $hummandb -buildver hg38 \
   -out ${OUT_FOLDER3}/Strelka2.indels.annovar \
   -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
   -nastring . -csvout -polish --otherinfo
```
