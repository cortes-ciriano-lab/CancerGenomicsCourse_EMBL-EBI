# Variant annotation
Variant annotation refers to the process of providing additional contextual information and functional interpretation for identified genetic variants. It involves assigning biological and clinical annotations to variants, such as their location in the genome, impact on gene function, population frequency and predicted pathogenicity.

In this section, we will learn to annotate somatic variant calls using [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to find pathogenic variants in cancer samples. 

This workflow continues with the results obtained from the [variant calling exercise](https://github.com/cortes-ciriano-lab/CancerGenomicsCourse_EMBL-EBI/blob/main/docs/VariantCalling.md#variant-calling), where we detected and filtered somatic variants in the COLO829 sample with two algorithms ([Strelka2](https://www.nature.com/articles/s41592-018-0051-x) and [MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-))

<div align="center">
<img src="/docs/VC_workflow.png" width="30%">
</div>

**Important**

In today's example, we will use this set of database annotations: *refGene, cytoBand, exac03, avsnp147, dbnsfp30a and gnomad_genome*. However, [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) is very flexible and can work with many others. Find in this [link](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) how to download them, and [here](https://annovar.openbioinformatics.org/en/latest/user-guide/download/#additional-databases) other additional ones that might be interesting for you.

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

## GATK - MuTect2 calls
- Path to vcf file with the PASS indels
```
VCF3=${FOLDER}/practical_results/variant_calling/mutect2_results/mutect2.filtered.pass.vcf
```
- Prepare vcf format to ANNOVAR format
```
perl $ANNOVAR/convert2annovar.pl -format vcf4old $VCF3 > ${OUT_FOLDER3}/Mutect2.pass.avinput
```
- Annotate variants
```
perl $ANNOVAR/table_annovar.pl ${OUT_FOLDER3}/Mutect2.pass.avinput \
   $hummandb -buildver hg38 \
   -out ${OUT_FOLDER3}/Mutect2.annovar \
   -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
   -nastring . -csvout -polish --otherinfo
```

### Looking for pathogenic cancer variants
This is an illustrative example of an ANNOVAR output looks like (`head -10 ${OUT_FOLDER3}/Mutect2.annovar.hg38_multianno.csv `):
```
Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,cytoBand,ExAC_ALL,ExAC_AFR,ExAC_AMR,ExAC_EAS,ExAC_FIN,ExAC_NFE,ExAC_OTH,ExAC_SAS,avsnp147,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,DANN_score,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,integrated_fitCons_score,integrated_confidence_value,GERP++_RS,phyloP7way_vertebrate,phyloP20way_mammalian,phastCons7way_vertebrate,phastCons20way_mammalian,SiPhy_29way_logOdds,gnomAD_genome_ALL,gnomAD_genome_AFR,gnomAD_genome_AMR,gnomAD_genome_ASJ,gnomAD_genome_EAS,gnomAD_genome_FIN,gnomAD_genome_NFE,gnomAD_genome_OTH,Otherinfo
chr7,22257,22257,A,G,"intergenic","NONE;LOC102723672","dist=NONE;dist=122188",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	54	93"
chr7,37501,37501,A,C,"intergenic","NONE;LOC102723672","dist=NONE;dist=106944",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,0.0004,0.0006,0,0,0,0.0009,0.0002,0,"hom	.	65	93"
chr7,51999,51999,A,T,"intergenic","NONE;LOC102723672","dist=NONE;dist=92446",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	142	93"
chr7,174700,174700,C,T,"upstream","LOC105375115","dist=220",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	154	93"
chr7,180027,180027,G,A,"intergenic","LOC105375115;FAM20C","dist=4014;dist=12942",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	124	93"
chr7,180109,180109,C,A,"intergenic","LOC105375115;FAM20C","dist=4096;dist=12860",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	128	93"
chr7,182440,182440,G,T,"intergenic","LOC105375115;FAM20C","dist=6427;dist=10529",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	84	78"
chr7,274049,274049,C,T,"intergenic","FAM20C;LOC442497","dist=13275;dist=105376",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	149	93"
chr7,313487,313487,G,A,"intergenic","FAM20C;LOC442497","dist=52713;dist=65938",.,.,"7p22.3",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"hom	.	150	93"
```

When looking for pathogenic variants, we have to check a few things:

- Filtering based on Impact Scores: check the section `LJB* (dbNSFP) non-synonymous variants annotation` [here](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/) for detecting the most deleterious mutations
- Functional Impact: focus on variants that have a high predicted functional impact on the gene or protein (e.g., stop-gain, missense, frameshift, splice site)
- Consider frame-shift mutations as highly deleterious
- Population Frequency: Consider the frequency of the variant in population databases like gnomAD or ExAC. Rare variants are more likely to be pathogenic, but this should be interpreted in the context of the specific disease and inheritance pattern.
- As we are working with cancer samples, if you have many mutations left after all the previous filters, start looking at cancer genes. For instance, check the [CENSUS genes](https://cancer.sanger.ac.uk/census) or external databases like [IntoGen](https://www.intogen.org)

#### **QUESTION**

Open the variant annotated file that you prefer (*Mutect2.annovar.hg38_multianno.csv* or *Strelka2.snvs.annovar.hg38_multianno.csv*) in your most desired viewer (excel, R...). Do you find any interesting somatic mutation in COLO829?
