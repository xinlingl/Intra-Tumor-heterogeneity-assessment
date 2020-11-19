# Tumor-heterogeneity-assessment

![alt tag](https://github.com/xinlingl/Tumor-heterogeneity-assessment/blob/main/workflow.jpg)

This tool use Bayesian clustering method to predict tumor heterogeneity from whole exsome sequencing data of patients with solid tumor 

step 1: Index the input whole exome sequencing files

Using samtools to index the whole exome sequencing files (.bam) of tumor and normal sample.



step 2: Extract location of mutations information from vcf file

Use bcftools to extract chromosome number and position of in the vcf file. Mutations on X and Y chromosomes and mitochondria are excluded. 



step 3: Calculate read counts

Use alleleCount to calculate count of A, G, C, T, and read depth at each locus in both normal and tumor sample. 



step 4: Calculate log ratio of total signal intensity between tumor and normal sample and B-allele frequency 


step 5: Generate rounded copy number of major and minor allele of each segment of genome


step 6: Extract allelic counts from mutation file 


step 7: Combine copy number information and allelic counts 


step 8: Generate cluster label and cellular prevalence of each somatic mutation using 
