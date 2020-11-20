# Tumor-heterogeneity-assessment

![alt tag](https://github.com/xinlingl/Tumor-heterogeneity-assessment/blob/main/workflow.jpg)

This tool use Bayesian clustering method to predict tumor heterogeneity from whole exsome sequencing data of patients with solid tumor 
<br />

**step 1: Index the input whole exome sequencing files**

Using samtools to index the whole exome sequencing files (.bam) of tumor and normal sample.

<br />

**step 2: Extract location of mutations information from vcf file**

Use bcftools to extract chromosome number and position of in the vcf file. Mutations on X and Y chromosomes and mitochondria are excluded.

<br />

**step 3: Calculate read counts**

Use alleleCount to calculate count of A, G, C, T, and read depth at each locus in both normal and tumor sample.

<br />

**step 4: Calculate log ratio of total signal intensity between tumor and normal sample and B-allele frequency**
<br />

**step 5: Generate rounded copy number of major and minor allele of each segment of genome**
<br />

**step 6: Extract allelic counts from mutation file**

Use bcftools to extract reference and variance allele counts of each somatic mutation in the mutation file. Somatic mutations on X and Y chromosome and mitochondria are excluded for simplicity.
<br />

**step 7: Combine copy number information and allelic counts**
<br />

**step 9: Generate cluster label and cellular prevalence of each somatic mutation**

Use PyClone with burin as 10,000 and default thin which is the number of samples to discard as burning for the MCMC chain. 
