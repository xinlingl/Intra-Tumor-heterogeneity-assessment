# Intra-tumor-heterogeneity-assessment

![alt tag](https://github.com/xinlingl/Tumor-heterogeneity-assessment/blob/main/workflow.jpg)

This tool uses Bayesian clustering method to predict intra-tumor heterogeneity from whole exsome sequencing data of patients with solid tumor 

<br />

**step 1: Index the input whole exome sequencing files**

Using samtools to index the whole exome sequencing files (.bam) of tumor and normal sample.

<br />

**step 2: Extract location of mutations information from VCF file**

Use bcftools to extract chromosome numbers and positions of somatic mutations in a VCF file. 

<br />

**step 3: Calculate read counts**

Use alleleCount to calculate count of A, G, C, T, and read depth at each locus in both normal and tumor sample.

<br />

**step 4: Calculate log ratio of total signal intensity between tumor and normal sample and B-allele frequency**

<br />

**step 5: Generate rounded copy number of major and minor allele of each segment of genome**

<br />

**step 6: Extract allelic counts from mutation file**

Use bcftools to extract reference and variance allele counts of each somatic mutation in the mutation file. 

<br />

**step 7: Combine copy number information and allelic counts**

<br />

**step 8: Generate cluster label and cellular prevalence of each somatic mutation**

Use PyClone with burin as 10,000 and default thin which is the number of samples to discard as burning for the MCMC chain. 
