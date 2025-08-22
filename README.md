# GWAS-End-to-End-Tutorial
This is a Step by Step Tutorial for GWAS (for Plink format), explaining what SNP data looks like, how to perform quality control and imputation procedures, and how to run GWAS.
Credits for the figures and explanations here go to the more in-depth tutorials: [GWASTutorial](https://cloufield.github.io/GWASTutorial) or [doi: 10.1002/mpr.1608](https://pmc.ncbi.nlm.nih.gov/articles/PMC6001694/). 

I aggregate all components of GWAS into one repository, to clear the confusion of what to download and when/how to run certain scripts. This means that you can replace my data with your .bed, .bim, .fam files and replicate the whole experiment. 

## Prerequisites
You need to have the following installed: ```plink/plink2, R ```

## Assumption: You have SNP data in PLINK format
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/plink.png" width="600">

In my case, I use Plink data (.bed, .bim, .fam files) downloaded from the ADNI1, ADNI2, and ADNIGO studies.

## What are the necessary steps in a GWAS?

### Quality Control --> Population Structure Modeling --> Lifting and Imputation --> Associative Analysis 

Quality Control is done at a sample-level (to remove bad individuals; e.g. contamination, swaps, relatedness, sex mismatches) and SNP-level (to remove bad variatns; e.g. missingness, low MAF, HWE failures).

### QC Steps:

Step 1: Handle missingness per individual and per SNP: Delete individuals with missingness >0.05.

Step 2: Handle sex discrepancy: Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

Step 3: Extract autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF <0.01).

Step 4: Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).

Step 5: Remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

Step 6: We exclude all individuals with a PI_HAT > 0.2 to remove cryptic relatedness, assuming a random population sample.

## Population Structure Modeling:

## Lifting and Imputation:

## Associative Analysis:

