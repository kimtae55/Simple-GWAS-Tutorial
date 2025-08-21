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
Step 1: Handle missingness per individual and per SNP: Delete individuals with missingness >0.05

Step 2: Check for sex discrepancy: Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

Step 3:

## Step 1:

