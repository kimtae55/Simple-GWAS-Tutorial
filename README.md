# Simple GWAS Tutorial
This is a Step by Step Tutorial for GWAS (for Plink format), explaining what SNP data looks like, how to perform quality control and imputation procedures, and how to run GWAS.
Credits for the figures and explanations here go to the more in-depth tutorials: [GWASTutorial1](https://cloufield.github.io/GWASTutorial), [GWASTutorial2](https://www.ncbi.nlm.nih.gov/pubmed/29484742), or [GWASTutorial3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6001694/). The codes are taken from [https://github.com/MareesAT/GWA_tutorial](https://github.com/MareesAT/GWA_tutorial).

I aggregate all components of GWAS into one repository for study purposes, and to clear the confusion of what to download and when/how to run certain scripts. This means that you can replace my data with your .bed, .bim, .fam files and replicate the whole experiment. 

The output is a genotype matrix of shape (n x p), where n is the number of subjects, and p is the number of SNPs extracted from GWAS. 

## Prerequisites
You need to have the following installed: ```plink/plink2, R ```

## Assumption: You have SNP data in PLINK format
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/plink.png" width="600">

In my case, I use Plink data (.bed, .bim, .fam files) downloaded from the ADNI1, ADNI2, and ADNIGO studies. For simplicity, all data, scripts, and outputs will be located in a single working directory. 

## What are the necessary steps in a GWAS?

### Quality Control --> Liftover and Imputation --> Population Structure Modeling --> Associative Analysis 

Quality Control is done at a sample-level (to remove bad individuals; e.g. contamination, swaps, relatedness, sex mismatches) and SNP-level (to remove bad variants; e.g. missingness, low MAF, HWE failures).
If you have multiple datasets from different studies/timepoints, you should run this part separately, then merge the results (see the later optional step for merging)

### QC Steps:

Step 1: Handle missingness per individual and per SNP: Delete individuals with missingness >0.05.
```
> plink --bfile ADNI_cluster_01_forward_757LONI --missing 
> Rscript --no-save hist_miss.R # Visualize missingness
> plink --bfile ADNI_cluster_01_forward_757LONI --geno 0.05 --make-bed --out ADNI_cluster_01_forward_757LONI_2
> plink --bfile ADNI_cluster_01_forward_757LONI_2 --mind 0.05 --make-bed --out ADNI_cluster_01_forward_757LONI_3
```

Step 2: Handle sex discrepancy: Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.
```
> plink --bfile ADNI_cluster_01_forward_757LONI_3 --check-sex 
> Rscript --no-save gender_check.R # Visualize sex dicrepancy check
> grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
> plink --bfile ADNI_cluster_01_forward_757LONI_3 --remove sex_discrepancy.txt --make-bed --out ADNI_cluster_01_forward_757LONI_4 
```

Step 3: Extract autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF <0.01).
```
> awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ADNI_cluster_01_forward_757LONI_4.bim > snp_1_22.txt
plink --bfile ADNI_cluster_01_forward_757LONI_4 --extract snp_1_22.txt --make-bed --out ADNI_cluster_01_forward_757LONI_5
> plink --bfile ADNI_cluster_01_forward_757LONI_5 --freq --out MAF_check
> Rscript --no-save MAF_check.R
> plink --bfile ADNI_cluster_01_forward_757LONI_5 --maf 0.01 --make-bed --out ADNI_cluster_01_forward_757LONI_6
```

Step 4: Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
```
> plink --bfile ADNI_cluster_01_forward_757LONI_6 --hardy
> awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
> Rscript --no-save hwe.R
> plink --bfile ADNI_cluster_01_forward_757LONI_6 --hwe 1e-6 --hwe-all --make-bed --out ADNI_cluster_01_forward_757LONI_7
```

Step 5: Remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.
```
> plink --bfile ADNI_cluster_01_forward_757LONI_7 --exclude range high-LD-regions-hg18-NCBI36.txt --indep-pairwise 50 5 0.2 --out indepSNP
> plink --bfile ADNI_cluster_01_forward_757LONI_7 --extract indepSNP.prune.in --het --out R_check
> Rscript --no-save check_heterozygosity_rate.R
> Rscript --no-save heterozygosity_outliers_list.R # # Outputs fail-het-qc.txt
> sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
> plink --bfile ADNI_cluster_01_forward_757LONI_7 --remove het_fail_ind.txt --make-bed --out ADNI_cluster_01_forward_757LONI_8 
```

Step 6: We exclude all individuals with a PI_HAT > 0.2 to remove cryptic relatedness, assuming a random population sample.
```
> plink --bfile ADNI_cluster_01_forward_757LONI_8 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
> awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
> Rscript --no-save Relatedness.R
> plink --bfile ADNI_cluster_01_forward_757LONI_8 --filter-founders --make-bed --out ADNI_cluster_01_forward_757LONI_9
> plink --bfile ADNI_cluster_01_forward_757LONI_9 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
> plink --bfile ADNI_cluster_01_forward_757LONI_9 --missing
> plink --bfile ADNI_cluster_01_forward_757LONI_9 --missing --out ADNI_cluster_01_forward_757LONI_9
> awk 'NR==FNR && FNR>1 {miss[$1" "$2]=$6; next} 
     FNR>1 {a=$1" "$2; b=$3" "$4; fa=(a in miss?miss[a]:1.0); fb=(b in miss?miss[b]:1.0); 
            if (fa>fb) print $1, $2; else print $3, $4}' 
    ADNI_cluster_01_forward_757LONI_9.imiss pihat_min0.2_in_founders.genome 
  | sort -u > remove_pihat0.2_lowcall.txt
> plink --bfile ADNI_cluster_01_forward_757LONI_9 --remove remove_pihat0.2_lowcall.txt --make-bed --out ADNI_cluster_01_forward_757LONI_10
```

### Lifting and Imputation:

Step 1: ADNI datasets are often on older genome builds (e.g., hg18/NCBI36). Before imputation, convert to hg19/GRCh37. This ensures all datasets use the same genome coordinates, resulting in .bed/.bim/.fam files aligned to GRCh37. 
```
# https://genviz.org/module-01-intro/0001/06/02/liftoverTools/
# https://genome.sph.umich.edu/wiki/LiftOver
# download the LiftOver from http://genome.ucsc.edu/cgi-bin/hgLiftOver
# download chain file hg18ToHg19.over.chain.gz from http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/
# Download files from: https://martha-labbook.netlify.app/posts/converting-snp-array-data/
> chmod a+x liftOver
> sudo mv ~/Downloads/liftOver /usr/local/bin/
> plink --bfile ADNI_cluster_01_forward_757LONI_10 --recode --tab --out ADNI_cluster_01_forward_757LONI_11
> export SNP_path="/Users/taehyo/Library/CloudStorage/Dropbox/NYU/Research/Research/Data/data_l0ipls/GWAS/taehyo/scripts"
> python $SNP_path/liftOverPlink-master/liftOverPlink.py --bin /usr/local/bin/liftOver --map ADNI_cluster_01_forward_757LONI_11.map --out ADNI_cluster_01_forward_757LONI_11_lifted --chain $SNP_path/hg18ToHg19.over.chain.gz
> python $SNP_path/liftoverPlink-master/rmBadLifts.py --map ADNI_cluster_01_forward_757LONI_11_lifted.map --out ADNI_cluster_01_forward_757LONI_11_good_lifted.map --log ADNI_cluster_01_forward_757LONI_11_bad_lifted.dat
> cut -f 2 ADNI_cluster_01_forward_757LONI_11_bad_lifted.dat > ADNI_cluster_01_forward_757LONI_11_snps_exclude.dat
> cut -f 4 ADNI_cluster_01_forward_757LONI_11_lifted.bed.unlifted | sed "/^#/d" >> ADNI_cluster_01_forward_757LONI_11_snps_exclude.dat 
> plink --file ADNI_cluster_01_forward_757LONI_11 --recode --out ADNI_cluster_01_forward_757LONI_12 --exclude ADNI_cluster_01_forward_757LONI_11_snps_exclude.dat
> plink --ped ADNI_cluster_01_forward_757LONI_12.ped --map ADNI_cluster_01_forward_757LONI_11_good_lifted.map --recode --out ADNI_cluster_01_forward_757LONI_13
> plink --file ADNI_cluster_01_forward_757LONI_13 --make-bed --out ADNI_cluster_01_forward_757LONI_14
> plink --bfile ADNI_cluster_01_forward_757LONI_14 --recode vcf --out ADNI_cluster_01_forward_757LONI_14_vcf
'''
# Python code to heck if it has been lifted properly
from snps import SNPs
s = SNPs("ADNI_cluster_01_forward_757LONI_14_vcf.vcf")
s.build
s.build_detected
s.assembly #'GRCh37'
'''
```

Step 2: Check for duplicates
```
plink --bfile ADNI_cluster_01_forward_757LONI_14 --list-duplicate-vars 
plink --bfile ADNI_cluster_01_forward_757LONI_14 --exclude duplicatedSNP.txt --make-bed --out ADNI_cluster_01_forward_757LONI_15
```

Step 3: Imputation via Michigan Imputation Server
We will convert the plink data to vcf format, then use the MIS to impute data. 
```
> wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
> wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
> unzip HRC-1000G-check-bim-v4.3.0.zip
> gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz 
> plink --freq --bfile ADNI_cluster_01_forward_757LONI_15 --out ADNI_cluster_01_forward_757LONI_15.freq
> perl $SNP_path/HRC-1000G-check-bim.pl -b ADNI_cluster_01_forward_757LONI_14.bim -f ADNI_cluster_01_forward_757LONI_14.freq.frq -r $SNP_path/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

## Create a sorted vcf.gz file using BCFtools:
# Download BCFtools (bcftools-1.14) from http://www.htslib.org/download/
> cd bcftools-1.14    # and similarly for bcftools and htslib
> ./configure --prefix=/Users/haishu/Applications/
> make
> make install
> export PATH=/Users/taehyo/Applications/bin:$PATH```

> for chr in {1..22}
do
	bcftools sort ADNI_cluster_01_forward_757LONI_15-updated-chr$chr.vcf -Oz -o ADNI_cluster_01_forward_757LONI_15-updated-chr$chr.vcf.gz
done

# Now upload to Michigan Imputation Server for imputation: https://imputationserver.sph.umich.edu/index.html#!pages/login

> for chr in $(seq 1 22)
do
	unzip -P 2kFzKqy7lKa1GR chr_$chr.zip
done

> bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou | 
bcftools view -Ou -i 'R2>0.3' |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o ADNI1_allchromosomes.converted.R2_0.3.vcf.gz

> plink --vcf ADNI1_allchromosomes.converted.R2_0.3.vcf.gz --double-id --make-bed --out ADNI1_allchromosomes.converted.R2_0.3
```

### Population Structure Modeling:

I repeat the above steps for ADNI1, ADNI2, and ADNIGO. If you only have a single dataset, you can skip the below merge process. 

(Optional Step): If you have multiple datasets, merge all of them now that they are aligned and imputed on the same genome build, followed by post-merge quality check. 


Step 1: Population stratification is corrected by extracting principal components (PCs) for each dataset separately using LD-pruned SNPs. The top PCs are used as covariates in GWAS to control for ancestry differences. 
```
Code will go here
```

### Associative Analysis:

Run GWAS on progression (e.g. ADAS-Cog, MMSE)
```
Code will go here
```

## Application: Sparse Canoninical Correlation Analysis using Imaging-Omics (FDG-PET, SNP)

Goal: To investigate how genetic variation (SNPs) and brain imaging features (FDG-PET ROIs) are related, and to assess how these associations change across disease stages (CN, MCI, AD). We use Sparse Canonical Correlation Analysis (SCCA) to identify low-dimensional, interpretable patterns linking high-dimensional SNP data with FDG-PET imaging phenotypes.

Input: 
- (n,p) SNP data matrix (split into CN, MCI, AD)
- (n,q) FDG-PET data matrix (split into CN, MCI, AD)
- Note that we did not use disease status (CN, MCI, AD) as outcomes for GWAS, because our goal is not to find SNPs associated with diagnosis, but rather SNPs linked to disease progression. Instead, we selected SNPs based on their association with continuous progression measures (e.g., cognitive decline, imaging biomarkers), so that the same SNP panel can be meaningfully compared across CN, MCI, and AD groups during SCCA.

Please read [this]() for the complete application details. 
