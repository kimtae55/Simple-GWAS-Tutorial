# Simple GWAS Tutorial
This is a Step by Step Tutorial for GWAS (for Plink format), explaining what SNP data looks like, how to perform quality control and imputation procedures, and how to run GWAS.
Credits for the figures and explanations here go to the more in-depth tutorials: [GWASTutorial1](https://cloufield.github.io/GWASTutorial), [GWASTutorial2](https://www.ncbi.nlm.nih.gov/pubmed/29484742), or [GWASTutorial3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6001694/). The codes are taken from [https://github.com/MareesAT/GWA_tutorial](https://github.com/MareesAT/GWA_tutorial), but modified to be plink2 compatible. 

I aggregate all components of GWAS into one repository for study purposes, and to clear the confusion of what to download and when/how to run certain scripts. This means that you can replace my data with your .bed, .bim, .fam files and replicate the whole experiment. 

The output is a genotype matrix of shape (n x p), where n is the number of subjects, and p is the number of SNPs extracted from GWAS. At the end of the tutorial is an application example of how one can utilize the processed SNP data.

## Prerequisites
You need to have the following installed: ```plink2, R```

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
> plink2 --bfile ADNI_cluster_01_forward_757LONI --missing 
> Rscript --no-save hist_miss.R # Visualize missingness
> plink2 --bfile ADNI_cluster_01_forward_757LONI --geno 0.05 --make-bed --out ADNI_cluster_01_forward_757LONI_2
> plink2 --bfile ADNI_cluster_01_forward_757LONI_2 --mind 0.05 --make-bed --out ADNI_cluster_01_forward_757LONI_3
```

Step 2: Handle sex discrepancy: Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.
```
> plink2 --bfile ADNI_cluster_01_forward_757LONI_3 --check-sex 
> Rscript --no-save gender_check.R # Visualize sex dicrepancy check
> grep "PROBLEM" plink2.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
> plink2 --bfile ADNI_cluster_01_forward_757LONI_3 --remove sex_discrepancy.txt --make-bed --out ADNI_cluster_01_forward_757LONI_4 
```

Step 3: Extract autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF <0.01).
```
> awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ADNI_cluster_01_forward_757LONI_4.bim > snp_1_22.txt
> plink2 --bfile ADNI_cluster_01_forward_757LONI_4 --extract snp_1_22.txt --make-bed --out ADNI_cluster_01_forward_757LONI_5
> plink2 --bfile ADNI_cluster_01_forward_757LONI_5 --freq --out MAF_check
> Rscript --no-save MAF_check.R
> plink2 --bfile ADNI_cluster_01_forward_757LONI_5 --maf 0.01 --make-bed --out ADNI_cluster_01_forward_757LONI_6
```

Step 4: Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
```
> plink2 --bfile ADNI_cluster_01_forward_757LONI_6 --hardy
> Rscript --no-save hwe.R
> plink2 --bfile ADNI_cluster_01_forward_757LONI_6 --hwe 1e-6 --make-bed --out ADNI_cluster_01_forward_757LONI_7
```

Step 5: Remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.
```
> plink2 --bfile ADNI_cluster_01_forward_757LONI_7 --exclude range high-LD-regions-hg18-NCBI36.txt --indep-pairwise 50 5 0.2 --out indepSNP
> plink2 --bfile ADNI_cluster_01_forward_757LONI_7 --extract indepSNP.prune.in --het --out R_check
> Rscript --no-save check_heterozygosity_rate.R
> plink2 --bfile ADNI_cluster_01_forward_757LONI_7 --remove het_fail_ind.txt --make-bed --out ADNI_cluster_01_forward_757LONI_8 
```

Step 6: We exclude all individuals with a PI_HAT > 0.2 to remove cryptic relatedness, assuming a random population sample.
```
# 1) LD-prune SNPs (used only for relatedness detection)
plink2 --bfile ADNI_cluster_01_forward_757LONI_8 --indep-pairwise 200 100 0.1 --out indepSNP

# 2) Use pruned SNPs to identify related individuals (temporary dataset)
plink2 --bfile ADNI_cluster_01_forward_757LONI_8 --extract indepSNP.prune.in --king-cutoff 0.10 --make-bed --out ADNI_relcheck_tmp

# 3) Save list of unrelated individuals
awk '{print $1, $2}' ADNI_relcheck_tmp.fam > unrelated.keep

# 4) Apply unrelated sample list to full dataset (keeps ~500k SNPs)
plink2 --bfile ADNI_cluster_01_forward_757LONI_8 --keep unrelated.keep --make-bed --out ADNI_cluster_01_forward_757LONI_10
```

### Liftover and Imputation:

Step 1: ADNI datasets are often on older genome builds (e.g., hg18/NCBI36). Before imputation, convert to hg19/GRCh37. This ensures all datasets use the same genome coordinates, resulting in .bed/.bim/.fam files aligned to GRCh37. 
```
# Check for original build
> plink2 --bfile ADNI_cluster_01_forward_757LONI_10 --chr 1 --recode vcf bgz --out ADNI_cluster_01_forward_757LONI_11_chr1
# Quick assembly check on the VCF
> pip install snps
> python - << 'PY'
from snps import SNPs
s = SNPs("ADNI_cluster_01_forward_757LONI_11_chr1.vcf.gz")
print("assembly:", s.assembly, "build_detected:", s.build_detected)
PY
# Convert plink data to txt .bed, since the tool takes in a different format:
> awk 'BEGIN{OFS="\t"} {print "chr"$1, $4-1, $4, $2}' ADNI_cluster_01_forward_757LONI_10.bim > ADNI_cluster_01_forward_757LONI_10_hg18.bed
```
Now go to http://genome.ucsc.edu/cgi-bin/hgLiftOver to convert from original build to GRCh37/hg19, download the lifted files.
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/liftover.png" width="600">

To convert back to PLINK format and check the build, follow the steps below:
```
# Get rsID of failed liftOvers
> grep -v '^#' hglft_genome_199168_dad920.err.txt | awk '{print $4}' > ADNI_cluster_01_forward_757LONI_10_exclude_snps.txt
# Convert UCSC mapped BED to plink map
> awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1); print $1, $4, 0, $3}' hglft_genome_199168_dad920.bed > ADNI_cluster_01_forward_757LONI_10_lifted.map
> plink2 --bfile ADNI_cluster_01_forward_757LONI_10 --exclude ADNI_cluster_01_forward_757LONI_10_exclude_snps.txt --make-bed --out ADNI_cluster_01_forward_757LONI_11
# --- INPUTS YOU ALREADY HAVE AT THIS POINT ---
# ADNI_cluster_01_forward_757LONI_11.{bed,bim,fam}        # post-exclude (failed liftover SNPs removed)
# hglft_genome_199168_dad920.bed                          # UCSC LiftOver "Converted" (mapped) BED (GRCh37)
# hglft_genome_199168_dad920.err.txt                      # UCSC LiftOver log (kept earlier for exclude list)
# Build PLINK update tables from the UCSC-mapped BED
#    lifted.chr : rsid <tab> newCHR
#    lifted.pos : rsid <tab> newBP  (use BED 'end' column as 1-based position)
> awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print $4, $1}' hglft_genome_199168_dad920.bed > lifted.chr
> awk 'BEGIN{OFS="\t"} {                    print $4, $3}' hglft_genome_199168_dad920.bed > lifted.pos
# Update CHR and BP on the _11 set (requires producing a sorted PGEN first), then convert back to BED/BIM/FAM for normal downstream use.
> plink2 --bfile ADNI_cluster_01_forward_757LONI_11 --update-chr lifted.chr 2 --update-map lifted.pos 2 --sort-vars --make-pgen --out ADNI_cluster_01_forward_757LONI_12_tmp
> plink2 --pfile ADNI_cluster_01_forward_757LONI_12_tmp --make-bed --out ADNI_cluster_01_forward_757LONI_12
# Check if the build is correct now
> plink2 --bfile ADNI_cluster_01_forward_757LONI_12 --recode vcf bgz --out ADNI_cluster_01_forward_757LONI_12_GRCh37
> python - << 'PY'
from snps import SNPs
s = SNPs("ADNI_cluster_01_forward_757LONI_12_GRCh37.vcf.gz")
print("Assembly:", s.assembly)
print("Build detected:", s.build_detected)
PY
```

Step 1.5: How to identify correct reference panel (HRC? 1000G? HapMap?)

We usually do a PCA step here to decide which reference panel our dataset is most aligned with, but since ADNI is >90% non-hispanic/hispanic caucasian, we use the HRC.r1-1.GRCh37 (src:[https://doi.org/10.1111/cns.14073]( https://doi.org/10.1111/cns.14073))

Step 2: Pre-imputation Correction given reference panel (HRC/1000G) 
```
> curl -L -O https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
> unzip -o HRC-1000G-check-bim-v4.3.0.zip
> curl -L -O https://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
> gunzip -f HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

> plink2 --bfile ADNI_cluster_01_forward_757LONI_12 --chr 1-22 --make-bed --out ADNI_cluster_01_forward_757LONI_12_auto
> plink2 --bfile ADNI_cluster_01_forward_757LONI_12_auto --freq --out ADNI_cluster_01_forward_757LONI_12_auto
> awk 'BEGIN{OFS=" "}
     NR==1 {next}
     {
       chr=$1; id=$2; ref=$3; alt=$4; afs=$6; n=$7;
       if (id=="." || id=="" || ref=="" || alt=="" || afs=="." || n=="." || afs=="" || n=="") next;
       if (alt ~ /,/) next;
       af = afs + 0.0; if (af<0) af=0; if (af>1) af=1;
       if (af <= 0.5) { a1=alt; a2=ref; maf=af } else { a1=ref; a2=alt; maf=1-af }
       nchr = 2 * (n + 0);
       printf "%s %s %s %s %.6f %d\n", chr, id, a1, a2, maf, nchr
     }' ADNI_cluster_01_forward_757LONI_12_auto_withIDs.afreq > ADNI_cluster_01_forward_757LONI_12_auto.frq
# Sanity check for correct conversion
# Expect large overlap (tens of thousands), not single digits.
> cut -d' ' -f2 ADNI_cluster_01_forward_757LONI_12_auto.frq       | sort -u > frq_ids.txt
> awk '{print $2}' ADNI_cluster_01_forward_757LONI_12_auto_withIDs.bim | sort -u > bim_ids.txt
> comm -12 frq_ids.txt bim_ids.txt | wc -l
> perl HRC-1000G-check-bim.pl -b ADNI_cluster_01_forward_757LONI_12_auto.bim -f ADNI_cluster_01_forward_757LONI_12_auto.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
> chmod +x Run-plink.sh
> bash Run-plink.sh
```
At this point, we have *chr1.vcf to *chr22.vcf

Step 4: Imputation via Michigan Imputation Server
We will convert the plink data to vcf format, then use the MIS to impute data. 
```
# Download BCFtools (bcftools-1.14) from http://www.htslib.org/download/
> cd bcftools-1.14    # and similarly for bcftools and htslib
> ./configure --prefix=/Users/taehyo/Applications/
> make
> make install
> export PATH=/Users/taehyo/Applications/bin:$PATH```

> for chr in {1..22}
do
	bcftools sort ADNI_cluster_01_forward_757LONI_12_auto-updated-chr$chr.vcf -Oz -o ADNI_cluster_01_forward_757LONI_12_auto-updated-chr$chr.vcf.gz
done

# Now upload to Michigan Imputation Server for imputation: https://imputationserver.sph.umich.edu/index.html#!pages/login

> for chr in $(seq 1 22)
do
	unzip -P 2kFzKqy7lKa1GR chr_$chr.zip
done

> bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou | 
bcftools view -Ou -i 'R2>0.3' |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o ADNI1_allchromosomes.converted.R2_0.3.vcf.gz

> plink2 --vcf ADNI1_allchromosomes.converted.R2_0.3.vcf.gz --double-id --make-bed --out ADNI1_allchromosomes.converted.R2_0.3
```

### Population Structure Modeling:

I repeat the above steps for ADNI1, ADNI2, and ADNIGO. If you only have a single dataset, you can skip the below merge process. 

(Optional Step): If you have multiple datasets, merge all of them now that they are aligned and imputed on the same genome build, followed by final post-merge quality check. 
```
# https://martha-labbook.netlify.app/posts/extracting-data-for-variants-common-in-both-file-sets/
> awk '{print $2}' ADNI1_allchromosomes.converted.R2_0.3.bim | sort > ADNI1_snp_sorted.txt
> awk '{print $2}' ADNIGO2_allchromosomes.converted.R2_0.3.bim | sort > ADNIGO2_snp_sorted.txt
> awk '{print $2}' ADNI3_allchromosomes.converted.R2_0.3.bim | sort > ADNI3_snp_sorted.txt

> comm -12 ADNI1_snp_sorted.txt ADNIGO2_snp_sorted.txt > intersect_ADNI1_GO2_snps.txt
> comm -12 intersect_ADNI1_GO2_snps.txt ADNI3_snp_sorted.txt > intersect_ADNI1_GO2_3_snps.txt
> wc -l intersect_ADNI1_GO2_3_snps.txt

> plink2 --bfile /scratch/hs120/ADNI_SNP/ADNI1_SNP_QC/ADNI1_allchromosomes.converted.R2_0.3 --extract intersect_ADNI1_GO2_3_snps.txt --make-bed --out ADNI1_intersect_snps
> plink2 --bfile /scratch/hs120/ADNI_SNP/ADNIGO2_SNP_QC/ADNIGO2_allchromosomes.converted.R2_0.3 --extract intersect_ADNI1_GO2_3_snps.txt --make-bed --out ADNIGO2_intersect_snps
> plink2 --bfile /scratch/hs120/ADNI_SNP/ADNI3_SNP_QC/ADNI3_allchromosomes.converted.R2_0.3 --extract intersect_ADNI1_GO2_3_snps.txt --make-bed --out ADNI3_intersect_snps

> plink2 --merge-list allfiles2merge.txt --make-bed --out ADNI1_GO2_3_merged_1
> awk '{print $0, $1":"$4}' ADNI1_GO2_3_merged_1.bim > post_imputation_updated_positions
> awk '{print $1":"$4}' ADNI1_GO2_3_merged_1.bim | sort | uniq -d > post_imputation_updated_duplicated_positions 
> grep -w -f  post_imputation_updated_duplicated_positions post_imputation_updated_positions | awk '{print $2}' > post_imputation_updated_duplicated_IDs
> plink2 --bfile ADNI1_GO2_3_merged_1 --exclude post_imputation_updated_duplicated_IDs --make-bed --out ADNI1_GO2_3_merged_2
> plink2 --bfile ADNI1_GO2_3_merged_2 --snps-only --make-bed --out ADNI1_GO2_3_merged_snps

# Post-QC
> plink2 --bfile ADNI1_GO2_3_merged_snps --maf 0.01 --hwe 1e-6 --hwe-all --geno 0.05 --make-bed --out ADNI1_GO2_3_merged_snps.R2_0.3.MAF_0.01.HWE_0.000001.geno_0.05

# Modify below file names and script to align with your data
> Rscript subjects_relatedness.R
> plink2 --bfile ADNI_cluster_01_forward_757LONI_14 --update-name ADNI1_snp_rename.txt --make-bed --out ADNI_cluster_01_forward_757LONI_14_renameID
> plink2 --bfile ADNI_cluster_01_forward_757LONI_14_renameID --extract ADNI_common_snp.txt --make-bed --out ADNI1_commonSNP_data
> plink2 --bfile ADNI_GO2_10 --update-name ADNIGO2_snp_rename.txt --make-bed --out ADNI_GO2_10_renameID
> plink2 --bfile ADNI_GO2_10_renameID --extract /Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI_common_snp.txt --make-bed --out ADNIGO2_commonSNP_data
> plink2 --bfile ADNI3_SNP_10 --update-name ADNI3_snp_rename.txt --make-bed --out ADNI3_SNP_10_renameID
> plink2 --bfile ADNI3_SNP_10_renameID --extract /Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI_common_snp.txt --make-bed --out ADNI3_commonSNP_data
> plink2 --merge-list allfiles2merge_CheckSubjectRelatedness.txt --make-bed --out ADNI1_GO2_3_merged_CheckSubjectRelatedness
> plink2 --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.2 --out indepSNP_ADNIall_CheckSubjectRelatedness
> plink2 --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness --filter-founders --make-bed --out ADNI1_GO2_3_merged_CheckSubjectRelatedness_2
> plink2 --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness_2 --extract indepSNP_ADNIall_CheckSubjectRelatedness.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders_ADNIall_CheckSubjectRelatedness
> plink2 --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness_2 --missing
# Generate a list of FID and IID of the individual(s) with a Pihat above 0.2, to check who had the lower call rate of the pair.
#        81   024_S_4084          Y       48    84007 0.0005714
#     ADNI3   024_S_6005          Y       51    84007 0.0006071


#       118   023_S_4035          Y       36    84007 0.0004285
#      591   023_S_0058          Y       15    84007 0.0001786


#       179   012_S_4094          Y       37    84007 0.0004404
#  ADNI3   002_S_6066          Y       67    84007 0.0007976


#      342   137_S_4466          Y       62    84007 0.000738
#     453   021_S_0159          Y       44    84007 0.0005238


#      396   011_S_4235          Y       65    84007 0.0007737
#    ADNI3   011_S_6303          Y       47    84007 0.0005595

> nano 0.2_low_call_rate_pihat_ADNIall_CheckSubjectRelatedness.txt
ADNI3   024_S_6005
118   023_S_4035
ADNI3   002_S_6066 
342   137_S_4466
396   011_S_4235 

> plink2 --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness_2 --remove 0.2_low_call_rate_pihat_ADNIall_CheckSubjectRelatedness.txt --make-bed --out ADNI1_GO2_3_merged_CheckSubjectRelatedness_3
> awk '{print $1"_"$2, $1"_"$2}' ADNI1_GO2_3_merged_CheckSubjectRelatedness_3.fam > ADNIall_CheckSubjectRelatedness_subj.txt 
> plink2 --bfile ADNI1_GO2_3_merged_snps.R2_0.3.MAF_0.01.HWE_0.000001.geno_0.05 --keep ADNIall_CheckSubjectRelatedness_subj.txt --make-bed --out ADNI1_GO2_3_merged_snps.R2_0.3.MAF_0.01.HWE_0.000001.geno_0.05.IBD_0.2
```

Step 1: Population stratification is corrected by extracting principal components (PCs) for each dataset separately using LD-pruned SNPs. The top PCs are used as covariates in GWAS to control for ancestry differences. 
```
> plink2 --bfile ADNI_cluster_01_forward_757LONI_14 --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.2 --out pruned_data
> plink2 --bfile ADNI_cluster_01_forward_757LONI_14 --extract pruned_data.prune.in --make-bed --out ADNI_cluster_01_forward_757LONI_14_pruned
> plink2 --bfile ADNI_cluster_01_forward_757LONI_14_pruned --pca 10 --out ADNI1_PCA

> plink2 --bfile ADNI_GO2_10 --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.2 --out pruned_data
> plink2 --bfile ADNI_GO2_10 --extract pruned_data.prune.in --make-bed --out ADNI_GO2_10_pruned
> plink2 --bfile ADNI_GO2_10_pruned --pca 10 --out ADNIGO2_PCA

> plink2 --bfile ADNI3_SNP_10 --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.2 --out pruned_data
> plink2 --bfile ADNI3_SNP_10 --extract pruned_data.prune.in --make-bed --out ADNI3_SNP_10_pruned
> plink2 --bfile ADNI3_SNP_10_pruned --pca 10 --out ADNI3_PCA
```

### Associative Analysis:

Run GWAS on progression (e.g. ADAS-Cog, MMSE)
- dummy variables of adni cohort (adni1, adni2, adnigo)
- gender, age
- race, ethnicity, education, marital status, 10 PC
- try to maintain 300-400 SNPs? make sure p >> n 
- response variable should be logistic (CN or AD), each regression containing 1 SNP, use pvalue threshold to retain reasonable # of SNP
- this is more straightforward to see rather than using other response variables
- then, before SCCA, do linear regression (regress out the same covariates, and then i
- check if our snp data contains SNP for APOE4 (double check)
- for FDG PET, take the average across each ROI, then regress 
```
# Modify below to fit your research question. 
> Rscript subjects_bl_covar.R
> plink2 --bfile ADNI1_GO2_3_merged_snps.R2_0.3.MAF_0.01.HWE_0.000001.geno_0.05.IBD_0.2 --keep ADNIall_ADandCN_subj.txt --make-bed --out ADNIall_ADandCN_1
> plink2 --bfile ADNIall_ADandCN_2 --logistic hide-covar --covar ADNIall_ADandCN_covar.txt --allow-no-sex --out ADNIall_ADandCN_test_result
> plink2 -bfile ADNIall_ADandCN_2 --extract sigSNPs.txt --make-bed --out ADNIall_ADandCN_3
> plink2 --bfile ADNIall_ADandCN_3 --recodeA --out ADNIall_ADandCN_sigSNPs
```

## Application: Sparse Canonical Correlation Analysis using Imaging-Omics (FDG-PET, SNP)

Goal: To investigate how genetic variation (SNPs) and brain imaging features (FDG-PET ROIs) are related, and to assess how these associations change across disease stages (CN, MCI, AD). We use Sparse Canonical Correlation Analysis (SCCA) to identify low-dimensional, interpretable patterns linking high-dimensional SNP data with FDG-PET imaging phenotypes.

Input: 
- (n,p) SNP data matrix (split into CN, MCI, AD)
- (n,q) FDG-PET data matrix (split into CN, MCI, AD)
- Note that we did not use disease status (CN, MCI, AD) as outcomes for GWAS, because our goal is not to find SNPs associated with diagnosis, but rather SNPs linked to disease progression. Instead, we selected SNPs based on their association with continuous progression measures (e.g., cognitive decline, imaging biomarkers), so that the same SNP panel can be meaningfully compared across CN, MCI, and AD groups during SCCA.

Please read [this]() for the complete application details. 
