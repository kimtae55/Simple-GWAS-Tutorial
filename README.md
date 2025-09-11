# Simple GWAS Tutorial for Plink2
This is a Step by Step Tutorial for GWAS (for Plink format), explaining what SNP data looks like, how to perform quality control and imputation procedures, and how to run GWAS.
Credits for the figures and explanations here go to the more in-depth tutorials: [GWASTutorial1](https://cloufield.github.io/GWASTutorial), [GWASTutorial2](https://www.ncbi.nlm.nih.gov/pubmed/29484742), or [GWASTutorial3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6001694/). Most of the codes are taken from [https://github.com/MareesAT/GWA_tutorial](https://github.com/MareesAT/GWA_tutorial), but modified to be compatible with plink2 (for macbook m1 and above users). 

You can replace my data with your own plink data (.bed, .bim, .fam) and replicate the whole experiment. 
- If you plan to use GRCh37/hg19 genome build, follow instructions [here](https://github.com/kimtae55/Simple-GWAS-Tutorial/tree/main/GRCh37_hg19).
- If you plan to use GRCh38/hg38 genome build, stay on this page.
- For research, GRCh38/hg38 is recommended as it is the most up-to-date human genome build (released 2013). 

At the end of the tutorial, I provide an application example of using the extracting SNPs and additional FDG-PET data to conduct an imaging-genetics sparse canonical correlation analysis. 

## Prerequisites
- Install ```plink2, R, python```
- Install BCFtools (bcftools-1.14) from [http://www.htslib.org/download/](http://www.htslib.org/download/), then: 
```
> cd bcftools-1.14    # and similarly for bcftools and htslib
> ./configure --prefix=/Users/taehyo/Applications/
> make
> make install
```
- Note: Add above executables to $PATH for convenience. 

## Assumption: You have SNP data in PLINK format
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/plink.png" width="600">

In my case, I use Plink data (.bed, .bim, .fam files) downloaded from the ADNI1, ADNI2, and ADNIGO studies. For simplicity, all data, scripts, and outputs will be located in a single working directory. 

In terminal, define your base file name (string that precedes your .bed, .bim, .fam filenames)
```
base="ADNI_cluster_01_forward_757LONI"
```

## What are the necessary steps in a GWAS?

### Quality Control --> Liftover and Imputation --> Population Structure Modeling --> Associative Analysis 

Quality Control is done at a sample-level (to remove bad individuals; e.g. contamination, swaps, relatedness, sex mismatches) and SNP-level (to remove bad variants; e.g. missingness, low MAF, HWE failures).
If you have multiple datasets from different studies/timepoints, you should run this part separately, then merge the results (see the later optional step for merging)

### QC Steps:

Step 1: Handle missingness per individual and per SNP: Delete individuals with missingness >0.05.
```
# Step 1. Visualize missingness
plink2 --bfile "$base" --missing
Rscript --no-save hist_miss.R

# Step 2. Filter SNPs with missingness > 5%
plink2 --bfile "$base" \
       --geno 0.05 \
       --make-bed \
       --out "${base}_2"

# Step 3. Filter individuals with missingness > 5%
plink2 --bfile "${base}_2" \
       --mind 0.05 \
       --make-bed \
       --out "${base}_3"
```

Step 2: Handle sex discrepancy: Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.
```
# Step 1. Check sex discrepancies
plink2 --bfile "${base}_3" \
       --check-sex max-female-xf=0.2 min-male-xf=0.8

# Step 2. Visualize sex discrepancy check
Rscript --no-save gender_check.R

# Step 3. Extract problematic samples
grep "PROBLEM" plink2.sexcheck | awk '{print $1, $2}' > sex_discrepancy.txt

# Step 4. Remove problematic samples
plink2 --bfile "${base}_3" \
       --remove sex_discrepancy.txt \
       --make-bed \
       --out "${base}_4"
```

Step 3: Extract autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF <0.01).
```
# Step 1. Extract autosomal SNPs (chr 1–22)
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' "${base}_4.bim" > snp_1_22.txt

# Step 2. Keep only autosomal SNPs
plink2 --bfile "${base}_4" \
       --extract snp_1_22.txt \
       --make-bed \
       --out "${base}_5"

# Step 3. Calculate allele frequencies
plink2 --bfile "${base}_5" \
       --freq \
       --out MAF_check

# Step 4. Visualize MAF distribution
Rscript --no-save MAF_check.R

# Step 5. Filter SNPs with MAF < 1%
plink2 --bfile "${base}_5" \
       --maf 0.01 \
       --make-bed \
       --out "${base}_6"
```

Step 4: Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
```
# Step 1. Compute Hardy-Weinberg Equilibrium (HWE) stats
plink2 --bfile "${base}_6" \
       --hardy

# Step 2. Visualize HWE results
Rscript --no-save hwe.R

# Step 3. Filter SNPs failing HWE (p < 1e-6)
plink2 --bfile "${base}_6" \
       --hwe 1e-6 \
       --make-bed \
       --out "${base}_7"
```

Step 5: Remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.
```
# Step 1. LD pruning — exclude long-range high-LD regions, keep independent SNPs
plink2 --bfile "${base}_7" \
       --exclude range high-LD-regions-hg18-NCBI36.txt \
       --indep-pairwise 50 5 0.2 \
       --out indepSNP

# Step 2. Check heterozygosity using pruned SNPs
plink2 --bfile "${base}_7" \
       --extract indepSNP.prune.in \
       --het \
       --out R_check

# Step 3. Visualize heterozygosity rates
Rscript --no-save check_heterozygosity_rate.R

# Step 4. Remove individuals failing heterozygosity QC
plink2 --bfile "${base}_7" \
       --remove het_fail_ind.txt \
       --make-bed \
       --out "${base}_8"
```

Step 6: We exclude all individuals with a PI_HAT > 0.2 to remove cryptic relatedness, assuming a random population sample.
```
# Step 1. LD-prune SNPs (used only for relatedness detection)
plink2 --bfile "${base}_8" \
       --indep-pairwise 200 100 0.1 \
       --out indepSNP

# Step 2. Use pruned SNPs to identify related individuals (temporary dataset)
plink2 --bfile "${base}_8" \
       --extract indepSNP.prune.in \
       --king-cutoff 0.10 \
       --make-bed \
       --out ADNI_relcheck_tmp

# Step 3. Save list of unrelated individuals
awk '{print $1, $2}' ADNI_relcheck_tmp.fam > unrelated.keep

# Step 4. Apply unrelated sample list to full dataset (keeps all SNPs, drops related individuals)
plink2 --bfile "${base}_8" \
       --keep unrelated.keep \
       --make-bed \
       --out "${base}_10"
```

### Liftover and Imputation:

Step 1: ADNI datasets are often on older genome builds (e.g., hg18/NCBI36). Before imputation, convert to hg19/GRCh37. This ensures all datasets use the same genome coordinates, resulting in .bed/.bim/.fam files aligned to GRCh37. You could also lift to hg38, which is the most current reference genome build, but I choose hg19 because it is easier to validate against existing works (e.g [doi:10.1016/j.neurobiolaging.2019.06.003](https://pmc.ncbi.nlm.nih.gov/articles/PMC6732252/), [doi: 10.1186/s12880-025-01782-2](https://pubmed.ncbi.nlm.nih.gov/40676546/))
```
# Step 1. Export chr1 VCF to check genome build
plink2 --bfile "${base}_10" \
       --chr 1 \
       --recode vcf bgz \
       --out "${base}_11_chr1"

# Step 2. Quick assembly check on VCF
pip install snps
python - << 'PY'
from snps import SNPs
s = SNPs("${base}_11_chr1.vcf.gz")
print("assembly:", s.assembly, "build_detected:", s.build_detected)
PY

# Step 3. Convert .bim to UCSC BED format for liftover (assumes original build = hg18)
awk 'BEGIN{OFS="\t"} {print "chr"$1, $4-1, $4, $2}' \
    "${base}_10.bim" > "${base}_10_hg18.bed"
```
Now go to http://genome.ucsc.edu/cgi-bin/hgLiftOver to convert from original build to GRCh38/hg38, download the lifted files.
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/liftover.png" width="600">

To convert back to PLINK format and check the build, follow the steps below:
```
# Files produced by UCSC LiftOver (you already have these)
#   - Converted BED (target build, e.g., GRCh38):   hglft_genome_199168_dad920.bed
#   - Error log with unmapped intervals:            hglft_genome_199168_dad920.err.txt

# 1) Get rsIDs of failed liftOvers (exclude list)
grep -v '^#' hglft_genome_199168_dad920.err.txt | awk '{print $4}' > "${base}_10_exclude_snps.txt"

# 2) (Optional) Create a PLINK .map from the converted BED (for reference/inspection)
#    Columns: chr (no 'chr' prefix), rsID, genetic_distance(=0), new_bp (BED end; 1-based)
awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1); print $1, $4, 0, $3}' \
    hglft_genome_199168_dad920.bed > "${base}_10_lifted.map"

# 3) Drop variants that failed to lift
plink2 --bfile "${base}_10" \
       --exclude "${base}_10_exclude_snps.txt" \
       --make-bed \
       --out "${base}_11"

# 4) Build PLINK update tables from the converted BED
#    lifted.chr : rsid <tab> newCHR(without 'chr')
#    lifted.pos : rsid <tab> newBP  (use BED 'end' col as 1-based position)
awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print $4, $1}' hglft_genome_199168_dad920.bed > lifted.chr
awk 'BEGIN{OFS="\t"} {                    print $4, $3}' hglft_genome_199168_dad920.bed > lifted.pos

# 5) Update CHR/BP on the post-exclude set; sort variants; write temporary PGEN; convert back to BED
plink2 --bfile "${base}_11" \
       --update-chr lifted.chr 2 \
       --update-map  lifted.pos 2 \
       --sort-vars \
       --make-pgen \
       --out "${base}_12_tmp"

plink2 --pfile "${base}_12_tmp" \
       --make-bed \
       --out "${base}_12"

# 6) Quick build sanity check by exporting a VCF and letting 'snps' detect assembly
plink2 --bfile "${base}_12" \
       --recode vcf bgz \
       --out "${base}_12_GRCh38"   # rename GRCh37 if that’s your target

python - <<PY
from snps import SNPs
s = SNPs("${base}_12_GRCh38.vcf.gz")
print("Assembly:", s.assembly, "Build detected:", s.build_detected)
PY
```

Step 2: Pre-imputation setup
```
# Step 1. Keep only autosomal chromosomes (1–22)
plink2 --bfile "${base}_12" \
       --chr 1-22 \
       --make-bed \
       --out "${base}_12_auto"

# Step 2. Export per-chromosome, bgzipped VCF files (ready for imputation)
for chr in {1..22}
do
  plink2 \
    --bfile "${base}_12_auto" \
    --chr "$chr" \
    --recode vcf bgz \
    --out "${base}_12_auto-updated-chr${chr}"
done

# Check SNPs exist
for chr in {1..22}; do
  echo "chr${chr}: $(bcftools view -H "${base}_12_auto-updated-chr${chr}.vcf.gz" | wc -l) variants"
done
```
At this point, we have *chr1.vcf to *chr22.vcf

Step 4: Imputation via Michigan Imputation Server
```
# Now upload to Michigan Imputation Server for imputation: https://imputationserver.sph.umich.edu/index.html#!pages/login

> for chr in $(seq 1 22)   
do
    unzip -P 'PASSWORD FROM EMAIL' chr_${chr}.zip
done

> bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou | 
bcftools view -Ou -i 'R2>0.3' |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o ADNI1_allchromosomes.converted.R2_0.3.vcf.gz

> plink2 --vcf ADNI1_allchromosomes.converted.R2_0.3.vcf.gz --double-id --make-bed --out ADNI1_allchromosomes.converted.R2_0.3
```

### Population Structure Modeling:

At this moment, I have ADNI1/GO/2 imputed with the HRC r1.1 (hg19) panel, which is European-centric. I assess ancestry and adjust for population structure by anchoring PCA to an external reference panel (HapMap3 or 1000G) which has more coverage beyond European ancestry. I compute allele-weighted PCs on the combined ADNI + external reference dataset so the axes reflect global ancestry, then project each ADNI cohort onto this shared PC space. Correct clustering is confirmed when ADNI samples group tightly with European reference populations (e.g., CEU/TSI), while outliers are identified if they drift toward non-European clusters. I will utilize the samples that align well with the HRC r1.1 (hg19) panel. The top 10 PCs are also saved to be included later as covariates for GWAS. 
```
```

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
