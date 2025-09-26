# Simple GWAS Tutorial for Plink2
This is a Step by Step Tutorial for GWAS (for Plink format), explaining what SNP data looks like, how to perform quality control and imputation procedures, and how to run GWAS.
Credits for the figures and explanations here go to the more in-depth tutorials: [GWASTutorial1](https://cloufield.github.io/GWASTutorial), [GWASTutorial2](https://www.ncbi.nlm.nih.gov/pubmed/29484742), or [GWASTutorial3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6001694/). The main pipeline is inspired by code from [https://github.com/MareesAT/GWA_tutorial](https://github.com/MareesAT/GWA_tutorial). 
 
Most GWAS tutorials out there only support plink and older genome build (e.g. hg18, hg19), so this tutorial provides plink2/hg38/TopMed reference panel compatible instructions. 

You can replace my data with your own plink data (.bed, .bim, .fam) and replicate the whole experiment. Note that GRCh38/hg38 is highly recommended as it is the most up-to-date human genome build (released 2013). 

At the end of the tutorial, I provide an application example of using the extracting SNPs and additional FDG-PET data to conduct an imaging-genetics sparse canonical correlation analysis. 

## Table of Contents
1. [Setup](#prerequisites)
2. [GWAS Step 1: Pre-imputation QC](#gwas-step-1-pre-imputation-quality-control-qc)
3. [GWAS Step 2: Liftover and Imputation with TopMed Imputation Server](#liftover-and-imputation)
4. [GWAS Step 3: Population Structure Modeling, Merging, and Post Imputation QC](#population-structure-modeling)
5. [GWAS Step 4: Associative Analysis](#associative-analysis)
6. [Application: Sparse Canonical Correlation Analysis using Imaging-Omics (FDG-PET, SNP)](#application-sparse-canonical-correlation-analysis-using-imaging-omics-fdg-pet-snp)


# Prerequisites

### 1. Install Required Software
- **plink & plink2** (merge functions are not fully implemented for plink2 as of 09.25.2025, but should be updated soon)
- **R**
- **Python**
- **BCFtools** (≥ v1.14) from [htslib.org](http://www.htslib.org/download/)

#### Example: Installing BCFtools
```bash
cd bcftools-1.14        # navigate into the bcftools source folder (repeat similarly for htslib if needed)
./configure --prefix=/Users/taehyo/Applications/
make
make install
```

👉 Add the installation directory (e.g., `/Users/taehyo/Applications/bin`) to your `$PATH` so the executables can be called directly in terminal.

### 2. Assumption: SNP Data in PLINK Format
Your input data should be in **PLINK binary format**:

- `.bed` → binary genotype data  
- `.bim` → SNP information  
- `.fam` → family/individual information  

<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/plink.png" width="600">

In this tutorial, I use PLINK data (`.bed/.bim/.fam`) downloaded from the **ADNI1, ADNI2, and ADNIGO** studies.  
For simplicity, all data, scripts, and outputs are stored in a single working directory.

#### Define your base file name
For all commands from now, `base` string is the path to your original `.bed/.bim/.fam` file.  
```bash
base="ADNI_cluster_01_forward_757LONI"
```

## GWAS Step 1: Pre-imputation Quality Control (QC)

Quality Control is done at a sample-level (to remove bad individuals; e.g. contamination, swaps, relatedness, sex mismatches) and SNP-level (to remove bad variants; e.g. missingness, low MAF, HWE failures).
If you have multiple datasets from different studies/timepoints, you should run this part separately, then merge the results (see the later optional step for merging)

### 1) Handle missingness per individual and per SNP: Delete individuals with missingness >0.05.
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

### 2) Handle sex discrepancy: Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.
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

Step 1: ADNI datasets are often on older genome builds (e.g., hg18/NCBI36). Before imputation, convert to GRCh38/hg38. This ensures all datasets use the same genome coordinates, resulting in .bed/.bim/.fam files aligned to GRCh38/hg38. If your data is already aligned to GRCh38/hg38, skip this step.
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
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/liftover.png" width="800">

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

# Step 2. Export per-chromosome, bgzipped VCF files with chr-prefix (hg38 style)
for chr in {1..22}; do
  plink2 \
    --bfile "${base}_12_auto" \
    --chr "$chr" \
    --output-chr chrM \
    --recode vcf bgz \
    --out "${base}_12_auto-chr${chr}"
  tabix -f "${base}_12_auto-chr${chr}.vcf.gz"
done

# Quick check: CHROM should be chr1..chr22
for chr in {1..22}; do
  echo -n "chr${chr}: "
  bcftools view -H "${base}_12_auto-chr${chr}.vcf.gz" | awk '{print $1}' | head -1
done

# modify make_topmed_ready.sh to set the correct $base, $REF_FASTA
# You can download the reference fasta file (hg38.analysisSet.fa.gz) from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/
# This will fix strand flip issues and make the vcf files ready for imputation. 
./make_topmed_ready.sh
```
At this point, you have *-chr1.topmed.vcf.gz, ..., *-chr22.topmed.vcf.gz ready for imputation. 

Step 4: Upload the above vcfs to TopMed Imputation Server ([https://imputation.biodatacatalyst.nhlbi.nih.gov/#!run/imputationserver%402.0.0-beta3](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!run/imputationserver%402.0.0-beta3)) for imputation. Once completed:
```
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

At this point, I have ADNI1/GO/2 imputed with the TopMed reference panel. I assess ancestry and adjust for population structure by anchoring PCA to an external reference panel (HapMap3 or 1000G). I compute allele-weighted PCs on the combined ADNI + external reference dataset so the axes reflect global ancestry, then project each ADNI cohort onto this shared PC space. The top 10 PCs are saved to be included later as covariates for GWAS. 
```
# Download 1000G reference panel 
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

# Per-chromosome VCFs + index
for chr in $(seq 1 22); do
  f="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
  curl -fL -C - -O "${BASE_URL}/${f}"
  curl -fL -C - -O "${BASE_URL}/${f}.tbi"
done

# Population panel (sample -> pop -> super-pop)
curl -fL -o integrated_call_samples_v3.20130502.ALL.panel \
  "https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel"


```

(Optional Step): If you have multiple datasets, merge all of them now that they are aligned and imputed on the same genome build, followed by final post-merge quality check. 
```
THREADS=1
MEM=8000

echo "[1/6] Build intersect of SNP IDs across all three (BIM col2)..."
awk '{print $2}' ADNI1_allchromosomes.converted.R2_0.3.bim     | sort -u > ids_1.txt
awk '{print $2}' ADNI2_SET1_allchromosomes.converted.R2_0.3.bim | sort -u > ids_2.txt
awk '{print $2}' ADNI2_SET2_allchromosomes.converted.R2_0.3.bim | sort -u > ids_3.txt

comm -12 ids_1.txt ids_2.txt > ids_12.txt
comm -12 ids_12.txt ids_3.txt > ids_intersect.txt
wc -l ids_intersect.txt
[ -s ids_intersect.txt ] || { echo "ERROR: no common SNP IDs."; exit 1; }

echo "[2/6] Subset each dataset to the intersect (plink2)..."
plink2 --threads $THREADS --memory $MEM --bfile ADNI1_allchromosomes.converted.R2_0.3     --extract ids_intersect.txt --make-bed --out ADNI1_intersect
plink2 --threads $THREADS --memory $MEM --bfile ADNI2_SET1_allchromosomes.converted.R2_0.3 --extract ids_intersect.txt --make-bed --out ADNI2_SET1_intersect
plink2 --threads $THREADS --memory $MEM --bfile ADNI2_SET2_allchromosomes.converted.R2_0.3 --extract ids_intersect.txt --make-bed --out ADNI2_SET2_intersect

echo "[3/6] Merge with plink 1.9..."
cat > merge_list.txt <<EOF
ADNI2_SET1_intersect.bed ADNI2_SET1_intersect.bim ADNI2_SET1_intersect.fam
ADNI2_SET2_intersect.bed ADNI2_SET2_intersect.bim ADNI2_SET2_intersect.fam
EOF

if ! plink --bfile ADNI1_intersect --merge-list merge_list.txt --make-bed --out ADNI_merged_1 ; then
  echo "Merge flagged problematic sites; excluding missnp and retrying..."
  plink --bfile ADNI1_intersect --exclude ADNI_merged_1-merge.missnp --make-bed --out ADNI1_intersect_fix
  plink --bfile ADNI2_SET1_intersect --exclude ADNI_merged_1-merge.missnp --make-bed --out ADNI2_SET1_intersect_fix
  plink --bfile ADNI2_SET2_intersect --exclude ADNI_merged_1-merge.missnp --make-bed --out ADNI2_SET2_intersect_fix
  cat > merge_list_fix.txt <<EOF
ADNI2_SET1_intersect_fix.bed ADNI2_SET1_intersect_fix.bim ADNI2_SET1_intersect_fix.fam
ADNI2_SET2_intersect_fix.bed ADNI2_SET2_intersect_fix.bim ADNI2_SET2_intersect_fix.fam
EOF
  plink --bfile ADNI1_intersect_fix --merge-list merge_list_fix.txt --make-bed --out ADNI_merged_1
fi

echo "[4/6] FAST de-dup by position (CHR:BP)..."
awk 'NR==FNR{ k=$1":"$4; c[k]++; next } { k=$1":"$4; if (c[k]>1) print $2 }' ADNI_merged_1.bim ADNI_merged_1.bim > ADNI_merged_1_dupIDs.txt
wc -l ADNI_merged_1_dupIDs.txt

if [ -s ADNI_merged_1_dupIDs.txt ]; then
  plink2 --threads $THREADS --memory $MEM --bfile ADNI_merged_1 --exclude ADNI_merged_1_dupIDs.txt --make-bed --out ADNI_merged_2
else
  cp ADNI_merged_1.bed ADNI_merged_2.bed
  cp ADNI_merged_1.bim ADNI_merged_2.bim
  cp ADNI_merged_1.fam ADNI_merged_2.fam
fi

echo "[5/6] Keep SNPs only..."
plink2 --threads $THREADS --memory $MEM --bfile ADNI_merged_2 --snps-only --make-bed --out ADNI_merged_snps

```

Step 1: Post Imputation/Merge QC
```
# ---------- restore sex from original fams (FID_IID -> SEX) ----------
awk '{if($5==1||$5==2) print $1"_"$2, $5}' ADNI_1/ADNI_cluster_01_forward_757LONI_11.fam > fid_iid_sex.map
awk '{if($5==1||$5==2) print $1"_"$2, $5}' ADNI_GO_2_SET_1/ADNI_GO_2_Forward_Bin_11.fam >> fid_iid_sex.map
awk '{if($5==1||$5==2) print $1"_"$2, $5}' ADNI_GO_2_SET_2/ADNI_GO2_GWAS_2nd_orig_BIN_11.fam >> fid_iid_sex.map
awk 'NR==FNR{m[$1]=$2; next} {iid=$2; if(iid in m) print $1, $2, m[iid]}' fid_iid_sex.map ADNI_merged_snps.fam > sex_update_resolved.txt

plink2 --threads 1 --memory 8000 --bfile ADNI_merged_snps --update-sex sex_update_resolved.txt --make-bed --out ADNI_qc1
awk '{print $5}' ADNI_qc1.fam | sort | uniq -c

# ---------- restore sex from original fams (FID_IID -> SEX) ----------
awk '{if($5==1||$5==2) print $1"_"$2, $5}' ADNI_1/ADNI_cluster_01_forward_757LONI_11.fam > fid_iid_sex.map
awk '{if($5==1||$5==2) print $1"_"$2, $5}' ADNI_GO_2_SET_1/ADNI_GO_2_Forward_Bin_11.fam >> fid_iid_sex.map
awk '{if($5==1||$5==2) print $1"_"$2, $5}' ADNI_GO_2_SET_2/ADNI_GO2_GWAS_2nd_orig_BIN_11.fam >> fid_iid_sex.map
awk 'NR==FNR{m[$1]=$2; next} {iid=$2; if(iid in m) print $1, $2, m[iid]}' fid_iid_sex.map ADNI_merged_snps.fam > sex_update_resolved.txt

plink2 --threads 1 --memory 8000 --bfile ADNI_merged_snps --update-sex sex_update_resolved.txt --make-bed --out ADNI_qc1
awk '{print $5}' ADNI_qc1.fam | sort | uniq -c

# ---------- missingness ----------
plink2 --threads 1 --memory 8000 --bfile ADNI_qc1 --missing --out ADNI_qc1_miss
plink2 --threads 1 --memory 8000 --bfile ADNI_qc1 --geno 0.05 --make-bed --out ADNI_qc2
plink2 --threads 1 --memory 8000 --bfile ADNI_qc2 --mind 0.05 --make-bed --out ADNI_qc3

# ---------- autosomes only + MAF ----------
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ADNI_qc3.bim > snp_autosomes.txt
plink2 --threads 1 --memory 8000 --bfile ADNI_qc3 --extract snp_autosomes.txt --make-bed --out ADNI_qc4
plink2 --threads 1 --memory 8000 --bfile ADNI_qc4 --freq --out ADNI_maf
plink2 --threads 1 --memory 8000 --bfile ADNI_qc4 --maf 0.01 --make-bed --out ADNI_qc5

# ---------- HWE ----------
plink2 --threads 1 --memory 8000 --bfile ADNI_qc5 --hardy --out ADNI_hwe
plink2 --threads 1 --memory 8000 --bfile ADNI_qc5 --hwe 1e-6 --make-bed --out ADNI_qc6

# ---------- heterozygosity outliers (hg38 long-range LD mask) ----------
plink2 --threads 1 --memory 8000 --bfile ADNI_qc6 --exclude range high-LD-regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.2 --out indepSNP
plink2 --threads 1 --memory 8000 --bfile ADNI_qc6 --extract indepSNP.prune.in --het --out ADNI_het
# Create het_fail_ind.txt from ADNI_het.het (your usual R/awk); if it exists, remove:
if [ -s het_fail_ind.txt ]; then
  plink2 --threads 1 --memory 8000 --bfile ADNI_qc6 --remove het_fail_ind.txt --make-bed --out ADNI_qc7
else
  cp ADNI_qc6.bed ADNI_qc7.bed; cp ADNI_qc6.bim ADNI_qc7.bim; cp ADNI_qc6.fam ADNI_qc7.fam
fi

# ---------- relatedness (automatic prune; no manual PI_HAT picks) ----------
plink2 --threads 1 --memory 8000 --bfile ADNI_qc7 --indep-pairwise 200 100 0.1 --out indepSNP_rel
plink2 --threads 1 --memory 8000 --bfile ADNI_qc7 --extract indepSNP_rel.prune.in --king-cutoff 0.10 --make-bed --out ADNI_rel_unrelated_tmp
awk '{print $1, $2}' ADNI_rel_unrelated_tmp.fam > unrelated.keep
plink2 --threads 1 --memory 8000 --bfile ADNI_qc7 --keep unrelated.keep --make-bed --out ADNI_qc_final
```
Now, we can finally proceed to population stratification and gwas using the ADNI_qc_final dataset

Step 2: Population stratification is corrected by extracting principal components (PCs) for each dataset separately using LD-pruned SNPs. The top PCs are used as covariates in GWAS to control for ancestry differences. 
```
# LD prune (exclude long-range LD regions on hg38 to avoid correlated SNPs)
plink2 --threads 1 --memory 8000 --bfile ADNI_qc_final --exclude range high-LD-regions-hg38-GRCh38.txt --indep-pairwise 200 100 0.1 --out indepSNP_pca

# PCA using pruned SNPs
plink2 --threads 1 --memory 8000 --bfile ADNI_qc_final --extract indepSNP_pca.prune.in --pca 20 approx --out ADNI_pca
```

Because the APOE4 allele is determined by a specific combination of two single nucleotide polymorphisms (SNPs): rs429358 and rs7412. For APOE4, the combination is defined by a C allele at position 19:44908684 (rs429358) and a C allele at position 19:44908822 (rs7412), I make sure my dataset contains these SNPs - I need to check again after running GWAS.
```
awk -F'\t' '$1==19 && ($4==44908684 || $4==44908822)' ADNI_qc_final.bim
# 19	chr19:44908684:T:C	0	44908684	C	T
# 19	chr19:44908822:C:T	0	44908822	T	C
```

### Associative Analysis:

I utilize the approporiate demographic data from ADNI studies to curate the covariates that I want to control for confouding, namely: study phase, gender, age, race, ethnicity, education, marital status, top 10 PC
```
Rscript gwas_processing.R

# Variable Description
**ID variables**
- `FID`, `IID`: Family and Individual IDs from PLINK `.fam`.  
  Used to align covariates with genotype data.

**Phenotype**
- Stored separately as `DIAG01` (CN = 1, AD = 2, others = -9).

**Demographics / Covariates**
- `AGE`: Age at baseline visit (years, numeric).
- `PTEDUCAT`: Years of education (numeric; e.g., 12 = high school, 16 = college).

**Sex (from `.fam`)**
- `PTGENDER_Female` = 1 if **female**, 0 if **male**.

**Race (from `PTRACCAT`)**  
Reference category = **White** (all race dummies = 0).
- `PTRACCAT_AmIndian`   = 1 if American Indian/Alaska Native, else 0  
- `PTRACCAT_Asian`      = 1 if Asian, else 0  
- `PTRACCAT_Black`      = 1 if Black, else 0  
- `PTRACCAT_HawaiianPI` = 1 if Native Hawaiian/Other Pacific Islander, else 0  
- `PTRACCAT_MoreOne`    = 1 if More than one race, else 0

**Ethnicity (from `PTETHCAT`)**  
Reference = **Not Hisp/Latino** (coded 0).
- `PTETHCAT_HispLatino` = 1 if Hisp/Latino, else 0

**Marital status (from `PTMARRY`)**  
Reference = **Married** (all marital dummies = 0).
- `PTMARRY_Divorced`     = 1 if Divorced (code 3), else 0  
- `PTMARRY_Widowed`      = 1 if Widowed (code 2), else 0  
- `PTMARRY_NeverMarried` = 1 if Never married (code 4), else 0

**Ancestry PCs**
- `PC1`–`PC10`: Top 10 principal components from genotype PCA (numeric).
```

Now, run GWAS. If you are interested in the model (glm) and the exact equation formulations, read [this](https://cloufield.github.io/GWASTutorial/06_Association_tests/#association-testing-basics)
```
plink2 \
  --bfile ADNI_qc_final \
  --keep ADNI_keep_ADvsCN.txt \
  --pheno ADNI_pheno_ADvsCN.txt \
  --pheno-name DIAG01 \
  --covar ADNI_covar_ADvsCN.txt \
  --covar-name \
    AGE,PTEDUCAT,PTGENDER_Female, \
    PTRACCAT_AmIndian,PTRACCAT_Asian,PTRACCAT_Black, \
    PTRACCAT_MoreOne,PTRACCAT_HawaiianPI, \
    PTETHCAT_HispLatino, \
    PTMARRY_Divorced,PTMARRY_Widowed,PTMARRY_NeverMarried, \
    PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --covar-variance-standardize \
  --glm hide-covar firth-fallback \
  --threads 1 \
  --memory 8000 \
  --out ADNI_GWAS_ADvsCN

# OUTPUT:
Start time: Thu Sep 25 00:06:22 2025
36864 MiB RAM detected; reserving 8000 MiB for main workspace.
Using 1 compute thread.
1477 samples (633 females, 844 males; 1477 founders) loaded from
ADNI_qc_final.fam.
7918996 variants loaded from ADNI_qc_final.bim.
1 binary phenotype loaded (294 cases, 441 controls).
--keep: 735 samples remaining.
21 covariates loaded from ADNI_covar_ADvsCN.txt.
735 samples (340 females, 395 males; 735 founders) remaining after main
filters.
294 cases and 441 controls remaining after main filters.
--covar-variance-standardize: 21 covariates transformed.
Calculating allele frequencies... done.
--glm logistic-Firth hybrid regression on phenotype 'DIAG01': 1%
```

We now visualize/select the significant SNPs using different thresholds:
```
# Refer to https://cloufield.github.io/GWASTutorial/Visualization/#qc-check for different visualization techniques
python gwaslab_plot.py 
```
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/gwas_mqq_plot.png" width="800">
<img src="https://github.com/kimtae55/GWAS-End-to-End-Tutorial/blob/main/figs/gwas_chr19_region.png" width="800">


(Optional: Preprocessing steps for Imaging Genetics SCCA application (see last section))
- try to maintain 300-400 SNPs? make sure p >> n 
- then, before SCCA, do linear regression (regress out the same covariates, and then i
- check if our snp data contains SNP for APOE4 (double check)
- for FDG PET, take the average across each ROI, then regress 


## Application: Sparse Canonical Correlation Analysis using Imaging-Omics (FDG-PET, SNP)

Goal: To investigate how genetic variation (SNPs) and brain imaging features (FDG-PET ROIs) are related, and to assess how these associations change across disease stages (CN, MCI, AD). We use Sparse Canonical Correlation Analysis (SCCA) to identify low-dimensional, interpretable patterns linking high-dimensional SNP data with FDG-PET imaging phenotypes.

Input: 
- (n,p) SNP data matrix (split into CN, MCI, AD)
- (n,q) FDG-PET data matrix (split into CN, MCI, AD)


Please read [this]() for the complete application details. 
