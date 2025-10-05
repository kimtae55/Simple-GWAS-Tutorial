# GWAS Tutorial for hg38/GRCh38 with PLINK1.9
Most GWAS tutorials available only support plink for older genome build (e.g. hg18, hg19), so this tutorial provides plink1.9/hg38/TopMed reference panel compatible instructions. 

Main quality control process had referenced the steps from [Andries Marees's GWS Tutorial](https://github.com/MareesAT/GWA_tutorial). [Liftover](https://martha-labbook.netlify.app/posts/converting-snp-array-data/) and [data merging](https://martha-labbook.netlify.app/posts/extracting-data-for-variants-common-in-both-file-sets/) were referenced from Martha Aquino. Post-imputation processing referenced [Hai Shu's steps](https://github.com/shu-hai/imputation).

This tutorial follows the following steps: Quality control, Liftover, Imputation, Data merging, GWAS test.

## Prerequisites

### Software Used
- **PLINK**  ≥ v1.90 (v1.90-b7.7 local, v1.90b6.21 hpc)
- **R 4.5.0**
- **Python 3.13.0**
- **BCFtools** ≥ v1.14 (v1.22 local, v1.14 hpc)
- **Perl 5.32.0** 

Parts of our process in imputation and merging are completed in HPC. 

### Input data in PLINK format
Your input data should be in PLINK binary format, with the following extensions: 
- `.bed` → binary genotype data
- `.bim` → SNP information
- `.fam` → family/individual information

Data of this tutorial can be downloaded from the ADNI1, ADNI2 + ADNIGO, and ADNI3 [project official site](https://adni.loni.usc.edu/data-samples/adni-data/). 

Clinical data from ADNI can be found from the same site. We used the file `ADNIMERGE.csv` for subject phenotypes and `DXSUM.csv` for diagosis information. 


## Quality Control
Quality control is done separatedly for each dataset in every cohort. Among cohorts, data for ADNI2 and ADNIGO cohort and ADNI3 cohort are separated to two sets for downloading from online database. We merge the datasets within each cohort to a single dataset for both ADNIGO2 and ADNI3 later after liftover to build hg38. 

Noted we did not include the population stratification in this step. Population information will be controled through PCA extraction for GWAS later. 

The following quality control is an example with ADNI1 (`ADNI_cluster_01_forward_757LONI`). 

1) Filter SNP & individuals with missingness > 5% 
```
# filter for snp
plink --bfile ADNI_cluster_01_forward_757LONI\
      --geno 0.05 \
      --make-bed --out ADNI_snpfilt 

# filter for individuals 
plink --bfile ADNI_snpfilt\
      --mind 0.05\
      --make-bed --out ADNI1_s1
```

2) Remove individuals with sex discrepancy. 
By default, the `--check-sex` command of PLINK filters out F estimates > 0.2 for females and < 0.8 for males, and flagged with "PROBLEM". Source command can be found https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex. We use this command to filter out all individuals failed to pass the threshold and remove from our data. 

```
plink --bfile ADNI1_s1 --check-sex
awk '$5 == "PROBLEM" {print $1, $2}' plink.sexcheck > sex_discrepancy.txt

plink --bfile ADNI1_s1\
      --remove sex_discrepancy.txt\
      --make-bed --out ADNI1_s2
```

3) Filter autosomal SNPs only and delete SNPs with minor allele frequency (MAF) < 0.01.
```
# Select chromosomes 1 to 22
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ADNI1_s2.bim > snp_1_22.txt
plink --bfile ADNI1_s2\
      --extract snp_1_22.txt\
      --make-bed --out ADNI1_s2_1

# Remove SNPs with a low MAF frequency
plink --bfile ADNI1_s2_1\
      --maf 0.01
      \--make-bed --out ADNI1_s3
```

4) Remove SNPs with Hardy-Weinberg equilibrium (HWE) > 1e^-6.
For our study, no control and sample group defined. 
```
plink --bfile ADNI1_s3\
      --hwe 1e-6\
      --hwe-all\
      --make-bed --out ADNI1_s4
```

5) Filter individuals with heterozygozygosity rate < 3SD from mean. The high linkage disequilibrium region lists had referenced from [Hannah Meyer's documentation on GWAS QC](https://github.com/cran/plinkQC/tree/master). File of the corresponding genome build should be used for each cohort (hg18 for ADNI1, hg19 for ADNIGO2 adn ADNI3). Noted that we will need the hg38 high-LD file for later quality control after liftover and data merging. 

```
# Remove high-LD regions
plink --bfile ADNI1_s4 \
      --exclude range high-LD-regions-hg18-NCBI36.txt \
      --indep-pairwise 50 5 0.2 \
      --out indepSNP

# Check heterozygosity 
plink --bfile ADNI1_s4 --extract indepSNP.prune.in --het --out R_check

Rscript --no-save heterozygosity_outliers_list.R
# Output of the command above: fail-het-qc.txt .
# Convert for PLINK compatible 
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.
plink --bfile ADNI1_s4\
      --remove het_fail_ind.txt\
      --make-bed --out ADNI1_s5
```

6) Check for relatedness, exclude individuals with pihat > 0.2

```
# check for relationships between individuals with a pihat > 0.2
plink --bfile ADNI1_s5\
      --extract indepSNP.prune.in\
      --genome\
      --min 0.2\
      --out pihat_min0.2

# Filter pairs with pihat > 0.9 (duplicates/MZ twins)
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome

# Only include founderes 
plink --bfile ADNI1_s5\
      --filter-founders\
      --make-bed --out ADNI1_s5_1

# Recompute relatedness among founders with a pihat > 0.2
plink --bfile ADNI1_s5_1\ 
      --extract indepSNP.prune.in\
      --genome\
      --min 0.2\
      --out pihat_min0.2_in_founders

# Compute per-individual missingness rate 
plink --bfile ADNI1_s5_1 --missing

# Drop the individaul with higher missingness in each pair to drop. This step is completed through Python. 
python3
import pandas as pd

genome = pd.read_csv("pihat_min0.2_in_founders.genome", delim_whitespace=True)
related = genome[genome["PI_HAT"] > 0.2]
imiss = pd.read_csv("plink.imiss", delim_whitespace=True, usecols=["FID", "IID", "F_MISS"])

to_remove = []
for _, row in related.iterrows():
    id1 = (row["FID1"], row["IID1"])
    id2 = (row["FID2"], row["IID2"])

    # retrieve F_MISS
    fmiss1 = imiss[(imiss["FID"] == id1[0]) & (imiss["IID"] == id1[1])]["F_MISS"].values[0]
    fmiss2 = imiss[(imiss["FID"] == id2[0]) & (imiss["IID"] == id2[1])]["F_MISS"].values[0]

    # Decide which one to remove
    to_remove.append(id1 if fmiss1 > fmiss2 else id2)

# output 
with open("0.2_low_call_rate_pihat.txt", "w") as f:
    for fid, iid in to_remove:
        f.write(f"{fid}\t{iid}\n")

quit()

# Remove the selected individuals 
plink --bfile ADNI1_s5_1\
      --remove 0.2_low_call_rate_pihat.txt\
      --make-bed --out ADNI1_s6
```

## Liftover to hg38/GRCh38
For the origin dataset, ADNI1 is in hg18/NCBI36 and the rest cohorts are in hg19/NCBI37. In this step, we lift all datasets to GRCh38/hg38; this step can be skipped if your data is already in build hg38/GRCh38. 

LiftOver and chain files can be found from [UCSC LiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver). Lifting files to lift from both hg18 and hg19 directly to hg38 available. 

Liftover is done separatedly for each dataset in every cohort. The following quality control is an continued example with ADNI1 (`ADNI1_s6`) from hg18 to hg38, which we used `hg18ToHg38.over.chain.gz` from UCSC. 

```
# Recode to .ped/.map 
plink --bfile ADNI1_s6 --recode --tab --out ADNI1_s6_1

# liftover on .map files
python3 liftOverPlink_backup.py \
  --bin /usr/local/bin/liftOver \
  --map ADNI1_s6_1.map \
  --out ADNI1_s6_1_lifted \
  --chain hg18ToHg38.over.chain.gz

# clean bad SNPs
python3 liftOverPlink_MarthaLab/rmBadLifts.py \
  --map ADNI1_s6_1_lifted.map \
  --out ADNI1_s6_1_good_lifted.map \
  --log ADNI1_s6_1_bad_lifted.dat

# create & clean exclusion list 
cut -f2 ADNI1_s6_1_bad_lifted.dat > ADNI1_s6_1_snps_exclude.dat
cut -f4 ADNI1_s6_1_lifted.bed.unlifted | sed "/^#/d" >> ADNI1_s6_1_snps_exclude.dat
plink --file ADNI1_s6_1\
      --recode\
      --out ADNI1_s6_s1\
      --exclude ADNI1_s6_1_snps_exclude.dat

# merge cleaned .ped with .map & convert back to plink binary
plink --ped ADNI1_s6_s1.ped\
      --map ADNI1_s6_1_good_lifted.map\
      --recode --out ADNI1_s6_s2
plink --file ADNI1_s6_s2 --make-bed --out ADNI1_s6_s3
```

After confirming lifting all datasets into hg38, we merge the two datasets in ADNIGO2 and ADNI3. Below is an example for merging the two sets in ADNIGO2.

```
plink --bfile ADNI2set1_s6_s3\
      --bmerge ADNI2set2_s6_s3.bed ADNI2set2_s6_s3.bim ADNI2set2_s6_s3.fam\
      --make-bed --out ADNI_GO2
```

## Imputation
Imputation is perform on the [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!). Compared to Michigan Imputation Server with reference panel HRC for hg38, the sample size and ancestry distribution for TOPMed dataset provides higher coverage on the TOPMed Imputation Server for hg38. ([Sengupta et al., 2023](https://www.sciencedirect.com/science/article/pii/S2666979X23001003))

Both Michigan Imputation Server and TOPMed Imputation Server requires [data preparation](https://topmedimpute.readthedocs.io/en/latest/). For GRCh38/hg38: 
- Both Server recommanded to use TOPMed reference panel with process from Will Rayner: http://www.well.ox.ac.uk/~wrayner/tools/, which we followed here in this tutorial. You will need `HRC-1000G-check-bim-v4.3.0.zip` from the Rayner's process and `ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz` from TOPMed for this step. 
- Chromosomes should be encoded with prefix 'chr'.

Imputation is done separatedly for each cohort; we merged datasets within ADNIGO2 and ADNI3, leading us to submit three jobs (for ADNI1, ADNIGO2, and ADNI3) to the Imputation Server. 

The following quality control is an continued example with ADNI1 (`ADNI1_s6_s3`). For ADNIGO2 and ADNI3, this preparation should performed on the merged data.

```
# prepre TOPMed dataset (Will Rayner)
perl CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
# change to /usr/bin/gunzip; output file: PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
gunzip PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
unzip HRC-1000G-check-bim-v4.3.0.zip

# get allele frequency
plink --freq --bfile ADNI1_s6_s3 --out ADNI1_s6_s3.freq

# cd /scratch/zy2412/ADNI_SNP_0916/ADNI1_QC
sbatch run_perl.sh  

# match the data with TOPMed. We completed this substep in HPC cluster.
perl HRC-1000G-check-bim.pl \
  -b ADNI1_s6_s3.bim \
  -f ADNI1_s6_s3.freq.frq \
  -r PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab \
  -h
sh Run-plink.sh

# sort for download from HPC
module load bcftools/intel/1.14

for chr in {1..22}
do
  bcftools sort ADNI1_s6_s3-updated-chr${chr}.vcf -Ou | \
  bcftools annotate --rename-chrs chr_rename_map.txt -Oz \
    -o ADNI1_s6_s3-updated-chr${chr}.vcf.gz
done
```

This should give a list of ADNI1_s6_s3-updated-chr*.vcf.gz for all available chromosome# ready to upload for imputation. 

On TOPMed Imputation Server, job for each cohort should be submit separately. We used inputs below for each job. 

```
Ref panel: TOPMed r3
rsq filter: OFF
Phasing: Eagle v2.4
Population: All
Mode: Quality control and imputation
```

Once completed, imputated data can be downloaded through cURL from the Server. We filter SNPs with R^2 > 3 after concatenating all SNPs files, for imputation quality control. Below is an example for ADNI1, same process need for ADNIGO2 and ADNI3. 

```
# The imputated files will have names of chr*.dose.vcf.gz after unzip

# Concatenate the files in HPC 
bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou | 

# Filter R^2 > 3 in HPC
bcftools view -Ou -i 'R2>0.3' |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o ADNI1_allchromosomes.converted.R2_0.3.vcf.gz

plink --vcf ADNI1_allchromosomes.converted.R2_0.3.vcf.gz --double-id --make-bed --out ADNI1_allchromosomes.converted.R2_0.3
```

## Merge cohorts into a single file
In this step, common variants of each cohorts were extracted for merging. Once merged, reperform quality control with same pipeline as in the previous step. 

1) Extract common SNP for merging 
```
# Extract SNP from .bim and sort numerically
awk '{print $2}' ADNI1_allchromosomes.converted.R2_0.3.bim | sort > ADNI1_snp_sorted.txt
awk '{print $2}' ADNI2_allchromosomes.converted.R2_0.3.bim | sort > ADNIGO2_snp_sorted.txt
awk '{print $2}' ADNI3_allchromosomes.converted.R2_0.3.bim | sort > ADNI3_snp_sorted.txt

# Extract common IDs in three cohorts
comm -12 ADNI1_snp_sorted.txt ADNIGO2_snp_sorted.txt > intersect_ADNI1_GO2_snps.txt
comm -12 intersect_ADNI1_GO2_snps.txt ADNI3_snp_sorted.txt > intersect_ADNI1_GO2_3_snps.txt
```

2) Merge data based on extracted common SNPs

```
# Extract data using intersect_ADNI1_GO2_3_snps.txt from imputed data
plink --bfile ADNI1_allchromosomes.converted.R2_0.3 --extract intersect_ADNI1_GO2_3_snps.txt --make-bed --out ADNI1_intersect_snps
plink --bfile ADNI2_allchromosomes.converted.R2_0.3 --extract intersect_ADNI1_GO2_3_snps.txt --make-bed --out ADNI2_intersect_snps
plink --bfile ADNI3_allchromosomes.converted.R2_0.3 --extract intersect_ADNI1_GO2_3_snps.txt --make-bed --out ADNI3_intersect_snps

## Merge the file sets to obtain a list of variants with multiple positions/chromosomes and multiple alleles
plink --merge-list allfiles2merge.txt --make-bed --out ADNI1_GO2_3_merged_1

# Remove SNPs duplicated by position
awk '{print $0, $1":"$4}' ADNI1_GO2_3_merged_1.bim > post_imputation_updated_positions
awk '{print $1":"$4}' ADNI1_GO2_3_merged_1.bim | sort | uniq -d > post_imputation_updated_duplicated_positions 
grep -w -f  post_imputation_updated_duplicated_positions post_imputation_updated_positions | awk '{print $2}' > post_imputation_updated_duplicated_IDs

plink \
--bfile ADNI1_GO2_3_merged_1 \
--exclude post_imputation_updated_duplicated_IDs \
--make-bed \
--out ADNI1_GO2_3_merged_2

# write out merged files
plink --bfile ADNI1_GO2_3_merged_2 --snps-only --make-bed --out ADNI1_GO2_3_merged_snps
```

3) Quality control on the merged file 

```
# the same thresholds used as previous quality control process
plink --bfile ADNI1_GO2_3_merged_snps \
      --maf 0.01 \
      --hwe 1e-6 \
      --hwe-all \
      --geno 0.05 \
      --make-bed \
      --out ADNI1_GO2_3_merged_snps_qc
```

Relatedness check for this step performed based on subject information extracted from the quality controled original cohort data before imputation. We performed this step in R. 

```
library("tidyverse")

ADNI1.snp <- read.table("ADNI1_s6_s3.bim",header=F)
ADNIGO2.snp <- read.table("ADNI_GO2.bim",header=F)
ADNI3.snp <- read.table("ADNI_3.bim",header=F)

ADNI1.snp$newID <- str_c(ADNI1.snp$V1, ADNI1.snp$V4, ADNI1.snp$V5, ADNI1.snp$V6,sep=":")
ADNIGO2.snp$newID <- str_c(ADNIGO2.snp$V1, ADNIGO2.snp$V4, ADNIGO2.snp$V5, ADNIGO2.snp$V6,sep=":")
ADNI3.snp$newID <- str_c(ADNI3.snp$V1, ADNI3.snp$V4, ADNI3.snp$V5, ADNI3.snp$V6,sep=":")

ADNI1.snp <- ADNI1.snp %>% select("V2", "newID")
ADNIGO2.snp <- ADNIGO2.snp %>% select("V2", "newID")
ADNI3.snp <- ADNI3.snp %>% select("V2", "newID")

write.table(ADNI1.snp, "ADNI1_snp_rename.txt", row.names = F, col.names = F, quote = F)
write.table(ADNIGO2.snp, "ADNI2_snp_rename.txt", row.names = F, col.names = F, quote = F)
write.table(ADNI3.snp, "ADNI3_snp_rename.txt", row.names = F, col.names = F, quote = F)

# Concatecate 
ADNI.common.snp <- intersect(intersect(ADNI1.snp$newID,ADNIGO2.snp$newID),ADNI3.snp$newID)
length(ADNI.common.snp) # 84524

ADNI1.snp.common <- ADNI1.snp$newID[ADNI1.snp$newID%in%ADNI.common.snp]
ADNIGO2.snp.common <- ADNIGO2.snp$newID[ADNIGO2.snp$newID%in%ADNI.common.snp]
ADNI3.snp.common <- ADNI3.snp$newID[ADNI3.snp$newID%in%ADNI.common.snp]

# remove duplicate SNPs
ADNI1.snp.com.duplicate <- ADNI1.snp.common[duplicated(ADNI1.snp.common)]
ADNIGO2.snp.com.duplicate <- ADNIGO2.snp.common[duplicated(ADNIGO2.snp.common)]
ADNI3.snp.com.duplicate <- ADNI3.snp.common[duplicated(ADNI3.snp.common)]

ADNI.common.snp.duplicate <-  c(ADNI1.snp.com.duplicate,ADNIGO2.snp.com.duplicate,ADNI3.snp.com.duplicate)
ADNI.common.snp <- setdiff(ADNI.common.snp, ADNI.common.snp.duplicate)

write.table(ADNI.common.snp, "ADNI_common_snp.txt", row.names = F, col.names = F, quote = F)
``` 

Based on the extracted `ADNI_common_snp.txt` we check for subject relatedness. 

```
# Extract SNPs based on the common SNP list 
plink --bfile ADNI1_s6_s3 --update-name ADNI1_snp_rename.txt --make-bed --out ADNI1_s6_s3_renameID
plink --bfile ADNI1_s6_s3_renameID --extract ADNI_common_snp.txt --make-bed --out ADNI1_commonSNP_data

plink --bfile ADNI_GO2 --update-name ADNI2_snp_rename.txt --make-bed --out ADNI_GO2_renameID
plink --bfile ADNI_GO2_renameID --extract ADNI_common_snp.txt --make-bed --out ADNI2_commonSNP_data

plink --bfile ADNI_3 --update-name ADNI3_snp_rename.txt --make-bed --out ADNI_3_renameID
plink --bfile ADNI_3_renameID --extract ADNI_common_snp.txt --make-bed --out ADNI3_commonSNP_data

# Merge extracted data for relatedness check 
plink --merge-list allfiles2merge_CheckSubjectRelatedness.txt --make-bed --out ADNI1_GO2_3_merged_CheckSubjectRelatedness

# Exclude high linkage disequilibrium regions, reference list must be for hg38. 
plink --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness --exclude range high-LD-regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.2 --out indepSNP_ADNIall_CheckSubjectRelatedness

# Continue with previous QC check for relatedness (same as previous QC)
plink --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness --filter-founders --make-bed --out ADNI1_GO2_3_merged_CheckSubjectRelatedness_2
plink --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness_2 --extract indepSNP_ADNIall_CheckSubjectRelatedness.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders_ADNIall_CheckSubjectRelatedness
plink --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness_2 --missing

# Check for missingness and remove the one with higher F_MISS (same as previous QC)
python3
import pandas as pd

genome = pd.read_csv("pihat_min0.2_in_founders_ADNIall_CheckSubjectRelatedness.genome", delim_whitespace=True)
related = genome[genome["PI_HAT"] > 0.2]
imiss = pd.read_csv("plink.imiss", delim_whitespace=True, usecols=["FID", "IID", "F_MISS"])

to_remove = []
for _, row in related.iterrows():
    id1 = (row["FID1"], row["IID1"])
    id2 = (row["FID2"], row["IID2"])

    # Get F_MISS values
    fmiss1 = imiss[(imiss["FID"] == id1[0]) & (imiss["IID"] == id1[1])]["F_MISS"].values[0]
    fmiss2 = imiss[(imiss["FID"] == id2[0]) & (imiss["IID"] == id2[1])]["F_MISS"].values[0]

    # Decide which one to remove
    to_remove.append(id1 if fmiss1 > fmiss2 else id2)

# save to file
with open("0.2_low_call_rate_pihat_ADNIall_CheckSubjectRelatedness.txt", "w") as f:
    for fid, iid in to_remove:
        f.write(f"{fid}\t{iid}\n")   # output has zero index

quit()


# Continue removal (same as previous QC)
plink --bfile ADNI1_GO2_3_merged_CheckSubjectRelatedness_2\
      --remove 0.2_low_call_rate_pihat_ADNIall_CheckSubjectRelatedness.txt\
      --make-bed --out ADNI1_GO2_3_merged_CheckSubjectRelatedness_3

awk '{print $1"_"$2, $1"_"$2}' ADNI1_GO2_3_merged_CheckSubjectRelatedness_3.fam > ADNIall_CheckSubjectRelatedness_subj.txt 

plink --bfile ADNI1_GO2_3_merged_snps_qc\
      --keep ADNIall_CheckSubjectRelatedness_subj.txt\
      --make-bed --out ADNI1_GO2_3_merged_snps.R2_0.3.MAF_0.01.HWE_0.000001.geno_0.05.IBD_0.2
```


## GWAS test
### Prepare covariant files 
Two clinical files, `ADNIMERGE.csv` and `DXSUM.csv`, (plus their corresponding dictionaries) are needed for this step. All files are available through ADNI online database. Within both files, the column `VISCODE` provides data recording phase; only the baseline `bl` observation should be included.

For covariants, we need to merge clinical data for below information:
- **subject list** from subject lists of pruned data for PCA.
- **covariants** we compiled sex (from datasets before imputation), six patient information (from `ADNIMERGE.csv`), and top 10 principal components from genotype PCA for ancestry controling. 
- **phenotype** from `DIAGNOSIS` of `DXSUM.csv`; in this tutorial, we only included CN and AD. 

```
# Prune data of each cohort & perform PCA
plink --bfile ADNI1_s6_s3 --exclude range high-LD-regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.2 --out pruned_data
plink --bfile ADNI1_s6_s3 --extract pruned_data.prune.in --make-bed --out ADNI1_s6_s3_pruned
plink --bfile ADNI1_s6_s3_pruned --pca 10 --out ADNI1_PCA

plink --bfile ADNI_GO2 --exclude range high-LD-regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.2 --out pruned_data
plink --bfile ADNI_GO2 --extract pruned_data.prune.in --make-bed --out ADNI_GO2_pruned
plink --bfile ADNI_GO2_pruned --pca 10 --out ADNIGO2_PCA

plink --bfile ADNI_3 --exclude range high-LD-regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.2 --out pruned_data
plink --bfile ADNI_3 --extract pruned_data.prune.in --make-bed --out ADNI_3_pruned
plink --bfile ADNI_3_pruned --pca 10 --out ADNI3_PCA

# Retreive covariants
Rscript subjects_bl_covar.R

# Recode to indicator variables for PLINK1.9
Rscript subjects_bl_covar_addindicator.R
```

The above will generate covariants of below: 
```
**ID variables**
- `FID`, `IID`: family and individual IDs from PLINK `.fam`.  

**Covariant: sex (from PLINK `.fam`)**
-  `SEX2`: 0 if **male**, 1 if **female**.

**Covariant: phase (COLPROT of ADNIMERGE)**
Reference category = **ADNI1** (coded 0)
- `COLPROT2`      = 1 if ADNI2, else 0  
- `COLPROT3`      = 1 if ADNIGO, else 0  
- `COLPROT4`      = 1 if ADNI3, else 0  

**Covariant: race (PTRACCAT of ADNIMERGE)**
Reference category = **White** (coded 0)
- `PTRACCAT2`      = 1 if Black, else 0  
- `PTRACCAT3`      = 1 if Asian, else 0  
- `PTRACCAT4`      = 1 if More than one race, else 0  
- `PTRACCAT5`      = 1 if American Indian / Alaskan Native, else 0  
- `PTRACCAT6`      = 1 if Native Hawaiian/Other Pacific Islander, else 0  

**Covariant: ethnicity (PTETHCAT of ADNIMERGE)**
Reference category = **Hisp/Latino** (coded 0)
- `PTETHCAT2`      = 1 if Not Hisp/Latino, else 0 

**Covariant: education (PTEDUCAT of ADNIMERGE)**
- `PTEDUCAT`: numeric

**Covariant: maritual status (PTMARRY of ADNIMERGE)**
Reference category = **Divorced** (coded 0)
- `PTMARRY2`      = 1 if Married, else 0  
- `PTMARRY3`      = 1 if Never married, else 0  
- `PTMARRY4`      = 1 if Widowed, else 0  

**Covariant: ancestry PCA**
- `PC1`–`PC10`: Top 10 principal components from genotype PCA (numeric).

**Phenotype: baseline status (DIAGNOSIS of DXSUM)**
- `status`: 1 if **CN**, 2 if **AD** 
```

### Association test
```
# Filter data matching the filtered subject list 
plink --bfile ADNI1_GO2_3_merged_snps.R2_0.3.MAF_0.01.HWE_0.000001.geno_0.05.IBD_0.2\
      --keep ADNIall_ADandCN_subj.txt\
      --make-bed --out ADNIall_ADandCN_1

# Attach phenotype information 
plink -bfile ADNIall_ADandCN_1 --pheno ADNIall_ADandCN_pheno.txt --make-bed --out ADNIall_ADandCN_2

# Associate test
plink --bfile ADNIall_ADandCN_2\
      --logistic hide-covar\
      --covar ADNIall_ADandCN_covar_plink.txt\
      --allow-no-sex\
      --out ADNIall_ADandCN_test_result

# get significant SNPs
Rscript getSigSNPs.R
plink --bfile ADNIall_ADandCN_2\
      --extract sigSNPs.txt\
      --make-bed --out ADNIall_ADandCN_3

# export results
plink --bfile ADNIall_ADandCN_3\
      --recodeA\
      --allow-no-sex\
      --out ADNIall_ADandCN_sigSNPs
```



