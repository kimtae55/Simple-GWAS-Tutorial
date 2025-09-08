library("tidyverse")
ADNI1.snp <- read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI_cluster_01_forward_757LONI_14.bim",header=F)
ADNIGO2.snp <- read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNIGO2_SNP_QC/ADNI_GO2_10.bim",header=F)
ADNI3.snp <- read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI3_SNP_QC/ADNI3_SNP_10.bim",header=F)

ADNI1.snp$newID <- str_c(ADNI1.snp$V1, ADNI1.snp$V4, ADNI1.snp$V5, ADNI1.snp$V6,sep=":")
ADNIGO2.snp$newID <- str_c(ADNIGO2.snp$V1, ADNIGO2.snp$V4, ADNIGO2.snp$V5, ADNIGO2.snp$V6,sep=":")
ADNI3.snp$newID <- str_c(ADNI3.snp$V1, ADNI3.snp$V4, ADNI3.snp$V5, ADNI3.snp$V6,sep=":")

ADNI1.snp <- ADNI1.snp %>% select("V2", "newID")
ADNIGO2.snp <- ADNIGO2.snp %>% select("V2", "newID")
ADNI3.snp <- ADNI3.snp %>% select("V2", "newID")

write.table(ADNI1.snp, "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI1_snp_rename.txt", row.names = F, col.names = F, quote = F)
write.table(ADNIGO2.snp, "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNIGO2_SNP_QC/ADNIGO2_snp_rename.txt", row.names = F, col.names = F, quote = F)
write.table(ADNI3.snp, "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI3_SNP_QC/ADNI3_snp_rename.txt", row.names = F, col.names = F, quote = F)

ADNI.common.snp <- intersect(intersect(ADNI1.snp$newID,ADNIGO2.snp$newID),ADNI3.snp$newID)
length(ADNI.common.snp) #84104

ADNI1.snp.common <- ADNI1.snp$newID[ADNI1.snp$newID%in%ADNI.common.snp]
ADNIGO2.snp.common <- ADNIGO2.snp$newID[ADNIGO2.snp$newID%in%ADNI.common.snp]
ADNI3.snp.common <- ADNI3.snp$newID[ADNI3.snp$newID%in%ADNI.common.snp]

#remove duplicate SNPs
ADNI1.snp.com.duplicate <- ADNI1.snp.common[duplicated(ADNI1.snp.common)]
ADNIGO2.snp.com.duplicate <- ADNIGO2.snp.common[duplicated(ADNIGO2.snp.common)]
ADNI3.snp.com.duplicate <- ADNI3.snp.common[duplicated(ADNI3.snp.common)]

ADNI.common.snp.duplicate <-  c(ADNI1.snp.com.duplicate,ADNIGO2.snp.com.duplicate,ADNI3.snp.com.duplicate)
ADNI.common.snp <- setdiff(ADNI.common.snp, ADNI.common.snp.duplicate)
length(ADNI.common.snp) #84007

write.table(ADNI.common.snp, "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI_common_snp.txt", row.names = F, col.names = F, quote = F)


