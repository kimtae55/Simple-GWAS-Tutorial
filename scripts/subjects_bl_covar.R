library("tidyverse")
library("intrval")
diag.info = read.csv("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_FDGPET/Clinica/CLINICAL_DATA/DXSUM_PDXCONV_ADNIALL.csv"
                    )[,c("Phase", "RID", "PTID","VISCODE2", "DXCHANGE", "DXCURREN", "DIAGNOSIS")]
diag.info <- diag.info %>% filter(VISCODE2=="bl")
diag.info$status[diag.info$Phase=="ADNI1"]=diag.info$DXCURREN[diag.info$Phase=="ADNI1"]
diag.info$status[diag.info$Phase%in%c("ADNIGO","ADNI2")]=diag.info$DXCHANGE[diag.info$Phase%in%c("ADNIGO","ADNI2")]
diag.info$status[diag.info$Phase=="ADNI3"]=diag.info$DIAGNOSIS[diag.info$Phase=="ADNI3"]

diag.info$status[!diag.info$status%in%c("1","2","3")]="Other"
diag.info$status[diag.info$status=="1"]="CN"
diag.info$status[diag.info$status=="2"]="MCI"
diag.info$status[diag.info$status=="3"]="AD"

table(diag.info$status, useNA="always")
#   AD    CN   MCI Other 
#  408   845  1060     6

diag.info <- diag.info %>% select("RID", "PTID", "status")
dim(diag.info)

diag.info=diag.info[diag.info$status%in%c("AD","CN"),]


ADNI1 = read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI_cluster_01_forward_757LONI_14_pruned.fam",header=F)
ADNIGO2 = read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNIGO2_SNP_QC/ADNI_GO2_10_pruned.fam",header=F)
ADNI3 = read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI3_SNP_QC/ADNI3_SNP_10_pruned.fam",header=F)

ADNIall = rbind(ADNI1,ADNIGO2,ADNI3)[,c(1,2,5)]
colnames(ADNIall)=c("FID", "IID", "SEX")

length(intersect(ADNIall$IID,diag.info$PTID)) #976

ADNIall = ADNIall[ADNIall$IID%in%diag.info$PTID,]

ADNIall <- ADNIall %>% left_join(diag.info, by=c("IID"="PTID"))
ADNIall$RID = NULL
dim(ADNIall)
length(unique(ADNIall$IID))

#get age
merge.info <- read.csv("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_FDGPET/Clinica/CLINICAL_DATA/ADNIMERGE.csv")[,c("RID", "PTID", "VISCODE","DX_bl","AGE")]
merge.info=merge.info[merge.info$VISCODE=="bl",c("PTID","AGE")]
ADNIall <- ADNIall %>% left_join(merge.info, by=c("IID"="PTID"))

# get top 10 PCs
ADNI1.pca = read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNI1_PCA.eigenvec",header=F)
ADNIGO2.pca = read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNIGO2_SNP_QC/ADNIGO2_PCA.eigenvec",header=F)
ADNI3.pca = read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI3_SNP_QC/ADNI3_PCA.eigenvec",header=F)

ADNIall.pca = rbind(ADNI1.pca,ADNIGO2.pca,ADNI3.pca)
pc = str_c("PC",1:10,sep="")
colnames(ADNIall.pca) = c("FID", "IID", pc)

ADNIall <- ADNIall %>% left_join(ADNIall.pca, by=c("FID","IID"))
dim(ADNIall) #976  15

ADNIall <- ADNIall %>% relocate(status,.after=PC10)
ADNIall$IID = str_c(ADNIall$FID,ADNIall$IID,sep="_")
ADNIall$FID = ADNIall$IID
ADNIall$status[ADNIall$status=="AD"]=2
ADNIall$status[ADNIall$status=="CN"]=1

ADNI.IBD <- read.table("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNIall_CheckSubjectRelatedness_subj.txt",header=F)

ADNIall <- ADNIall[ADNIall$IID%in%ADNI.IBD[,1],]

dim(ADNIall) # 973  15

write.table(ADNIall[,1:14], "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNIall_ADandCN_covar.txt", row.names = F, quote = F)
write.table(ADNIall[,1:2], "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNIall_ADandCN_subj.txt", row.names = F, quote = F)
write.table(ADNIall[,c(1:2,15)], "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNIall_ADandCN_pheno.txt", row.names = F, quote = F)







