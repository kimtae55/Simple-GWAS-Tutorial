library("tidyverse")
library("interval")

### import diagnosis 
diag.info = read.csv("DXSUM_11Aug2025.csv")[,c("PHASE", "RID", "PTID","VISCODE2", "DIAGNOSIS")]

diag.info <- diag.info %>%
  filter(VISCODE2 == "bl") %>%
  mutate(
    status = as.character(DIAGNOSIS),
    status = ifelse(!status %in% c("1", "2", "3"), "Other", status),
    status = recode(status, "1" = "CN", "2" = "MCI", "3" = "AD")
  )

diag.info <- diag.info %>% select("RID", "PTID", "status")
diag.info=diag.info[diag.info$status%in%c("AD","CN"),]

### get age, sex, batch, race, ethnicity, education, marital status
merge.info <- read.csv("ADNIMERGE_18Aug2025.csv") %>%
  select(RID, PTID, VISCODE, DX_bl, AGE, COLPROT, PTRACCAT, PTETHCAT, PTEDUCAT, PTMARRY) %>%
  filter(VISCODE == "bl")
  
# remove unknown or blank
merge.info <- merge.info %>%
  filter(PTRACCAT != "Unknown", PTETHCAT != "Unknown", !PTMARRY %in% c("Unknown", ""))

### import pruned data for subject list
ADNI1 = read.table("ADNI1_s6_s3.fam",header=F)
ADNIGO2 = read.table("ADNI_GO2.fam",header=F)
ADNI3 = read.table("ADNI_3.fam",header=F)
ADNIall = rbind(ADNI1, ADNIGO2, ADNI3)[,c(1, 2, 5)]
colnames(ADNIall)=c("FID", "IID", "SEX")

merge.info <- merge.info[merge.info$VISCODE=="bl",c("PTID","AGE","COLPROT","PTRACCAT","PTETHCAT","PTEDUCAT","PTMARRY")]
ADNIall <- ADNIall %>% left_join(merge.info, by=c("IID"="PTID"))

length(intersect(ADNIall$IID,diag.info$PTID))
ADNIall = ADNIall[ADNIall$IID%in%diag.info$PTID,]
ADNIall <- ADNIall %>% left_join(diag.info, by=c("IID"="PTID"))
ADNIall$RID = NULL

### get top 10 PCs
ADNI1.pca = read.table("ADNI1_PCA.eigenvec",header=F)
ADNIGO2.pca = read.table("ADNIGO2_PCA.eigenvec",header=F)
ADNI3.pca = read.table("ADNI3_PCA.eigenvec",header=F)

ADNIall.pca = rbind(ADNI1.pca,ADNIGO2.pca,ADNI3.pca)
pc = str_c("PC",1:10,sep="")
colnames(ADNIall.pca) = c("FID", "IID", pc)

ADNIall <- ADNIall %>% left_join(ADNIall.pca, by=c("FID","IID"))
ADNIall <- ADNIall %>% relocate(status,.after=PC10)






# clean up clinical data
##### recode for category
# diagnosis
ADNIall$status[ADNIall$status=="CN"]=1
ADNIall$status[ADNIall$status=="AD"]=2

# cohort
ADNIall$COLPROT[ADNIall$COLPROT=="ADNI1"]=1
ADNIall$COLPROT[ADNIall$COLPROT=="ADNI2"]=2
ADNIall$COLPROT[ADNIall$COLPROT=="ADNIGO"]=3
ADNIall$COLPROT[ADNIall$COLPROT=="ADNI3"]=4

# race
ADNIall$PTRACCAT[ADNIall$PTRACCAT=="White"]=1
ADNIall$PTRACCAT[ADNIall$PTRACCAT=="Black"]=2
ADNIall$PTRACCAT[ADNIall$PTRACCAT=="Asian"]=3
ADNIall$PTRACCAT[ADNIall$PTRACCAT=="More than one"]=4
ADNIall$PTRACCAT[ADNIall$PTRACCAT=="Am Indian/Alaskan"]=5
ADNIall$PTRACCAT[ADNIall$PTRACCAT=="Hawaiian/Other Pl"]=6

# ethnicity
ADNIall$PTETHCAT[ADNIall$PTETHCAT=="Hisp/Latino"]=1
ADNIall$PTETHCAT[ADNIall$PTETHCAT=="Not Hisp/Latino"]=2

# marital
ADNIall$PTMARRY[ADNIall$PTMARRY=="Divorced"]=1
ADNIall$PTMARRY[ADNIall$PTMARRY=="Married"]=2
ADNIall$PTMARRY[ADNIall$PTMARRY=="Never married"]=3
ADNIall$PTMARRY[ADNIall$PTMARRY=="Widowed"]=4





###### remove unneeded columns & remove blanks
ADNIall <- ADNIall[!apply(ADNIall == "" | is.na(ADNIall), 1, any), ]
ADNIall <- ADNIall[ , c(setdiff(names(ADNIall), "status"), "status") ]

### remove related
ADNIall$IID = str_c(ADNIall$FID,ADNIall$IID,sep="_")
ADNIall$FID = ADNIall$IID

ADNI.IBD <- read.table("ADNIall_CheckSubjectRelatedness_subj.txt",header=F)
ADNIall <- ADNIall[ADNIall$IID%in%ADNI.IBD[,1],]
dim(ADNIall) # 1114   20




# write files
write.table(ADNIall[,1:19], "ADNIall_ADandCN_covar.txt", row.names = F, quote = F)
write.table(ADNIall[,1:2], "ADNIall_ADandCN_subj.txt", row.names = F, quote = F)
write.table(ADNIall[,c(1:2,20)], "ADNIall_ADandCN_pheno.txt", row.names = F, quote = F)







