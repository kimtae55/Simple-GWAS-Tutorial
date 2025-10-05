library(dplyr)

covar <- read.table("ADNIall_ADandCN_covar.txt", header=TRUE)

# convert all categorical vars to factors
covar <- covar %>%
  mutate(across(c(SEX, PTETHCAT, COLPROT, PTRACCAT, PTMARRY), as.factor))

X <- model.matrix(~ AGE + SEX + COLPROT + PTRACCAT + PTETHCAT + PTEDUCAT + PTMARRY + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
                  data = covar)


# bind back
covar_dum <- cbind(covar[,c("FID","IID")], X[,-1])

write.table(covar_dum, "ADNIall_ADandCN_covar_plink.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)

