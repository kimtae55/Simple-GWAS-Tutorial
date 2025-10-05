library(data.table)

test.rlt <- fread("ADNIall_ADandCN_test_result.assoc.logistic")

# threshold 5e-5
sum(na.omit(test.rlt$P)<5e-5)
test.rlt$P[is.na(test.rlt$P)]=1
sigSNPs <- test.rlt$SNP[test.rlt$P<5e-5]
length(sigSNPs)

write.table(sigSNPs, "sigSNPs.txt", row.names = F, quote = F, col.names = F)