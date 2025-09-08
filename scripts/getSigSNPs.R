library(data.table)
#library(GWASTools)
#library(qvalue)
#library(IHWpaper)#CLfdr
test.rlt <- fread("/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/ADNIall_ADandCN_test_result.assoc.logistic")
sum(na.omit(test.rlt$P)<5e-5)
test.rlt$P[is.na(test.rlt$P)]=1
sigSNPs <- test.rlt$SNP[test.rlt$P<5e-5]
length(sigSNPs)

write.table(sigSNPs, "/Users/haishu/Desktop/genetic_imaging/ADNI_data/ADNI_SNP/ADNI1_SNP_QC/sigSNPs.txt", row.names = F, quote = F, col.names = F)


#qqPlot(is.na(test.rlt$P))

#chisq <- qchisq(1-test.rlt$P,1)
#median(chisq,na.rm=T)/qchisq(0.5,1)

#test.rlt[is.na(test.rlt$P),]

#sum(p.adjust(na.omit(test.rlt$P),"BY")<0.1)

#hist(na.omit(test.rlt$P),nclass=20)

#qobj <- qvalue(na.omit(test.rlt$P),fdr.level =0.2)
#sum(qobj$significant)

#obj.clfdr <- clfdr(test.rlt$P[!is.na(test.rlt$P)],test.rlt$CHR[!is.na(test.rlt$P)],0.25)
#sum(rejected_hypotheses(obj.clfdr))
