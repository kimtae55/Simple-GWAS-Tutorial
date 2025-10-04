import gwaslab as gl
import numpy as np

# 1. Load your PLINK2 GWAS results
# The key is to tell gwaslab which column has the p-values, effect allele, etc.
# For plink2 --glm files, use fmt="plink2"
gwas_file = "ADNI_GWAS_ADvsCN1.DIAG01.glm.logistic.hybrid"

# Load PLINK2 GWAS results
ss = gl.Sumstats(gwas_file, fmt="plink2")

df = ss.data.copy()
df = df.sort_values("P")
df[df["P"] < 5e-8 ].to_csv("GWAS_hits_suggestive_5e-8.csv",    index=False)
df[df["P"] < 5e-5 ].to_csv("GWAS_hits_suggestive_5e-5.csv",    index=False)
df[df["P"] < 1e-4 ].to_csv("GWAS_hits_suggestive_1e-4.csv",    index=False)
df[df["P"] < 2e-4 ].to_csv("GWAS_hits_suggestive_2e-4.csv",    index=False)
df[df["P"] < 3e-4 ].to_csv("GWAS_hits_suggestive_3e-4.csv",    index=False)
df[df["P"] < 4e-4 ].to_csv("GWAS_hits_suggestive_4e-4.csv",    index=False)
df[df["P"] < 5e-4 ].to_csv("GWAS_hits_suggestive_5e-4.csv",    index=False)
df[df["P"] < 5e-2 ].to_csv("GWAS_hits_suggestive_5e-2.csv",    index=False)

print('top sig snps saved')

# Manhattan plot and QQ plot
ss.plot_mqq(anno=True,
            save="gwas_mqq_plot.png", 
            save_args={"dpi":300, "facecolor":"white"}, 
            build='38')
ss.plot_mqq(mode="r", 
            anno=True,
            region=(19,44000000,46000000),   # 44–46 Mb region
            region_grid=True,
            build="38",
            save="gwas_chr19_region.png", 
            save_args={"dpi":300, "facecolor":"white"})

