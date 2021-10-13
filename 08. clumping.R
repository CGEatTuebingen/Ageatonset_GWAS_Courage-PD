# Clumping
## A. Load relevant libraries
library(tidyverse)
library(data.table)

## B. Create chromosome specific input file for running clumping in Plink (R platform)
### Use the meta-analysis output generated in 04.linear_regression step (meta_analysis_courage_05052021.tbl)
sumstat <- fread("meta_analysis_courage_05052021.tbl")
sumstat_snpid <- left_join(sumstat, hrc_data); rm(sumstat)
sumstat_snpid  <- sumstat_snpid %>% rename(P = `P-value`)
sumstat_snpid  <- sumstat_snpid  %>% mutate(Allele1 = toupper(Allele1))
sumstat_snpid  <- sumstat_snpid  %>% mutate(Allele2 = toupper(Allele2))
sumstat_snpid  <- sumstat_snpid %>% separate("MarkerName", c("chr", "bp"), ":")
sumstat_snpid$chr <- as.numeric(as.character(sumstat_snpid$chr)) 
sumstat_snpid$bp <- as.numeric(as.character(sumstat_snpid$bp))
sumstat_snpid <- sumstat_snpid %>% filter(P < 0.00001)
uniq <- unique(unlist(sumstat_snpid$CHR))
dataset <- "courage_all"
for (i in 1:length(uniq)){
  data_1 <- subset(sumstat_snpid, CHR == uniq[i])
  myfile <- paste0("clump_input19df_mychr","_",uniq[i],"_", dataset,".txt" )
  write.table(data_1, myfile,row.names=FALSE,sep="\t", quote = FALSE) 
}

### Content of the sample input file "clump_input19df_mychr_5_courage_all.txt"
# CHR	SNP	BP	A1	SE	P
# 5	rs1155605	102799858	T	0.294	9.777e-07
# 5	rs983119	102761802	A	0.2916	6.905e-06
# 5	rs2099077	102764048	A	0.2902	5.255e-06
# 5	rs12515712	102781210	T	0.2961	2.615e-06
# 5	rs12518613	102788222	T	0.2935	8.719e-07
# 5	rs1014772	102788681	A	0.294	1.236e-06
# 5	rs1428490	102798969	A	0.294	9.777e-07
# 5	rs1428489	102810436	T	0.3097	7.434e-06
# 5	rs4142813	102765227	T	0.2916	6.597e-06
# 5	rs17155310	102757059	A	0.3029	8.176e-06
# 5	rs62365561	102795941	A	0.2934	6.906e-07
# 5	rs1507768	102817024	T	0.3118	6.899e-06
# 5	rs10065975	102768968	T	0.2897	1.603e-06
# 5	rs12110237	102796072	A	0.294	1.263e-06
# 5	rs76833279	102783574	A	0.2943	1.505e-06
# 5	rs10515345	102786159	A	0.2932	9.005e-07
# 5	rs10479252	102776916	T	0.2953	1.913e-06
# 5	rs1864084	102791137	A	0.2949	3.098e-06


## C.  Run for each chromosome on Plink platform
plink --bfile ./HRC_ld/chr1 --clump clump_input19df_mychr_1_courage_all.txt --clump-field P --clump-kb 250 --clump-p1 1 --clump-r2 0.1 --clump-snp-field SNP --memory 120000 --out chr1_clump_courage19df_all 

## D. Merge all the chromosome specific output files on R platform
masterlist =list.files(pattern = "19df_all\\.clumped$")
l <- lapply(masterlist, fread)
dt <- rbindlist(l)
dt <- dt %>% select(CHR,SNP, BP, P, SP2) %>%  mutate(SP2 = gsub("\\(1\\)", "", SP2))
write.csv(dt, "courage_clumping19df_all_results.csv")

### Content of the sample output file "courage_clumping19df_all_results.csv"
# CHR	SNP	BP	P	SP2
# 1	1	rs116768866	30030233	8.08E-06	NONE
# 2	10	rs113744107	63326120	4.68E-06	rs111364217
# 3	14	rs1046099	93649501	3.09E-06	NONE
# 4	16	rs6498855	62458351	6.85E-06	rs11859572,rs9938742
# 5	2	rs138663448	204719819	4.58E-06	NONE
# 6	2	rs6761931	136150512	7.29E-06	NONE
# 7	2	rs834139	158073096	8.98E-06	rs834138
# 8	21	rs9978043	22664753	9.44E-06	rs8126533,rs232475
# 9	3	rs62245756	72054924	5.47E-06	rs62245757
# 10	4	rs61621377	105504104	4.24E-06	rs62329658,rs10488867,rs10488868,rs17211918
# 11	4	rs147125739	94226785	8.41E-06	NONE
# 12	4	rs979511	143360851	9.03E-06	NONE
# 13	5	rs62365561	102795941	6.91E-07	rs17155310,rs983119,rs2099077,rs4142813,rs10065975,rs10479252,rs12515712,rs76833279,rs10515345,rs12518613,rs1014772,rs1864084,rs12110237,rs1428490,rs1155605,rs1428489,rs1507768
# 14	7	rs2041640	151224309	4.96E-06	rs2041639
# 15	8	rs13267794	91889302	2.10E-06	rs34649228,rs13248156

### In addition, various SNP annotation tools were used for functional annotation and and identification of overlapping and nearby genes
# http://db.systemsbiology.net/kaviar/
# https://snipa.helmholtz-muenchen.de/snipa3/
# https://www.snp-nexus.org/v4/
