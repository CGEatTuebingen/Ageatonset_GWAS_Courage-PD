## Generate Locuszoom input file
### Filter chromosome 4 SNPs
metanalysis_courage_chr4 <-
  fread("meta_analysis_courage_05052021.tbl") %>% separate(MarkerName, into = c("CHR", "POS"), sep = ":") %>% filter(HetDf >
                                                                                                                      18) %>% filter(HetISq < 50) %>% rename(
                                                                                                                        EA = Allele1,
                                                                                                                        OA = Allele2,
                                                                                                                        beta = Effect,
                                                                                                                        se = StdErr,
                                                                                                                        Isq = HetISq
                                                                                                                      ) %>% mutate(CHR = as.numeric(CHR)) %>% filter(CHR == "4")


### Filter SNPs on 2kb on either side of the specific loci e.g. BST1 loci and create output file specific to locuszoom 
metanalysis_consortium_chr4_2kb_BST1 <-
  metanalysis_consortium_chr4 %>% mutate(POS = as.numeric(POS)) %>%
  filter(POS %in% c(15537348:15937348)) %>% mutate(SNP = paste(CHR, POS, sep =
                                                                 ":")) %>% rename("P-value" = "PVALUE") %>% select(c(SNP, EA, OA, beta, se, "P-value", Isq)) %>% mutate(EA = toupper(EA))  %>% mutate(OA = toupper(OA))
metanalysis_consortium_chr4_2kb_BST1$CHR <-
  paste("chr", metanalysis_consortium_chr4_2kb_BST1$CHR, sep = "")
metanalysis_consortium_chr4_2kb_BST1$SNP <-
  paste(metanalysis_consortium_chr4_2kb_BST1$CHR ,
        metanalysis_consortium_chr4_2kb_BST1$SNP,
        sep = "")
metanalysis_consortium_chr4_2kb_BST1  %>%  select(c(SNP, "P-value")) %>%  rename(MarkerName = SNP , P.value = "P-value") %>%
  write.table(
    "locuszoom_consortium_BST1.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

## Header of the input file for locuszoom plot                                                                                                                      ) %>% mutate(CHR = as.numeric(CHR)) %>% filter(CHR == "4")
# MarkerName	P.value
# chr4:15840444	0.05687
# chr4:15697925	0.2247
# chr4:15704663	0.3483
# chr4:15733454	0.04854
# chr4:15734652	0.04854
# chr4:15681719	0.7246
# chr4:15823964	0.06463
# chr4:15707001	0.211
# chr4:15647537	0.6079
# chr4:15728806	0.05818
# chr4:15728079	0.02503
# chr4:15765183	0.08403
# chr4:15811762	0.8286
# chr4:15715877	0.5744
# chr4:15837382	0.6321


### Go to http://locuszoom.org/genform.php?type=yourdata
# Upload txt file to Path to your file
# Give the name of SNP of interest e.g. chr4:15737348 in the section specific region to display along with a flanking size of 200kb
# Click on plot data 
# A pdf is generated. Use edit feature in pdf to rename the loci names with SNPids.
# Extract images from pdf as TIF file

