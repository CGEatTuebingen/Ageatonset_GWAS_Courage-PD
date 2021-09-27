# Generate Locuszoom input file
## Filter chromosome 4 SNPs
metanalysis_consortium_chr4 <-
  fread("consortium_summarystatistics.txt") %>% mutate(MarkerName = str_sub(MarkerName, 4)) %>% separate(MarkerName, into = c("CHR", "POS"), sep = ":") %>%
  rename(
    EA = REF,
    OA = ALT,
    beta = EFFECT,
    se = SE,
    Isq = HetISq
  ) %>%
  mutate(CHR = as.numeric(CHR)) %>%
  filter(CHR == "4")


## Filter SNPs on 2kb on either side of the specific loci e.g. BST1 loci and create output file specific to locuszoom 
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


# Go to http://locuszoom.org/genform.php?type=yourdata
# Upload txt file to Path to your file
# Give the name of SNP of interest e.g. chr4:15737348 in the section specific region to display along with a flanking size of 200kb
# Click on plot data 
# A pdf is generated. Use edit feature in pdf to rename the loci names with SNPids.
# Extract images from pdf as TIF file




