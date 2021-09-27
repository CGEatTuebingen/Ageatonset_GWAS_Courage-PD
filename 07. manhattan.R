#Manhattan plot

## A. Load relevant libraries
library(tidyverse)
library(data.table)
library("reshape2")
library("ggrepel")

## B. Provide metanalysis output file as input
sumstat <- fread("consoritum_metal_output.tbl")

## C. Annotation of the file with SNP names using HRC dataset
hrc_data <- fread("hrc_annotations.txt")
hrc_data <- hrc_data %>% select(MarkerName, SNP)
sumstat_snpid <- left_join(sumstat, hrc_data) %>% sumstat_snpid %>% rename(P = `P-value`) %>% separate("MarkerName", c("CHR", "BP"), ":")

## D. Filter for relevant df and Isq; A HetDf>18 corresponds to df of 19 or minimum of 20 cohorts 
sumstat_snpid <- sumstat_snpid %>% filter(HetDf>18 & HetISq<=50)
sumstat_manhattan_input <- sumstat_snpid %>% select(SNP, CHR, BP, P) 
sumstat_manhattan_input$CHR <- as.numeric(sumstat_manhattan_input$CHR) 
sumstat_manhattan_input$BP <- as.numeric(sumstat_manhattan_input$BP) 
sumstat_manhattan_input$P <- as.numeric(sumstat_manhattan_input$P) 

## E. Provide list of SNPids and genes for annotation
consorium_genelist_manhattan <- fread("genelist_cosortium.txt")

## F. Fomatting the dataset
gwas <- sumstat_manhattan_input
hits <- consorium_genelist_manhattan
gwas$log10Praw <- -1*log(gwas$P, base = 10)
gwas$log10P <- ifelse(gwas$log10Praw > 40, 40, gwas$log10Praw)
gwas$Plevel <- NA
gwas$Plevel[gwas$P < 5E-08] <- "possible"
gwas$Plevel[gwas$P > 5E-08] <- "likely"
gwasFiltered <- gwas
snpsOfInterest <- hits$SNP
gwasToPlotUnsorted <- merge(gwasFiltered, hits, by = "SNP", all.x = T)
gwasToPlot <- gwasToPlotUnsorted[order(gwasToPlotUnsorted$CHR,gwasToPlotUnsorted$BP),]
plotting <- gwasToPlot %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(gwasToPlot, ., by=c("CHR"="CHR")) %>%
  arrange(ordered(CHR), BP) %>%
  mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(log10P>6, "yes", "no")) 
axisdf <- plotting %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


## G. Plotting
thisManhattan <- ggplot(plotting, aes(x=BPcum, y=log10P)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey60", "grey30"), 22 )) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(limits = c(0, 9) ) +
  geom_point(data=subset(plotting, is_highlight=="yes" & Plevel == "likely"), color = "black", size=2) +
  geom_point(data=subset(plotting, is_highlight=="yes" & Plevel == "possible"), color = "red", size=2) +
  geom_label_repel( data=subset(plotting, is_annotate=="yes"), aes(label=GENE, fill=factor(STATUS)), alpha = 0.5,  size=2) +
  scale_fill_manual(values = c("white", "black")) +
  geom_hline(yintercept=6, linetype="dashed", color = "black") +
  geom_hline(yintercept=7.301, linetype="dashed", color = "red") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text=element_text(size=16,  family="serif")
  ) +
  xlab("Chromosome") +
  ylab("-log10P")

## H. Saving the plot (Time-consuming: 5-15 minutes)
ggsave(plot = thisManhattan, file = "manhattanPlot_consortium_date.png", units = c("in"), dpi = 600, height = 5.5, width = 12)

