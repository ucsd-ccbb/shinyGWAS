library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(R.utils)

getwd()
# setwd("/Users/dchilinfuentes/CCBB_projects/shinyGWAS-main")
fuma_snps_df = fread("../FUMA_ASD_job58887/snps.txt", sep="\t")
head(fuma_snps_df)
row.names(fuma_snps_df) <- fuma_snps_df$rsID
fuma_snps_df$negLogP <- -log10(fuma_snps_df$gwasP)

# ggplot(fuma_snps_df, aes(factor(GenomicLocus), fill=func)) + geom_bar(position="stack") + theme_classic() + xlab("Genomic Loci") + theme(axis.text.x = element_text(size=10, angle=90))
# ggplot(fuma_snps_df, aes(func, CADD, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(face="bold", size=14, angle=90))
# ggplot(fuma_snps_df, aes(func, minChrState, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(face="bold", size=14, angle=90))
# ggplot(fuma_snps_df, aes(func, negLogP, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(face="bold", size=14, angle=90))
# Reproduce for individual genomic locus (drop down selection)

# eQTL
fuma_eqtl <- fread("../FUMA_ASD_job58887/eqtl.txt", sep="\t")
fuma_eqtl <- fuma_eqtl %>%
  group_by(uniqID) %>%
  filter(p==min(p)) %>%
  select(c("uniqID", "p", "tissue", "symbol"))
fuma_eqtl <- fuma_eqtl[!duplicated(fuma_eqtl$uniqID), ]
fuma_eqtl <- merge(fuma_snps_df, fuma_eqtl, by="uniqID")
#extract unique tissues for  eqtl tissue igv track
eQTL_tissues <- unique(fuma_eqtl$tissue)



# GWAS Catalog
fuma_gwas_cat <- fread("../FUMA_ASD_job58887/gwascatalog.txt", sep="\t")
fuma_gwas_cat <- fuma_gwas_cat %>%
  group_by(snp) %>%
  filter(P==min(P)) %>%
  select(c("snp", "P", "Trait", "PMID", "Context"))
fuma_gwas_cat <- fuma_gwas_cat[!duplicated(fuma_gwas_cat$snp), ]
fuma_gwas_cat <- merge(fuma_snps_df, fuma_gwas_cat, by.x="rsID", by.y="snp")


#GWAS Catalog used in Manhattan Tab 
gwasCatalog <- fread("../data/GWAS_catalog_v1.0.2_signif_only_filtered_reordered_renamed.txt")
#make sure col1 = chr, col2=start, and col3=end
gwasCatalog$chr <- as.character(gwasCatalog$chr) #make sure this col are "character"
gwasCatalog$start <- as.numeric(gwasCatalog$start) #numeric
gwasCatalog$end <- as.numeric(gwasCatalog$end) #numeric



                       
# ggplot(fuma_gwas_cat, aes(factor(GenomicLocus), fill=Trait)) + geom_bar(position="stack") + theme_classic() + xlab("Genomic Loci")

# Genes per locus
fuma_genes <- fread("../FUMA_ASD_job58887/genes.txt", sep="\t")
fuma_genes <- fuma_genes[!duplicated(fuma_genes$symbol), ]
# ggplot(fuma_genes, aes(x=pLI, y=ncRVIS)) + geom_point() + theme_classic()
# ggplot(fuma_genes, aes(fct_infreq(GenomicLocus))) + geom_bar() + xlab("Genomic Locus") + ylab("genes mapped to locus") + theme_classic()

###  Manhattan #####
gwas <- fread("../data/iPSYCH-PGC_ASD_Nov2017.gz", sep="\t", select = c(1,2,3,9)) %>% filter(P < 0.01)
mypalette <- c("black", "gray")
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line
# merge fuma df to add CADD score, and RDB to gwas df
subset_fuma_snps_df <- subset(fuma_snps_df, select=c(rsID, chr, pos, gwasP, func, CADD, RDB))
colnames(subset_fuma_snps_df) <-c("SNP", "CHR", "BP", "P", "func", "CADD", "RDB")
new_gwas <- merge(gwas,subset_fuma_snps_df,by=c("CHR", "SNP","BP","P"),  all.x=TRUE)
#add GTEx links for eQTLs
url = "https://gtexportal.org/home/snp/"
new_gwas$GTEX <- paste0(url, new_gwas$SNP)

