library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(R.utils)

setwd("/Users/adammark/projects/shiny/shinyGWAS")
fuma_snps_df = fread("FUMA_ASD_job58887/snps.txt", sep="\t")
row.names(fuma_snps_df) <- fuma_snps_df$rsID
fuma_snps_df$negLogP <- -log10(fuma_snps_df$gwasP)

# ggplot(fuma_snps_df, aes(factor(GenomicLocus), fill=func)) + geom_bar(position="stack") + theme_classic() + xlab("Genomic Loci") + theme(axis.text.x = element_text(size=10, angle=90))
# ggplot(fuma_snps_df, aes(func, CADD, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(face="bold", size=14, angle=90))
# ggplot(fuma_snps_df, aes(func, minChrState, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(face="bold", size=14, angle=90))
# ggplot(fuma_snps_df, aes(func, negLogP, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(face="bold", size=14, angle=90))
# Reproduce for individual genomic locus (drop down selection)

# eQTL
fuma_eqtl <- fread("FUMA_ASD_job58887/eqtl.txt", sep="\t")
fuma_eqtl <- fuma_eqtl %>%
  group_by(uniqID) %>%
  filter(p==min(p)) %>%
  select(c("uniqID", "p", "tissue", "symbol"))
fuma_eqtl <- fuma_eqtl[!duplicated(fuma_eqtl$uniqID), ]
fuma_eqtl <- merge(fuma_snps_df, fuma_eqtl, by="uniqID")

# GWAS Catalog
fuma_gwas_cat <- fread("FUMA_ASD_job58887/gwascatalog.txt", sep="\t")
fuma_gwas_cat <- fuma_gwas_cat %>%
  group_by(snp) %>%
  filter(P==min(P)) %>%
  select(c("snp", "P", "Trait", "PMID", "Context"))
fuma_gwas_cat <- fuma_gwas_cat[!duplicated(fuma_gwas_cat$snp), ]
fuma_gwas_cat <- merge(fuma_snps_df, fuma_gwas_cat, by.x="rsID", by.y="snp")
                       
# ggplot(fuma_gwas_cat, aes(factor(GenomicLocus), fill=Trait)) + geom_bar(position="stack") + theme_classic() + xlab("Genomic Loci")

# Genes per locus
fuma_genes <- fread("FUMA_ASD_job58887/genes.txt", sep="\t")
fuma_genes <- fuma_genes[!duplicated(fuma_genes$symbol), ]
# ggplot(fuma_genes, aes(x=pLI, y=ncRVIS)) + geom_point() + theme_classic()
# ggplot(fuma_genes, aes(fct_infreq(GenomicLocus))) + geom_bar() + xlab("Genomic Locus") + ylab("genes mapped to locus") + theme_classic()

# Manhattan
gwas <- fread("iPSYCH-PGC_ASD_Nov2017.gz", sep="\t", select = c(1,2,3,9))# %>% filter(P < 0.05)
mypalette <- c("black", "gray")
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

