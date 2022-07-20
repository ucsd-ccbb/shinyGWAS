library(data.table)
library(ggplot2)
library(forcats)
library(R.utils)
library(VariantAnnotation)
library(circlize)
library(RCircos)
library(dplyr)
options(stringsAsFactors=FALSE)



#input directories
FUMA_Dir = "../data/FUMA_Results"
GWAS_Catalog_Dir = "../data/GWAS_Catalog"
Summary_Stats_Dir = "../data/Summary_Stats"

#snps file
fuma_snps_df <- fread(file.path(FUMA_Dir, 'snps.txt'), sep='\t')
row.names(fuma_snps_df) <- fuma_snps_df$rsID
fuma_snps_df$negLogP <- -log10(fuma_snps_df$gwasP)

# Reproduce for individual genomic locus (drop down selection)

# eQTL
fuma_eqtl_df <- fread(file.path(FUMA_Dir, 'eqtl.txt'), sep='\t')
#parse & filter
fuma_eqtl <- fuma_eqtl_df %>%
  filter(p<5E-8) %>%
  dplyr::select(c("uniqID", "p", "tissue", "symbol"))
fuma_eqtl <- merge(fuma_snps_df, fuma_eqtl, by="uniqID")
#extract unique tissues for  eqtl tissue igv track
eQTL_tissues <- unique(fuma_eqtl$tissue)

#### absolute code ####
# GWAS Catalog 
# fuma_gwas_cat <- fread("../data/FUMA_Results/gwascatalog.txt", sep="\t")
# fuma_gwas_cat <- fuma_gwas_cat %>%
#   group_by(snp) %>%
#   filter(P==min(P)) %>%
#   dplyr::select(c("snp", "P", "Trait", "PMID", "Context"))
# fuma_gwas_cat <- fuma_gwas_cat[!duplicated(fuma_gwas_cat$snp), ]
# fuma_gwas_cat <- merge(fuma_snps_df, fuma_gwas_cat, by.x="rsID", by.y="snp")

### Tracks###
#to create tracks make sure your dataframe has col1=chr, col2=start, col3=end. Otherwise igvShiny wont read it

#GWAS Catalog used in Manhattan Tab 
gwasCatalog <- fread(file.path(GWAS_Catalog_Dir, "GWAS_catalog_v1.0.2_signif_only_filtered_reordered_renamed.txt"))
gwasCatalog$chr <- as.character(gwasCatalog$chrom) #make sure this col are "character"
gwasCatalog$start <- as.numeric(gwasCatalog$start) #numeric
gwasCatalog$end <- as.numeric(gwasCatalog$end) #numeric

#RBD scores
RDB_score <- fuma_snps_df %>%
  dplyr::select(c("chr", "pos", "RDB"))
RDB_score$chr <- as.character(RDB_score$chr) #make sure this col are "character"
RDB_score$start <- as.numeric(RDB_score$pos) #numeric
RDB_score$end <- RDB_score$start + 1 #numeric
# reorder the columns using select
RDB_score <- RDB_score %>%
  dplyr::select(c("chr", "start", "end", "RDB"))


#CADD scores
CADD_scores_df <- fuma_snps_df %>%
  dplyr::select(c("chr", "pos", "CADD"))
CADD_scores_df$chr <- as.character(CADD_scores_df$chr) #make sure this col are "character"
CADD_scores_df$start <- as.numeric(CADD_scores_df$pos) #numeric
CADD_scores_df$end <- CADD_scores_df$start + 1 #numeric
# reorder the columns using select
CADD_scores_df <- CADD_scores_df %>%
  dplyr::select(c("chr", "start", "end", "CADD"))

# Genes per locus
fuma_genes <- fread(file.path(FUMA_Dir, "genes.txt"), sep="\t")
fuma_genes <- fuma_genes[!duplicated(fuma_genes$symbol), ]

###  Manhattan #####
gwas <- fread(file.path(Summary_Stats_Dir, "cd-meta.txt.gz"), sep="\t", select = c(1,2,3,5)) %>% 
  rename(P = 'SCAN-P') %>%
  filter(P < 0.01)

mypalette <- c("black", "gray")
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line
# merge fuma df to add CADD score, and RDB to gwas df
subset_fuma_snps_df <- subset(fuma_snps_df, select=c(rsID, chr, pos, gwasP, func, CADD, RDB))
colnames(subset_fuma_snps_df) <-c("SNP", "CHR", "POS", "P", "func", "CADD", "RDB")
new_gwas <- merge(gwas, subset_fuma_snps_df,by=c("CHR", "SNP","POS","P"),  all.x=TRUE)
#add GTEx links for eQTLs
url = "https://gtexportal.org/home/snp/"
new_gwas$GTEX <- paste0(url, new_gwas$SNP)

# Circos
l <- list.files(FUMA_Dir)
if(length(grep("GenomicRiskLoci.txt|snps.txt|genes.txt|ci.txt|eqtl.txt", l)) < 5) {
  stop("Missing Input Files!")
}

# Load GenomicRiskLoci data
loci_all <- as.data.frame(fread(file.path(FUMA_Dir, "GenomicRiskLoci.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE))
loci_all_filt <- loci_all[,c(1,3,4,5,7,8)] 
loci_all_filt[,3] <- paste0("chr", loci_all_filt[,3])
rm(loci_all)

# Load snps data
snps_all <- as.data.frame(fread("../data/FUMA_Results/snps.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
snps_all_filt <- cbind.data.frame(snps_all[,c(3,4,4,8)], snps_all[["r2"]])
snps_all_filt <- snps_all_filt[!is.na(snps_all_filt[,4]),]
rm(snps_all)

# Load genes data
gns_all <- as.data.frame(fread(file="../data/FUMA_Results/genes.txt", sep="\t", header=TRUE))
tmp_loci <- unlist(lapply(strsplit(gns_all[["GenomicLocus"]], ":"), function(x) x[[1]]))
gns_all_filt <- cbind.data.frame(gns_all[,c(3,4,5,2)], Loci=tmp_loci, ENS=as.character(gns_all[,1]))
gns_all_filt[,1] <- paste0("chr", gns_all_filt[,1])
gns_all_filt_inp <- gns_all_filt[,-5]
gns_all_filt_inp$ENS <- as.character(gns_all_filt_inp$ENS)
gns_all_filt_inp$V5 <- "id=0"
gns_all_filt_inp$V5[which(gns_all$ciMap == "Yes")] <- "id=1"
if(sum(colnames(gns_all) %in% "eqtlMapSNPs") > 0) {
  gns_all_filt_inp$V5[which(gns_all$eqtlMapSNPs > 0 & gns_all$ciMap == "Yes")] <- "id=2"
}
rm(gns_all)

# Load ci data
ci_all <- as.data.frame(fread(file="../data/FUMA_Results/ci.txt", sep="\t", header=TRUE))
if(ncol(ci_all) >= 11) {
  ci_all <- ci_all[ci_all[,11] == 1,]
}
ci_all_filt <- ci_all[ci_all[,8] == "intra",]
rm(ci_all)

# Load eqtl data
eqtl_all <- as.data.frame(fread(file="../data/FUMA_Results/eqtl.txt", sep="\t", header=TRUE))
if(ncol(eqtl_all) >= 14) {
  eqtl_all <- eqtl_all[eqtl_all[,14] == 1,]
}
eqtl_all_filt <- eqtl_all[,c(4,6,11,12)]
eqtl_all_filt_gns <- eqtl_all_filt[eqtl_all_filt[,1] %in% gns_all_filt$ENS,]
rm(eqtl_all)

# Process SNPs
snps_all_filt[,3] <- snps_all_filt[,3]+1
snps_all_filt$logP <- -log(snps_all_filt[,4]+0.0000,10)
snps_all_filt$color <- sapply(snps_all_filt[,5], function(x) {
  if(x >= 0.8) {"red"} else if(x >= 0.6) {"dark orange"} else if(x >= 0.4) {"forest green"} else if(x >= 0.2) {
    "blue"} else {"grey30"}
})
snps_all_split <- lapply(split(snps_all_filt, snps_all_filt$chr), function(x) {
  xx <- x[order(x[,3],decreasing = FALSE),]
  xx[,1] <- paste0("chr", xx[,1])
  if(nrow(xx) > 15000) {
    return(xx[1:15000,])
  } else {
    return(xx)
  }
})
names(snps_all_split) <- paste0("chr", names(snps_all_split))

## Process chromatin interactions
ci_all_filt <- ci_all_filt[order(ci_all_filt[,1]),]
ci_all_filt_pos <- cbind.data.frame(ci_all_filt[,1], do.call(rbind.data.frame,
                                                             lapply(strsplit(ci_all_filt[,2], split =":"), function(x) {
                                                               x1 <- paste0("chr", x[[1]])
                                                               c(x1, strsplit(x[[2]], split="-")[[1]])
                                                             })),
                                    do.call(rbind.data.frame, lapply(strsplit(ci_all_filt[,3], split =":"), function(x) {
                                      x1 <- paste0("chr", x[[1]])
                                      c(x1, strsplit(x[[2]], split="-")[[1]])
                                    })), ci_all_filt[,4])
ci_all_filt_pos <- cbind.data.frame(ci_all_filt_pos[,1:2],
                                    apply(ci_all_filt_pos[,3:4],2, as.numeric),
                                    ci_all_filt_pos[,5],
                                    apply(ci_all_filt_pos[,3:4],2, as.numeric),
                                    rep("dark orange", nrow(ci_all_filt_pos)))
colnames(ci_all_filt_pos) <- paste0("V", seq(1:ncol(ci_all_filt_pos)))
ci_all_filt_pos <- ci_all_filt_pos[order(ci_all_filt_pos[,8]),]
ci_all_split <- lapply(split(ci_all_filt_pos, ci_all_filt_pos[,2]), function(x) {
  if(nrow(x) > 10000) {
    return(x[1:10000,])
  } else {
    return(x)
  }
})


## Process eqtls

eqtl_all_filt_gns_inp <- eqtl_all_filt_gns[,c(3,4,4,3,4,4,1,1,2)]
eqtl_all_filt_gns_inp[,3] <- as.numeric(eqtl_all_filt_gns_inp[,2]) + 1
eqtl_all_filt_gns_inp[,5] <- gns_all_filt[,2][match(eqtl_all_filt_gns_inp[,8], gns_all_filt[,6])]
eqtl_all_filt_gns_inp[,6] <- gns_all_filt[,3][match(eqtl_all_filt_gns_inp[,8], gns_all_filt[,6])]
eqtl_all_filt_gns_inp[,7] <- "forest green"
colnames(eqtl_all_filt_gns) <- paste0("V", seq(1:ncol(eqtl_all_filt_gns)))
eqtl_all_split <- lapply(split(eqtl_all_filt_gns_inp, eqtl_all_filt_gns_inp[,1]), function(x) {
  xx <- x[order(x[,9],decreasing = FALSE),]
  xx[,1] <- xx[,4] <- paste0("chr", xx[,1])
  if(nrow(xx) > 100000) {
    return(xx[1:100000,-c(8,9)])
  } else {
    return(xx[,-c(8,9)])
  }
})
names(eqtl_all_split) <- paste0("chr", names(eqtl_all_split))

## Generate regions
loci_all_split <- split(loci_all_filt, loci_all_filt[,3])
chrFileList <- lapply(names(snps_all_split), function(c) {
  if(nrow(snps_all_split[[c]] > 0)) {
    generateChromosomeFiles(c, loci_all_split[[c]], ci_all_split[[c]], snps_all_split[[c]],
                            gns_all_filt[gns_all_filt[,1] %in% c,])#c, loci, ci, snp, gene
  }
})
names(chrFileList) <- names(snps_all_split)

## Create circos input based on regions

gns_all_filt_inp$col <- "dark orange"
gns_all_filt_inp$col[as.character(gns_all_filt_inp[,5]) == "id=2"] <- "red"
gns_all_filt_inp$col[as.character(gns_all_filt_inp[,5]) == "id=0"] <- "forest green"

gns_all_filt_inp2 <- processSubChr(gns_all_filt_inp, chrFileList)
gns_all_filt_inp2 <- gns_all_filt_inp2[!duplicated(gns_all_filt_inp2),-5]
gns_all_filt_inp2 <- gns_all_filt_inp2[grep("_", gns_all_filt_inp2[,1]),]

snps_all_split_inp <- lapply(1:length(snps_all_split), function(i) {
  processSubChr(snps_all_split[[i]], chrFileList, chr=names(snps_all_split)[i])[,-c(4,5)]
})
names(snps_all_split_inp) <- names(snps_all_split)

ci_all_split_inp <- lapply(1:length(ci_all_split), function(i) {
  processSubChr(ci_all_split[[i]][-c(1)], chrFileList, chr=names(ci_all_split)[i], link=TRUE)
})
names(ci_all_split_inp) <- names(ci_all_split)

eqtl_all_split_inp <- lapply(1:length(eqtl_all_split), function(i) {
  processSubChr(eqtl_all_split[[i]], chrFileList, chr=names(eqtl_all_split)[i], link=TRUE)
})
names(eqtl_all_split_inp) <- names(eqtl_all_split)

## Generate circos for each chromosome

dataList <- list(gns=gns_all_filt_inp2, snps=snps_all_split_inp, chrFileList=chrFileList, ci=ci_all_split_inp, eqtl=eqtl_all_split_inp)
