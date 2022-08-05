library(data.table)
library(ggplot2)
library(forcats)
library(R.utils)
library(VariantAnnotation)
library(circlize)
library(RCircos)
library(dplyr)
options(stringsAsFactors=FALSE)

getwd()

source("scripts/Functions.R")


#input directories
FUMA_Dir = "data/FUMA_Results"
GWAS_Catalog_Dir = "data/GWAS_Catalog"
Summary_Stats_Dir = "data/Summary_Stats"



####################################################################
############## LOAD FILES HERE #####################################
####################################################################
####### FUMA FILES ####### 
# SNPS 
print("Loading snps.txt File")
fuma_snps_df <- loadInputFile(FUMA_Dir, "snps.txt")
row.names(fuma_snps_df) <- fuma_snps_df$rsID
fuma_snps_df$negLogP <- -log10(fuma_snps_df$gwasP)
print("----- snps.txt loaded -----")

# eQTL
print("Loading eqtl.txt File")
fuma_eqtl_df <- loadInputFile(FUMA_Dir, "eqtl.txt")
print("----- eqtl.txt loaded -----")

# Genes per locus
print("Loading genes.txt File")
fuma_genes_df <- loadInputFile(FUMA_Dir, "genes.txt")
fuma_genes_df <- fuma_genes_df[!duplicated(fuma_genes_df$symbol), ]
print("----- genes.txt Loaded -----")

# Genes per locus
print("Loading GenomicRiskLoci.txt File")
fuma_loci_df <- loadInputFile(FUMA_Dir, "GenomicRiskLoci.txt")
print("----- Genomic Risck Loci file Loaded -----")

# Load ci data
print("Loading ci.txt File")
fuma_ci_df <- as.data.frame(loadInputFile(FUMA_Dir, "ci.txt"))
print("----- ci.txt file Loaded -----")



####### GWAS Catalog ####### 
print("Loading GWAS_catalog_v1.0.2 File")
gwasCatalog <- loadInputFile(GWAS_Catalog_Dir, "GWAS_catalog_v1.0.2_signif_only_filtered_reordered_renamed.txt")
print("----- GWAS Catalog Loaded -----")

####### GWAS Summary ####### 
print("Loading GWAS Summary File")
gwasSummary <- loadInputFile(Summary_Stats_Dir, "assoc_results.LEWYX8.glm-FILTERED-TESTADD.linear", colNames = c("CHR","BP","P")) %>%
  filter(P<0.01)
print("----- GWAS Summary Loaded -----")


##### Circos #####
print("Checking Circos input files")
l <- list.files(FUMA_Dir)
if(length(grep("GenomicRiskLoci.txt|snps.txt|genes.txt|ci.txt|eqtl.txt", l)) < 5) {
  stop("Missing Input Files!")
}
print("All input files exist for Circos plots")
##### END OF LOADING FILES ######


####################################################################
################ PARSE FILES HERE ################################
####################################################################
### MANAHATTAN TAB ####
# eQTL
print("Pasring eQTL file for eQTL Tissue Tracks")
fuma_eqtl <- fuma_eqtl_df %>%
  filter(p<5E-8) %>%
  dplyr::select(c("uniqID", "p", "tissue", "symbol"))
fuma_eqtl <- merge(fuma_snps_df, fuma_eqtl, by="uniqID")
eQTL_tissues <- unique(fuma_eqtl$tissue) #extract unique tissues for  eqtl tissue igv track
print("----- eQTL Tissues ready -----")

#GWAS Catalog 
print("Preparing GWAS Catalog Track")
gwasCatalog <- prepDfToIvgTrack(gwasCatalog, colNames = c('chr', 'start', 'end'), reArrange=FALSE)
print("----- GWAS Catalog ready -----")

#RBD scores (uses snps file)
print("Preparing RDB Score Track")
RDB_score <- prepDfToIvgTrack(fuma_snps_df, colNames=c("chr", "pos", "RDB"), reArrange=TRUE)
RDB_score <- getRDBSupportingInfo(RDB_score)
print("----- RDB Score ready -----")

#CADD scores (uses snps file)
print("Preparing CADD Score Track")
CADD_scores <- prepDfToIvgTrack(fuma_snps_df, colNames=c("chr", "pos", "CADD"), reArrange=TRUE)
print("----- CADD Scores ready -----")


#GWAS manhattan plot (uses snps file)
print("Preparing GWAS Track")
gwas <- gwasSummary %>%
  dplyr::select(c("CHR","BP","P")) %>%
  filter(P<0.01)
mypalette <- c("black", "gray")
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line
# merge fuma df to add CADD score, and RDB to gwas df
subset_fuma_snps_df <- subset(fuma_snps_df, select=c(rsID, chr, pos, gwasP, func, CADD, RDB))
colnames(subset_fuma_snps_df) <-c("SNP", "CHR", "BP", "P", "func", "CADD", "RDB")
new_gwas <- merge(gwas,subset_fuma_snps_df,by=c("CHR", "BP","P"),  all.x=TRUE)
#add GTEx links for eQTLs
url = "https://gtexportal.org/home/snp/"
new_gwas$GTEX <- paste0(url, new_gwas$SNP)
rm(gwasSummary)
rm(gwas)
print("----- GWAS Track ready -----")

####### End of Manhattan Tab ########
####################################
### Circos TAB #### (absolute code. Code has been moved to CircosPlots.R)

# # Genes per locus
# # fuma_genes <- as.data.frame(fuma_genes_df)
# # fuma_genes <- fuma_genes[!duplicated(fuma_genes$symbol), ]
# # print("----- Genes ready -----")
# 
# # Genomic Risk Loci 
# # loci_all <- as.data.frame(fuma_loci_df)
# # loci_all_filt <- loci_all[,c(1,3,4,5,7,8)]
# # loci_all_filt[,3] <- paste0("chr", loci_all_filt[,3])
# # rm(fuma_loci_df)
# # rm(loci_all)
# # print("----- Genomic Risck Loci ready -----")
# # 
# # # filter snps data
# # snps_all <- as.data.frame(fuma_snps_df)
# # snps_all_filt <- cbind.data.frame(snps_all[,c(3,4,4,8)], snps_all[["r2"]])
# # snps_all_filt <- snps_all_filt[!is.na(snps_all_filt[,4]),]
# # rm(snps_all)
# # print("----- SNPs filtred and ready -----")
# 
# 
# 
# # filter genes data
# gns_all <- as.data.frame(fuma_genes_df)
# gns_all[,3] <- paste0("chr", gns_all[,3])
# 
# filtered_index <- which(gns_all$chr %in% unique(loci_all_filt$chr)) #filter out chr that are not found in the genomic risk loci file
# gns_all_filt_chr <- gns_all[filtered_index,]
# gns_all_filt_chr$GenomicLocus<-as.character(gns_all_filt_chr$GenomicLocus)
# tmp_loci <- unlist(lapply(strsplit(gns_all_filt_chr[["GenomicLocus"]], ":"), function(x) x[[1]]))
# gns_all_filt <- cbind.data.frame(gns_all_filt_chr[,c(3,4,5,2)], Loci=tmp_loci, ENS=as.character(gns_all_filt_chr[,1]))
# # gns_all_filt[,1] <- paste0("chr", gns_all_filt[,1])
# gns_all_filt_inp <- gns_all_filt[,-5]
# gns_all_filt_inp$ENS <- as.character(gns_all_filt_inp$ENS)
# gns_all_filt_inp$V5 <- "id=0"
# gns_all_filt_inp$V5[which(gns_all$ciMap == "Yes")] <- "id=1"
# if(sum(colnames(gns_all) %in% "eqtlMapSNPs") > 0) {
#   gns_all_filt_inp$V5[which(gns_all$eqtlMapSNPs > 0 & gns_all$ciMap == "Yes")] <- "id=2"
# }
# rm(fuma_genes_df)
# rm(gns_all)
# rm(gns_all_filt_chr)
# rm(filtered_index)
# print("----- Genes File  filtered  -----")
# 
# 
# # Load ci data
# ci_all <- as.data.frame(fuma_ci_df)
# if(ncol(ci_all) >= 11) {
#   ci_all <- ci_all[ci_all[,11] == 1,]
# }
# ci_all_filt <- ci_all[ci_all[,8] == "intra",]
# # rm(fuma_ci_df)
# rm(ci_all)
# print("----- ci File  Loaded  -----")
# 
# 
# # Load eqtl data
# # eqtl_all <- as.data.frame(loadInputFile(FUMA_Dir, "eqtl.txt"))
# eqtl_all <- as.data.frame(fuma_eqtl_df)
# if(ncol(eqtl_all) >= 14) {
#   eqtl_all <- eqtl_all[eqtl_all[,14] == 1,]
# }
# eqtl_all_filt <- eqtl_all[,c(4,6,11,12)]
# eqtl_all_filt_gns <- eqtl_all_filt[eqtl_all_filt[,1] %in% gns_all_filt$ENS,]
# rm(fuma_eqtl_df) #no longer needed
# rm(eqtl_all)
# print("----- eqtl File  loaded and filtered  -----")
# 
# 
# # Process SNPs
# print("Processing and filtering snps....")
# snps_all_filt[,3] <- snps_all_filt[,3]+1
# snps_all_filt$logP <- -log(snps_all_filt[,4]+0.0000,10)
# snps_all_filt$color <- sapply(snps_all_filt[,5], function(x) {
#   if(x >= 0.8) {"red"} else if(x >= 0.6) {"dark orange"} else if(x >= 0.4) {"forest green"} else if(x >= 0.2) {
#     "blue"} else {"grey30"}
# })
# print("....")
# snps_all_split <- lapply(split(snps_all_filt, snps_all_filt$chr), function(x) {
#   xx <- x[order(x[,3],decreasing = FALSE),]
#   xx[,1] <- paste0("chr", xx[,1])
#   if(nrow(xx) > 15000) {
#     return(xx[1:15000,])
#   } else {
#     return(xx)
#   }
# })
# names(snps_all_split) <- paste0("chr", names(snps_all_split))
# print("...processing and filtering for SNPs completed")
# 
# 
# ## Process chromatin interactions
# print("Processing and filtering crhomatin interactions")
# ci_all_filt <- ci_all_filt[order(ci_all_filt[,1]),]
# ci_all_filt_pos <- cbind.data.frame(ci_all_filt[,1], do.call(rbind.data.frame,
#                                                              lapply(strsplit(ci_all_filt[,2], split =":"), function(x) {
#                                                                x1 <- paste0("chr", x[[1]])
#                                                                c(x1, strsplit(x[[2]], split="-")[[1]])
#                                                              })),
#                                     do.call(rbind.data.frame, lapply(strsplit(ci_all_filt[,3], split =":"), function(x) {
#                                       x1 <- paste0("chr", x[[1]])
#                                       c(x1, strsplit(x[[2]], split="-")[[1]])
#                                     })), ci_all_filt[,4])
# ci_all_filt_pos <- cbind.data.frame(ci_all_filt_pos[,1:2],
#                                     apply(ci_all_filt_pos[,3:4],2, as.numeric),
#                                     ci_all_filt_pos[,5],
#                                     apply(ci_all_filt_pos[,3:4],2, as.numeric),
#                                     rep("dark orange", nrow(ci_all_filt_pos)))
# print("...")
# colnames(ci_all_filt_pos) <- paste0("V", seq(1:ncol(ci_all_filt_pos)))
# ci_all_filt_pos <- ci_all_filt_pos[order(ci_all_filt_pos[,8]),]
# ci_all_split <- lapply(split(ci_all_filt_pos, ci_all_filt_pos[,2]), function(x) {
#   if(nrow(x) > 10000) {
#     return(x[1:10000,])
#   } else {
#     return(x)
#   }
# })
# print("... chromatin processing and filtering completed")
# 
# 
# ## Process eqtls
# print("Starting eqtl processing ....")
# eqtl_all_filt_gns_inp <- eqtl_all_filt_gns[,c(3,4,4,3,4,4,1,1,2)]
# eqtl_all_filt_gns_inp[,3] <- as.numeric(eqtl_all_filt_gns_inp[,2]) + 1
# eqtl_all_filt_gns_inp[,5] <- gns_all_filt[,2][match(eqtl_all_filt_gns_inp[,8], gns_all_filt[,6])]
# eqtl_all_filt_gns_inp[,6] <- gns_all_filt[,3][match(eqtl_all_filt_gns_inp[,8], gns_all_filt[,6])]
# eqtl_all_filt_gns_inp[,7] <- "forest green"
# colnames(eqtl_all_filt_gns) <- paste0("V", seq(1:ncol(eqtl_all_filt_gns)))
# eqtl_all_split <- lapply(split(eqtl_all_filt_gns_inp, eqtl_all_filt_gns_inp[,1]), function(x) {
#   xx <- x[order(x[,9],decreasing = FALSE),]
#   xx[,1] <- xx[,4] <- paste0("chr", xx[,1])
#   if(nrow(xx) > 100000) {
#     return(xx[1:100000,-c(8,9)])
#   } else {
#     return(xx[,-c(8,9)])
#   }
# })
# names(eqtl_all_split) <- paste0("chr", names(eqtl_all_split))
# print("... eqtl processing completed")
# 
# # Generate regions
# print("Genearting regions ...")
# loci_all_split <- split(loci_all_filt, loci_all_filt[,3])
# chrFileList <- lapply(names(snps_all_split), function(c) {
#   if(nrow(snps_all_split[[c]] > 0)) {
#     generateChromosomeFiles(c, loci_all_split[[c]], ci_all_split[[c]], snps_all_split[[c]],
#                             gns_all_filt[gns_all_filt[,1] %in% c,])#c, loci, ci, snp, gene
#   }
# })
# names(chrFileList) <- names(snps_all_split)
# print("... Regions generated completed")
# 
# # Create circos input based on regions
# print("creating circos input based on regions")
# gns_all_filt_inp$col <- "dark orange"
# gns_all_filt_inp$col[as.character(gns_all_filt_inp[,5]) == "id=2"] <- "red"
# gns_all_filt_inp$col[as.character(gns_all_filt_inp[,5]) == "id=0"] <- "forest green"
# 
# gns_all_filt_inp2 <- processSubChr(gns_all_filt_inp, chrFileList)
# gns_all_filt_inp2 <- gns_all_filt_inp2[!duplicated(gns_all_filt_inp2),-5]
# gns_all_filt_inp2 <- gns_all_filt_inp2[grep("_", gns_all_filt_inp2[,1]),]
# 
# snps_all_split_inp <- lapply(1:length(snps_all_split), function(i) {
#   processSubChr(snps_all_split[[i]], chrFileList, chr=names(snps_all_split)[i])[,-c(4,5)]
# })
# names(snps_all_split_inp) <- names(snps_all_split)
# 
# ci_all_split_inp <- lapply(1:length(ci_all_split), function(i) {
#   processSubChr(ci_all_split[[i]][-c(1)], chrFileList, chr=names(ci_all_split)[i], link=TRUE)
# })
# names(ci_all_split_inp) <- names(ci_all_split)
# 
# eqtl_all_split_inp <- lapply(1:length(eqtl_all_split), function(i) {
#   processSubChr(eqtl_all_split[[i]], chrFileList, chr=names(eqtl_all_split)[i], link=TRUE)
# })
# names(eqtl_all_split_inp) <- names(eqtl_all_split)
# 
# ## Generate circos for each chromosome
# print('Generating circos for each chromosome')
# dataList <- list(gns=gns_all_filt_inp2, snps=snps_all_split_inp, chrFileList=chrFileList, ci=ci_all_split_inp, eqtl=eqtl_all_split_inp)



print("--- data_prep.R loaded ---")