library(dplyr)
library(tidyr)


#parse data
source("scripts/Functions.R")
getwd()
#input directories
inputdir = "/Users/dchilinfuentes/CCBB_projects/shinyGWAS_Repo copy/CCBB_shinyGWAS_Prototype/app/data/FUMA_Results"

##snps
FUMA_SNPs_df <- as.data.frame(loadInputFile(inputdir, "snps.txt")) %>%
  dplyr::select(c('uniqID','chr','pos','rsID','non_effect_allele','effect_allele','MAF','gwasP','GenomicLocus','nearestGene',
                  'dist','func','CADD','RDB','minChrState','commonChrState')) 
row.names(FUMA_SNPs_df) <- FUMA_SNPs_df$rsID
FUMA_SNPs_df$negLogP <- -log10(FUMA_SNPs_df$gwasP)


#eqt
fuma_eqtl_df <- as.data.frame(loadInputFile(inputdir, "eqtl.txt"))
fuma_eqtl_df <- fuma_eqtl_df[order(fuma_eqtl_df$uniqID, decreasing=TRUE),]
fuma_eqtl_df <- fuma_eqtl_df[!duplicated(fuma_eqtl_df$uniqID),]
row.names(fuma_eqtl_df) <- fuma_eqtl_df$uniqID

#gwas cat
fuma_gwas_cat <- as.data.frame(loadInputFile(inputdir, "gwascatalog.txt"))
fuma_gwas_cat <- fuma_gwas_cat[order(fuma_gwas_cat$snp, decreasing=TRUE),]
fuma_gwas_cat <- fuma_gwas_cat[!duplicated(fuma_gwas_cat$snp),]
row.names(fuma_gwas_cat) <- fuma_gwas_cat$snp

FUMA_annot_df = data.frame()

for (focal_loc in unique(FUMA_SNPs_df$GenomicLocus)){ # loop over genomic loci
  # focal_loc = 13
  print(focal_loc)
  focal_df <- subset(FUMA_SNPs_df, GenomicLocus==focal_loc)
  print(dim(focal_df))
  # # check if there are any eQTLs
  focal_fuma_id=focal_df$uniqID
  focal_fuma_id
  index <- which(focal_fuma_id %in% fuma_eqtl_df$uniqID)
  intersecting_ids <- focal_df$uniqID[index]
  intersecting_ids
  if (length(intersecting_ids)>0){
    row.names(focal_df) = focal_fuma_id
    focal_eqtl <- fuma_eqtl_df[intersecting_ids,]
    colnames(focal_eqtl)
    #sort by pval
    focal_eqtl <- focal_eqtl[order(focal_eqtl$p),]
    # focal_eqtl <- focal_eqtl[!duplicated(focal_eqtl$snp),]
    
    #reshape focal df
    focal_df_intIDs <- focal_df[rownames(focal_df) %in% intersecting_ids, ]  # Extract rows from data
    #add columns to focal df
    focal_df_intIDs$top_eQTL_pval <- focal_eqtl$p
    focal_df_intIDs$top_eQTL_tissue <- focal_eqtl$tissue
    focal_df_intIDs$top_eQTL_gene <- focal_eqtl$symbol
    focal_df_intIDs$top_gwas_cat_pval<- NA
    focal_df_intIDs$top_gwas_cat_trait<-NA
    focal_df_intIDs$top_gwas_cat_PMID<-NA
    focal_df_intIDs$top_gwas_cat_context <-NA
    #add NAs to the remainder of that focal df
    focal_df_noIntIDs <- focal_df[!rownames(focal_df) %in% intersecting_ids, ]  # Extract rows from data
    # focal_df_noIntIDs <- focal_df_noIntIDs[!duplicated(focal_df_noIntIDs$snp),]
    
    focal_df_noIntIDs$top_eQTL_pval<- NA
    focal_df_noIntIDs$top_eQTL_tissue<-NA
    focal_df_noIntIDs$top_eQTL_gene<-NA   
    focal_df_noIntIDs$top_gwas_cat_pval<- NA
    focal_df_noIntIDs$top_gwas_cat_trait<-NA
    focal_df_noIntIDs$top_gwas_cat_PMID<-NA
    focal_df_noIntIDs$top_gwas_cat_context <-NA
    focal_df <- rbind(focal_df_intIDs, focal_df_noIntIDs) #merge
  }else{
    focal_df$top_eQTL_pval<- NA
    focal_df$top_eQTL_tissue<-NA
    focal_df$top_eQTL_gene<-NA
    focal_df$top_gwas_cat_pval<- NA
    focal_df$top_gwas_cat_trait<-NA
    focal_df$top_gwas_cat_PMID<-NA
    focal_df$top_gwas_cat_context <-NA
    
  }
  # # check if there are any prior hits in GWAS catalog
  focal_fuma_id=focal_df$rsID
  index <- which(focal_df$rsID %in% fuma_gwas_cat$snp)
  intersecting_ids <- focal_df$rsID[index]
  intersecting_ids
  if (length(intersecting_ids)>0){
    row.names(focal_df) <- focal_df$rsID
    focal_gwas_cat <- fuma_gwas_cat[intersecting_ids,]
    #sort
    focal_gwas_cat <- focal_gwas_cat[order(focal_gwas_cat$P),]
    #reshape focal df
    focal_df_intIDs_GWAS <- focal_df[rownames(focal_df) %in% intersecting_ids, ]  # Extract rows from data
    focal_df_intIDs_GWAS$top_eQTL_pval<- NA
    focal_df_intIDs_GWAS$top_eQTL_tissue<-NA
    focal_df_intIDs_GWAS$top_eQTL_gene<-NA
    # focal_df_intIDs <- rbind(focal_df_intIDs, focal_df_intIDs_GWAS)
    #add columns to focal df
    focal_df_intIDs_GWAS$top_gwas_cat_pval <-  focal_gwas_cat$P
    focal_df_intIDs_GWAS$top_gwas_cat_trait <-  focal_gwas_cat$Trait
    focal_df_intIDs_GWAS$top_gwas_cat_PMID <-  focal_gwas_cat$PMID
    focal_df_intIDs_GWAS$top_gwas_cat_context <-  focal_gwas_cat$Context
    #add NAs to the remainder of that focal df
    focal_df_noIntIDs_GWAS <- focal_df[!rownames(focal_df) %in% intersecting_ids, ]  # Extract rows from data
    if (dim(focal_df_noIntIDs_GWAS)[1] != 0) {
      # focal_df_intIDs <- rbind(focal_df_intIDs, focal_df_noIntIDs_GWAS)
      focal_df_noIntIDs_GWAS$top_gwas_cat_pval<- NA
      focal_df_noIntIDs_GWAS$top_gwas_cat_trait<-NA
      focal_df_noIntIDs_GWAS$top_gwas_cat_PMID<-NA
      focal_df_noIntIDs_GWAS$top_gwas_cat_context <-NA
    }
    
    focal_df_gwas <- rbind(focal_df_intIDs_GWAS, focal_df_noIntIDs_GWAS)
    
    focal_df <- rbind(focal_df, focal_df_gwas)
    
  }else{
    focal_df$top_gwas_cat_pval<- NA
    focal_df$top_gwas_cat_trait<-NA
    focal_df$top_gwas_cat_PMID<-NA
    focal_df$top_gwas_cat_context <-NA
    
  }
  FUMA_annot_df <- rbind(FUMA_annot_df, focal_df)
  
}
#PLOT
FUMA_annot_plot <- FUMA_annot_df %>%
  dplyr::select(c(GenomicLocus, top_gwas_cat_trait))
test <- na.omit(FUMA_annot_plot)

#plot
gwasBarplot <- ggplot(test, aes(x=factor(top_gwas_cat_trait), fill=factor(GenomicLocus))) +
  geom_bar(position="stack") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(legend.text = element_text( size = 7))+
  # ylim(0, 20) +
  xlab("Top GWAS Catalog Trait") +
  labs(fill = "Genomic Locus") +
  scale_x_discrete(guide = guide_axis(angle = 90)) 
  
  # theme(axis.text.x = element_text(size=10, angle=90))
# gwasBarplot

#### plot the other way 
#create bar plot
# FUMA_annot_plot <- FUMA_annot_df %>%
#   dplyr::select(c(GenomicLocus, top_gwas_cat_trait))
# test <- na.omit(FUMA_annot_plot)
# test$values <-1
# 
# test<- test %>%
#   group_by(top_gwas_cat_trait, GenomicLocus) %>%
#   # filter(Capture != "") %>%  # filter out captured ones (handling)
#   summarise(Count = n())   #get the count for each fish type (long format)
#   #spread(top_gwas_cat_trait, Count) #%>%# Use the spread() function from tidyr package to convert the data from long to wide format
#   # dplyr::select(test , -c(GenomicLocus))
# 
# test <- as.data.frame(t(test))
# colnames(test) <- test[1,] 
# test <- test[-1,]
# 
# 
# test[] <- lapply(test, function(x) as.numeric(as.character(x)))
# sapply(test, class)
# 
# 
# test <- test[rowSums(test[], na.rm=TRUE)>1,]
# 
# traits <- rownames(test)
# test$top_gwas_cat_trait <- traits
# locuses <- colnames(test)
# 
# 
# gwasBarplot <- ggplot(test, aes(factor(locuses), fill=c(top_gwas_cat_trait))) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(legend.position=c(0.5,1)) +
#   theme(legend.text = element_text( size = 7))+
#   theme(axis.text.x = element_text(size=10, angle=450))
# 
# gwasBarplot



