source("scripts/CircosFunctions.R")
# setwd("/Users/adammark/projects/shiny/shinyGWAS/FUMA_ASD_job58887/")
# setwd("/Users/artnasamran/Documents/shinyGWAS/FUMA_ASD_job58887/")
# setwd("data/FUMA_Results")

l <- list.files("data/FUMA_Results")
if(length(grep("GenomicRiskLoci.txt|snps.txt|genes.txt|ci.txt|eqtl.txt", l)) < 5) {
  stop("Missing Input Files!")
}

# Load GenomicRiskLoci data
loci_all <- as.data.frame(fread("data/FUMA_Results/GenomicRiskLoci.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
loci_all_filt <- loci_all[,c(1,3,4,5,7,8)] 
loci_all_filt[,3] <- paste0("chr", loci_all_filt[,3])
rm(loci_all)

# Load snps data
snps_all <- as.data.frame(fread("data/FUMA_Results/snps.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
snps_all_filt <- cbind.data.frame(snps_all[,c(3,4,4,8)], snps_all[["r2"]], rsID=snps_all[,2])
snps_all_filt <- snps_all_filt[!is.na(snps_all_filt[,4]),]
rm(snps_all)

# Load genes data
gns_all <- as.data.frame(fread(file="data/FUMA_Results/genes.txt", sep="\t", header=TRUE))
tmp_loci <- as.numeric(unlist(lapply(strsplit(as.character(gns_all[["GenomicLocus"]]), ":"), function(x) x[[1]])))
gns_all_filt <- cbind.data.frame(gns_all[,c(3,4,5,2)], Loci=tmp_loci, ENS=as.character(gns_all[,1]))
gns_all_filt[,1] <- paste0("chr", gns_all_filt[,1])
gns_all_filt$V5 <- "id=0"
gns_all_filt$V5[which(gns_all$ciMap == "Yes")] <- "id=1"
if(sum(colnames(gns_all) %in% "eqtlMapSNPs") > 0) {
  gns_all_filt$V5[which(gns_all$eqtlMapSNPs > 0 & gns_all$ciMap == "Yes")] <- "id=2" 
}
gns_all_filt_inp <- do.call(rbind.data.frame,as.list(apply(gns_all_filt, 1, function(x) {
  if(x[1] %in% loci_all_filt$chr) {
    tmp.loci <- loci_all_filt[loci_all_filt$chr %in% x[1],,drop=FALSE]
    if(as.numeric(x[5]) %in% as.numeric(tmp.loci$GenomicLocus)) {
      return(x[-5])
    }
  }
})))    
colnames(gns_all_filt_inp) <- colnames(gns_all_filt)[-5]
gns_all_filt_inp$ENS <- as.character(gns_all_filt_inp$ENS)
rm(gns_all)

# Load ci data
ci_all <- as.data.frame(fread(file="data/FUMA_Results/ci.txt", sep="\t", header=TRUE))
if(ncol(ci_all) >= 11) {
  ci_all <- ci_all[ci_all[,11] == 1,]
}
ci_all_filt <- ci_all[ci_all[,8] == "intra",]
rm(ci_all)

# Load eqtl data
eqtl_all <- as.data.frame(fread(file="data/FUMA_Results/eqtl.txt", sep="\t", header=TRUE))
if(ncol(eqtl_all) >= 14) {
  eqtl_all <- eqtl_all[eqtl_all[,14] == 1,]
}
eqtl_all_filt <- eqtl_all[,c(4,6,11,12)]
eqtl_all_filt_gns <- eqtl_all_filt[eqtl_all_filt[,1] %in% gns_all_filt$ENS,]
rm(eqtl_all)

## Process SNPs 
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
  processSubChr(snps_all_split[[i]], chrFileList, chr=names(snps_all_split)[i])[,-c(4,5)][,c(1:3,5,4,6)]
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

# for(chr in names(snps_all_split_inp)) {
#   plotCircosByChr(chr, dataList) 
# }


print("--- CircosPlots.R loaded ---")
