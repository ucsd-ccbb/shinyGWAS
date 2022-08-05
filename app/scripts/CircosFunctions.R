## functions and packages

#library(IRdisplay)
## Processing functions
generateChromosomeFiles <- function(c, loci, ci, snp, gene) {
  loci2 <- cbind.data.frame(loci, loci[,5:6])
  if(!is.null(ci)) {
    for (i in 1:nrow(ci)) {
      ii <- ci[i,]
      if(min(as.numeric(ii[3]), as.numeric(ii[6])) < loci2[loci2[,1] == as.numeric(ii[1]), 7]) {
        loci2[loci2[,1] == as.numeric(ii[1]),7] <- min(as.numeric(ii[3]), as.numeric(ii[6]))
      }
      if(max(as.numeric(ii[4]), as.numeric(ii[7])) > loci2[loci2[,1] == as.numeric(ii[1]), 8]) {
        loci2[loci2[,1] == as.numeric(ii[1]),8] <- max(as.numeric(ii[4]), as.numeric(ii[7]))
      }
    }}
  if(nrow(gene) > 0) {
    for (i in 1:nrow(gene)) {
      g <- gene[i,]
      if(as.numeric(g[1,2]) < loci2[loci2[,1] == as.numeric(g[1,5]), 7]) {
        loci2[loci2[,1] == as.numeric(g[1,5]),7] <- as.numeric(g[1,2])
      }
      if(as.numeric(g[1,3]) > loci2[loci2[,1] == as.numeric(g[1,5]), 8]) {
        loci2[loci2[,1] == as.numeric(g[1,5]),8] <- as.numeric(g[1,3])
      }
    }}
  cur_pos <- 0
  tmp_end <- c()
  for (i in 1:nrow(loci2)) {
    l <- loci2[i,]
    if (cur_pos == 0) {
      if (((as.numeric(l[7])-1000)/1000000) <= 0) {
        tmp_start <- 0
      } else {
        tmp_start <- (((as.numeric(l[7])-1000)/1000000)-1)*1000000
        cur_pos = as.numeric(l[8])
      }
    } else if ((((as.numeric(l[7])-1000)/1000000)-1)-(((cur_pos+1000)/1000000)+1) <= 1) {
      cur_pos = max(cur_pos, as.numeric(l[8]))
    } else {
      tmp_end <- c(tmp_end, ((((cur_pos+1000)/1000000)+1)*1000000))
      tmp_start <- c(tmp_start, ((((as.numeric(l[7])-1000)/1000000)-1)*1000000))
      cur_pos = as.numeric(l[8])
    }
  }
  tmp_end <- c(tmp_end, ((((cur_pos+1000)/1000000)+1)*1000000))
  reg = data.frame(chr=rep(c, length(tmp_start)), start=tmp_start, end=tmp_end)
  
  reg$value <- "#59bfff"
  reg$value2 <- "#cde7f7"
  reg[,1] <- paste0(unique(as.character(reg[,1])), "_",seq(1:nrow(reg)))
  highlights <- loci[,c(3,5,6)]
  highlights$value <- highlights$value2 <- "#1520A6"
  reg_highlights <- rbind.data.frame(reg, highlights)
  reg_highlights <- reg_highlights[order(reg_highlights[,2]),]
  j <- 0
  for (i in 1:nrow(reg_highlights)) {
    if(reg_highlights[i,1] == c) {
      reg_highlights[i,1] <- paste0(reg_highlights[i,1], "_", j)
    } else {
      j <- j + 1
    }
  }
  reg_highlights 
}

processSubChr <- function(inp, chrList=chrFileList, chr=NULL, link=FALSE) {
  if(is.null(chr)) {
    tmp_df <- as.data.frame(t(apply(inp, 1, function(x) {
      chr <- x[1]
      as.character(assignSubChr(x, chrList, chr))
    })))
  } else {
    tmp_df <- as.data.frame(t(apply(inp, 1, function(x) as.character(assignSubChr(x, chrList, chr, link)))))
  }
  tmp_df[,2] <- as.numeric(tmp_df[,2])     
  tmp_df[,3] <- as.numeric(tmp_df[,3])     
  if(link) {
    tmp_df[,5] <- as.numeric(tmp_df[,5])     
    tmp_df[,6] <- as.numeric(tmp_df[,6])
  }
  colnames(tmp_df) <- paste0("V", seq(1:ncol(tmp_df)))
  tmp_df[grep("_", tmp_df[,1]),]
}

assignSubChr <- function(x, chrList, chr, link=FALSE) {
  tmp_chr <- chrList[[chr]]
  if(nrow(tmp_chr) > 0) {
    for (i in 1:nrow(tmp_chr)) {
      if (x[1] == chr) {
        if(as.numeric(x[2]) > tmp_chr[i,2] & as.numeric(x[2]) < tmp_chr[i,3]) {
          if (link) { 
            x[1] <- x[4] <- tmp_chr[i,1]
          } else {
            x[1] <- tmp_chr[i,1] 
          }
          return(as.character(x))
        }
      }
    }
  }     
}

## Plot function by chr
plotCircosByChr <- function(chr, dataList) {
  gns <- dataList[["gns"]]
  snps <- dataList[["snps"]]
  chrFile <- dataList[["chrFileList"]]
  ci <- dataList[["ci"]]
  eqtl <- dataList[["eqtl"]]
  if(nrow(gns[grep(paste0(chr, "_"), gns[,1]),]) < 1) {
    message(paste0("No genes for ", chr, "!"))
  } else {
    snps[[chr]][,4] <- as.numeric(snps[[chr]][,4])
    diff <- max(snps[[chr]][,4])
    if(diff/10 > 1) {
      if(diff/100 > 1) {
        incr <- seq(25, diff, by=25)
      } else {
        incr <- seq(2.5, diff, by=2.5)
      }
    } else  {
      incr <- seq(1, diff, by=1)
    }
    col_text <- "grey40"
    circos.clear()
    #pdf(paste0("CircosPlot_", chr), width=5, height=5)
    circos.par("start.degree" = 90)
    circos.genomicInitialize(chrFile[[chr]], plotType = NULL)
    circos.genomicTrackPlotRegion(snps[[chr]], bg.border=F,
                                  panel.fun = function(region, value, ...) {
                                    circos.genomicPoints(region, value, cex=0.1, col=value[[3]], border=NA, bg.border=F, ...)
                                    tmp.xlim <- get.cell.meta.data("xlim")[1]
                                    tmp.xlim.adj <- 0.001 * tmp.xlim
                                    circos.text(tmp.xlim+tmp.xlim.adj, get.cell.meta.data("ycenter"), 
                                                labels = "-log(gwasP)", col="grey40", cex=0.7,
                                                facing = "clockwise")
                                    for (i in incr) {
                                      circos.yaxis(at=i, labels.cex=0.5, lwd=0, tick.length=0, 
                                                   labels.col=col_text, col="grey40")
                                      circos.segments(x0=0, x1=max(snps[[chr]]$V3), y0=incr, y1=incr, 
                                                      lwd=0.6, lty="11", col="grey90")
                                    }

                                    
                                  },cell.padding = c(0, 0, 0, 0), track.margin = c(0.07,0.05)
    )
    for (i in unique(snps[[chr]][,1])) {
      tmp.snps <- snps[[chr]][snps[[chr]][,1] %in% i,]
      max.snp <- tmp.snps[as.numeric(tmp.snps[,4]) == max(as.numeric(tmp.snps[,4])),,drop=FALSE]
      max.snp[,4] <- as.numeric(max.snp[,4])
      circos.genomicText(region=max.snp[,1:3,drop=FALSE], value=max.snp[,-c(1:3),drop=FALSE], cex=0.7, 
                         labels.column=2, numeric.column=1, track.index=1, sector.index=max.snp[1,1],
                         col="grey40", adj=c(0.5,-0.3))
    }
    circos.genomicTrackPlotRegion(chrFile[[chr]], ylim = c(0, 1), bg.border = NA, track.height = 0.05,
                                  panel.fun = function(region, value, ...) {
                                    tmp.xlim = get.cell.meta.data("xlim", get.current.sector.index(), 
                                                                  track.index = get.current.track.index()) 
                                    tmp.major.at = NULL
                                    tmp.major.at = seq(floor(tmp.xlim[1]/1000000)*1000000, tmp.xlim[2], by = 1000000)
                                    tmp.major.at = c(tmp.major.at, tmp.major.at[length(tmp.major.at)] + 1000000)
                                    tmp.major.tick.labels = paste(round(tmp.major.at/1000000, digits=2),
                                                                  "Mb", sep = "")
                                    
                                    circos.genomicAxis(h="top",tickLabelsStartFromZero = FALSE, major.by = 1000000,
                                                       major.at=tmp.major.at, labels=tmp.major.tick.labels,
                                                       track.index = get.current.track.index(), labels.cex=0.5,
                                                       labels.facing="clockwise", col=col_text, labels.col=col_text)
                                    
                                    col = value[[1]]
                                    circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = col, 
                                                       border = F, ...)
                                    xlim = get.cell.meta.data("xlim")
                                    circos.rect(xlim[1], 0, xlim[2], 1, border = "grey60")
                                  }, cell.padding = c(0, 0, 0, 0), track.margin = c(0.01,0.1)
    )
    circos.genomicLabels(gns_all_filt_inp2[grep(paste0(chr, "_"), gns[,1]),], #nslabels.column = 4, 
                         labels.column=4, cex=0.55, font=8, side = "outside", padding=0.03, track.margin=c(0.01, 0.02),
                         col = gns[[6]], line_col = "grey70", connection_height = 0.02
    )
    
    circos.genomicTrackPlotRegion(chrFile[[chr]], ylim = c(0, 1), bg.border = NA, track.height = 0.05,
                                  panel.fun = function(region, value, ...) {
                                    col = value[[2]]
                                    circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = col, 
                                                       border = F, ...)
                                    xlim = get.cell.meta.data("xlim")
                                    circos.rect(xlim[1], 0, xlim[2], 1, border = "grey60")
                                  }, cell.padding = c(0, 0, 0, 0), track.margin = c(-0.001,0.001)
    )
    if (!is.null(ci_all_split_inp[[chr]])) {
      circos.genomicLink(region1 = ci[[chr]][,1:3], ci[[chr]][,4:6], 
                         col=ci[[chr]][,7], lwd=0.01,border=NA, h.ratio=0.5)
    }
    if (!is.null(eqtl_all_split_inp[[chr]])) {
      circos.genomicLink(region1 = eqtl[[chr]][,1:3], eqtl[[chr]][,4:6], 
                         col=eqtl[[chr]][,7], lwd=0.01,border=NA, h.ratio=0.5)
    }
    #dev.off()
  }}


print("--- CircosFunctions.R loaded ---")
