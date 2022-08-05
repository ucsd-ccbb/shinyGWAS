#load input file
loadInputFile = function(fileDir, fileName, colNames='all'){
  df <- fread(file.path(fileDir, fileName), sep='\t', header=TRUE)
  if (tolower(colNames[1]) != 'all'){
    print(paste0("Loading file with specific columns"))
    print(colNames)
    return(df)
  }
  else{
    print('Loading file with all columns')
    return(df)
  }
}

#prepares df for an IGV track
prepDfToIvgTrack = function(df, colNames, reArrange=FALSE){
  df$chr <- as.character(df[[colNames[1]]]) #make sure this col are "character"
  df$start <- as.numeric(df[[colNames[2]]]) #numeric
  df$end <- as.numeric(df[[colNames[2]]] + 1) #numeric
  #rearrange columns 
  if (reArrange=="TRUE"){
    df <- df %>%
      dplyr::select(c("chr", "start", "end", colNames[3]))
    return(df)
  }
  else{
    return(df)
  }
}


getRDBSupportingInfo = function(RDB_score){
  RDB_score <-  na.omit(RDB_score)
  for (snp in 1:length(RDB_score$RDB)){
    if (RDB_score$RDB[snp] == "1a"){
      RDB_score$Info[snp] <- "1a: eQTL + TF binding + matched TF motif + matched DNase Footprint + DNase peak"
      RDB_score$RDB[snp] <- 1
    }
    if (RDB_score$RDB[snp] == "1b"){
      RDB_score$Info[snp] <- "1b: eQTL + TF binding + any motif + DNase Footprint + DNase peak"
      RDB_score$RDB[snp] <- 1
    }
    if (RDB_score$RDB[snp] == "1c"){
      RDB_score$Info[snp] <- "1c: eQTL + TF binding + matched TF motif ++ DNase peak"
      RDB_score$RDB[snp] <- 1
    }
    if (RDB_score$RDB[snp] == "1d"){
      RDB_score$Info[snp] <- "1da: eQTL + TF binding + any motif + DNase peak"
      RDB_score$RDB[snp] <- 1
    }
    if (RDB_score$RDB[snp] == "1e"){
      RDB_score$Info[snp] <- "1e: eQTL + TF binding + matched TF motif "
      RDB_score$RDB[snp] <- 1
    }
    if (RDB_score$RDB[snp] == "1f"){
      RDB_score$Info[snp] <- "1f: eQTL + TF binding/DNase peak"
      RDB_score$RDB[snp] <- 1
    }
    if (RDB_score$RDB[snp] == "2a"){
      RDB_score$Info[snp] <- "2a: TF binding + matched TF motif + matched DNase Footprint + DNase peak"
      RDB_score$RDB[snp] <- 2
    }
    if (RDB_score$RDB[snp] == "2b"){
      RDB_score$Info[snp] <- "2b: TF binding + any motif + DNase Footprint + DNase peak"
      RDB_score$RDB[snp] <- 2
    }
    if (RDB_score$RDB[snp] == "2c"){
      RDB_score$Info[snp] <- "2c: TF binding + matched TF motif + DNase peak"
      RDB_score$RDB[snp] <- 2
    }
    if (RDB_score$RDB[snp] == "3a"){
      RDB_score$Info[snp] <- "3a: TF binding + any motif + DNase peak"
      RDB_score$RDB[snp] <- 3
    }
    if (RDB_score$RDB[snp] == "3b"){
      RDB_score$Info[snp] <- "3b: TF binding + matched TF motif"
      RDB_score$RDB[snp] <- 3
    }
    if (RDB_score$RDB[snp] == "4"){
      RDB_score$Info[snp] <- "4: TF binding + DNase peak"
      RDB_score$RDB[snp] <- 4
    }
    if (RDB_score$RDB[snp] == "5"){
      RDB_score$Info[snp] <- "5: TF binding OR DNase peak"
      RDB_score$RDB[snp] <- 5
    }
    if (RDB_score$RDB[snp] == "6"){
      RDB_score$Info[snp] <- "6: Motif hit"
      RDB_score$RDB[snp] <- 6
    }
    if (RDB_score$RDB[snp] == "7"){
      RDB_score$Info[snp] <- "7: Other"
      RDB_score$RDB[snp] <- 7
    }
  }
  RDB_score$RDB <- as.numeric(RDB_score$RDB)
  return(RDB_score)
}


print("--- Functions.R loaded ---")
