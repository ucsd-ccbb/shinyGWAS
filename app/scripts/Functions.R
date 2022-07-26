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
