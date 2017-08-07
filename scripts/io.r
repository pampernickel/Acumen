#::: functions for processing Acumen data (Nexus)

source("./scripts/vis.r")

`%ni%` <- Negate(`%in%`)

readData <- function(dir, field="", control.plates=""){
  list.files(dir, pattern=".csv", full.names=T) -> files
  lapply(files, function(x) suppressWarnings(readLines(x))) -> fc
  
  # select relevant portion of files
  if (field != ""){
    # find row in which field is and extract the plate immediately below it
    lapply(fc, function(x){
      d <- NA
      grep(field, x, ignore.case = T) -> rind
      grep("Statistic", x, ignore.case = T) -> sind
      if (length(sind) > 0){
        sind[c(which(sind %in% rind), which(sind %in% rind)+1)] -> lim
        strsplit(x[(lim[1]+1):(lim[2]-1)], ",") -> mat
        
        # check mat lengths
        length(unique(sapply(mat, length))) -> l
        if (l == 1){
          d <- matrix(0, nrow=length(mat)-1, ncol=length(mat[[1]])-1)
          for (i in 2:length(mat)){
            d[i-1,]<- as.numeric(mat[[i]][2:length(mat[[i]])])
          }
          colnames(d) <- mat[[1]][2:length(mat[[1]])]
          rownames(d) <- sapply(mat[2:length(mat)], function(y) y[1])
        }
      }
      return(d)
    }) -> mat
    names(mat) <- files
    mat[which(sapply(mat, function(x) is.null(nrow(x))) %in% F)] -> mat
  }
  
  # if there are control plates, remove them and renumber plates accordingly
  if (control.plates != ""){
    unlist(lapply(control.plates, function(x) grep(x, names(mat)))) -> rm.ind
    if (length(rm.ind) > 0){
      mat[-rm.ind] -> mat
    }
  }
  
  return(mat) 
}

# map mat to different drugs
procPlates <- function(mat, drug.map, control.map, pos, neg, fin.map=""){
  # control.map is a .csv file that contains information on positive and
  # negative controls for an experiment
  # pos and neg are positive and negative controls, to be chosen among the drugs
  # in the control.map
  if (fin.map != ""){
    read.csv(fin.map, stringsAsFactors = F) -> fm
    sapply(names(mat), function(x) fm$code[which(fm$file %in% x)]) -> code.map
  }
  
  # process control map
  cm <- read.csv(control.map, stringsAsFactors = F)
  
  if (pos %ni% cm$Name || neg %ni% cm$Name){
    stop("Your control names must match content under the `Name' column of control.map.")
  } else {
    #visualize controls
    normalize(control.map, mat, pos, neg) -> matn
  }
}

normalize <- function(control.map, mat, pos, neg){
  cm <- read.csv(control.map, stringsAsFactors = F)
  dnorm <- list()
  
  strsplit(names(mat), "\\/") -> nn
  gsub("_plate.csv", "", sapply(nn, function(x) x[length(x)])) -> nn
  
  for (i in 1:length(mat)){
    mat[[i]] -> d
    apply(cm, 1, function(x){
      # get collection of coords
      which(rownames(d) %in% x[3]) -> start
      nrow(d) -> end
      seq(start, end, by=as.numeric(x[5])) -> int
      as.numeric(d[int,as.numeric(x[4])])
    }) -> vals
    names(vals) <- cm$Name
    
    if (length(unique(cm$Name)) != length(vals)){
      # collapse vals with the same name
      lapply(unique(cm$Name), function(x)
        as.numeric(unlist(vals[which(names(vals) %in% x)]))) -> vals
      names(vals) <- unique(cm$Name)
    }
    
    # extract the part of d that does not fall under these column/row combos
    # and normalize based on the specified pos and neg;
    # normalize based on the ff. formula:
    # data-pos.control/neg.control-pos.control
    mean(as.numeric(vals[[which(names(vals) %in% pos)]])) -> pc
    mean(as.numeric(vals[[which(names(vals) %in% neg)]])) -> nc
    (d-pc)/(nc-pc) -> dnorm[[i]]
  }
  names(dnorm) <- nn
  return(dnorm)
}