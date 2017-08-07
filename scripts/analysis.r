getDelta <- function(matn, p1="KO", p2="WT", drug.map){
  # rank drugs based on differences between wt and ko
  # p1 = prefix set 1; p2 = prefix set 2;
  # the assumption is that the files are automatically named,
  # and have tags as WT/KO
  if (!is.loaded("gdata")){
    library(gdata)    
  }
  read.xls(drug.map) -> dm
  
  # convert dm$position..384. in row, column coordinates of DM
  grep(p1, names(matn)) -> p1
  grep(p2, names(matn)) -> p2
  
  if (length(p1) == length(p2)){
    diffn <- list()
    for (i in 1:length(p1)){
      matn[[p1[i]]]-matn[[p2[i]]] -> diffn[[i]]
    }
    names(diffn) <- names(matn)[p1]
    
    # append diff measurements to dms
    dm[order(as.character(dm$plate.ID..384.)),c(3,4,6)] -> dms
    as.character(dms$plate.ID..384.) -> dms$plate.ID..384.
    as.character(unique(dms$plate.ID..384.)) -> ids
    if (length(ids) == length(diffn)){
      all.res <- matrix(0, nrow=0, ncol=3)
      for (i in 1:length(ids)){
        dms[which(dms$plate.ID..384. %in% ids[i]),] -> dmss
        apply(dmss, 1, function(x){
          as.character(as.matrix(x[2])) -> coords
          substr(coords, 1, 1) -> row
          as.numeric(substr(coords, 2, nchar(coords))) -> col
          diffn[[i]][which(rownames(diffn[[i]]) %in% row), which(colnames(diffn[[i]]) %in% col)] -> res
          matn[[p1[i]]][which(rownames(matn[[p1[i]]]) %in% row), 
                        which(colnames(matn[[p1[i]]]) %in% col)] -> r1
          matn[[p2[i]]][which(rownames(matn[[p2[i]]]) %in% row), 
                        which(colnames(matn[[p2[i]]]) %in% col)] -> r2
          c(r1, r2, res) -> res
          return(res)
        }) -> res
        rbind(all.res, t(res)) -> all.res
      }
      cbind(dms, all.res) -> dms
      as.character(dms$Drug.name) -> dms$Drug.name
      colnames(dms)[(ncol(dms)-2):ncol(dms)] <- c("KO", "WT", "diff")
      dms[rev(order(dms$diff)),] -> dms
    } else {
      stop("Number of plates specified on drug map does not match the number of files.")
    }
  } else {
    stop("Not all samples have a pair. Please check your data folder to ensure that none of the files are missing.")
  }
  return(dms)
}

getTop <- function(dms, perc=.05, dir=c("pos", "neg", "all")){
  # get top x% of extremes in dms
  ceiling(perc*nrow(dms)) -> fin
  if (dir %in% "pos"){
    dms[1:fin,] -> tops
  } else if (dir %in% "neg"){
    dms[(nrow(dms)-fin-1):nrow(dms),] -> tops
  } else {
    # get top x% and bottom x%
    dms[c(c(1:fin), 
          c((nrow(dms)-fin-1):nrow(dms))),] -> tops
  }
  return(tops)
}