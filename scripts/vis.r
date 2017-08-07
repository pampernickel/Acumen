
visControls <- function(control.map, mat, ret.mode="image"){
  # control.map is a .csv file that contains information on positive and
  # negative controls for an experiment
  # mat are a series of matrices generated using readData(); note that the
  # assumption is that the controls are in the same order across the plates
  if (!is.loaded("ggplot2")){
    library(ggplot2)
  }
  
  cm <- read.csv(control.map, stringsAsFactors = F)
  vis.mat <- matrix(0, nrow=0, ncol=3)
  colnames(vis.mat) <- c("drug", "value", "plate")
  for (i in 1:length(mat)){
    names(mat)[i] -> nn
    unlist(strsplit(nn, "\\/")) -> nn
    gsub("_plate.csv", "", nn[length(nn)]) -> nn
    
    apply(cm, 1, function(x){
      # get collection of coords
      mat[[i]] -> d
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
    
    # collapse list to values that can be used for ggplot histogram
    for (j in 1:length(vals)){
      cbind(vals[[j]], rep(names(vals)[j], length(vals[[j]])), 
            rep(nn, length(vals[[j]]))) -> t
      colnames(t) <- colnames(vis.mat)
      rbind(vis.mat, t) -> vis.mat
    }
  }
  
  as.data.frame(vis.mat) -> vis.mat
  as.numeric(as.character(vis.mat$drug)) -> vis.mat$drug
  as.character(vis.mat$plate) -> vis.mat$plate
  sapply(strsplit(vis.mat$plate, "\\."), function(x) x[2]) -> vis.mat$plate
  if (ret.mode == "image"){
    ggplot(vis.mat, aes(x=drug, fill=value))+geom_histogram(bins=20, alpha=0.5)+facet_wrap(~plate)+
      theme(strip.text = element_text(size=15,face="bold"),
            axis.text.x=element_text(size=14),
            axis.text.y=element_text(size=14),
            axis.title=element_text(size=14)) -> p
  } else {
    vis.mat -> p #return matrix instead
  }
    
  return(p)
}

visTops <- function(dms, tops){
  if (!is.loaded("ggplot2")){
    library(ggplot2)
  }
  
  labs <- rep("o", nrow(dms))
  labs[which(dms$Drug.name %in% tops$Drug.name)] <- "top"
  cbind(dms, labs) -> vis
  
  # sort x based on current order on vis
  vis[order(vis$diff),] -> vis
  vis$Drug.name <- factor(vis$Drug.name, levels = vis$Drug.name)
  if (length(which(duplicated(as.character(vis$Drug.name)))) > 0){
    vis[-which(duplicated(as.character(vis$Drug.name))),] -> vis
  }
  
  ggplot(vis, aes(x=Drug.name, y=diff, fill = labs))+
    geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4))+
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(face="bold",angle=90)) -> p
  return(p)
}

visPlateWControls <- function(vis.mat, mat, control.map){
  if (!is.loaded("ggplot2")){
    library(ggplot2)
  }
  
  # add contents of mat to t1, which is generated using visControls
  cm <- read.csv(control.map, stringsAsFactors = F)
  unique(cm$Column) -> rm.cols
  for (i in 1:length(mat)){
    names(mat)[i] -> nn
    unlist(strsplit(nn, "\\/")) -> nn
    gsub("_plate.csv", "", nn[length(nn)]) -> nn
    sapply(strsplit(nn, "\\."), function(x) x[2]) -> nn
    as.vector(mat[[i]][,-rm.cols]) -> vals
    cbind(vals, rep("other", rep(length(vals))), rep(nn, length(vals))) -> t
    colnames(t) <- colnames(vis.mat)
    rbind(vis.mat, t) -> vis.mat
  }
  
  as.data.frame(vis.mat) -> vis.mat
  as.numeric(as.character(vis.mat$drug)) -> vis.mat$drug
  as.character(vis.mat$plate) -> vis.mat$plate
  ggplot(vis.mat, aes(x=drug, fill=value))+geom_histogram(bins=20, alpha=0.5)+facet_wrap(~plate)+
    theme(strip.text = element_text(size=15,face="bold"),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title=element_text(size=14)) -> p
  return(p)
}

subsetDF <- function(df, colname){
  # generic function to extract df subset based on elements
  # of a column name
  unique(df[,which(colnames(df) %in% colname)]) -> elems
  which(colnames(df) %in% colname) -> cind
  lapply(elems, function(x) 
    df[which(df[,cind] %in% x),]) -> dfs
  names(dfs) <- elems
  return(dfs)
}

makePlate <- function(pf=384){
  # extract position limits and create plate with these dimensions
  if (pf == 384){
    rep(LETTERS[1:16], 24) -> rl
    rep(1:16, 24) -> rn
    cbind(rn[order(rn)], rep(1:24, 16)) -> mat
    paste("0", mat[which(nchar(mat[,2]) == 1),2], sep="") -> mat[which(nchar(mat[,2]) == 1),2]
    paste(rl[order(rl)], mat[,2], sep="") -> nn
    cbind(nn, mat) -> mat
    as.data.frame(mat) -> mat
    as.numeric(as.character(mat[,2])) -> mat[,2]
    as.numeric(as.character(mat[,3])) -> mat[,3]
    colnames(mat) <- c("well", "row", "column")
  }
  
  return(mat)
}

plotActiveDrugs <- function(tops, plate="all", pf=384){
  if (!is.loaded("ggplot2")){
    library(ggplot2)
  }
  
  if (plate == "all"){
    subsetDF(tops, "plate.ID..384.") -> dfs
    vis <- matrix(0, nrow=0, ncol=5)
    colnames(vis) <- c("well", "row", "col", "tag", "plate")
    for (i in 1:length(dfs)){
      makePlate() -> p
      lab <- rep("n", nrow(p))
      lab[which(p[,1] %in% dfs[[i]]$position..384.)] <- "y"
      cbind(p, lab, rep(names(dfs)[i], nrow(p))) -> t
      colnames(t) <- colnames(vis)
      rbind(vis, t) -> vis
    }
    as.data.frame(vis) -> vis
    as.numeric(as.character(vis$row)) -> vis$row
    as.numeric(as.character(vis$col)) -> vis$col
    scale_color_manual(values=c("lightgrey", "red")) -> man.col
    as.character(unique(vis$plate)) -> plates
    vis$plate <- factor(vis$plate, 
                        levels = plates[order(plates)])
    
    # invert row order such that 1,1 starts from top to bottom rather than
    # from the bottom
    vis$row <- factor(vis$row, 
                      levels = unique(vis$row)[rev(order(unique(vis$row)))])
    as.numeric(as.character(vis$row))
    ggplot(vis, aes(x=col, y=row, color=tag))+
      geom_point()+facet_wrap(~plate)+man.col+
      #scale_y_continuous(minor_breaks = seq(0, 16, 4))+
      theme(panel.background = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(face="bold",colour="black"),
        legend.position = "bottom",
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> p
  } else {
    
  }
  return(p)
}