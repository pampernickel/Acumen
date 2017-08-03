
visControls <- function(control.map, mat){
  # control.map is a .csv file that contains information on positive and
  # negative controls for an experiment
  # mat are a series of matrices generated using readData(); note that the
  # assumption is that the controls are in the same order across the plates
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
  
  pdf("./results/testVis.pdf",w=30,h=20)
  ggplot(vis.mat, aes(x=drug, fill=value))+geom_histogram(bins=20, alpha=0.5)+facet_wrap(~plate)+
    theme(strip.text = element_text(size=15,face="bold"),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title=element_text(size=14))
  dev.off()
}