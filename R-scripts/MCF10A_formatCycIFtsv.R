files <- list.files(path="C:/Users/maayanlab/MCF10A/MCF10A_app/app/static/data/cycif/values", pattern="*.tsv", full.names=T, recursive=FALSE)
lapply(files, function(x) {
  t <- read.table(x, header=T,check.names = FALSE) # load file
  t <- as.data.frame(t)
  t$agent <- tolower(t$agent)
  name <- paste("CycIF_",t$agent[1],"_",t$time[1],"h.tsv",sep = "")
  t[,4] <- paste(tolower(t$agent),t$concentration,t$time,t$compartment,sep = "_")
  t <- t[,-c(1:3)]
  colnames(t)[1] <- ""
  out <- t(t)
  write.table(out, paste("C:/Users/maayanlab/MCF10A/MCF10A_app/app/static/clustergrammer/tsv/",name,sep = ""), sep="\t", quote=F, row.names=T, col.names=F)
    
})

# write.table(out, "C:/Users/maayanlab/MCF10A/MCF10A_app/app/static/clustergrammer/tsv", sep="\t", quote=F, row.names=T, col.names=T)