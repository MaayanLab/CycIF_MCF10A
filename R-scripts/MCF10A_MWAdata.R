drugs <- c("Alpelisib","Palbociclib","Trametinib","Vorinostat")
lapply(drugs,filegen)

filegen <- function(x) {
  temp <- read_delim(paste("C:/Users/maayanlab/Desktop/Workspace/MWA_data/",x,"_ratios.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
  times <- unique(temp$`time(hr)`)
  lapply(times,function(y) write.table(temp[temp$`time(hr)` == y,],file = paste("MWA_",x,"_",y,"hr.tsv",sep=""),quote=FALSE, row.names = FALSE, sep = '\t'))
}