library(readxl)

## Extract relevant rows and columns from raw dataset
raw <- read_excel("C:/Users/maayanlab/Desktop/Workspace/LINCS_P100_DIA_Plate43_annotated_minimized_2017-02-15_11-53-50.processed.xlsx")
data <- raw[-c(1:11,14,15,18:22),-c(2,4:13)]
genes <- data[,c(1,2)]
genes <- genes[-c(1:4),]
genes[,1] <- lapply(genes[,1],function(x) sub('^.*?\\_','',x))
colnames(genes) <- c("","pr_gene_symbol")
data <- data[,-c(1,2)]
data[data=="Taxol"] <- "Paclitaxel"
data[data=="BYL719"] <- "Alpelisib"

drugs <- unique(c(data[3,]))
for (i in 1:length(drugs)) {
  drugName <- drugs[[i]]
  tempData <- data[,c(data[3,]==drugName)]
  colnames(tempData) <- paste(round(as.numeric(tempData[2,]),3),tempData[4,],tempData[1,],sep = "_")
  tempData <- tempData[-c(1:4),]
  tempData <- data.matrix(tempData)
  tempFile <- data.frame(matrix(NA,nrow = nrow(tempData),ncol = 0))
  row.names(tempFile) <- c(genes[,1][[1]])
  j <- 1
  while(j < ncol(tempData)) {
    if (strsplit(colnames(tempData)[j+1],"_")[[1]][3]!=2) {
      tempFile[,ncol(tempFile)+1] <- tempData[,j]
      j = j+1
    }
    else {
      tempFile[,ncol(tempFile)+1] <- rowMeans(tempData[,j:(j+1)],na.rm = TRUE)
      j = j+2
    }
    colnames(tempFile)[ncol(tempFile)] <- substring(colnames(tempData)[j-1],1,nchar(colnames(tempData)[j-1])-2)
  }
  write.table(tempFile, paste("C:/Users/maayanlab/MCF10A/MCF10A_app/app/static/clustergrammer/tsv/P100_",tolower(drugName),".tsv",sep = ""), sep="\t", quote=F, row.names=T, col.names=NA)
  tempFile <- cbind(genes[,2],tempFile)
  write.table(tempFile, paste("C:/Users/maayanlab/Desktop/P100_",tolower(drugName),".tsv",sep = ""), sep="\t", quote=F, row.names=T, col.names=NA)
}