##See for MCF 10A first dataset details: http://lincs.hms.harvard.edu/db/datasets/20303/main
source('http://bioconductor.org/biocLite.R')
biocLite('preprocessCore')

library(readxl)
library(plotly)
library(preprocessCore)



## Generate new HMS dataset data frames and labels
##Note: Preprocessed datasets manually in excel to speed up upload time into R
for(i in 3:8) {
  ##Retrieve HMS datasets
  name=paste("HMS_Dataset_2030",i,sep = "")
  assign(name,read_excel(paste("C:/Users/maayanlab/Desktop/Workspace/HMS_data/",name,".xlsx",sep="")))
  
  ##Extract drug name, dose concentration, combination for plot
  assign(paste(name,"labels",sep = "_"),get(name)[,1:3])
  assign(name,get(name)[,-c(1:3)])
}

##Cytosolic data compilation
HMS_Dataset_20303_labels[,3] <- paste(HMS_Dataset_20303_labels$`DrugName+Conc`,"24hr",sep = "+")
HMS_Dataset_20303_complete <- cbind(HMS_Dataset_20303_labels[,3],HMS_Dataset_20303)
HMS_Dataset_20304_labels[,3] <- paste(HMS_Dataset_20304_labels$`DrugName+Conc`,"48hr",sep = "+")
HMS_Dataset_20304_complete <- cbind(HMS_Dataset_20304_labels[,3],HMS_Dataset_20304)
HMS_Dataset_20305_labels[,3] <- paste(HMS_Dataset_20305_labels$`DrugName+Conc`,"72hr",sep = "+")
HMS_Dataset_20305_complete <- cbind(HMS_Dataset_20305_labels[,3],HMS_Dataset_20305)
HMS_Dataset_cytosolic <- rbind(HMS_Dataset_20303_complete,HMS_Dataset_20304_complete,HMS_Dataset_20305_complete)
HMS_Dataset_cytosolic[,1] <- gsub("+","_",HMS_Dataset_cytosolic[,1],fixed = TRUE)
colnames(HMS_Dataset_cytosolic) <- c("DrugName+Conc","RP-Ap32;ratIgG(H+L)","Rb(pS807;pS811);Rb(pS807);Rb(pS811)", "H3(pS10)", "E2F-1", "FOXO1a", "c-Myc", "Actin", "p53;mouseIgG(H+L)/goat", "p21", "p27Kip1","S6(pS235;pS236)","H2A.X(pS139)","Cyclin-E1","HCSCellMaskDeepRedReagent","FOXO3a;rabbitIgG(H+L)/goat/highlycross-adsorbed","ERK-1(pT202;pY204);ERK-1(pT202);ERK-2(pT185;pY187);ERK-2(pT185)","Cyclin-D1","KI-67","EGFR","PCNA","MitoTrackerGreenReagent")

##Cytosolic data mean and zscore
cytosol_mean <- data.frame()
cytosol_zscore <- data.frame()
numDrugs <- length(unique(HMS_Dataset_cytosolic$`DrugName+Conc`))
for (j in 2:length(colnames(HMS_Dataset_cytosolic))) {
  for (i in 1:numDrugs) {
    temp <- HMS_Dataset_cytosolic[HMS_Dataset_cytosolic$`DrugName+Conc`==unique(HMS_Dataset_cytosolic$`DrugName+Conc`)[i],j]
    cytosol_mean[i,j-1] <- mean(temp)
    temp <- (temp - mean(HMS_Dataset_cytosolic[,j]))/sd(HMS_Dataset_cytosolic[,j])
    cytosol_zscore[i,j-1] <- mean(temp)
  }
}
row.names(cytosol_mean) <- unique(HMS_Dataset_cytosolic$`DrugName+Conc`)
colnames(cytosol_mean) <- colnames(HMS_Dataset_cytosolic)[2:22]
row.names(cytosol_zscore) <- unique(HMS_Dataset_cytosolic$`DrugName+Conc`)
colnames(cytosol_zscore) <- colnames(HMS_Dataset_cytosolic)[2:22]

##Nuclear data compilation
HMS_Dataset_20306_labels[,3] <- paste(HMS_Dataset_20306_labels$`DrugName+Conc`,"24hr",sep = "+")
HMS_Dataset_20306_complete <- cbind(HMS_Dataset_20306_labels[,3],HMS_Dataset_20306)
HMS_Dataset_20307_labels[,3] <- paste(HMS_Dataset_20307_labels$`DrugName+Conc`,"48hr",sep = "+")
HMS_Dataset_20307_complete <- cbind(HMS_Dataset_20307_labels[,3],HMS_Dataset_20307)
HMS_Dataset_20308_labels[,3] <- paste(HMS_Dataset_20308_labels$`DrugName+Conc`,"72hr",sep = "+")
HMS_Dataset_20308_complete <- cbind(HMS_Dataset_20308_labels[,3],HMS_Dataset_20308)
HMS_Dataset_nuclear <- rbind(HMS_Dataset_20306_complete,HMS_Dataset_20307_complete,HMS_Dataset_20308_complete)
HMS_Dataset_nuclear[,1] <- gsub("+","_",HMS_Dataset_nuclear[,1],fixed = TRUE)
colnames(HMS_Dataset_nuclear) <- c("DrugName+Conc","RP-Ap32;ratIgG(H+L)","Rb(pS807;pS811);Rb(pS807);Rb(pS811)", "H3(pS10)", "E2F-1", "FOXO1a", "c-Myc", "Actin", "p53;mouseIgG(H+L)/goat", "p21", "p27Kip1","S6(pS235;pS236)","H2A.X(pS139)","Cyclin-E1","HCSCellMaskDeepRedReagent","FOXO3a;rabbitIgG(H+L)/goat/highlycross-adsorbed","ERK-1(pT202;pY204);ERK-1(pT202);ERK-2(pT185;pY187);ERK-2(pT185)","Cyclin-D1","KI-67","EGFR","PCNA","MitoTrackerGreenReagent")

##Cytosolic data mean and zscore
nucleus_mean <- data.frame()
nucleus_zscore <- data.frame()
numDrugs <- length(unique(HMS_Dataset_nuclear$`DrugName+Conc`))
for (j in 2:length(colnames(HMS_Dataset_nuclear))) {
  for (i in 1:numDrugs) {
    temp <- HMS_Dataset_nuclear[HMS_Dataset_nuclear$`DrugName+Conc`==unique(HMS_Dataset_nuclear$`DrugName+Conc`)[i],j]
    nucleus_mean[i,j-1] <- mean(temp)
    temp <- (temp - mean(HMS_Dataset_nuclear[,j]))/sd(HMS_Dataset_nuclear[,j])
    nucleus_zscore[i,j-1] <- mean(temp)
  }
}
row.names(nucleus_mean) <- unique(HMS_Dataset_nuclear$`DrugName+Conc`)
colnames(nucleus_mean) <- colnames(HMS_Dataset_nuclear)[2:22]
row.names(nucleus_zscore) <- unique(HMS_Dataset_nuclear$`DrugName+Conc`)
colnames(nucleus_zscore) <- colnames(HMS_Dataset_nuclear)[2:22]

##Averaged quantile normalization
nucleus_quantnorm <- nucleus_mean
temp <- as.data.frame(normalize.quantiles(as.matrix(HMS_Dataset_nuclear[,2:22])))
temp <- cbind(HMS_Dataset_nuclear[,1], temp)
colnames(temp) <- colnames(HMS_Dataset_nuclear)
numDrugs <- length(unique(temp$`DrugName+Conc`))
for (j in 2:length(colnames(temp))) {
  for (i in 1:numDrugs) {
    nucleus_quantnorm[i,j-1]  <- mean(temp[temp$`DrugName+Conc`==unique(temp$`DrugName+Conc`)[i],j])
  }
}

cytosol_quantnorm <- cytosol_mean
temp <- as.data.frame(normalize.quantiles(as.matrix(HMS_Dataset_cytosolic[,2:22])))
temp <- cbind(HMS_Dataset_cytosolic[,1], temp)
colnames(temp) <- colnames(HMS_Dataset_cytosolic)
numDrugs <- length(unique(temp$`DrugName+Conc`))
for (j in 2:length(colnames(temp))) {
  for (i in 1:numDrugs) {
    cytosol_quantnorm[i,j-1]  <- mean(temp[temp$`DrugName+Conc`==unique(temp$`DrugName+Conc`)[i],j])
  }
}

##Save cytosolic mean and zscore
##NOTE: Remove 'DrugName+Conc' from first line of file after saving it
write.table(cytosol_mean, file = 'cytosol_mean.tsv', quote=FALSE, sep = '\t')
write.table(cytosol_zscore, file = 'cytosol_zscore.tsv', quote=FALSE, sep = '\t')

##Save nuclear mean and zscore
write.table(nucleus_mean, file = 'nucleus_mean.tsv', quote=FALSE, sep = '\t')
write.table(nucleus_zscore, file = 'nucleus_zscore.tsv', quote=FALSE, sep = '\t')

##Save single cell data
write.table(HMS_Dataset_cytosolic, file = 'cytosol_allcells.tsv', quote=FALSE, row.names = FALSE, sep = '\t')
write.table(HMS_Dataset_nuclear, file = 'nucleus_allcells.tsv', quote=FALSE, row.names = FALSE, sep = '\t')

##Save average quantile normalization
write.table(cytosol_quantnorm, file = 'cytosol_quantnorm.tsv', quote=FALSE, row.names = FALSE, sep = '\t')
write.table(nucleus_quantnorm, file = 'nucleus_quantnorm.tsv', quote=FALSE, row.names = FALSE, sep = '\t')

##Heatmap plotting
plot_ly(x = colnames(cytosol_mean), y = row.names(cytosol_mean), z = as.matrix(cytosol_mean), colorscale = "Greys", type = "heatmap") %>% layout(title = "cytosol_mean")
plot_ly(x = colnames(cytosol_zscore), y = row.names(cytosol_zscore), z = as.matrix(cytosol_zscore), colorscale = "Greys", type = "heatmap") %>% layout(title = "cytosol_zscore")
plot_ly(x = colnames(nucleus_mean), y = row.names(nucleus_mean), z = as.matrix(nucleus_mean), colorscale = "Greys", type = "heatmap") %>% layout(title = "nucleus_mean")
plot_ly(x = colnames(nucleus_zscore), y = row.names(nucleus_zscore), z = as.matrix(nucleus_zscore), colorscale = "Greys", type = "heatmap") %>% layout(title = "nucleus_zscore")