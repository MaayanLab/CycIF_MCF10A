library(readxl)

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
HMS_Dataset_20303_labels[,3] <- "24"
HMS_Dataset_20304_labels[,3] <- "48"
HMS_Dataset_20305_labels[,3] <- "72"
HMS_Dataset_cytlabels <- rbind(HMS_Dataset_20303_labels,HMS_Dataset_20304_labels,HMS_Dataset_20305_labels)
HMS_Dataset_cytlabels[,4] <- "cytosol"
HMS_Dataset_cytdata <- rbind(HMS_Dataset_20303,HMS_Dataset_20304,HMS_Dataset_20305)
HMS_Dataset_cytdata <- cbind(HMS_Dataset_cytlabels,HMS_Dataset_cytdata)

##Nuclear data compilation
HMS_Dataset_20306_labels[,3] <- "24"
HMS_Dataset_20307_labels[,3] <- "48"
HMS_Dataset_20308_labels[,3] <- "72"
HMS_Dataset_nuclabels <- rbind(HMS_Dataset_20306_labels,HMS_Dataset_20307_labels,HMS_Dataset_20308_labels)
HMS_Dataset_nuclabels[,4] <- "nucleus"
HMS_Dataset_nucdata <- rbind(HMS_Dataset_20306,HMS_Dataset_20307,HMS_Dataset_20308)
HMS_Dataset_nucdata <- cbind(HMS_Dataset_nuclabels,HMS_Dataset_nucdata)

HMS_Dataset_all <- rbind(HMS_Dataset_cytdata,HMS_Dataset_nucdata)
colnames(HMS_Dataset_all) <- c("agent","concentration","time","compartment","RP-Ap32","Rb(pS807;pS811)", "H3(pS10)", "E2F-1", "FOXO1a", "c-Myc", "Actin", "p53", "p21", "p27Kip1","S6(pS235;pS236)","H2A.X(pS139)","Cyclin-E1","HCSCellMaskDeepRed","FOXO3a","ERK-1/ERK-2*","Cyclin-D1","KI-67","EGFR","PCNA","MitoTrackerGreen")
# colnames(HMS_Dataset_all) <- c("agent","concentration","time","compartment","RP-Ap32;ratIgG(H+L)","Rb(pS807;pS811);Rb(pS807);Rb(pS811)", "H3(pS10)", "E2F-1", "FOXO1a", "c-Myc", "Actin", "p53;mouseIgG(H+L)/goat", "p21", "p27Kip1","S6(pS235;pS236)","H2A.X(pS139)","Cyclin-E1","HCSCellMaskDeepRedReagent","FOXO3a;rabbitIgG(H+L)/goat/highlycross-adsorbed","ERK-1(pT202;pY204);ERK-1(pT202);ERK-2(pT185;pY187);ERK-2(pT185)","Cyclin-D1","KI-67","EGFR","PCNA","MitoTrackerGreenReagent")

write.table(HMS_Dataset_all,file = "MCF10A_CycIF_allcells.tsv",quote=FALSE, row.names = FALSE, sep = '\t')

##Calculate mean zscore for each drug, concentration, time, compartment combination
drugName <- ""
drugTime <- 0

#allzscores <- data.frame(matrix(NA,nrow = 0,ncol = 25))

zscore <- function(x,y) mean((x-mean(y))/sd(y))

numDrugs <- length(unique(HMS_Dataset_all$agent))
for (i in 1:numDrugs) {
  drugName = unique(HMS_Dataset_all$agent)[i]
  temp_drug <- HMS_Dataset_all[HMS_Dataset_all$agent==drugName,]
  numTime <- length(unique(temp_drug$time))
  for (j in 1:numTime) {
    drugTime <- unique(temp_drug$time)[j]
    temp_time <- temp_drug[temp_drug$time==drugTime,]
    temp_time[,26] <- paste(temp_time$agent,temp_time$concentration,temp_time$time,temp_time$compartment,sep = "_") ##V26
    combinations <- expand.grid(unique(temp_time$concentration),unique(temp_time$compartment))
    combinations[,3] <- paste(temp_time$agent[1],combinations$Var1,temp_time$time[1],combinations$Var2,sep = "_")
    numCombinations <- nrow(combinations)
    tempFile <- data.frame(matrix(NA,nrow = 0,ncol = 25))
    for (k in 1:numCombinations) {
      temp <- temp_time[temp_time$V26==combinations$V3[k],]
      #tempFile[k,5:25] <- t(as.data.frame(apply(temp[,5:25],2,mean)))
      if (temp$compartment[1] =="nucleus") tempFile[k,5:25] <- t(as.data.frame(mapply(zscore,temp[,5:25],HMS_Dataset_nucdata[,5:25])))
      else tempFile[k,5:25] <- t(as.data.frame(mapply(zscore,temp[,5:25],HMS_Dataset_cytdata[,5:25])))
      tempFile[k,1:4] <- t(as.data.frame(strsplit(combinations[k,3],split = "[_]")))
      colnames(tempFile) <- colnames(temp)[1:25]
    }
    #allzscores <- rbind(allzscores,tempFile)
    write.table(tempFile,file = paste("CycIF_",drugName,"_",drugTime,"h.tsv",sep = ""),quote=FALSE, row.names = FALSE, sep = '\t')
  }
}
#colnames(allzscores) <- colnames(HMS_Dataset_all)
#write.table(allzscores,file = "CycIF_allmeanzscores.tsv",quote=FALSE, row.names = FALSE, sep = '\t')
