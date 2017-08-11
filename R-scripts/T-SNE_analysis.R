## Run once at start of session
library(devtools)
install_version("plotly", version = "4.5.6", repos = "http://cran.us.r-project.org")

library(readxl)
library(Rtsne)
library(plotly)

for(i in 3:8) {
  ptm <- proc.time()
  ##Retrieve HMS datasets
  name=paste("HMS_Dataset_2030",i,sep = "")
  assign(name,read_excel(paste("C:/Users/maayanlab/Desktop/Workspace/HMS_data/",name,".xlsx",sep="")))

  ##Extract drug name, dose concentration, combination for plot
  assign(paste(name,"labels",sep = "_"),get(name)[,1:3])
  assign(name,get(name)[,-c(1:3)])

  ##Run TSNE (4D)
  ##https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf
  tsne_name=paste(name,"tsne",sep="_")
  assign(tsne_name,Rtsne(get(name),dims=4))
  proc.time() - ptm
  
  ##Save results (and retrieve) Note default save location is in My Documents
  saveRDS(tsne_name,file = paste(tsne_name,".rds",sep = ""))
  # tsne_name <- readRDS(paste(tsne_name,".rds",sep = ""))
}

##3D Plot
##Note column name in color can be varied to get differently colorized points
##https://plot.ly/r/reference/#scatter3d
#plot_ly(data = as.data.frame(get(tsne_name)$Y), type="scatter3d", mode="markers", x = get(tsne_name)$Y[,1], y = get(tsne_name)$Y[,2], z = get(tsne_name)$Y[,3], color=get(paste(name,"labels",sep = "_"))$'DrugName+Conc')
name=paste("HMS_Dataset_2030",8,sep = "")
tsne_name=paste(name,"tsne",sep="_")
plot_ly(data = as.data.frame(get(tsne_name)$Y), type="scatter3d", mode="markers", x = get(tsne_name)$Y[,1], y = get(tsne_name)$Y[,2], z = get(tsne_name)$Y[,3], color=get(paste(name,"labels",sep = "_"))$DrugName)

##2D Plot
##plot_ly(data = as.data.frame(tsne_out$Y), type="scatter3d", x = tsne_out$Y[,1], y = tsne_out$Y[,2], z=tsne_out$Y[,3], visible = TRUE, mode='markers', color=get(paste(name,"drug+dose",sep = "_"))$'DrugName+Conc')