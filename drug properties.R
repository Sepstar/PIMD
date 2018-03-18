setwd("F:/R Workspace/")
load("F:/R Workspace/Input files/affni_fusion.Rdata")
C = 32 # number of clusters
cluster_group = spectralClustering(affni_fusion, C)
druglist = as.data.frame(rownames(affni_fusion))
## chemical feature fluctuation
DrugBank_CFmean = t(as.matrix(read.csv("F:/R Workspace/Input files/DrugBank_CF_mean.csv", header = T)))
CF = read.csv("F:/R Workspace/Input files/Chemical Feature.csv", header = T)
CF_cluster_mean = c()
pvalue = c()
special = c()
for(i in 1:C){
  cluster_index = which(cluster_group == i)
  cluster_drug = druglist[cluster_index,]
  cluster_drug = as.character(cluster_drug)
  CF_cluster = c()
  CF_cluster_repurpose = c()
  for(drug in cluster_drug){
    CFindex = which(CF[,1] == drug)
    CF_cluster_temp = as.matrix(CF[CFindex,c(1,5:18)])
    CF_cluster = rbind(CF_cluster, CF_cluster_temp)
  }
  rm(CF_cluster_temp)
  CF_cluster_mean_temp = as.matrix(apply(CF_cluster[,-1],2,mean,na.rm = TRUE))
  CF_cluster_mean = cbind(CF_cluster_mean, CF_cluster_mean_temp)
  pvalue_temp = matrix(data = NA,nrow = nrow(CF_cluster_mean),ncol = 1)
  for(j in 2:15){
    if(all(is.na(CF_cluster[,j])) == FALSE){
      if(length(which(mean(CF_cluster[,j]) == CF_cluster[,j])) == length(CF_cluster[,j])){
        special_temp = c(i,colnames(CF_cluster)[j])
        special = rbind(special, special_temp)
      }else{
        ttest = t.test(CF_cluster[,j], mu = DrugBank_CFmean[j-1,])
        pvalue_temp[j-1,1] = ttest$p.value
      }
    }
  }
  pvalue = cbind(pvalue, pvalue_temp)
  rm(CF_cluster_mean_temp)
}
for(i in 1:nrow(pvalue)){
  for(j in 1:ncol(pvalue)){
    if(pvalue[i,j]>=0.05 | is.na(pvalue[i,j]) == TRUE){
      CF_cluster_mean[i,j] = 0
    }
  }
}
CF_cluster_meanandpvlaue = rbind(CF_cluster_mean,pvalue)
write.csv(CF_cluster_meanandpvlaue,file = "F:/R Workspace/Output files/CF_cluster_meanandpvalue.csv")
for(i in 1:nrow(pvalue)){
  for(j in 1:ncol(pvalue)){
    if(pvalue[i,j]>=0.05 | is.na(pvalue[i,j]) == TRUE){
      CF_cluster_mean[i,j] = 0
    }else{
      CF_cluster_mean[i,j] = (CF_cluster_mean[i,j][[1]] - DrugBank_CFmean[i,1][[1]])/abs(DrugBank_CFmean[i,1][[1]])
    }
  }
}
CF_cluster_meanandpvlaue = rbind(CF_cluster_mean,pvalue)
write.csv(CF_cluster_meanandpvlaue,file = "F:/R Workspace/Output files/CF_cluster_meanandpvalue(fluctuation).csv")

## CP_JD_OD fluctuation
DrugBank_CPJDODmean = t(as.matrix(read.csv("F:/R Workspace/Input files/DrugBank_ODJDCP_mean.csv", header = T)))
CPJDOD = read.csv("F:/R Workspace/Input files/CP_JD_OD.csv", header = T)
CPJDOD_cluster_mean = c()
pvalue = c()
special = c()
for(i in 1:C){
  cluster_index = which(cluster_group == i)
  cluster_drug = druglist[cluster_index,]
  cluster_drug = as.character(cluster_drug)
  CPJDOD_cluster = c()
  CPJDOD_cluster_repurpose = c()
  for(drug in cluster_drug){
    CPJDODindex = which(CPJDOD[,1] == drug)
    CPJDOD_cluster_temp = as.matrix(CPJDOD[CPJDODindex,c(1,3:83)])
    CPJDOD_cluster = rbind(CPJDOD_cluster, CPJDOD_cluster_temp)
  }
  rm(CPJDOD_cluster_temp)
  CPJDOD_cluster_mean_temp = as.matrix(apply(CPJDOD_cluster[,-1],2,mean,na.rm = TRUE))
  CPJDOD_cluster_mean = cbind(CPJDOD_cluster_mean, CPJDOD_cluster_mean_temp)
  pvalue_temp = matrix(data = NA,nrow = nrow(CPJDOD_cluster_mean),ncol = 1)
  for(j in 2:82){
    if(all(is.na(CPJDOD_cluster[,j])) == FALSE){
      if(length(which(mean(CPJDOD_cluster[,j]) == CPJDOD_cluster[,j])) == length(CPJDOD_cluster[,j])){
        special_temp = c(i,colnames(CPJDOD_cluster)[j])
        special = rbind(special, special_temp)
      }else{
        ttest = t.test(CPJDOD_cluster[,j], mu = DrugBank_CPJDODmean[j-1,])
        pvalue_temp[j-1,1] = ttest$p.value
      }
    }
  }
  pvalue = cbind(pvalue, pvalue_temp)
  rm(CPJDOD_cluster_mean_temp)
}
for(i in 1:nrow(pvalue)){
  for(j in 1:ncol(pvalue)){
    if(pvalue[i,j]>=0.05 | is.na(pvalue[i,j]) == TRUE){
      CPJDOD_cluster_mean[i,j] = 0
    }
  }
}
CPJDOD_cluster_meanandpvlaue = rbind(CPJDOD_cluster_mean,pvalue)
write.csv(CPJDOD_cluster_meanandpvlaue,file = "F:/R Workspace/Output files/CPJDOD_cluster_meanandpvalue.csv")
for(i in 1:nrow(pvalue)){
  for(j in 1:ncol(pvalue)){
    if(pvalue[i,j]>=0.05 | is.na(pvalue[i,j]) == TRUE){
      CPJDOD_cluster_mean[i,j] = 0
    }else{
      CPJDOD_cluster_mean[i,j] = (CPJDOD_cluster_mean[i,j][[1]] - DrugBank_CPJDODmean[i,1][[1]])/abs(DrugBank_CPJDODmean[i,1][[1]])
    }
  }
}
CPJDOD_cluster_meanandpvlaue = rbind(CPJDOD_cluster_mean,pvalue)
write.csv(CPJDOD_cluster_meanandpvlaue,file = "F:/R Workspace/Output files/CPJDOD_cluster_meanandpvalue(fluctuation).csv")

