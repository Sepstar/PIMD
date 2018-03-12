setwd("F:/R Workspace/")
library(SNFtool)
## first, input 3 drug similarity networks
side = as.matrix(read.csv("F:/R Workspace/Input files/1_side effect.csv", header = T, row.names = 1))
chem = as.matrix(read.csv("F:/R Workspace/Input files/2_chem.csv", header = T, row.names = 1))
target = as.matrix(read.csv("F:/R Workspace/Input files/3_target.csv", header = T, row.names = 1))
## set the parameters
K = 20; # number of neighbors, usually (10~30)
alpha = 0.5; # hyperparameter, usually (0.3~0.8)
T = 20; # Number of Iterations, usually (10~20)
affini_side = affinityMatrix(1-side, K, alpha)
affini_chem = affinityMatrix(1-chem, K, alpha)
affini_target = affinityMatrix(1-target, K, alpha)
affini_fusion = SNF(list(affini_side,affini_chem,affini_target), K, T)
rownames(affini_fusion) = rownames(side)
colnames(affini_fusion) = rownames(side)
## ATC enrichment analysis
C = 32 # number of clusters
group = spectralClustering(affini_fusion, C)
CID = rownames(chem)
CIDandgroup = cbind(CID,group)
N = 593 # number of all drugs
CIDandatc = read.csv("F:/R Workspace/Input files/drugs-ATC.csv", header = F)
atc = c("0", "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")
odds = matrix(data = NA, nrow = C, ncol = 15)
pvalues = matrix(data = NA, nrow = C, ncol = 15)
for(i in 1:C){
  temp = which(CIDandgroup[,2] == i)
  index = as.numeric(CIDandgroup[temp, 1])
  n = length(index)# number of drugs in the ith community
  club = data.frame()
  for(j in index){
    for(k in 1:length(CIDandatc[,1])){
      if(j == CIDandatc[k,1]){
        club = rbind(club, CIDandatc[k,])
      }
    }
  }
  atcindex = 0
  for(l in atc){
    atcindex = atcindex + 1
    ktemp = c()
    for(q in 1:length(club[,2])){
      if(l == club[q,2]){
        ktemp = c(ktemp, club[q,1])
      }
    }
    mtemp = c()
    k = length(unique(ktemp))# number of drugs with specific ATC code in the ith community
    for(p in 1:length(CIDandatc[,1])){
      if(l == CIDandatc[p,2]){
        mtemp = c(mtemp, CIDandatc[p,1])
      }
    }
    m = length(unique(mtemp))# number of drugs with specific ATC code in the whole network
    odds[i,atcindex] = (k*N)/(n*m)
    pvalues[i,atcindex] = phyper(k-1, m, N-m, n, lower.tail=FALSE)
  }
}
rm(temp, ktemp, mtemp, index, atcindex, club, n, k, m)
## filter odds with p-value smaller than 0.05
for(i in 1:nrow(pvalues)){
  for(j in 1:ncol(pvalues)){
    if(pvalues[i,j]>=0.05){
      odds[i,j] = 0
    }
  }
}
oddsandpvlaue = rbind(odds,pvalues)
colnames(oddsandpvlaue) = c("0", "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")
write.csv(oddsandpvlaue,file = "F:/R Workspace/Output files/ATC enrichment.csv")

## NMIs among drug similarity networks
affini_fusion_sidechem = SNF(list(affini_side, affini_chem), K, T)
affini_fusion_sidetarget = SNF(list(affini_side, affini_target), K, T)
affini_fusion_chemtarget = SNF(list(affini_chem, affini_target), K, T)
combine = list(affini_fusion, affini_fusion_sidechem, affini_fusion_sidetarget, affini_fusion_chemtarget, side, chem, target)
Comcordance_matrix = concordanceNetworkNMI(combine, 32)
rownames(Comcordance_matrix) = c("affini_fusion", "affini_fusion_sidechem", "affini_fusion_sidetarget", "affini_fusion_chemtarget", "side", "chem", "target")
colnames(Comcordance_matrix) = c("affini_fusion", "affini_fusion_sidechem", "affini_fusion_sidetarget", "affini_fusion_chemtarget", "side", "chem", "target")
write.csv(Comcordance_matrix, file = "F:/R Workspace/Output files/Comcordance matrix.csv")

## data type contribution
affini_list = NULL
affini_list = c(affini_list, list(affini_side))
affini_list = c(affini_list, list(affini_chem))
affini_list = c(affini_list, list(affini_target))
contribution_matrix = affini_fusion
contribution_matrix[which(contribution_matrix != 0)] = 0
for (i in 1:nrow(contribution_matrix))
{
  for (j in i:nrow(contribution_matrix))
  {
    if (j != i)
    {
      temp = matrix(data = 0,nrow = 3,ncol = 1)
      temp[1] = affini_list[[1]][i,j]
      temp[2] = affini_list[[2]][i,j]
      temp[3] = affini_list[[3]][i,j]
      rank_temp = order(-temp)
      temp = temp[rank_temp]
      diff = (temp[1]-temp[2])/temp[2]
      diff = c(diff,(temp[2]-temp[3])/temp[3])
      if (diff[1]>0.1)# if the first was more than 10% higher than the second, the edge was attributed to the first data type.
      {
        contribution_matrix[i,j] = rank_temp[1]
      }else# if the difference between the first and the second is less than 10%, then look at the difference between the second and the third
      {
        if (diff[2]>0.1)# if the difference between the second and the third is greater than 10%, the edge was attributed to both the first and the second data type
        {
          a = c(rank_temp[1],rank_temp[2])
          if (length(intersect(a,c(1,2)))==2)
          {
            contribution_matrix[i,j] = 4
          }
          if (length(intersect(a,c(1,3)))==2)
          {
            contribution_matrix[i,j] = 5
          }
          if (length(intersect(a,c(2,3)))==2)
          {
            contribution_matrix[i,j] = 6
          }
        }else# if the difference between the second and the third is less than 10%, the edge was attributed to all three data types
        {
          contribution_matrix[i,j] = 7
        }
      }
    }
  }
}
contribution_sta = table(contribution_matrix[which(contribution_matrix!=0)])
pie_col=c("skyblue","green","pink","yellow","orange","red","bisque","burlywood","cadetblue","gray","chocolate","cornsilk","firebrick","khaki","black")
pie(contribution_sta,border = NA,col = pie_col,clockwise = T,cex=0.8)
