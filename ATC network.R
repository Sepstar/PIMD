setwd("F:/R Workspace/")
edge_all = read.csv("F:/R Workspace/Input files/all edges.csv", header = F)
edge_all = edge_all[order(edge_all[,3], decreasing = T),]
top = list()
for(i in 1:15){
  temp = 175528*i/100
  top[[i]] = edge_all[1:temp,]
}
## mean value
edge = matrix(data = NA, nrow = 15*14/2, ncol = 3)
atc = c("0", "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")
datasource = top[[5]]
k = 1
for(i in 1:(length(atc)-1)){
  for(j in (i+1):length(atc)){
    temp_fir = datasource[which(datasource[,4]==atc[i]),]
    temp_fir_sec = temp_fir[which(temp_fir[,5]==atc[j]),]
    temp_sec = datasource[which(datasource[,4]==atc[j]),]
    temp_sec_fir = temp_sec[which(temp_sec[,5]==atc[i]),]
    temp_all = rbind(temp_fir_sec, temp_sec_fir)
    edge[k, 1] = atc[i]
    edge[k, 2] = atc[j]
    edge[k, 3] = mean(temp_all[,3])
    k = k + 1 
  }
}
edge[which(edge[,3]=="NaN"),3] = 0
write.csv(edge, file = "F:/R Workspace/Output files/edge(mean).csv")

## number of shared drugs
drug_atc = read.csv("F:/R Workspace/Input files/drugs-ATC.csv", header = F)
drug_atc = unique(drug_atc)
edge_mat = matrix(data = NA, nrow = 15, ncol = 15)
atc = c("0", "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")
temp = list()
for(i in 1:length(atc)){
  temp[[i]] = drug_atc[which(drug_atc[,2]==atc[i]),1]
}
for(i in 1:15){
  for(j in 1:15){
    edge_mat[i,j] = length(intersect(temp[[i]],temp[[j]]))
  }
}
rownames(edge_mat) = c("0", "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")
colnames(edge_mat) = c("0", "A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")

edge_list = matrix(data = NA, nrow = nrow(edge_mat)*(nrow(edge_mat)-1)/2, ncol = 3)
k=1
for (i in 1:15)
{
  for (j in i:15)
  {
    if (j!=i){
      edge_list[k,1] = rownames(edge_mat)[i]
      edge_list[k,2] = rownames(edge_mat)[j]
      edge_list[k,3] = edge_mat[i,j]
      k = k + 1
    }
  }
}
write.table(edge_list, file = "F:/R Workspace/Output files/edge(shared drugs).txt", quote = F, sep = "\t", row.names = F, col.names = F, na = "")


