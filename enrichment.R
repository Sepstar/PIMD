setwd("F:/R Workspace/")
load("F:/R Workspace/Input files/affni_fusion.Rdata")
C = 32 # number of clusters
cluster_group = spectralClustering(affni_fusion, C)
druglist = as.data.frame(rownames(affni_fusion))
## superclass enrichment
CC = read.csv("F:/R Workspace/Input files/superclass.csv", header = T)
N = 593
superclass = c("Organic Acids and Derivatives","Benzenoids","Organonitrogen Compounds",
               "Homogeneous Non-metal Compounds","Heterocyclic Compounds","Organophosphorus Compounds",
               "Phenylpropanoids and Polyketides","Organooxygen Compounds","Mixed Metal/Non-metal Compounds",
               "Lipids","Alkaloids and Derivatives","Lignans and Norlignans","Not Available")
pvalues = matrix(data = NA, nrow = 13, ncol = C)
odds = matrix(data = NA, nrow = 13, ncol = C)
for(i in 1:C){
  cluster_index = which(cluster_group == i)
  cluster_drug = druglist[cluster_index,]
  index = as.character(cluster_drug)
  n = length(index)
  club = data.frame()
  for(j in index){
    for(k in 1:length(CC[,1])){
      if(j == CC[k,1]){
        club = rbind(club, CC[k,])
      }
    }
  }
  CCindex = 0
  for(l in superclass){
    CCindex = CCindex + 1
    ktemp = c()
    for(q in 1:length(club[,8])){
      if(l == club[q,8]){
        ktemp = c(ktemp, club[q,1])
      }
    }
    mtemp = c()
    k = length(unique(ktemp))
    for(p in 1:length(CC[,1])){
      if(l == CC[p,8]){
        mtemp = c(mtemp, CC[p,1])
      }
    }
    m = length(unique(mtemp))
    odds[CCindex,i] = (k*N)/(n*m)
    pvalues[CCindex,i] = phyper(k-1, m, N-m, n, lower.tail=FALSE)
  }
}
for(i in 1:nrow(pvalues)){
  for(j in 1:ncol(pvalues)){
    if(pvalues[i,j]>=0.05){
      # pvalues[i,j] = 0
      odds[i,j] = 0
    }
  }
}
rownames(pvalues) = c("Organic Acids and Derivatives","Benzenoids","Organonitrogen Compounds",
                      "Homogeneous Non-metal Compounds","Heterocyclic Compounds","Organophosphorus Compounds",
                      "Phenylpropanoids and Polyketides","Organooxygen Compounds","Mixed Metal/Non-metal Compounds",
                      "Lipids","Alkaloids and Derivatives","Lignans and Norlignans","Not Available")
rownames(odds) = c("Organic Acids and Derivatives","Benzenoids","Organonitrogen Compounds",
                   "Homogeneous Non-metal Compounds","Heterocyclic Compounds","Organophosphorus Compounds",
                   "Phenylpropanoids and Polyketides","Organooxygen Compounds","Mixed Metal/Non-metal Compounds",
                   "Lipids","Alkaloids and Derivatives","Lignans and Norlignans","Not Available")
write.csv(pvalues,file = "F:/R Workspace/Output files/superclass enrichment pvalues.csv")
write.csv(odds,file = "F:/R Workspace/Output files/superclass enrichment odds.csv")

## ADMET enrichment
ADMET = read.csv("F:/R Workspace/Input files/ADMET.csv", header = T)
N = 593
HIA = c("+", "-", "Not Available")
BBB = c("+", "-", "Not Available")
Caco = c("+", "-", "Not Available")
Psub = c("Non-substrate", "Substrate", "Not Available")
Pinhi1 = c("Non-inhibitor", "inhibitor", "Not Available")
Pinhi2 = c("Non-inhibitor", "inhibitor", "Not Available")
Renal = c("Non-inhibitor", "inhibitor", "Not Available")
CYP4502C9 = c("Non-inhibitor", "inhibitor", "Not Available")
CYP4502D6 = c("Non-inhibitor", "inhibitor", "Not Available")
CYP4503A4 = c("Non-inhibitor", "inhibitor", "Not Available")
CYP4501A2 = c("Non-inhibitor", "inhibitor", "Not Available")
CYP4502C1 = c("Non-inhibitor", "inhibitor", "Not Available")
CYP450inhi = c("Low CYP Inhibitory Promiscuity", "High CYP Inhibitory Promiscuity", "Not Available")
Amestest = c("Non AMES toxic", "AMES toxic", "Not Available")
Carcinogenicity = c("Non-carcinogens", "Carcinogens", "Not Available")
Biodegradation = c("Ready biodegradable", "Not ready biodegradable", "Not Available")
hERGinhi1 = c("Weak inhibitor", "Strong inhibitor", "Not Available")
hERGinhi2 = c("Non-inhibitor", "Inhibitor", "Not Available")
pvalues = matrix(data = NA, nrow = 54, ncol = C)
odds = matrix(data = NA, nrow = 54, ncol = C)
for(i in 1:C){
  cluster_index = which(cluster_group == i)
  cluster_drug = druglist[cluster_index,]
  index = as.character(cluster_drug)
  n = length(index)
  club = data.frame()
  for(j in index){
    for(k in 1:length(ADMET[,1])){
      if(j == ADMET[k,1]){
        club = rbind(club, ADMET[k,])
      }
    }
  }
  termcount = 5
  admetindex = 0
  termsall =  rbind(HIA,BBB,Caco,Psub,Pinhi1,Pinhi2,Renal,CYP4502C9,CYP4502D6,CYP4503A4,CYP4501A2,CYP4502C1,CYP450inhi,Amestest,Carcinogenicity,Biodegradation,hERGinhi1,hERGinhi2)
  for(term in 1:nrow(termsall)){
    for(l in termsall[term,]){
      admetindex = admetindex + 1
      ktemp = c()
      for(q in 1:length(club[,termcount])){
        if(l == club[q,termcount]){
          ktemp = c(ktemp, club[q,1])
        }
      }
      mtemp = c()
      k = length(unique(ktemp))
      for(p in 1:length(ADMET[,1])){
        if(l == ADMET[p,termcount]){
          mtemp = c(mtemp, ADMET[p,1])
        }
      }
      m = length(unique(mtemp))
      odds[admetindex,i] = (k*N)/(n*m)
      pvalues[admetindex,i] = phyper(k-1, m, N-m, n, lower.tail=FALSE)
    }
    termcount = termcount + 1
  }
}
for(i in 1:nrow(pvalues)){
  for(j in 1:ncol(pvalues)){
    if(pvalues[i,j]>=0.05){
      # pvalues[i,j] = 0
      odds[i,j] = 0
    }
  }
}
rownames(pvalues) = c(HIA,BBB,Caco,Psub,Pinhi1,Pinhi2,Renal,CYP4502C9,CYP4502D6,CYP4503A4,CYP4501A2,CYP4502C1,CYP450inhi,Amestest,Carcinogenicity,Biodegradation,hERGinhi1,hERGinhi2)
rownames(odds) = c(HIA,BBB,Caco,Psub,Pinhi1,Pinhi2,Renal,CYP4502C9,CYP4502D6,CYP4503A4,CYP4501A2,CYP4502C1,CYP450inhi,Amestest,Carcinogenicity,Biodegradation,hERGinhi1,hERGinhi2)
write.csv(pvalues,file = "F:/R Workspace/Output files/ADMET enrichment pvalues.csv")
write.csv(odds,file = "F:/R Workspace/Output files/ADMET enrichment odds.csv")

## target-based enrichment analysis
library(mygene)
library(clusterProfiler)
library(org.Hs.eg.db)
genesymbol = as.data.frame(read.csv("F:/R Workspace/Input files/genesymbol.csv", header = F))
druglist = as.data.frame(rownames(affni_fusion))
gene_group_list=NULL
for(i in 1:C){
  cluster_index = which(cluster_group == i)
  cluster_drug = druglist[cluster_index,]
  cluster_drug = as.character(cluster_drug)
  genesymbol_cluster = c()
  for(drug in cluster_drug){
    genesymbolindex = which(genesymbol[,1] == drug)
    genesymbol_cluster_temp = as.matrix(genesymbol[genesymbolindex,2])
    genesymbol_cluster = rbind(genesymbol_cluster, genesymbol_cluster_temp)
  }
  genesymbol_cluster = unique(genesymbol_cluster)
  temp=queryMany(genesymbol_cluster, scopes="symbol", fields="entrezgene", species="human")[["entrezgene"]]
  temp=temp[which(!is.na(temp))]
  gene_group_list=c(gene_group_list,list(temp))
}

gene_group_GOMF=NULL
gene_group_GOBP=NULL
gene_group_GOCC=NULL
gene_group_KEGG=NULL
for (i in 1:C)
{
  print(i)
  temp_MF=enrichGO(gene_group_list[[i]], 'org.Hs.eg.db', ont="MF", pvalueCutoff=0.05)@result[,c(1,2,7)]
  temp_CC=enrichGO(gene_group_list[[i]], 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.05)@result[,c(1,2,7)]
  temp_BP=enrichGO(gene_group_list[[i]], 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05)@result[,c(1,2,7)]
  temp_KEGG=enrichKEGG(gene_group_list[[i]], organism="human", pvalueCutoff=0.05)@result[,c(1,2,7)]
  gene_group_GOMF=c(gene_group_GOMF,list(temp_MF))
  gene_group_GOCC=c(gene_group_GOCC,list(temp_CC))
  gene_group_GOBP=c(gene_group_GOBP,list(temp_BP))
  gene_group_KEGG=c(gene_group_KEGG,list(temp_KEGG))
  rm(temp_MF)
  rm(temp_BP)
  rm(temp_CC)
  rm(temp_KEGG)
}

gene_all_GOBP=NULL
gene_all_GOCC=NULL
gene_all_GOMF=NULL
gene_all_KEGG=NULL
for (i in 1:C)
{
  gene_all_GOBP=rbind(gene_all_GOBP,gene_group_GOBP[[i]][,1:2])
  gene_all_GOCC=rbind(gene_all_GOCC,gene_group_GOCC[[i]][,1:2])
  gene_all_GOMF=rbind(gene_all_GOMF,gene_group_GOMF[[i]][,1:2])
  gene_all_KEGG=rbind(gene_all_KEGG,gene_group_KEGG[[i]][,1:2])
}
gene_all_GOBP=as.matrix(unique(gene_all_GOBP))
gene_all_GOCC=as.matrix(unique(gene_all_GOCC))
gene_all_GOMF=as.matrix(unique(gene_all_GOMF))
gene_all_KEGG=as.matrix(unique(gene_all_KEGG))

gene_GOBP_matrix=matrix(data = 1,nrow = nrow(gene_all_GOBP),ncol=C)
rownames(gene_GOBP_matrix)=gene_all_GOBP[,2]
gene_GOCC_matrix=matrix(data = 1,nrow = nrow(gene_all_GOCC),ncol=C)
rownames(gene_GOCC_matrix)=gene_all_GOCC[,2]
gene_GOMF_matrix=matrix(data = 1,nrow = nrow(gene_all_GOMF),ncol=C)
rownames(gene_GOMF_matrix)=gene_all_GOMF[,2]
gene_KEGG_matrix=matrix(data = 1,nrow = nrow(gene_all_KEGG),ncol=C)
rownames(gene_KEGG_matrix)=gene_all_KEGG[,2]
colnames(gene_GOBP_matrix)=paste("cluster_",c(1:C),sep = "")
colnames(gene_GOCC_matrix)=paste("cluster_",c(1:C),sep = "")
colnames(gene_GOMF_matrix)=paste("cluster_",c(1:C),sep = "")
colnames(gene_KEGG_matrix)=paste("cluster_",c(1:C),sep = "")

for (i in 1:C)
{
  print(i)
  temp=gene_group_GOBP[[i]]
  if (nrow(temp)!=0)
  {
    for (j in 1:nrow(temp))
    {
      gene_GOBP_matrix[which(temp[j,2]==gene_all_GOBP[,2]),i]=temp[j,3]
    }
  }
  temp=gene_group_GOCC[[i]]
  if (nrow(temp)!=0)
  {
    for (j in 1:nrow(temp))
    {
      gene_GOCC_matrix[which(temp[j,2]==gene_all_GOCC[,2]),i]=temp[j,3]
    }
  }
  temp=gene_group_GOMF[[i]]
  if (nrow(temp)!=0)
  {
    for (j in 1:nrow(temp))
    {
      gene_GOMF_matrix[which(temp[j,2]==gene_all_GOMF[,2]),i]=temp[j,3]
    }
  }
  temp=gene_group_KEGG[[i]]
  if (nrow(temp)!=0)
  {
    for (j in 1:nrow(temp))
    {
      gene_KEGG_matrix[which(temp[j,2]==gene_all_KEGG[,2]),i]=temp[j,3]
    }
  }
  rm(temp)
}

gene_GOBP_matrix=-log10(gene_GOBP_matrix)
gene_GOCC_matrix=-log10(gene_GOCC_matrix)
gene_GOMF_matrix=-log10(gene_GOMF_matrix)
gene_KEGG_matrix=-log10(gene_KEGG_matrix)

C = 32
colnames(gene_GOBP_matrix)=paste("cluster ",c(1:C),sep = "")
colnames(gene_GOCC_matrix)=paste("cluster ",c(1:C),sep = "")
colnames(gene_GOMF_matrix)=paste("cluster ",c(1:C),sep = "")
colnames(gene_KEGG_matrix)=paste("cluster ",c(1:C),sep = "")

gene_GOBP_matrix_4=gene_GOBP_matrix
gene_GOCC_matrix_4=gene_GOCC_matrix
gene_GOMF_matrix_4=gene_GOMF_matrix
gene_KEGG_matrix_4=gene_KEGG_matrix

gene_GOBP_matrix_4[which(gene_GOBP_matrix_4>4)]=4
gene_GOCC_matrix_4[which(gene_GOCC_matrix_4>4)]=4
gene_GOMF_matrix_4[which(gene_GOMF_matrix_4>4)]=4
gene_KEGG_matrix_4[which(gene_KEGG_matrix_4>4)]=4
gene_GOBP_matrix_4=gene_GOBP_matrix_4[,which(colSums(gene_GOBP_matrix_4)!=0)]
gene_GOCC_matrix_4=gene_GOCC_matrix_4[,which(colSums(gene_GOCC_matrix_4)!=0)]
gene_GOMF_matrix_4=gene_GOMF_matrix_4[,which(colSums(gene_GOMF_matrix_4)!=0)]
gene_KEGG_matrix_4=gene_KEGG_matrix_4[,which(colSums(gene_KEGG_matrix_4)!=0)]

pheatmap(gene_GOBP_matrix_4,color=colorRampPalette(c("GhostWhite","firebrick3"))(1000),cluster_cols=T,cluster_rows=T,border_color=NA,fontsize=10,show_rownames = F,show_colnames = T
         ,fontsize_col=10)
pheatmap(gene_GOCC_matrix_4,color=colorRampPalette(c("GhostWhite","firebrick3"))(1000),cluster_cols=T,cluster_rows=T,border_color=NA,fontsize=10,show_rownames = F,show_colnames = T
         ,fontsize_col=10)
pheatmap(gene_GOMF_matrix_4,color=colorRampPalette(c("GhostWhite","firebrick3"))(1000),cluster_cols=T,cluster_rows=T,border_color=NA,fontsize=10,show_rownames = F,show_colnames = T
         ,fontsize_col=10)
pheatmap(gene_KEGG_matrix_4,color=colorRampPalette(c("GhostWhite","firebrick3"))(1000),cluster_cols=T,cluster_rows=T,border_color=NA,fontsize=10,show_rownames = T,show_colnames = T
         ,fontsize_row=3,fontsize_col=10)