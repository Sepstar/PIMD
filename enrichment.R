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
