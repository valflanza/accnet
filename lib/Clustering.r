library(dplyr)
library(tidyr)
library(mclust)
library(cluster)


setwd(runPath)
Net = read.table(TableFile, header = TRUE)
Net.matrix = Net %>% group_by(Source, Target) %>% summarise(value = 1) %>% spread(Source, value, fill = 0)
tmp = Net.matrix$Target
Net.matrix = Net.matrix[, -1]
rownames(Net.matrix) = tmp
rm(tmp)


### Performing GU clustering


Net.dist = dist(Net.matrix, method = "binary")

## Bayesian clustering
Net.mclust = Mclust(as.matrix(Net.dist))
Net.mclust.cluster = as.data.frame(Net.mclust$classification)
Net.mclust.cluster$ID = rownames(Net.mclust.cluster)
colnames(Net.mclust.cluster) = c("MclustCluster", "ID")

## Hierarchical clustering

Net.GU.dend = hclust(Net.dist, method = "average")

qnames = c(0.75, 0.85, 0.90, 0.95, 0.99)
q = quantile(Net.GU.dend$height, probs = qnames)
j = 1
Net.GU.tmp = list()
for (i in q)
{
  Net.GU.tmp[[j]] = as.data.frame(cutree(Net.GU.dend, h = i))
  
  colnames(Net.GU.tmp[[j]]) = c("Ctmp")
  Net.GU.tmp[[j]]$tmp = "GU"
  Net.GU.tmp[[j]] = Net.GU.tmp[[j]] %>% unite(Cluster, Ctmp, tmp, sep =
                                                "")
  Net.GU.tmp[[j]]$ID = rownames(Net.GU.tmp[[j]])
  colnames(Net.GU.tmp[[j]]) = c(paste(c("GU", qnames[j]), collapse = ""), "ID")
  
  
  if (j > 1)
  {
    Net.GU.cluster = full_join(Net.GU.cluster, Net.GU.tmp[[j]])
  } else{
    Net.GU.cluster = Net.GU.tmp[[j]]
  }
  j = j + 1
  
}
rm(Net.GU.tmp)
Net.GU.cluster = Net.GU.cluster[sort(colnames(Net.GU.cluster))]


### Performing Protein clustering

Net.prot.matrix = Net %>% group_by(Source) %>% mutate(rep = n()) %>% filter(rep >
                                                                              1) %>% group_by(Source, Target) %>% summarise(value = 1) %>% spread(Target, value, fill = 0)
tmp = Net.prot.matrix$Source
Net.prot.matrix = Net.prot.matrix[, -1]
rownames(Net.prot.matrix) = tmp
rm(tmp)

Net.prot.dist = dist(Net.prot.matrix, method = "binary")
Net.prot.dend = hclust(Net.prot.dist)
qnames = c(0.75, 0.85, 0.90, 0.95, 0.99)
q = quantile(Net.prot.dend$height, probs = qnames)
j = 1
Net.prot.tmp = list()
for (i in q)
{
  Net.prot.tmp[[j]] = as.data.frame(cutree(Net.prot.dend, h = i))
  
  colnames(Net.prot.tmp[[j]]) = c("Ctmp")
  Net.prot.tmp[[j]]$tmp = "P"
  Net.prot.tmp[[j]] = Net.prot.tmp[[j]] %>% unite(Cluster, Ctmp, tmp, sep =
                                                    "")
  Net.prot.tmp[[j]]$ID = rownames(Net.prot.tmp[[j]])
  colnames(Net.prot.tmp[[j]]) = c(paste(c("PT", qnames[j]), collapse = ""), "ID")
  
  
  if (j > 1)
  {
    Net.prot.cluster = full_join(Net.prot.cluster, Net.prot.tmp[[j]])
  } else{
    Net.prot.cluster = Net.prot.tmp[[j]]
  }
  j = j + 1
  
}
rm(Net.prot.tmp)
Net.prot.cluster = Net.prot.cluster[sort(colnames(Net.prot.cluster))]

### Summarising outputs

Net.Clusters = full_join(Net.prot.cluster, Net.GU.cluster)
Net.Clusters = full_join(Net.Clusters, Net.mclust.cluster)
Net.Clusters[is.na(Net.Clusters)] = ""
write.table(
  Net.Clusters,
  "Cluster.csv",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

