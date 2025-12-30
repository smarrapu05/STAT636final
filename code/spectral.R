rm(list=ls())
library(irlba)

Amb <- as.matrix(read.csv("C:/Users/Saketh/Documents/STAT636/Project/adjMb.csv",header=FALSE))
Ahb <- as.matrix(read.csv("C:/Users/Saketh/Documents/STAT636/Project/adjHb.csv",header=FALSE))
sum(Amb != Ahb)


#Spectral Embedding
d <- 3

eigMb <- eigen(Amb)
Xmb <- eigMb$vectors[, 1:d] %*% diag(sqrt(abs(eigMb$values[1:d])))

eigHb <- eigen(Ahb)
Xhb <- eigHb$vectors[,1:d] %*% diag(sqrt(abs(eigHb$values[1:d])))

M <- cbind(Xmb,Xhb)
svdRes <- svd(M)
Vcommon <- svdRes$u[,1:d]

Rmb <- Amb %*% Vcommon
Rhb <- Ahb %*% Vcommon

diffs <- sqrt(rowSums((Rmb - Rhb)^2))

plot(diffs, type="l", col="darkblue", lwd=2,
     main="Difference in Spectral Embedding: Hb vs Mb",
     xlab="Residue Index", 
     ylab="Euclidean Distance in Latent Space")

mostdiff <- order(diffs, decreasing=TRUE)[1:5]
cat("The most structurally divergent residues are at indices:", mostdiff, "/n")

mapinfo <- read.csv("C:/Users/Saketh/Documents/STAT636/Project/nodeMapping.csv")
mapinfo[mostdiff,]
mapinfo$diff <- diffs

#spectral clustering
set.seed(1)
k <- 5
kmeansRes <- kmeans(Rhb,centers=k,nstart=25)
plot(Rhb[,1],Rhb[,2], col=kmeansRes$cluster, pch=19,
     main="Spectral Clusters of Hemoglobin",
     xlab = "Latent Vector 1", ylab = "Latent Vector 2")
plot(kmeansRes$cluster, type="p", pch="|", cex=1.5, col=kmeansRes$cluster,
     main="Clustering along the Sequence Hemoglobin",
     xlab="Residue Index", ylab="Cluster ID", yaxt="n")

print(data.frame(
  Index = mostdiff,
  Cluster = kmeansRes$cluster[mostdiff]
))

mostOutlierCluster <- names(which.max(table(kmeansRes$cluster[mostdiff])))
cat(paste("The most different cluster for Hemoglobin appears to be cluster", mostOutlierCluster, "/n"))
outliermembers <- which(kmeansRes$cluster == mostOutlierCluster)
print(outliermembers)

mapinfo$kmeansclusterhb <- kmeansRes$cluster

kmeansResmb <- kmeans(Rmb,centers=k,nstart=25)
plot(Rmb[,1],Rmb[,2], col=kmeansResmb$cluster, pch=19,
     main="Spectral Clusters of Myoglobin",
     xlab = "Latent Vector 1", ylab = "Latent Vector 2")
plot(kmeansResmb$cluster, type="p", pch="|", cex=1.5, col=kmeansResmb$cluster,
     main="Clustering along the Sequence Myoglobin",
     xlab="Residue Index", ylab="Cluster ID", yaxt="n")

print(data.frame(
  Index = mostdiff,
  Cluster = kmeansResmb$cluster[mostdiff]
))

mostOutlierCluster <- names(which.max(table(kmeansResmb$cluster[mostdiff])))
cat(paste("The most different cluster for Myoglobin appears to be cluster", mostOutlierCluster, "/n"))
outliermembers <- which(kmeansResmb$cluster == mostOutlierCluster)
print(outliermembers)

mapinfo$kmeansclustermb <- kmeansResmb$cluster
write.csv(mapinfo, "C:/Users/Saketh/Documents/STAT636/Project/nodeMappingreal.csv")

#Visualizations
library(plotly)
library(igraph)
dfMb <- data.frame(
  x = Rmb[,1],y = Rmb[,2], z = Rmb[,3],
  Cluster = as.factor(mapinfo$kmeansclusterhb),
  Residue = paste0(mapinfo$Mb_Residue,mapinfo$Mb_PDB_ID)
)
dfHb <- data.frame(
  x = Rhb[,1],y = Rhb[,2], z = Rhb[,3],
  Cluster = as.factor(mapinfo$kmeansclusterhb),
  Residue = paste0(mapinfo$Hb_Residue,mapinfo$Hb_PDB_ID)
)

gMb <- graph_from_adjacency_matrix(Amb, mode="undirected",diag=FALSE)
layoutMb <- layout_with_fr(gMb)
png("C:/Users/Saketh/Documents/STAT636/Project/MyoglobinNetwork.png",width=800,height=800)
plot(gMb, layout=layoutMb, 
     vertex.size=4, 
     vertex.label=NA, 
     vertex.color="skyblue",
     edge.color=adjustcolor("darkgrey", alpha.f=0.6),
     main="Myoglobin Network")
dev.off()

gHb <- graph_from_adjacency_matrix(Ahb, mode="undirected",diag=FALSE)
layoutHb <- layout_with_fr(gHb)
png("C:/Users/Saketh/Documents/STAT636/Project/HemoglobinNetwork.png",width=800,height=800)
plot(gHb, layout=layoutHb, 
     vertex.size=4, 
     vertex.label=NA, 
     vertex.color="salmon",
     edge.color=adjustcolor("darkgrey", alpha.f=0.6),
     main="Hemoglobin Network")
dev.off()

fig_joint <- plot_ly() %>%
  add_markers(data = dfMb, x = ~x, y = ~y, z = ~z, 
              color = I("skyblue"), name = "Myoglobin",
              text = ~Residue, opacity = 0.5, 
              marker = list(size=5)) %>%
  add_markers(data = dfHb, x = ~x, y = ~y, z = ~z, 
              color = I("salmon"), name = "Hemoglobin",
              text = ~Residue, opacity = 0.6, 
              marker = list(size=5)) %>%
  layout(title = "Joint Spectral Embedding: Latent Space Alignment",
         scene = list(xaxis = list(title = "Latent Dim 1"),
                      yaxis = list(title = "Latent Dim 2"),
                      zaxis = list(title = "Latent Dim 3")))
print(fig_joint)


mbClust <- plot_ly(data = dfMb, x = ~x, y = ~y, z = ~z,
                        color = ~Cluster, colors = "viridis",
                        text = ~Residue, type = "scatter3d", mode = "markers",
                        marker = list(size=5)) %>%
  layout(title = "Myoglobin: Spectral Clusters",
         scene = list(aspectmode='cube'))
print(mbClust)
hbClust <- plot_ly(data = dfHb, x = ~x, y = ~y, z = ~z,
                        color = ~Cluster, colors = "viridis",
                        text = ~Residue, type = "scatter3d", mode = "markers",
                        marker = list(size=5)) %>%
  layout(title = "Hemoglobin: Spectral Clusters",
         scene = list(aspectmode='cube'))
print(hbClust)

#Descriptive Statistics
getNetworkStats <- function(g,name){
  numNodes <- vcount(g)
  numEdges <- ecount(g)
  dens <- edge_density(g)
  avgDeg <- mean(degree(g))
  trans <- transitivity(g,type="global")
  avgpath <- mean_distance(g,directed = FALSE)
  diam <- diameter(g, directed=FALSE)
  return(data.frame(
    Network = name,
    Nodes = numNodes,
    Edges = numEdges,
    Density <- dens,
    Average_Degree <- avgDeg,
    Transitivity <- trans,
    Average_Path <- avgpath,
    Diameter <- diam))
}
statsMb <- getNetworkStats(gMb,"Myoglobin")
statsHb <- getNetworkStats(gHb,"Hemoglobin")
finalTable <- rbind(statsMb,statsHb)
write.csv(finalTable,"C:/Users/Saketh/Documents/STAT636/Project/NetworkStats.csv",row.names=FALSE)
