# Microbiome-network-analysis-SparCC
This is a R script to analysis micrbiome OTU/ASV dataset in order to build a network.
We used SparCC to do this analsysis


getwd()
# Loading required packages
library(igraph)
#library(Hmisc) (for correlation)
library(Matrix)
# Install SpiecEasi package
install.packages("devtools")
library(devtools)
install_github("zdk123/SpiecEasi", force = TRUE)
devtools::install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(sparcc)


# Load the data with the OTU table: otudata.csv
otu.table<-read.csv(file.choose(), header=T, row.names = 1)

# Read taxonomy file associated with OTU table into new object: otu_taxonomy.csv
tax<-read.csv(file.choose(),header=T, row.names = 1, check.name = T)

# Check how many OTUs we have
dim(otu.table)

# Keep the OTUs with more than 10 counts
otu.table.filter<-otu.table[,colSums(otu.table)>10]

# Check for the decrease in the number of OTUs
dim(otu.table.filter)
head(otu.table.filter)
# SparCC network (This will take some time to run depending on size of dataset: minutes to hours)
sparcc.matrix <- sparcc(otu.table.filter)



head(sparcc.matrix)
cor=sparcc.matrix$Cor # select all the correlation values
cor
sparcc.cutoff <- 0.5
sparcc.adj <- ifelse(abs(sparcc.matrix$Cor) >= sparcc.cutoff, 1, 0)

# Add OTU names to rows and columns
rownames(sparcc.adj) <- colnames(otu.table.filter)
colnames(sparcc.adj) <- colnames(otu.table.filter)

# Create an adjacency matrix in igraph format
net.grph=graph.adjacency(sparcc.adj, mode = "undirected",weighted=TRUE,
                         diag=FALSE)

# Calculate edge weight == level of correlation
edgew<-E(net.grph)$weight

# Identify isolated nodes
bad.vs<-V(net.grph)[degree(net.grph) == 0]

# Remove isolated nodes
net.grph <-delete.vertices(net.grph, bad.vs)

# Plot the graph object
plot(net.grph,
     vertex.size=4,
     vertex.frame.color="black",
     edge.curved=F,
     edge.width=1.5,
     layout=layout.auto,
     edge.color=ifelse(edgew<0,"red","blue"),
     vertex.label=NA,
     vertex.label.color="black",
     vertex.label.family="Times New Roman",
     vertex.label.font=2)

# Plot the graph object
plot(net.grph,
     vertex.size=5,
     vertex.frame.color="black",
     edge.curved=F,
     edge.width=1.5,
     edge.color=ifelse(edgew<0,"red","blue"),
     vertex.label.color="blue",
     vertex.label.family="Times New Roman",
     vertex.label.font=0.5)

library(igraph)

# Assume net.grph is your network graph object
net <- net.grph

# Perform Walktrap community detection
wt <- walktrap.community(net)

# Get cluster memberships
membership <- membership(wt)

# Get the size of each cluster
cluster_sizes <- sizes(wt)

# Identify the 5 largest clusters
top_clusters <- order(cluster_sizes, decreasing = TRUE)[1:5]

# Assign colors to the 5 largest clusters
cluster_colors <- rainbow(5)

# Color nodes based on cluster membership
V(net)$color <- "gray"  # Default color for smaller clusters
for (i in 1:5) {
  V(net)$color[membership == top_clusters[i]] <- cluster_colors[i]
}

# Calculate hub scores for keystone taxa identification
net_hs <- hub_score(net)$vector

# Highlight keystone taxa by increasing their size
V(net)$size <- 5  # Default size
for (i in 1:5) {
  cluster_nodes <- V(net)[membership == top_clusters[i]]
  keystone_node <- cluster_nodes[which.max(net_hs[cluster_nodes])]
  V(net)$size[keystone_node] <- 10  # Double the size for keystone taxa
}

# Plot the graph
plot(net,
     vertex.frame.color = "black",
     edge.curved = FALSE,
     edge.width = 1.5,
     edge.color = ifelse(edgew < 0, "red", "blue"),
     vertex.label.color = "blue",
     vertex.label.family = "Times New Roman",
     vertex.label.font = 0.5,
     main = "Network with Clusters and Keystone Taxa")

# Add legend
legend("topright", legend = paste("Cluster", 1:5), col = cluster_colors, pch = 19, pt.cex = 1.5)

tax <- read.csv(file.choose(), header = TRUE, check.names = TRUE)
library(dplyr)
taxonomy_data <- tax %>%
  mutate(Names = paste(Family, Genus, sep = "_"))


View(taxonomy_data)
# Assuming the taxonomy data contains two columns: "OTU_ID" and "Taxonomy"
# "OTU_ID" should match with the node names in your network

names(taxonomy_data)[1] <- "OTU"
# Create a mapping from OTU to Taxonomy
taxonomy_mapping <- setNames(taxonomy_data$Names, taxonomy_data$OTU)

# Replace node names with taxonomy names
node_names <- V(net)$name
# Replace node names with taxonomy names using the mapping
# Replace node names with taxonomy names using the mapping
taxonomy_names <- sapply(node_names, function(node) {
  if (node %in% names(taxonomy_mapping)) {
    return(taxonomy_mapping[[node]])
  } else {
    return(node)  # Keep the original node name if not found in the mapping
  }
})

# Assign taxonomy names to the nodes in the network
V(net)$name <- taxonomy_names

# Now, your network nodes have been replaced with taxonomy names

plot(net)
     vertex.frame.color = "black",
     edge.curved = FALSE,
     edge.width = 1.5,
     edge.color = ifelse(edgew < 0, "red", "blue"),
     vertex.label.color = "blue",
     vertex.label.family = "Times New Roman",
     vertex.label.font = 0.5,
     main = "Network with Clusters and Keystone Taxa")

plot(net)
