# ---- Libraries ----
library(BiocGenerics) 
library(DESeq2)
library(psych) 
library(NetworkToolbox)
library(ggplot2)
library(GGally)
library(sna)
library(network)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(ggplot2)
library(Matrix)
library(igraph)
library(DESeq2)
library(maftools)
library(vegan)
library(enrichR,quietly = T)
library(aricode)
library(SNFtool)
library(biomaRt)
library(dplyr)

#1: --- DATA FROM TGCA ---

# Define project and create project directory
proj <- "TCGA-LUAD" 
dir.create(file.path(proj)) 

# Query for RNA-Seq data for primary tumor samples
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Primary Tumor")

# Uncomment the following line to download data
# GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api") 

# Prepare and extract data for primary tumor samples
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata") 
rna.expr.data.C <- assay(rna.data.C)

# View gene information for the primary tumor data
View(BiocGenerics::as.data.frame(rowRanges(rna.data.C))) 
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

# Query for RNA-Seq data for normal tissue samples
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling",  
                        data.type = "Gene Expression Quantification",  
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

# Uncomment the following line to download normal tissue data
# GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api") 

# Prepare and extract data for normal tissue samples
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata") 
rna.expr.data.N <- assay(rna.data.N) 
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.N)) 

# Check if gene information is consistent between primary tumor and normal tissue
all(na.omit(genes.info2) == na.omit(genes.info)) #TRUE

# Query clinical data for the project
clinical.query <- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)

# View clinical data
View(clinical.query)

# Check distribution of AJCC pathologic stage
table(clinical.query$ajcc_pathologic_stage)

# Boxplot of age at index by AJCC pathologic stage
boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query, 
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )

# View RNA-Seq data for normal samples
View(rna.expr.data.N)

# Print dimensions of the datasets
dim(rna.expr.data.C) 
dim(rna.expr.data.N) 

# Check the number of unique patients in clinical data
length(unique(clinical.query$submitter_id))

#1.2: --- Data cleaning ---

# Check the number of columns and patient IDs in the cancer and normal samples
ncol(rna.expr.data.C) 
head(colnames(rna.expr.data.C)) 
head(substr(colnames(rna.expr.data.N), 1,12))

# Check the dimensions and uniqueness of patient IDs in the normal dataset
dim(rna.expr.data.N) 
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #no duplicates 

# Check the dimensions and uniqueness of patient IDs in the cancer dataset
dim(rna.expr.data.C) 
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #duplicates! 539-516

# Extract unique patient IDs from the cancer data
patients.C <- substr(colnames(rna.expr.data.C), 1,12) 
sort(table(patients.C))

# Identify unique patients with only one sample
unique.patients.C <- names(which(table(patients.C) == 1))

# Get the index of unique patients in the dataset
idx.unique.pats <- match(unique.patients.C, substr(colnames(rna.expr.data.C), 1,12))

# Keep only the samples from unique patients
expr.C <- as.data.frame(rna.expr.data.C[, idx.unique.pats]) 
expr.N <- as.data.frame(rna.expr.data.N)

# Rename patients with a shorter ID
colnames(expr.C) <- substr(colnames(expr.C), 1,12) 
unique(colnames(expr.C))

colnames(expr.N) <- substr(colnames(expr.N), 1,12) 
unique(colnames(expr.N))

# Check common and unique patients between cancer and normal samples
intersect(colnames(expr.N), colnames(expr.C)) #52 instead of 59
setdiff(colnames(expr.N), colnames(expr.C)) 
# 7 normal samples do not have a cancerous sample to compare to. Let's remove them
match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N)) #idx to remove
expr.N <- expr.N[,-c(6, 20, 21, 32, 38, 45, 51)] 

# Check if normal samples are a subset of cancer samples
length(intersect(colnames(expr.N), colnames(expr.C))) #52

# Check the actual counts for cancer and normal data
typeof(expr.C[1,1]) #ok - integer
any(is.na(expr.C)) #ok - false
any(is.nan(as.matrix(expr.C))) #ok - false

typeof(expr.N[1,1]) #ok - integer
any(is.na(expr.N)) #ok - false
any(is.nan(as.matrix(expr.N))) #ok - false

# Consider only patients with both normal and cancer samples
expr.C <- expr.C[, colnames(expr.N)]

#1.3: --- Normalizing data with Deseq2 --- 

# Check if the row names of cancer and normal datasets are the same
all(rownames(expr.C) == rownames(expr.N)) 
full.data <- cbind(expr.N, expr.C)

# Print dimensions of the full data
dim(full.data) #104
full.data <- data.frame(full.data)

# Create metadata for cancer/normal condition
metad <- rep("cancer", 104) 
metad[1:52] <- "normal" 
metad 
metad <- data.frame(metad) 
rownames(metad) <- colnames(full.data) 
colnames(metad)[1] <- "condition" 
metad[,1] <- as.factor(metad[,1]) 
full.data <- cbind(rownames(full.data), full.data)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData=full.data,  
                              colData=metad,  
                              design= ~condition, 
                              tidy=TRUE)

# View the DESeq2 counts
View(counts(dds)) 
dim(counts(dds))

# Filter genes with at least ten counts in 90% of the samples
(104*90)/100 
keep <- rowSums(counts(dds) >= 10) >= 93 
dds <- dds[keep,] 
dim(counts(dds)) #15634 genes on 104 samples

# Normalize the counts
dds <- estimateSizeFactors(dds) 
normalized_counts <- counts(dds, normalized=TRUE) 
sum(rowSums(normalized_counts == 0) == 104) #no null rows

# Separate normalized counts into cancer and normal data
filtr.expr.n <- as.data.frame(normalized_counts[, 1:52]) 
filtr.expr.c <- as.data.frame(normalized_counts[, 53:104]) 
# Cancerous sample names were added a ".1" in full.data because  
# they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,52)

# 2: --- DIFFERENTIALLY EXPRESSED GENES (DEGs) ---

# Calculates the log2 fold change (FC) between cancer and normal conditions for each gene.
fc <- log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n))
names(fc) <- rownames(filtr.expr.c)

# Performs paired t-tests comparing cancer and normal conditions for each gene.
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) {
  t.test(as.numeric(filtr.expr.c[i, ]), as.numeric(filtr.expr.n[i, ]), paired = TRUE)$p.value
})

# Adjusts the p-values for multiple comparisons using the False Discovery Rate (FDR) method.
pval.fc.fdr <- p.adjust(pval.fc, method = "fdr")

# Combines fold change values and adjusted p-values into a single data frame.
expr.table <- data.frame(fc = round(fc, 2), pval.fc.fdr = pval.fc.fdr)

# Identifies differentially expressed genes (DEGs) that meet the specified thresholds.
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.95 & expr.table$pval.fc.fdr <=  0.01, ])

# Output the number of DEGs
cat("Number of DEGs:", length(deg.genes), "\n") #898

# The difference between overexpressed and underexpressed genes is based on their expression levels,
# i.e., how much a gene is activated to produce RNA and proteins in a given condition.
expr.table$diffexpressed <- "NO"
# A gene is overexpressed when its expression level is higher in one condition than in another,
# meaning the gene produces more RNA or proteins in cancer compared to healthy tissue.
expr.table$diffexpressed[expr.table$fc >= 1.95 & expr.table$pval.fc.fdr <= 0.01] <- "UP"
# A gene is underexpressed when its expression level is lower in one condition than in another,
# meaning the gene produces less RNA or proteins in cancer compared to healthy tissue.
expr.table$diffexpressed[expr.table$fc <= -1.95 & expr.table$pval.fc.fdr <= 0.01] <- "DOWN"

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
summary(expr.table$diffexpressed) 
# DOWN    NO    UP 
#  357 14736   541  

# Generates a volcano plot to visualize the results of differential expression analysis.
ggplot(data = expr.table, aes(x = fc, y = -log10(pval.fc.fdr), col = diffexpressed)) +  
  geom_point() +
  xlab("Fold Change (log2)") + 
  ylab("-log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.01), col = "red") +
  geom_vline(xintercept = 1.95, col = "red") +
  geom_vline(xintercept = -1.95, col = "red") +
  theme_minimal()

# A higher prevalence of overexpressed genes compared to underexpressed genes suggests a preferential
# activation of specific biological processes (e.g., cell proliferation, survival) in cancer.

# SUBSET normalized expression data for DEGs
# Extracts expression data for differentially expressed genes (DEGs) in cancer and normal conditions.
deg_expr_c <- filtr.expr.c[deg.genes, ]  # Cancer condition
deg_expr_n <- filtr.expr.n[deg.genes, ]  # Normal condition

# 2.1: --- Distinct subtypes based on gene expression data ---
# Connects to the Ensembl database.
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Converts ENSG without removing version numbers
# Maps Ensembl gene IDs (with versions removed) to external gene names.
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gsub("\\.\\d+$", "", rownames(deg_expr_c)),  # Remove versions just for the match
  mart = ensembl
)

# Combines original ENSG IDs with gene symbols
gene_mapping_with_versions <- data.frame(
  ENSG_with_versions = rownames(deg_expr_c),
  ENSG_without_versions = gsub("\\.\\d+$", "", rownames(deg_expr_c)
  ))

deg_symbols <- merge(
  gene_mapping_with_versions,
  gene_mapping,
  by.x = "ENSG_without_versions",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
head(deg_symbols)

# Creates a table of gene-patient-expression data with mapped gene symbols.
# Removes duplicates based on external_gene_name.
valid_genes_unique <- valid_genes[!duplicated(valid_genes$external_gene_name), ]
# Filters the expression matrix to retain only unique genes.
deg_expr_c_unique <- deg_expr_c[valid_genes_unique$ENSG_with_versions, , drop = FALSE]
# Replaces row names with gene symbols.
rownames(deg_expr_c_unique) <- valid_genes_unique$external_gene_name
head(deg_expr_c_unique)

# Proximal proliferative, proximal inflammatory, and terminal respiratory unit
# Marker genes for subtypes.
proliferative_markers <- c("KRAS", "EGFR", "MYC")
inflammatory_markers <- c("IL6", "CXCL10", "TNF")
tru_markers <- c("NKX2-1", "SFTPC", "NAPSIN")

# Finds the markers present in DEGs
deg_proliferative <- intersect(deg_symbols$external_gene_name, proliferative_markers)
deg_inflammatory <- intersect(deg_symbols$external_gene_name, inflammatory_markers)
deg_tru <- intersect(deg_symbols$external_gene_name, tru_markers)

# Outputs the identified markers for each category.
cat("Markers found in DEGs:\n")
cat("Proliferative:", paste(deg_proliferative, collapse = ", "), "\n")
cat("Inflammatory:", paste(deg_inflammatory, collapse = ", "), "\n")
cat("TRU:", paste(deg_tru, collapse = ", "), "\n")

# 3: --- CO-EXPRESSION NETWORKS ---
# Applies log2 transformation to normalize expression data for both cancer and normal conditions.
log_deg_expr_c <- log2(deg_expr_c + 1)
log_deg_expr_n <- log2(deg_expr_n + 1)

# Computes Pearson correlation coefficients for gene expression in cancer and normal datasets.
cor_cancer <- cor(t(log_deg_expr_c), method = "pearson")
cor_normal <- cor(t(log_deg_expr_n), method = "pearson")

# Thresholds correlation values to create binary adjacency matrices for cancer and normal networks.
threshold <- 0.7
adj_cancer <- ifelse(abs(cor_cancer) >= threshold, 1, 0)
adj_normal <- ifelse(abs(cor_normal) >= threshold, 1, 0)

# 3.1: --- Create Graph Objects ---
# Generates graph representations from the adjacency matrices for both networks.
graph_cancer <- graph_from_adjacency_matrix(adj_cancer, mode = "undirected", diag = FALSE)
graph_normal <- graph_from_adjacency_matrix(adj_normal, mode = "undirected", diag = FALSE)

# Output the adjacency matrix for inspection 
# Uncomment the following lines to inspect the adjacency matrices.
# View(adj_cancer)
# View(adj_normal)

# 3.2: --- Analyze Network Properties (Scale-Free Topology) ---
# Calculates the degree (number of connections) of each node in the networks.
deg_cancer <- igraph::degree(graph_cancer)
deg_normal <- igraph::degree(graph_normal)

# Check for zero-degree nodes
# Counts the number of isolated nodes (zero-degree) in both networks.
cat("Zero-degree nodes - Cancer Network:", sum(deg_cancer == 0), "\n")
cat("Zero-degree nodes - Normal Network:", sum(deg_normal == 0), "\n")

# Exclude zero-degree nodes
# Filters out nodes with zero-degree to focus on connected components.
deg_cancer_nonzero <- deg_cancer[deg_cancer > 0]
deg_normal_nonzero <- deg_normal[deg_normal > 0]
length(deg_normal_nonzero)
length(deg_cancer_nonzero)

# Plots histograms of node degree distributions for cancer and normal networks.
hist(deg_cancer_nonzero, breaks = 30, main = "Degree Distribution (Cancer Network)", 
     xlab = "Degree", col = "lightblue", border = "black")
hist(deg_normal_nonzero, breaks = 30, main = "Degree Distribution (Normal Network)", 
     xlab = "Degree", col = "lightgreen", border = "black")

# Compute R-squared values to assess linearity in the log-log plot
deg_cancer_freq <- table(deg_cancer_nonzero) / length(deg_cancer_nonzero)
deg_cancer_log <- log10(as.numeric(names(deg_cancer_freq)))
freq_cancer_log <- log10(as.numeric(deg_cancer_freq))

deg_normal_freq <- table(deg_normal_nonzero) / length(deg_normal_nonzero)
deg_normal_log <- log10(as.numeric(names(deg_normal_freq)))
freq_normal_log <- log10(as.numeric(deg_normal_freq))

# Fits a linear model to the log-log data and computes the R-squared value for the cancer network.
cancer_lm <- lm(freq_cancer_log ~ deg_cancer_log)
cancer_r_squared <- summary(cancer_lm)$r.squared
cat("Cancer Network - R-squared (Log-Log Fit):", cancer_r_squared, "\n")

# Fits a linear model to the log-log data and computes the R-squared value for the normal network.
normal_lm <- lm(freq_normal_log ~ deg_normal_log)
normal_r_squared <- summary(normal_lm)$r.squared
cat("Normal Network - R-squared (Log-Log Fit):", normal_r_squared, "\n")

# Interpret Results
# Evaluates whether the networks exhibit scale-free topology based on R-squared values.
if (cancer_r_squared > 0.8) {
  cat("Cancer Network: Likely follows a scale-free topology.\n")
} else {
  cat("Cancer Network: Does not strongly exhibit scale-free topology.\n")
}

if (normal_r_squared > 0.8) {
  cat("Normal Network: Likely follows a scale-free topology.\n")
} else {
  cat("Normal Network: Does not strongly exhibit scale-free topology.\n")
}

# Apply power-law fitting to verify scale-free properties.
fit_power_law(deg_cancer)
fit_power_law(deg_normal)

# Note: Our network does not follow a scale-free topology, so identifying hubs is unnecessary (see guidelines).
# The following hub-related code is retained for soem considerations.

# 3.3 --- Identify Hubs ---
# Identifies hubs in cancer and normal networks using the 95th percentile as a threshold.
hub_threshold_cancer <- quantile(deg_cancer, 0.95)
hub_threshold_normal <- quantile(deg_normal, 0.95)

hubs_cancer <- names(deg_cancer[deg_cancer >= hub_threshold_cancer])
hubs_normal <- names(deg_normal[deg_normal >= hub_threshold_normal])

cat("Cancer Network hubs:", length(hubs_cancer), "\n")
cat("Normal Network hubs:", length(hubs_normal), "\n")

# Compare Hubs
# Identifies hubs unique to each network.
unique_hubs_cancer <- setdiff(hubs_cancer, hubs_normal)
unique_hubs_normal <- setdiff(hubs_normal, hubs_cancer)

cat("Hubs unique to Cancer Network:", unique_hubs_cancer, "\n")
cat("Hubs unique to Normal Network:", unique_hubs_normal, "\n")

# Visualizes the hub subnetwork for the cancer network.
# Identify hub IDs in the cancer network
hubs_cancer_ids <- match(hubs_cancer, rownames(adj_cancer))  # Finds the indices of hubs in the adjacency matrix.

# Identify the neighborhood of hubs in the cancer network
hubs_cancer_neigh <- c()
for (f in hubs_cancer_ids) {
  hubs_cancer_neigh <- append(hubs_cancer_neigh, neighbors(graph_cancer, f))  # Finds neighbors.
}
hubs_cancer_neigh <- unique(hubs_cancer_neigh)  # Removes duplicates.

# Get the names of the neighboring nodes
hubs_cancer_neigh_names <- rownames(adj_cancer)[hubs_cancer_neigh]  # Names of neighbors.

# Create a subgraph with hubs and their neighbors
subnet <- unique(c(hubs_cancer, hubs_cancer_neigh_names))  # Combines hubs and neighbors.
hub_cancer_adj <- adj_cancer[subnet, subnet]  # Adjacency matrix for the subgraph.

# Create a network for the subgraph
net_hub <- network(hub_cancer_adj, matrix.type = "adjacency", ignore.eval = FALSE, names.eval = "weights")

# Analyze the density of the subgraph
cat("Network Density (Hub Subnetwork):", network.density(net_hub), "\n")

# Count positive and negative edges
cat("Positive Edges:", sum(hub_cancer_adj > 0) / 2, "\n")
cat("Negative Edges:", sum(hub_cancer_adj < 0) / 2, "\n")

# Assign attributes to nodes (hub or non-hub with distinct colors)
net_hub %v% "type" <- ifelse(network.vertex.names(net_hub) %in% hubs_cancer, "hub", "non-hub")
net_hub %v% "color" <- ifelse(net_hub %v% "type" == "hub", "tomato", "skyblue")

# Extract edges as a list
edges <- as.matrix(as.edgelist(net_hub))  # Converts edge list to a matrix.
edge_start <- network.vertex.names(net_hub)[edges[, 1]]  # Start node (names).
edge_end <- network.vertex.names(net_hub)[edges[, 2]]    # End node (names).

# Determine edge colors
edge_color <- ifelse(
  edge_start %in% hubs_cancer & edge_end %in% hubs_cancer, "green",  # Hub-to-hub.
  ifelse(
    edge_start %in% hubs_cancer | edge_end %in% hubs_cancer, 
    "orange",  # Hub-to-non-hub.
    "gray"  # Non-hub-to-non-hub.
  )
)

# Set edge colors directly in the `network` object
net_hub %e% "ecolor" <- edge_color

# Visualize the subgraph using ggnet2
ggnet2(net_hub,  
       color = "color", 
       alpha = 0.9, 
       size = 2, 
       edge.color = "ecolor", 
       edge.alpha = 0.9,  
       edge.size = 0.15) +
  guides(size = "none") +
  ggtitle("Hub Subnetwork (Cancer Network)")

# Visualizes the hub subnetwork for the normal network.
# Identify hub IDs in the normal network
hubs_normal_ids <- match(hubs_normal, rownames(adj_normal))  # Finds the indices of hubs in the adjacency matrix.

# Identify the neighborhood of hubs in the normal network
hubs_normal_neigh <- c()
for (f in hubs_normal_ids) {
  hubs_normal_neigh <- append(hubs_normal_neigh, neighbors(graph_normal, f))  # Finds neighbors.
}
hubs_normal_neigh <- unique(hubs_normal_neigh)  # Removes duplicates.

# Get the names of the neighboring nodes
hubs_normal_neigh_names <- rownames(adj_normal)[hubs_normal_neigh]  # Names of neighbors.

# Create a subgraph with hubs and their neighbors
subnet_normal <- unique(c(hubs_normal, hubs_normal_neigh_names))  # Combines hubs and neighbors.
hub_normal_adj <- adj_normal[subnet_normal, subnet_normal]  # Adjacency matrix for the subgraph.

# Create a network for the subgraph
net_hub_normal <- network(hub_normal_adj, matrix.type = "adjacency", ignore.eval = FALSE, names.eval = "weights")

# Analyze the density of the subgraph
cat("Network Density (Hub Subnetwork - Normal):", network.density(net_hub_normal), "\n")

# Count positive and negative edges
cat("Positive Edges:", sum(hub_normal_adj > 0) / 2, "\n")
cat("Negative Edges:", sum(hub_normal_adj < 0) / 2, "\n")

# Assign attributes to nodes (hub or non-hub with distinct colors)
net_hub_normal %v% "type" <- ifelse(network.vertex.names(net_hub_normal) %in% hubs_normal, "hub", "non-hub")
net_hub_normal %v% "color" <- ifelse(net_hub_normal %v% "type" == "hub", "tomato", "skyblue")

# Extract edges as a list
edges_normal <- as.matrix(as.edgelist(net_hub_normal))  # Converts edge list to a matrix.
edge_start_normal <- network.vertex.names(net_hub_normal)[edges_normal[, 1]]  # Start node (names).
edge_end_normal <- network.vertex.names(net_hub_normal)[edges_normal[, 2]]    # End node (names).

# Determine edge colors
edge_color_normal <- ifelse(
  edge_start_normal %in% hubs_normal & edge_end_normal %in% hubs_normal, "green",  # Hub-to-hub.
  ifelse(
    edge_start_normal %in% hubs_normal | edge_end_normal %in% hubs_normal, 
    "orange",  # Hub-to-non-hub.
    "gray"  # Non-hub-to-non-hub.
  )
)

# Set edge colors directly in the `network` object
net_hub_normal %e% "ecolor" <- edge_color_normal

# Visualize the subgraph using ggnet2
ggnet2(net_hub_normal,  
       color = "color", 
       alpha = 0.9, 
       size = 2, 
       edge.color = "ecolor", 
       edge.alpha = 0.9,  
       edge.size = 0.15) +
  guides(size = "none") +
  ggtitle("Hub Subnetwork (Normal Network)")

# 4: --- DIFFERENTIAL CO-EXPRESSED NETWORKS ---
# 4.1 --- Compute Differential Co-expression Network (DCN) ---
# Fisher Z-transformation function
fisher_z <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}
# Input: Correlation matrices for cancer and normal conditions
# Using `cor_cancer` and `cor_normal` from previous sections
n_cancer <- ncol(deg_expr_c)  # Number of samples in cancer condition
n_normal <- ncol(deg_expr_n)  # Number of samples in normal condition

# Fisher z-transformation
z_cancer <- fisher_z(cor_cancer)
z_normal <- fisher_z(cor_normal)

# Differential z-score computation
z_diff <- (z_cancer - z_normal) / sqrt((1 / (n_cancer - 3)) + (1 / (n_normal - 3)))

# Thresholding to create binary adjacency matrix
threshold <- 3  # Threshold for |Z| < 3
adj_dcn <- abs(z_diff) > threshold
diag(adj_dcn) <- 0  # Remove self-loops

graph_dcn <- graph_from_adjacency_matrix(adj_dcn, mode = "undirected", diag = FALSE)

# 4.2 --- Analyze the Differential Co-expression Network ---
# Identify hubs as nodes in the top 5% by degree
hub_threshold_dcn <- quantile(deg_dcn, 0.95)
hubs_dcn <- names(deg_dcn[deg_dcn >= hub_threshold_dcn])
cat("Hubs identified in the network:\n", paste(hubs_dcn, collapse = ", "), "\n")

# Assign node type (hub or non-hub) in the network
net_diff %v% "type" <- ifelse(network.vertex.names(net_diff) %in% hubs_dcn, "hub", "non-hub")

# Assign colors to nodes (hub = purple, non-hub = cyan)
net_diff %v% "color" <- ifelse(net_diff %v% "type" == "hub", "red", "blue")

# Check if either of the nodes connected by an edge is a hub and set the edge colors
node_types_diff <- net_diff %v% "type"

# Use apply to iterate over edges and color edges connecting at least one hub
edge_colors_diff <- apply(as.matrix(net_diff, matrix.type = "edgelist"), 1, function(edge) {
  node1 <- as.numeric(edge[1])  # Convert to numeric
  node2 <- as.numeric(edge[2])
  
  if (node_types_diff[node1] == "hub" || node_types_diff[node2] == "hub") {
    "orange"  # Color the edge orange if at least one node is a hub
  } else {
    "grey"  # Otherwise, color the edge light green
  }
})

# Set edge colors in the network
network::set.edge.attribute(net_diff, "edgecolor", edge_colors_diff)

# Check the distribution of edge colors
table(net_diff %e% "edgecolor")

# Visualize the network with ggnet2
ggnet2(
  net_diff, 
  color = "color", 
  alpha = 0.7, 
  size = 2, 
  edge.color = "edgecolor", 
  edge.alpha = 1, 
  edge.size = 0.15
) + 
  guides(size = "none") + 
  theme(panel.background = element_rect(fill = "white")) + 
  ggtitle("Subnetwork with Hubs and Non-Hubs")

# 5: --- PATIENT SIMILARITY NETWORK ---
# 5.1: --- Compute the Patient Similarity Network using cancer gene expression profile (consider only the list of DEGs) ---

# Select DEGs expression data (Cancer condition)
patient_similarity <- cor(log2(deg_expr_c + 1), method = "pearson")
dim(patient_similarity)
# Threshold the correlation matrix to create a binary adjacency matrix
# Keep only strong correlations above a defined threshold
threshold <- 0.7
adj_psn <- ifelse(abs(patient_similarity) >= threshold, 1, 0)
# Remove self-loops (no patient is "similar" to itself in the adjacency matrix)
diag(adj_psn) <- 0
# Create the Patient Similarity Network (PSN) as an igraph object
graph_psn <- graph_from_adjacency_matrix(adj_psn, mode = "undirected", diag = FALSE)
# Visualize the PSN
plot(graph_psn, vertex.size = 5, vertex.label = NA,
     main = "Patient Similarity Network (Cancer DEGs)",
     vertex.color = "lightblue", edge.color = "gray")
# Analyze the PSN
# Compute basic network properties
cat("Number of nodes (patients):", vcount(graph_psn), "\n")
cat("Number of edges (similarity connections):", ecount(graph_psn), "\n")
cat("Network density:", edge_density(graph_psn), "\n")
# Compute degree distribution
deg_psn <- degree(graph_psn)
hist(deg_psn, breaks = 30, main = "Degree Distribution (Patient Similarity Network)",
     xlab = "Degree", col = "lightblue", border = "black")
# Identify patient hubs (top 5% by degree)
hub_threshold_psn <- quantile(deg_psn, 0.95)
hubs_psn <- names(deg_psn[deg_psn >= hub_threshold_psn])
cat("Patient hubs:", hubs_psn, "\n")

# 5.2: --- Perform the community detection (Louvain algorithm to the PSN) and represent the community structure ---

# Apply the Louvain algorithm to the PSN
louvain_communities <- cluster_louvain(graph_psn)  # `graph_psn` is the graph of the PSN

# Add community membership as a vertex attribute
V(graph_psn)$community <- membership(louvain_communities)

# Print the number of communities detected
cat("Number of communities detected:", length(unique(V(graph_psn)$community)), "\n")

# Generate a visually distinct color palette for the communities
community_colors <- rainbow(length(unique(V(graph_psn)$community)), alpha = 0.8)

# Assign colors to nodes based on their community
V(graph_psn)$color <- community_colors[V(graph_psn)$community]

# Assign a black border to nodes for better contrast
V(graph_psn)$frame.color <- "black"

# Adjust edge color and transparency for better readability
E(graph_psn)$color <- "gray70"
E(graph_psn)$width <- 0.5

# Use a layout that minimizes node overlap
layout_louvain <- layout_with_kk(graph_psn) # Kamada-Kawai layout for better spacing

# Plot with enhanced aesthetics
plot(
  graph_psn, 
  layout = layout_louvain, 
  vertex.label = NA,                 # No vertex labels for clarity
  vertex.size = 8,                   # Increased node size for better visibility
  vertex.color = V(graph_psn)$color, # Node color based on community
  vertex.frame.color = V(graph_psn)$frame.color, # Node border color
  edge.color = E(graph_psn)$color,   # Edge color
  edge.width = E(graph_psn)$width,   # Edge width
  main = "PSN Louvain Community Structure"
)

# Analyze community sizes
community_sizes <- sizes(louvain_communities)
cat("Community sizes:", community_sizes, "\n")

# 5.3: --- Similarity Network Fusion ---
# 5.3.1 Apply the algorithm to the 2-layer PSN obtained using gene expression profile and mutational profile

# Mutational data
quetab2query.mut <- GDCquery(
  project = proj, 
  data.category = "Simple Nucleotide Variation", 
  data.type = "Masked Somatic Mutation")
# GDCdownload(quetab2query.mut)
mut.data <- GDCprepare(quetab2query.mut) 
#mut.data_aux <- mut.data
mut.data$Tumor_Sample_Barcode <-substr(mut.data$Tumor_Sample_Barcode, 1,12)
length(unique(mut.data$Tumor_Sample_Barcode) )
colnames(filtr.expr.c) <- gsub("-1$", "", colnames(filtr.expr.c))
mut.data <-  mut.data[mut.data$Tumor_Sample_Barcode %in% colnames(filtr.expr.c),]
maf1 <- mut.data %>% read.maf

# Useful plots 
datatable(getSampleSummary(maf1), filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10),  rownames = FALSE)
plotmafSummary(maf = maf1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
oncoplot(maf = maf1, top = 10, removeNonMutated = TRUE)
lollipopPlot( maf = maf1,gene = 'TTN',  showMutationRate = F)
lollipopPlot( maf = maf1,gene = 'MUC16',  showMutationRate = F)
somaticInteractions(maf = maf1, top = 20, pvalue = c(0.05, 0.01))

# Binary mutation matrix cols = patiens, rows = genes
# cells to indicate if a specific gene is muated in a patient
mut_bin_mat <- table(mut.data$Tumor_Sample_Barcode, mut.data$Hugo_Symbol) > 0 
mut_bin_mat <- as.matrix(mut_bin_mat)

similarity_mat_mut <- 1 - as.matrix(vegdist(mut_bin_mat, method = "jaccard", binary = TRUE))
diag(similarity_mat_mut) <- 0
head(similarity_mat_mut)

# Normalize both
diag(patient_similarity) = 0
norm_exp <- standardNormalization(patient_similarity)
norm_mut <- standardNormalization(similarity_mat_mut)

common_ids <- intersect(rownames(norm_exp), rownames(norm_mut))
norm_exp <- norm_exp[common_ids, common_ids]
norm_mut <- norm_mut[common_ids, common_ids]

range(norm_exp) # -5.509984  2.571265
range(norm_mut) # -1.870727  6.399997

all.equal(dim(norm_exp), dim(norm_mut)) # Should return TRUE
all.equal(rownames(norm_exp), rownames(norm_mut)) # Should return TRUE
all.equal(colnames(norm_exp), colnames(norm_mut)) # Should return TRUE

similarity_list <- list(norm_exp, norm_mut)
# Apply Similarity Network Fusion (SNF)
fused_network <- SNF(similarity_list, K = 30, t = 50)  # K = neighbors, t = iterations
range(fused_network)
# Threshold the correlation matrix to create a binary adjacency matrix
# Keep only strong correlations above a defined threshold
quantile(fused_network)
threshold <- 0.01
adj_fnet <- ifelse(abs(fused_network) >= threshold, 1, 0)
# Remove self-loops (no patient is "similar" to itself in the adjacency matrix)
diag(adj_fnet) <- 0
# Convert fused network to igraph object
fused_graph <- graph_from_adjacency_matrix(adj_fnet, mode = "undirected", diag = FALSE)

# Apply Louvain algorithm
community <- cluster_louvain(fused_graph)

# Set a better layout for the graph
layout <- layout_with_fr(fused_graph) # Fruchterman-Reingold layout

# Classify edges based on whether they connect nodes in the same or different communities
edge_colors <- ifelse(
  membership(community)[ends(fused_graph, E(fused_graph))[, 1]] ==
    membership(community)[ends(fused_graph, E(fused_graph))[, 2]],
  "gray",  # Same community
  "lightblue"      # Different communities
)

# Plot with distinct edge colors
plot(
  community, fused_graph, 
  layout = layout,
  vertex.size = 6,                        # Adjust size of the nodes
  vertex.color = membership(community),   # Color nodes by community
  vertex.frame.color = "black",           # Add a border to the nodes
  vertex.label = NA,                      # Remove vertex labels
  edge.color = edge_colors,               # Color edges by intra- or inter-community
  edge.width = E(fused_graph)$weight * 5, # Scale edge width by weight
  edge.curved = 0.2,                      # Add slight curvature to edges
  main = "Community Detection on Fused Network",
  sub = paste("Number of Communities:", length(unique(membership(community))))
)

clusters <- membership(community)
table(clusters)  # Number of patients in each cluster

# 5.3.2 Compare the community structure with that obtained in 5.1
# For the Patient Similarity Network (PSN)
# Get community membership from the PSN
psn_community_membership <- membership(louvain_communities)  # louvain_communities is the PSN community object
psn_communities <- split(names(psn_community_membership), psn_community_membership)

# Print nodes for each PSN community
cat("Patient Similarity Network (PSN) - Nodes by Community:\n")
for (community_id in names(psn_communities)) {
  cat("Community", community_id, ":", paste(psn_communities[[community_id]], collapse = ", "), "\n")
}

# For the Fused Network
# Get cluster membership from the Fused Network
fused_community_membership <- membership(community)  # community is the Fused Network community object
fused_communities <- split(names(fused_community_membership), fused_community_membership)

# Print nodes for each cluster of the Fused Network
cat("\nFused Network - Nodes by Cluster:\n")
for (cluster_id in names(fused_communities)) {
  cat("Cluster", cluster_id, ":", paste(fused_communities[[cluster_id]], collapse = ", "), "\n")
}

# Comparison between the two networks
cat("\nComparison of communities between the two networks:\n")
common_nodes <- intersect(names(psn_community_membership), names(fused_community_membership))

for (node in common_nodes) {
  psn_comm <- psn_community_membership[node]
  fused_comm <- fused_community_membership[node]
  if (psn_comm != fused_comm) {
    cat("Node", node, "is in PSN Community", psn_comm, "and in Fused Cluster", fused_comm, "\n")
  }
}

# 6: --- ADDITIONAL ANALYSIS ---
# 6.1: --- Identify patients with high expression levels of IL6 and SFTPC ---
genes_of_interest <- c("IL6", "SFTPC")

# Check if the genes of interest are present in the expression matrix
missing_genes <- setdiff(genes_of_interest, rownames(deg_expr_c_unique))
if (length(missing_genes) > 0) {
  stop(paste("The following genes are not present in the matrix:", paste(missing_genes, collapse = ", ")))
}

# Subset the expression matrix for the genes of interest
expression_subset <- deg_expr_c_unique[genes_of_interest, , drop = FALSE]

# Define a threshold for high expression levels (e.g., 75th percentile)
thresholds <- apply(expression_subset, 1, quantile, probs = 0.75, na.rm = TRUE)

# Identify patients with expression levels above the threshold for each gene
high_expression_patients <- lapply(rownames(expression_subset), function(gene) {
  colnames(expression_subset)[expression_subset[gene, ] >= thresholds[gene]]
})
names(high_expression_patients) <- genes_of_interest

# Print patients with high expression levels for each gene
cat("Patients with high expression levels:\n")
for (gene in genes_of_interest) {
  cat(gene, ":", paste(high_expression_patients[[gene]], collapse = ", "), "\n")
}

# Analyze community membership in the PSN
psn_membership <- membership(louvain_communities)  # Community membership in the PSN
psn_community_patients <- lapply(high_expression_patients, function(patients) {
  communities <- psn_membership[patients]
  return(table(communities))
})
cat("Distribution across PSN communities:\n")
print(psn_community_patients)

# Analyze cluster membership in the Fused Network
fused_membership <- membership(community)  # Cluster membership in the Fused Network
fused_cluster_patients <- lapply(high_expression_patients, function(patients) {
  clusters <- fused_membership[patients]
  return(table(clusters))
})
cat("Distribution across Fused Network clusters:\n")
print(fused_cluster_patients)

# Barplot for PSN
barplot(
  unlist(psn_community_patients),
  main = "Distribution across PSN communities",
  xlab = "Communities",
  ylab = "Number of Patients",
  col = "lightblue"
)

# Barplot for Fused Network
barplot(
  unlist(fused_cluster_patients),
  main = "Distribution across Fused Network clusters",
  xlab = "Clusters",
  ylab = "Number of Patients",
  col = "lightgreen"
)

# Patients with high expression levels and their respective PSN communities and Fused Network clusters
cat("Patients with high expression levels and their PSN community and Fused Network cluster membership:\n")

for (gene in names(high_expression_patients)) {
  cat("\nGene:", gene, "\n")
  
  # Patients with high expression levels for the gene
  patients <- high_expression_patients[[gene]]
  
  # PSN communities of the patients
  psn_patient_communities <- psn_community_membership[patients]
  psn_patient_communities <- split(patients, psn_patient_communities)
  
  cat("  PSN Communities:\n")
  for (community_id in names(psn_patient_communities)) {
    cat("    Community", community_id, ":", paste(psn_patient_communities[[community_id]], collapse = ", "), "\n")
  }
  
  # Fused Network clusters of the patients
  fused_patient_clusters <- fused_community_membership[patients]
  fused_patient_clusters <- split(patients, fused_patient_clusters)
  
  cat("  Fused Network Clusters:\n")
  for (cluster_id in names(fused_patient_clusters)) {
    cat("    Cluster", cluster_id, ":", paste(fused_patient_clusters[[cluster_id]], collapse = ", "), "\n")
  }
}

# 6.2: --- Gene Enrichment Analysis in PSN Communities and Fused Network Clusters ---
# Check available databases
dbs <- listEnrichrDbs()
print(dbs)

# Select databases of interest (e.g., GO_Biological_Process_2021, KEGG_2021_Human)
selected_dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")

# Function to perform enrichment analysis with enrichR
perform_enrichment_enrichr <- function(genes_list) {
  enrichment_results <- lapply(genes_list, function(genes) {
    enrichr(genes, databases = selected_dbs)
  })
  return(enrichment_results)
}

# Enrichment for PSN Communities
# Extract genes associated with each PSN community
# Filter patients belonging to both the expression matrix and PSN communities
common_patients <- intersect(colnames(deg_expr_c_unique), names(psn_community_membership))

# Filter the expression matrix for common patients
deg_expr_c_filtered <- deg_expr_c_unique[, common_patients, drop = FALSE]

# For each PSN community, find expressed genes
psn_community_genes <- lapply(unique(psn_community_membership[common_patients]), function(community) {
  # Find patients belonging to this community
  patients_in_community <- names(psn_community_membership[common_patients]
                                 [psn_community_membership[common_patients] == community])
  
  # Find genes expressed (value > 0) in at least one patient in the community
  genes_expressed <- rownames(deg_expr_c_filtered)[rowSums(
    deg_expr_c_filtered[, patients_in_community, drop = FALSE]) > 0]
  return(genes_expressed)
})

# Assign community names to the list
names(psn_community_genes) <- unique(psn_community_membership[common_patients])

# Display expressed genes for a specific PSN community as an example
cat("Expressed genes in PSN Community 1:\n")
print(psn_community_genes[["1"]])

# Compute enrichment results for PSN communities
psn_enrichment <- perform_enrichment_enrichr(psn_community_genes)

# Enrichment for Fused Network Clusters
# Extract genes associated with each Fused Network cluster
# Filter patients belonging to both the expression matrix and Fused Network clusters
common_patients_fused <- intersect(colnames(deg_expr_c_unique), names(fused_membership))

# Filter the expression matrix for common patients
deg_expr_c_filtered_fused <- deg_expr_c_unique[, common_patients_fused, drop = FALSE]

# For each Fused Network cluster, find expressed genes
fused_cluster_genes <- lapply(unique(fused_membership[common_patients_fused]), function(cluster) {
  # Find patients belonging to this cluster
  patients_in_cluster <- names(
    fused_membership[common_patients_fused][fused_membership[common_patients_fused] == cluster])
  
  # Find genes expressed (value > 0) in at least one patient in the cluster
  genes_expressed <- rownames(deg_expr_c_filtered_fused)[rowSums(
    deg_expr_c_filtered_fused[, patients_in_cluster, drop = FALSE]) > 0]
  return(genes_expressed)
})

# Assign cluster names to the list
names(fused_cluster_genes) <- unique(fused_membership[common_patients_fused])

# Display expressed genes for a specific cluster as an example
cat("Expressed genes in Fused Network Cluster 1:\n")
print(fused_cluster_genes[["1"]])

# Compute enrichment results for Fused Network clusters
fused_enrichment <- perform_enrichment_enrichr(fused_cluster_genes)

# Examine Results
# Display results for a specific PSN community (e.g., Community 1)
print(psn_enrichment[["1"]][["GO_Biological_Process_2021"]])

# Display results for a specific Fused Network cluster (e.g., Cluster 1)
print(fused_enrichment[["1"]][["GO_Biological_Process_2021"]])

# 6.3: --- Pathway Enrichment Analysis for Communities and Clusters ----

# Function to find enriched pathways for a specific community or cluster
find_enriched_pathways <- function(patients, enrichment_results, gene_name) {
  enriched_pathways <- lapply(names(enrichment_results), function(group) {
    pathways <- enrichment_results[[group]][["GO_Biological_Process_2021"]]
    significant_pathways <- pathways[pathways$Adjusted.P.value <= 0.05, ]
    if (nrow(significant_pathways) > 0) {
      significant_pathways$Community_Cluster <- group
      significant_pathways$Gene <- gene_name
    }
    return(significant_pathways)
  })
  enriched_pathways <- do.call(rbind, enriched_pathways)
  return(enriched_pathways)
}

# Enrichment analysis for PSN communities
psn_il6_pathways <- find_enriched_pathways(high_expression_patients[["IL6"]], psn_enrichment, "IL6")
psn_sftpc_pathways <- find_enriched_pathways(high_expression_patients[["SFTPC"]], psn_enrichment, "SFTPC")

# Enrichment analysis for Fused Network clusters
fused_il6_pathways <- find_enriched_pathways(high_expression_patients[["IL6"]], fused_enrichment, "IL6")
fused_sftpc_pathways <- find_enriched_pathways(high_expression_patients[["SFTPC"]], fused_enrichment, "SFTPC")

# Display Results
# Combine results into a single data frame for visualization
all_pathways <- rbind(psn_il6_pathways, psn_sftpc_pathways, fused_il6_pathways, fused_sftpc_pathways)

# Ensure the data frame has the necessary columns
# Check for "Term", "Combined.Score", "Gene", and "Community_Cluster" columns
if (!all(c("Term", "Combined.Score", "Gene", "Community_Cluster") %in% colnames(all_pathways))) {
  stop("The all_pathways data frame must contain Term, Combined.Score, Gene, and Community_Cluster columns")
}

# Filter for the most relevant pathways (e.g., top 10 per Community_Cluster)
top_pathways <- all_pathways %>%
  group_by(Community_Cluster) %>%
  slice_max(order_by = Combined.Score, n = 10) %>%
  ungroup()

# Barplot of enriched pathways
ggplot(top_pathways, aes(x = reorder(Term, -Combined.Score), y = Combined.Score, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  facet_wrap(~ Community_Cluster, scales = "free") +
  labs(
    title = "Enriched Pathways by Gene and Group",
    x = "Biological Pathway",
    y = "Combined Score",
    fill = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

# Compare Results with LUAD Subtypes
# Direct comparison for pathways associated with LUAD subtypes
subtype_pathways <- c("cell proliferation", "inflammatory response", "terminal respiratory unit development")

linked_pathways <- all_pathways[grep(paste(subtype_pathways, collapse = "|"), all_pathways$Term, ignore.case = TRUE), ]
cat("Enriched Pathways associated with LUAD subtypes:\n")
print(linked_pathways)

# Analysis of Patients with High Expression Levels
# Group patients with high expression levels for each gene by community/cluster
analyze_patients_groups <- function(gene_name, community_results, cluster_results) {
  cat("\nResults for", gene_name, ":\n")
  cat("PSN Communities:\n")
  for (community in names(community_results[[gene_name]])) {
    patients <- community_results[[gene_name]][[community]]
    cat("  Community", community, ":", paste(patients, collapse = ", "), "\n")
  }
  
  cat("Fused Network Clusters:\n")
  for (cluster in names(cluster_results[[gene_name]])) {
    patients <- cluster_results[[gene_name]][[cluster]]
    cat("  Cluster", cluster, ":", paste(patients, collapse = ", "), "\n")
  }
}

# Print results for IL6 and SFTPC
analyze_patients_groups("IL6", psn_community_patients, fused_cluster_patients)
analyze_patients_groups("SFTPC", psn_community_patients, fused_cluster_patients)

# 6.2: --- Compute a different centrality index (CI) and check the overlap between the 5% of the nodes with 
# highest CI values and the degree-based hubs ---

# Compute betweenness centrality for Cancer and Normal networks
betweenness_cancer <- betweenness(graph_cancer, directed = FALSE)
betweenness_normal <- betweenness(graph_normal, directed = FALSE)

# Threshold: top 5% nodes by betweenness centrality
hub_threshold_betweenness_cancer <- quantile(betweenness_cancer, 0.95)
hub_threshold_betweenness_normal <- quantile(betweenness_normal, 0.95)

hubs_betweenness_cancer <- names(betweenness_cancer[betweenness_cancer >= hub_threshold_betweenness_cancer])
hubs_betweenness_normal <- names(betweenness_normal[betweenness_normal >= hub_threshold_betweenness_normal])

cat("Betweenness hubs - Cancer Network:", length(hubs_betweenness_cancer), "\n")
cat("Betweenness hubs - Normal Network:", length(hubs_betweenness_normal), "\n")

# Compare hubs from degree and betweenness centrality
overlap_cancer <- intersect(hubs_cancer, hubs_betweenness_cancer)
overlap_normal <- intersect(hubs_normal, hubs_betweenness_normal)

cat("Overlap between Degree-based and Betweenness-based hubs (Cancer):", length(overlap_cancer), "\n")
cat("Overlap between Degree-based and Betweenness-based hubs (Normal):", length(overlap_normal), "\n")

# Print nodes in common
cat("Common hubs (Cancer):", overlap_cancer, "\n")
cat("Common hubs (Normal):", overlap_normal, "\n")

# 6.3: --- Differential Co-expressed Network Analysis: Extracting Positive and Negative Subnetworks and 
# Identifying Key Hubs ---

# Create subnetworks
# Subnetwork with positive connections
adj_positive <- ifelse(z_diff > 0, 1, 0)
diag(adj_positive) <- 0  # Remove self-loops

# Subnetwork with negative connections
adj_negative <- ifelse(z_diff < 0, 1, 0)
diag(adj_negative) <- 0  # Remove self-loops

# Convert to igraph objects
graph_positive <- graph_from_adjacency_matrix(adj_positive, mode = "undirected", diag = FALSE)
graph_negative <- graph_from_adjacency_matrix(adj_negative, mode = "undirected", diag = FALSE)

# Calculate degrees for both subnetworks
degree_positive <- degree(graph_positive)
degree_negative <- degree(graph_negative)

# Identify positive and negative hubs separately
# Calculate the threshold for the highest degree hubs (positive and negative)
threshold_positive <- quantile(degree_positive, 0.95)
threshold_negative <- quantile(degree_negative, 0.95)

# Hubs with the highest positive degrees
hubs_positive <- names(degree_positive[degree_positive >= threshold_positive])

# Hubs with the highest negative degrees
hubs_negative <- names(degree_negative[degree_negative >= threshold_negative])

# Print hubs and their degrees
cat("Hubs with highest positive degree values:\n")
positive_hub_degrees <- degree_positive[hubs_positive]
cat(paste(hubs_positive, ":", positive_hub_degrees), "\n")

cat("Hubs with highest negative degree values:\n")
negative_hub_degrees <- degree_negative[hubs_negative]
cat(paste(hubs_negative, ":", negative_hub_degrees), "\n")

# Sort by degree (descending)
sorted_positive_hubs <- sort(positive_hub_degrees, decreasing = TRUE)
sorted_negative_hubs <- sort(negative_hub_degrees, decreasing = TRUE)

# Print hubs sorted by degree
cat("\nSorted Hubs with Highest Positive Degree Values:\n")
cat(paste(names(sorted_positive_hubs), ":", sorted_positive_hubs), "\n")

cat("\nSorted Hubs with Highest Negative Degree Values:\n")
cat(paste(names(sorted_negative_hubs), ":", sorted_negative_hubs), "\n")

