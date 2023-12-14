library("tximportData")
library("tximport")
library("readr")
library("DESeq2")
library("dplyr")
library("GenomicFeatures")
library("ggplot2")
library("apeglm")
library(umap)
library(tidyr)
library(pheatmap)

#creates transcript database, gets transcript Ids with corresponding geneIDs, gets rid of rows with N/A
txdb <- makeTxDbFromGFF("~/biol6150/ProjectSubmissions/Group23/Project56/GRCh38_latest_genomic.gff.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene = tx2gene %>% filter(!is.na(TXNAME))
head(tx2gene)
write.table(tx2gene, "~/biol6150/ProjectSubmissions/Group23/Project56/tx2gene.RefSeq.All.tsv", 
            sep = ",", 
            row.names = FALSE)

# Read tx2gene if already created before.
tx2gene <- read_csv("~/biol6150/ProjectSubmissions/Group23/Project56/tx2gene.RefSeq.All.tsv") %>% as.data.frame()
head(tx2gene)

# Read the samples map file.
sample_info <- read.table("~/biol6150/ProjectSubmissions/Group23/Project56//RNAseqanalysis/quant/SampleInfo.txt", header=TRUE) 
samples = sample_info %>% pull(SampleID)
print(samples)

file_paths = paste0("~/biol6150/ProjectSubmissions/Group23/Project56//RNAseqanalysis/quant/", samples, "/quant.sf")
print(file_paths)

# Create a txi object. Read Salmon SF files.
txi.salmon <- tximport(file_paths, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
print(head(txi.salmon$counts))


# Deseq needs specific format for the map file.
sampleTable <- data.frame(condition = factor(c(rep("UC",15), rep("CD",15))
))
rownames(sampleTable) <- colnames(txi.salmon$counts)
sampleTable

#Create the dds object. This prepares the data for DESeq2 to run diffrential expression.
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)


## Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
head(res)

# Check how many transcripts pass the filter by adjusted p value.
table(res$padj<0.05)

# Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)

#Visualize the results.
res = res %>% as.data.frame()

ggplot(res, aes(x = log2FoldChange,
                y = -log10(padj))) + 
  geom_point() +
  theme_bw(base_size = 24)

#Get significant results.
res_sig = res %>% filter(padj < 0.05) 
dim(res_sig)
head(res_sig)

write.csv(as.data.frame(res_sig), file = "~/biol6150/ProjectSubmissions/Group23/Project56/deseq2_p0.05_results.csv")



# Volcano plot
# Filter for specific genes in volcano plot
# Mapping of gene IDs to custom labels
annotations <- data.frame(
  Gene = c("XM_011545822.3", "XM_006711210.3"),
  Label = c("AQP8", "ETV3")
)

# Create a new column for gene names from row names
res_sig$Gene <- rownames(res_sig)
res_sig_annotated <- merge(res_sig, annotations, by = "Gene", all.x = TRUE)

# Extract the subset for the genes to annotate
annotated_genes <- subset(res_sig_annotated, Gene %in% annotations$Gene)

# Base plot
p <- ggplot(res_sig_annotated, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(pch = 21, fill = "magenta") +  # Default points
  ylim(c(0, 35)) +
  xlim(c(-30, 30)) +
  theme_bw(base_size = 24)

# Add the annotated genes with different properties
p <- p + geom_point(data = annotated_genes,
                    aes(x = log2FoldChange, y = -log10(padj)),
                    pch = 24,  
                    size = 5,  
                    fill = "green") 

# Add text labels for annotated genes using custom labels, positioned above the triangles
p <- p + geom_text(data = annotated_genes,
                   aes(x = log2FoldChange, y = -log10(padj), label = Label),
                   vjust = -0.5,  
                   hjust = 0.5,
                   size = 6)  

# Print the plot
p

## Log fold change shrinkage. This is a different package which performs shrinkage of log2fold.
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_UC_vs_CD", type="apeglm")
head(resLFC)

resLFC_sif = resLFC %>% as.data.frame() %>% filter(padj < 0.05)
dim(resLFC_sif)
head(resLFC_sif)

#Volcano plot for just the significant results.
ggplot(resLFC_sif, aes(x = log2FoldChange,
                       y = -log10(padj))) + 
  geom_point(pch = 21, fill = "magenta") +
  ylim(c(0,100)) +
  xlim(c(-50,50)) +
  theme_bw(base_size = 24)

#PCA plot
normalized_counts <- counts(dds, normalized=TRUE)
pca_data <- prcomp(t(normalized_counts))
pca_df <- as.data.frame(pca_data$x)
pca_df$condition <- sampleTable$condition

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2", color = "Condition")


#Box plot
#Plot gene of interests between UC and CD
normalized_counts <- counts(dds, normalized=TRUE)

# Specify your genes of interest
genes_of_interest <- c("XM_054330817.1", "NM_001257397.2")

# Subset for the genes of interest
normalized_counts_subset <- normalized_counts[rownames(normalized_counts) %in% genes_of_interest, ]
print(normalized_counts_subset)

# Transpose the data frame
data_t <- t(normalized_counts_subset)

# Convert to data frame and add row names as a new column
data_long <- as.data.frame(data_t)
data_long$Sample <- rownames(data_long)

# Use gather from tidyr to convert to long format
data_melted <- gather(data_long, key = "Gene", value = "Expression", -Sample)

# Add group information
data_melted$Group <- ifelse(as.numeric(sub("X", "", data_melted$Sample)) <= 15, "UC", "CD")

ggplot(data_melted, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free") +
  labs(title = "Gene Expression in UC and CD Groups", x = "Gene", y = "Expression Level") +
  scale_fill_manual(values = c("UC" = "blue", "CD" = "red")) +
  theme_minimal()

#Heatmap
# Extract normalized counts
#normalized_counts <- counts(dds, normalized = TRUE)

# Subset for the 857 genes in res_sig
#selected_genes <- rownames(res_sig)
#normalized_counts_subset <- normalized_counts[selected_genes, ]

#log_counts_subset <- log2(normalized_counts_subset + 1)

#sample_groups <- sampleTable # Replace 'group_column' with the actual column name

#pheatmap(log_counts_subset,
#         scale = "row",
#         cluster_rows = TRUE,
#         cluster_cols = TRUE,
#         annotation_col = data.frame(Group = sample_groups),
#         show_rownames = TRUE,
#         show_colnames = TRUE,
#         fontsize_row = 2, # Adjust the font size if necessary
#         fontsize_col = 2,
#         width = 857, 
#         height = 857)







#UMAP
# Perform PCA on normalized counts
pca_result <- prcomp(t(normalized_counts))

# Choose the number of principal components
# For instance, selecting the first 50 principal components
num_pc_to_use <- 20
pca_data <- pca_result$x[, 1:num_pc_to_use]

# Run UMAP on the PCA-reduced data
umap_result <- umap(pca_data)

# Prepare the UMAP data frame for plotting
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$condition <- sampleTable$condition

# Plot the UMAP results
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Plot", x = "UMAP Dimension 1", y = "UMAP Dimension 2", color = "Condition")
