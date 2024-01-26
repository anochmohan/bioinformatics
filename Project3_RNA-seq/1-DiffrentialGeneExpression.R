#Libraries Needed
library("tximportData")
library("tximport")
library("readr")
library("DESeq2")
library("dplyr")
library("GenomicFeatures")
library("ggplot2")
library("apeglm")
library("stringr")

# Define the directory for reference files. TximportData already has these for us for Gencode reference files
# We will make specific reference files using RefSeq gff3 files and GenomicFeatures package.
txdb <- makeTxDbFromGFF(file="/storage/ice-shared/biol6150/Data/DiffrentialExpression/Reference/GRCh38_latest_genomic.gff.gz")

# Read tx2gene if already created before.
tx2gene <- read_csv("/storage/ice-shared/biol6150/Data/DiffrentialExpression/Reference/tx2gene.RefSeq.All.tsv") %>% as.data.frame()

# Read the samples map file. # ADJUST SAMPLE MAP PATH HERE
sample_info <- read.table("/home/hice1/cwijeyesekera3/scratch/Data/Project5/SampleInfo.txt", header=TRUE) 
samples = sample_info %>% pull(SampleID)

#ADJUST FILE PATH HERE
file_paths = paste0("/home/hice1/cwijeyesekera3/scratch/Data/Project5/renamed_Quant/", samples, ".sf")

# Create a txi object. Read Salmon SF files.
txi.salmon <- tximport(file_paths, type = "salmon", txOut = TRUE, tx2gene = tx2gene)

# Deseq needs specific format for the map file.
sampleTable <- data.frame(condition = factor(c(rep("Normal",10), rep("Prostate_Cancer",10))))
rownames(sampleTable) <- colnames(txi.salmon$counts)

#Create the dds object. This prepares the data for DESeq2 to run diffrential expression.
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)

## Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
head(res)

write.csv(res, file="/home/hice1/cwijeyesekera3/scratch/Data/Project5/DE_results_filtered.csv")

res <- res[order(res$padj), ]
head(res)

#Visualize the results.
res = res %>% as.data.frame()

volc_plot <- ggplot(res, aes(x = log2FoldChange,
                             y = -log10(padj))) + 
  geom_point() +
  ggtitle(str_wrap("Volcano Plot of Differentially Expressed Genes", width=30)) +
  theme_bw(base_size = 20)
ggsave("/home/hice1/cwijeyesekera3/scratch/Data/Project5/DE_Results_Volcano.png",plot=volc_plot)

#Get significant results.
res_sig = res %>% filter(padj < 0.05) 
dim(res_sig)

#Volcano plot for just the significant results.
volc_sig_plot <- ggplot(res_sig, aes(x = log2FoldChange,
                                     y = -log10(padj))) + 
  geom_point(pch = 21, fill = "magenta") +
  ggtitle(str_wrap("Volcano Plot of Significantly Differentially Expressed Genes",width=30)) +
  theme_bw(base_size = 20)
ggsave("/home/hice1/cwijeyesekera3/scratch/Data/Project5/DE_SigResults_Volcano.png",plot=volc_sig_plot)




data <- read.csv("~/Charith/GaTech R/DE_results_filtered.csv")
selected_transcripts <- data %>%
  arrange(padj) %>%
  head(2)



#Gets the top 2 Differentially Expressed Genes
data <- read.csv("/home/hice1/cwijeyesekera3/scratch/Data/Project5/DE_results_filtered.csv")
selected_transcripts <- data %>%
  arrange(padj) %>%
  head(2)
transcript_names <- selected_transcripts[, 1]

#Gets Normalized count data
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
ncts <- counts(dds,normalized=True)
# Gets normalized count data for two genes 
selected_genes_data <- ncts[ncts$transcript %in% transcript_names, ]
Control_data <- selected_genes_data[, 1:10]
Disease_data <- selected_genes_data[, 11:20]

#Combining Data
combined_data <- bind_rows(
  mutate(control_data, Condition = "Control"),
  mutate(disease_data, Condition = "Disease")
  
  # Create box plot
  ggplot(combined_data, aes(x = Condition, y = value)) +
    geom_boxplot() +
    geom_jitter(position = position_jitter(0.2), color = "blue") +
    labs(title = "Comparison of Count Data between Control and Prostate Cancer Data",
         x = "Condition",
         y = "Normalized Counts")
