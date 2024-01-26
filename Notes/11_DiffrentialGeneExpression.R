library("tximportData")
library("tximport")
library("readr")
library("DESeq2")
library("dplyr")
library("GenomicFeatures")
library("ggplot2")
library("apeglm")


# Define the directory for reference files. TximportData already has these for us for Gencode reference files
# We will make specific reference files using RefSeq gff3 files and GenomicFeatures package.
txdb <- makeTxDbFromGFF(file="/storage/ice-shared/biol6150/Data/DiffrentialExpression/Reference/GRCh38_latest_genomic.gff.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene = tx2gene %>% filter(!is.na(TXNAME))
head(tx2gene)
write.table(tx2gene, "/storage/ice-shared/biol6150/Data/DiffrentialExpression/Reference/tx2gene.RefSeq.All.tsv", 
            sep = ",", 
            row.names = FALSE)

# Read tx2gene if already created before.
tx2gene <- read_csv("/storage/ice-shared/biol6150/Data/DiffrentialExpression/Reference/tx2gene.RefSeq.All.tsv") %>% as.data.frame()
head(tx2gene)

# Read the samples map file.
sample_info <- read.table("/storage/ice-shared/biol6150/Data/DiffrentialExpression/Quant/SampleInfo.txt", header=TRUE) 
samples = sample_info %>% pull(SampleID)
print(samples)

file_paths = paste0("/storage/ice-shared/biol6150/Data/DiffrentialExpression/Quant/", samples, ".sf")
print(file_paths)

# Create a txi object. Read Salmon SF files.
txi.salmon <- tximport(file_paths, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
print(head(txi.salmon$counts))


# Deseq needs specific format for the map file.
sampleTable <- data.frame(condition = factor(c(rep("Treated",3), rep("Untreated",3))
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

#Volcano plot for just the significant results.
ggplot(res_sig, aes(x = log2FoldChange,
                y = -log10(padj))) + 
  geom_point(pch = 21, fill = "magenta") +
  ylim(c(0,100)) +
  xlim(c(-15,15)) +
  theme_bw(base_size = 24)


## Log fold change shrinkage. This is a different package which performs shrinkage of log2fold.
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_Untreated_vs_Treated", type="apeglm")
head(resLFC)

resLFC_sif = resLFC %>% as.data.frame() %>% filter(padj < 0.05)
dim(resLFC_sif)
head(resLFC_sif)

#Volcano plot for just the significant results.
ggplot(resLFC_sif, aes(x = log2FoldChange,
                    y = -log10(padj))) + 
  geom_point(pch = 21, fill = "magenta") +
  ylim(c(0,100)) +
  xlim(c(-15,15)) +
  theme_bw(base_size = 24)

