getwd()
setwd("/Users/Irfan/Downloads/Worker_R_test/")

#Step 1: Read the datasets to R
file1 = "SampleTypeData.txt"
sampleTypeData = read.csv(file1, header = TRUE, sep = "\t", dec = ".")
dim(sampleTypeData)

file2 = "PairData_Trimmo_ENSEMBL_HTSeq_Edited.table"
PairData = read.csv(file2, header = TRUE, sep = "\t", dec = ".")
dim(PairData)

#Step 2: Select relevant data columns and data rows
row.names(PairData) <- PairData$DB_ID
PairData_remove2col <- PairData[, -c(1:2)] 
dim(PairData_remove2col)

n<-dim(PairData_remove2col)[1]
PairData_remove2col4rows<-PairData_remove2col[1:(n-4),] 
dim(PairData_remove2col4rows)

last4rowsPairData <- tail(PairData, n=4)
row.names(last4rowsPairData) # control data row names
# "suougibma__"            "lauQa_wol_oot__"        "dengila_ton__"          "euqinu_ton_tnemngila__"

#Step 3: Clean the column names from data matrix
library(dplyr)
library(stringr)
PairData_remove2col4rows_cleanColumns  <- PairData_remove2col4rows %>%
                                          rename_all(~ str_remove(., "\\_S.*"))

names(PairData_remove2col4rows_cleanColumns) <- gsub(x = names(PairData_remove2col4rows_cleanColumns), pattern = "\\.", replacement = "-")  

#Step 4: Reorder the sample types data
Ordered_sampleTypeData <- sampleTypeData[ order(match(sampleTypeData$Sample.ID.Corrected, colnames(PairData_remove2col4rows_cleanColumns))), ]
dim(Ordered_sampleTypeData)
row.names(Ordered_sampleTypeData) <- Ordered_sampleTypeData$Sample.ID.Corrected

#check to find out if two arrays have equal values, i-e same order in count data columns and sample data rows 
all(row.names(Ordered_sampleTypeData)==colnames(PairData_remove2col4rows_cleanColumns))

#Step 5: Load data to DESeq2 object 
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=PairData_remove2col4rows_cleanColumns, 
                              colData=Ordered_sampleTypeData[, c(4:8,15)], 
                              design= ~ Sample.type.3, tidy = FALSE)

dds

featureData <- data.frame(gene=rownames(PairData_remove2col4rows_cleanColumns))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds) #it carries metadata with information on the meaning of the columns


# Step 6: Relevel the design

dds$condition <- relevel(dds$Sample.type.3, ref = "CTRL_gut_CTRL")

# Step 7: Analyze Differential expression of genes

#Differential expression analysis
dds <- DESeq(dds)
png("qc-dispersions.png", 800, 600, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))
resultsNames(dds)


library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(dds$condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
#install.packages('gplots')
library(gplots)
png("qc-heatmap-samples.png", w=800, h=600, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(80, "black", "white"),
          ColSideColors=mycols[dds$condition], RowSideColors=mycols[dds$condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
png("qc-pca.png", 600, 600, pointsize=30)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

# Step 8: Extract results from DESeq2 Object, what sample groups are compared?

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
#Diagnostic plot
png("maplot.png", 800, 600, pointsize=25)
plotMA(res, ylim=c(-2,2))
dev.off()

# Gene clustering
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
png("gene_clustering_plot.png", 800, 600, pointsize=20)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column", margin=c(7, 7),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()



## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

# Step 9: Extract the overall view of the results

output <- as.data.frame(res) 
summary(output)


# baseMean         log2FoldChange        lfcSE            stat             pvalue     
# Min.   :      0.0   Min.   :-15.489   Min.   :0.000   Min.   :-11.649   Min.   :0.000  
# 1st Qu.:      0.0   1st Qu.:  0.125   1st Qu.:0.345   1st Qu.:  0.046   1st Qu.:0.016  
# Median :      0.0   Median :  0.657   Median :0.994   Median :  0.339   Median :0.376  
# Mean   :    219.5   Mean   :  0.835   Mean   :2.922   Mean   :  0.851   Mean   :0.439  
# 3rd Qu.:      8.2   3rd Qu.:  1.392   3rd Qu.:7.809   3rd Qu.:  1.998   3rd Qu.:0.875  
# Max.   :1000588.9   Max.   : 26.161   Max.   :7.845   Max.   : 10.206   Max.   :1.000  
# NA's   :28232     NA's   :28232   NA's   :28232     NA's   :28280  
# padj      
# Min.   :0.00   
# 1st Qu.:0.01   
# Median :0.11   
# Mean   :0.26   
# 3rd Qu.:0.46   
# Max.   :1.00   
# NA's   :39572  

#Reference: https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf

