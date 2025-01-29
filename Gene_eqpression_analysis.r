#setwd("~/Videos/mahm_microarray")
#######Directory to perform analysis.

# micro array analysis of  prostate cancer(PCa) and bone metastasis data.
## micro array machine GPL8469
# GSE26964
# Mahima Bhardwaj
#dt. 21/01/2023
################################################################
#  Differential expression analysis with limma
library(GEOquery)
library(limma)
#install.packages("umap")
library(umap)
library(tidyverse)
library(ggplot2)
library(affy)
library(dplyr)
library(DESeq2)
# load series and platform data from GEO

gset <- getGEO("GSE43332", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("00110011110000")
sml <- strsplit(gsms, split="")[[1]]

#another way
#gsms <- "111000"
#sml <- strsplit(gsms, split="")[[1]]
#eleminate sample marked as x
sel<-which(sml!="X")
sml<-sml[sel]
gset<-gset[,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
ex = exprs(gset)
head(ex)
#write.csv(exprs(gset),"D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/14776_matrix.csv")

#head(fData(gset))
#head

features <- fData(gset)
colnames(features)
head(features$mrna_assignment)
dim(features)
write.table(head(features))

features <-  select(features,"SPOT_ID")

full_output <- cbind(features,exprs(gset))
head(full_output)
write.table(head(full_output))
write.csv(full_output,"D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/43332_matrix.csv")

#normalisation
expression_data1 <- read.csv("32269_matrix.csv", header = TRUE)
data_unique39494 <- distinct(expression_data, Gene, .keep_all = TRUE)
write.csv(data_unique,"D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/32269_unique.csv")
expression_39494 <- data_unique39494 %>% column_to_rownames("Gene")
normalized_39494 <- normalizeQuantiles(expression_39494)

#write.table(head(expression_data))

write.csv(normalized_39494,"D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/normalized_32269.csv")

expression_data2 <- read.csv("55715_matrix.csv", header = TRUE)
data_unique55715 <- distinct(expression_data, Gene, .keep_all = TRUE)
expression_55715 <- data_unique55715 %>% column_to_rownames("Gene")
normalized_55715 <- normalizeQuantiles(expression_55715)

write.csv(normalized_55715,"D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/normalized_55715.csv")


#write.table(head(expression_data))
expression_data <- expression_data %>% column_to_rownames("X")
normalized_data1 <- normalizeQuantiles(expression_data)
write.csv(normalized_data,"D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/14776norm_matrix.csv")


# Merge the quantile normalized data
read_2034 <- read.csv("D:/Vignan_work/P.hD_2022/PhD_Bone/microarray_analysis/microarray_main/2034_matrix1.csv", header = TRUE)
read_14776 <- read.csv("14776_matrix1.csv", header = TRUE)
read_103357 <- read.csv("103357_matrix1.csv", header = TRUE)
read_55715 <- read.csv("55715_matrix1.csv", header = TRUE)
read_137842 <- read.csv("137842_matrix1.csv", header = TRUE)


merged_data_r <- merge(read_2034,read_14776,read_103357,read_55715,read_137842)



norm_39494 <- read.csv("normalized_39494.csv", header = TRUE)
norm_55715 <- read.csv("normalized_55715.csv", header = TRUE)
merged_data <- merge(norm_39494, norm_55715, by = "Gene_name", method='COMBAT')
head(merged_data)
 #  Batch Effect Correction using ComBat
# Assume 'mergedData' contains the merged and quantile normalized data

# Create a batch vector to specify the batches (platforms) for each sample
batch <- c(rep("GPL6480", nrow(norm_39494)), rep("GPL6947", nrow(norm_55715)))

# Apply ComBat for batch effect correction
install.packages("BiocManager")
BiocManager::install("sva")
library(sva)
combatData <- ComBat(merged_data, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
 

#annotation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys=rownames(gset),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")
dim(anno_palmieri)
head(anno_palmieri)

anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL)) #(filter all with NA values)
dim(anno_palmieri)

# assign samples to groups and set up design matrix

gs <- factor(sml)
groups <- make.names(c("BCa", "BM"),unique = TRUE)
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model. for linear model analysis.

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)  # 0.01 refers to the "assumed proportion of genes which are differentially expressed"
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
#IF you want to extrect all gene put
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(ex))
#write.csv(tT,"~/Videos/mahm_microarray/Output_data/Top_diffrentially_exp_gene.csv")
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","SEQUENCE","miRNA_ID_LIST","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")



# # alternatively (preferably? :))
# design_alt <- model.matrix(~ description, gset)
# fit_alt <- lmFit(gset, design_alt)
# cont.matrix.alt <- makeContrasts(descriptionG1, levels=design_alt)
# fit2_alt <- contrasts.fit(fit_alt, cont.matrix.alt)
# fit2_alt <- eBayes(fit2_alt, 0.01)
# tT_alt <- topTable(fit2_alt, adjust="fdr", sort.by="B", number=250)

#######################################################################

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "pink", border = "grey", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 2      # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))



# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE26964", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE26964", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE26964")

