library(limma)
library(edgeR)
library(tximport)
library(tidyverse)

##read in transcript ID to gene ID table
tx2gene <- read.table("directory/transcript_id_gene_id.tsv", header = FALSE)

##read in counts for samples and convert transcript IDs to gene IDs - number of samples can be changed based on comparisons
dir <- "directory/name/here"
files <- file.path(dir, 1:9, "quant.sf")
names(files) <- paste0("sample", 1:9)
txi <- tximport(files, type = "salmon", countsFromAbundance = "scaledTPM", tx2gene = tx2gene, geneIdCol = V2, txIdCol = V1)

##filter for genes which have more than 0.1 cpm for at least 2 of the samples
txi_subset_counts <- data.frame(txi$counts)
dgList <- DGEList(counts=txi_subset_counts, genes=rownames(txi_subset_counts))
countsPerMillion <- cpm(dgList)
##keep only genes with counts per million greater than 0.1 for at least 3 replicates (or 2 if using 2 replicates)
countCheck <- countsPerMillion > 0.1
keep <- which(rowSums(countCheck) >= 3)
dgList <- dgList[keep,]

##normalize using TMM
dgList <- calcNormFactors(dgList, method="TMM")

##plot PCA with sample names or points
plotMDS(dgList, gene.selection = "common")
plotMDS(dgList, gene.selection = "common", pch = 19)

##calculate variance of principal components
PCA <- prcomp(dgList)
summary(PCA)

#Input of design and model with target file
targets <-readTargets("directory/targets.txt")
type<-factor(targets$Type)
batch<-factor(targets$Batch)
##for no batch correction
design<-model.matrix(~0+type)

#voom normalization
y<-voom(dgList,design,plot=TRUE, normalize.method="none")
fit<-lmFit(y,design)
fit2<-eBayes(fit)
out<-topTable(fit2,number=nrow(fit2), sort.by = "F", genelist=fit2$genes)
write.table(out,file="limma_voom_0.1cpm.txt",row.names = FALSE, sep="\t",quote=F)
#coefficients in fit2 correspond to log2cpm values for each sample

#Contrasts for comparing any two conditions - these will give you adjusted p-values
#Wt vs Mut
cont_1<-makeContrasts("typeMut-typeWt",levels=design)
fit_1<-contrasts.fit(fit,cont_1)
fit_1<-eBayes(fit_1)
output_1<-topTable(fit_1,adjust="BH",number=nrow(fit_1))
write.table(output_1,file="contrast_mut_vs_wt_0.1cpm.txt",row.names = FALSE, sep="\t",quote=F)
