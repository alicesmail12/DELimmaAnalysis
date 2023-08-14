install.packages("BiocManager")
install.packages("tidyverse")
install.packages("dplyr")

BiocManager::install(version = "3.17")
BiocManager::install("Biobase")
BiocManager::install("limma")
BiocManager::install("GO.db")
install.packages("org.Hs.eg.db", repos="http://bioconductor.org/packages/3.17/data/annotation")

library(tidyverse)
library(dplyr)
library(Biobase)
library(limma)

setwd("/Users/alicesmail/Desktop/Programming/DE Analysis")

# Import data
x_df = read.csv("expression_matrix.csv")
x = data.matrix(x_df[,-1])
rownames(x) = x_df[[1]]

f = read.csv("feature_data.csv")
f <- f %>%
  column_to_rownames('probe')

p = read.csv("phenotype_data.csv")
p <- p %>%
  column_to_rownames('sample')

# DE boxplot for 1 gene
boxplot(x[27,] ~ p[, "disease"], main=f[27, "symbol"])

# Create ExpressionSet
eset <- ExpressionSet(assayData = x, 
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))
dim(eset)

# DE boxplot from eset
boxplot(exprs(eset)[27, ] ~ pData(eset)[, "disease"],
        main = fData(eset)[27, "symbol"])

# Create linear model
design <- model.matrix(~0 + disease, data=pData(eset))
head(design)
colSums(design)
table(pData(eset)[, "disease"])

# Testing with limma
cm <- makeContrasts(status = diseaseyes - diseaseno,
                    levels=design)
cm

# Fit model
fit <- lmFit(eset, design)
head(fit$coefficients)

fit2 <- contrasts.fit(fit, contrasts = cm)
head(fit2$coefficients)

# Calculate t-stat
fit2 <- eBayes(fit2)
fit2

# Results
results <- decideTests(fit2)
summary(results)

# Plot distribution
plotDensities(eset)

# Log transform
exprs(eset) <- log(exprs(eset))
plotDensities(eset, legend=FALSE)

# Quantile normalise
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, legend=FALSE)

# PCA
plotMDS(eset, labels=pData(eset)[, "disease"],
        gene.selection="common")

# Top 3 DE genes
topTable(fit2, number=3)

# Obtain all results
stats <- topTable(fit2, number=nrow(fit2), sort.by="none")
dim(stats)

# Histogram of p-values
hist(stats[, "P.Value"])

# Volcano plot
volcanoplot(fit2, highlight=10,
            names=fit2$genes[, "symbol"],
            xlab=)

# GO enrichment
enrich_go <- goana(fit2, geneid=entrez, species="Hs")
topGO(enrich_go, ontology="BP", number=10)






