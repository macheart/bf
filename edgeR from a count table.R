# Here, I access a count table and a pheno data from the R-Bioinformatics-Cookbook (2019 Packt Publishing) github reference to:
# 1: Create a count matrix with 9 samples and 14k+ Flybase genes && group the samples based on fly larvae stage
# 2: Collect the grouping variables and the count matrix into a DGEList object required for edgeR
# 3: Apply filterByExpr(), and calcNormFactors() to filter and normalize the dataset respectively.
# 4: Create a model `design` vector for application (w/ the DGEList object) in the i) glmQLFit() "model fitting" and ii) glmQLFTest() "model testing" 
# iii) Save the `result` "model test" vector into a statistical summary `res` vector
# 5: Append ENTREZID, gene symbols, and names to `res` object

# The focus of this assignment was to establish code for conventional differential expression assessments using edgeR 
# Code applicable for subsequent gene ontology or pathway analysis on gene lists is accessible in a different repository.

#### 1 Load the count data
library(readr)
count_dataframe <- readr::read_tsv(file.path(getwd(),"datasets","Chapter01","modencodefly_count_table.txt"))
genes <- count_dataframe[['gene']]
count_dataframe[['gene']] <- NULL
count_matrix <- as.matrix(count_dataframe)
rownames(count_matrix) <- genes
pheno_data <- readr::read_table2(file.path(getwd(),"datasets","Chapter01","modencodefly_phenodata.txt"))

## Specify experiments of interest:
experiments_of_interest <- c("L1Larvae","L2Larvae")
columns_of_interest <- which(pheno_data[['stage']] %in% experiments_of_interest)

## Form the grouping factor:
library(magrittr)
grouping <- pheno_data[['stage']][columns_of_interest] %>% forcats::as_factor()

## Form the subset of count data:
counts_of_interest <- count_matrix[,columns_of_interest] ###### matrix 9 samples, 14k+ FBGn

#### 2 Create the DGE object:
library(edgeR)
count_dge <- edgeR::DGEList(counts=counts_of_interest,group=grouping) # 14869 9

#### 3 preliminary filtering and normalization
keep=filterByExpr(count_dge) # remove gene counts below the count threshold
count_dge_filt=count_dge[keep,,keep.lib.sizes=FALSE] # filter the count matrix 10111 9
count_dge_filt_norm=calcNormFactors(count_dge_filt) # Normalization of samples

#### 4 Perform differential expression analysis
design <- model.matrix(~grouping) # sample classification
eset_dge <- edgeR::estimateDisp(count_dge_filt_norm,design) # provides AveLogCPM and dispersion estimates 

fit <- glmQLFit(eset_dge,design) # Fit a quasi-likelihood negative binomial generalized log-linear model to count data, option: robust=TRUE
result<-edgeR::glmQLFTest(fit,coef=2) # test of model
res=as.data.frame(topTags(result)) # result: logFC, logCPM, F, PValue, FDR

#### 5 Append ENTREZID, gene symbols, and names to `res` 
library(limma)
library(GO.db)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(annotation)
res$symbol = mapIds(org.Dm.eg.db,
                    keys=row.names(res),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

res$entrez = mapIds(org.Dm.eg.db,
                    keys=row.names(res),
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name= mapIds(org.Dm.eg.db,
                 keys=row.names(res),
                 column="GENENAME",
                 keytype="ENSEMBL",
                 multiVals="first")
write.csv(res,file="res_GO.csv")

