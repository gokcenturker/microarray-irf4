library(GEOquery)
# loading the package so that R can use its functions
# after every session we should load the packages again

d <- getGEO("GSE22036", GSEMatrix = TRUE) [[1]]

# getGEO function is used to import the data we will use from GEO so we type its accession number
# GSEMatrix = TRUE means that we want the processed data instead of the raw data derived from the scanner
# [[1]] : GEO might return several datasets. This happens even there is only one dataset. So we ask for
# the first item of the box that might contain mulitple items/datasets

# At this point d carries all the information

expr <- exprs(d)
# this function takes the expression matrix from d and stores it in expr. 
# Each row represents a probe
# Each column represents a sample
# Each number tells how strongly the gene was detected in that sample

meta <- pData(d)
# this function is used to extract the "phenotype data"
# each row is a sample

dim(expr)
# function gives the dimensions of the matrix. It tells us about the number of probes and samples.
# We use this to check if we extracted the correct data

head(meta)
# gives sample information
colnames(expr)
# brings the names of the columns in expression matrix, which shows the sample names

