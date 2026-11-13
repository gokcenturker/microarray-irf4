install.packages("BiocManager") 
# this allows me to install the tool called BiocManager,
# which is used for installing R packages for biology

BiocManager::install(c("GEOquery", "limma"))
# Here I use the tool I installed in the previous line to install GEOquery, which 
# will download datasets from GEO and limma for statistical analysis

install.packages("ggplot2")
# this plotting package will be used for volcano plot

install.packages("pheatmap")
# this package will be used to draw heatmaps

# in order to check whether all the packages are installed, we can run the following lines in the console:
# library(GEOquery)
# library(limma)
# library(ggplot2)
# library(pheatmap)