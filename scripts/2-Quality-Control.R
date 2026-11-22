# Before starting the statistical analysis, I will check the data quality

library(GEOquery)
library(pheatmap)
d <- getGEO("GSE22036", GSEMatrix = TRUE) [[1]]
expr <- exprs(d)
meta <- pData(d)

summary(as.vector(expr))
# this line gives a statistical summary of the dataset
# microarray data is normally on a log2 scale
# however we see a maximum much larger than log2 scale
# this shows me that the data hasn't been log-transformed

expr <- log2(expr)
summary(expr)
# raw microarray values are skewed, yet statistical tests we wil apply mainly requires 
# symmetrical datasets
# fold changes make the distribution more symmetric, therefore the variance is less dependent on the mean.

meta[, c("title", "geo_accession")]
# this gives me the all rows and combination of "title" and "geo_accession" columns of the meta table
# used this do decide on the groups

# grepl means "Does this text exist appear inside?"
# ifelse is a decision maker with the pattern: ifelse(question, answer if TRUE, answer if FALSE).

has_irf4 <- grepl("IRF4", meta$title) 
labels <- ifelse(has_irf4, "knockdown", "control")
# for each group, this code checks whether it is knockdown or control
# so far the result is text

group <- factor(labels, levels = c("control", "knockdown"))
# factor stores catgeroies. it basically know that there are only two possible values, and it knows their order.
# levels shows the categories and the order. 

table(group)

# After we see 2 and 2 for knockdown and control, we continue controls with checking if four samples look comparable.
# For this we produce two plots.

# BARPLOT
    cols <- ifelse(group == "knockdown", "firebrick", "steelblue")
    # we created a list of colors. If the answer to KD is true, it becomes red, if answer is no, it becomes blue
    
    boxplot(expr, las = 2, col = cols, main ="Expression per sample (log2)")
    # las rotates the sample names so they are readable
    # main is the title
    # col uses the colors we created above
    
    # we want the four boxes to look about the same height and position. 
    # If there is a box clearly higher or lower than others, this can be a warning

# DENSITY PLOT
    plot(density(expr[, 1]), col = cols[1], ylim = c(0, 0.35),
         main = "Density per sample")                         # estimates the shape of distribution of the first sample
    
    # we want to repeat this for all 4 samples, so create a loop
    
    for (i in 2:ncol(expr)) lines(density(expr[, i]), col =cols[i])
    legend("topright", legend = colnames(expr), col = cols, lty = 1, cex = 1)
    # ncol(expr) counts the columns of expr, so it returns 4
    # 2:4 means the numbers 2,3,4. we start at 2 because we already plotted the first sample
    # the instruciton starting with lines.. is repeated
    # density() estimates the shape of the sample distribution
    # lines() adds that curve to the existing plot, so plot() starts a new one while lines() draws to the existing plot
    # col = cols[i], from the color vector we created, it assigns the color corresponding to i-th item. red for KD and blue for control
    
    # How to read the density plot? four curves should mostly overlap
    # first bump on the left is the probes with very little or no expression, meaning the background noise
    # we should note where the first bump ends, because we will use it for further control
    # in our data, it is around 4

# FILTER LOW-EXPRESSION PROBES
    
    threshold <- 4
    keep <- rowSums(expr > threshold) >= 2
    table(keep)
    expr_filt <- expr[keep, ]
    dim(expr_filt)
    
    # expr > threshold gives TRUE/FALSE for every value
    # rowSums() counts in how many samples it is above the threshold for each probe
    # >= 2 keeps probes that are above the threshold in at least 2 samples. So that a gene expressed in only one group is not lost
    # expr[keep, ] keeps only those row.
    # dim() gives the new number of rows and columns
   
    # as a result 13442 probes were removed (FALSE) and 41234 were kept (TRUE)
  
# PCA
    # Principle Component Analysis takes the 41.234 probes and squeezes them into a few summary axes. This enables us to see in a
    # 2D plot, whether the samples group the way we expect
    pca <- prcomp(t(expr_filt), scale. = FALSE) # transposes the matrix so that samples are rows, and probes are columns. 
                                                  #prcomp needs samples as rows
                                                  #scale. = FALSE: do not rescale each probe. The data are already log2.
    
    var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)  
        # pca$sdev^2 is the variance of each PC. Dividing by the total and multiplying by 100 gives the percanteg of variation of each PC.
        # round(..., 1) keeps one decimal.
    var_exp

    plot(pca$x[, 1], pca$x[, 2], pch=19, cex=2, col=cols,   # are the positions of the samples on PC1 and PC2, pch =19 is solid circle in plot
         xlab = paste0("PC1 (", var_exp[1], "%)"),          # paste0 means sticking text together, to organize the axis labels
         ylab = paste0("PC2 (", var_exp[2], "%)"),
         main = "PCA of filtered data")
    text(pca$x[, 1], pca$x[, 2], labels = meta$title, pos =3, cex=1)  # writes each sample's title next to its point, pos =3 puts it above

# IRF KNOCKDOWN CHECK
    # we will investigate whether the experiment worked at the gene level
    colnames(fData(d))  # d was our GEO object.  
                          # fData is the annotation table.Each row is a probe and each column is an information about it
                          # colnames prints the column names
    fdat <- fData(d) # stored the annotation table in a variable called fdat
    irf4_probes <- rownames(fdat)[grepl("^IRF4$", fdat$'Gene Symbol')] # fdat$'Gene Symbol' takes the gene symbol column
                                                                      # grepl("^IRF4$",...) asks does the gene symbol match this pattern in each probe
                                                                      # ^ means text starts here, $ means text ends here. So we look for exact IRF4 match not sth like IRF4-GT11
                                                                      # rownames(fdat)[] takes probe ids, only for rows where the answer was TRUE
    irf4_probes
    round(expr[irf4_probes, , drop = FALSE], 2) # drop = FALSE keeps the result as a small table. Without it, if there were only one probe
                                              # R would turn it into a plain list and lose the sample names.