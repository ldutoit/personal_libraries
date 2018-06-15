#This library is intended for future used with the package DESEQ2 it contains a few functions that are useful
#It can be loaded by source("DESEQ2_shortcuts.R")

     library(DESeq2)
     library(matrixStats)
    

def PCAwithScreePlot(x){}

     bds <- estimateSizeFactors(dds)
     cds=estimateDispersions(bds)
     vsd=varianceStabilizingTransformation(cds)
     plotPCA(vsd,intgroup=c("batch"))
     
     #How to get PCA scree plot?
     
     ## calculate the variance for each gene
     rv <- rowVars(assay(vsd))

     ## select the ntop genes by variance
     select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
     
     ## perform a PCA on the data in assay(x) for the selected genes
     pca <- prcomp(t(assay(vsd)[select,]))
     
     ## the contribution to the total variance for each component
     percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
     
     ##plot the "percentVar"
     scree_plot=data.frame(percentVar)
     scree_plot[,2]<- c(1:dim(cds)[2])

 	 colnames(scree_plot)<-c("variance","component_number")
     ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")

