#!/usr/bin/Rscript

##' @name DE.analysis
##' @rdname DE.analysis
##' @title Automated Differential Expression Pipeline for RNA-Seq Counts using DESeq2 or edgeR
#' @param M Character. Choose method, either "DESeq2" or "edgeR".
#' @param IN A numeric matrix or data frame of raw RNA-Seq read counts (file path).
#' @param S Numeric pair. The number of samples in the dataset (x:y).
#' @param OUT Character. Prefix for the exported text and PDF files.
#' @param LF Numeric. Log fold change threshold
#' @param FD Numeric. adjusted p-value threshold.
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics abline hist legend
#' @importFrom stats dist prcomp
#' @importFrom utils head read.delim write.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
#' @importFrom genefilter rowVars
#' @importFrom SummarizedExperiment colData assay
#' @importFrom limma decideTests plotMDS
#' @importFrom data.table data.table
#' @import DESeq2
#' @import edgeR
#' @import dplyr
#' @import magrittr
#' @import calibrate
#' @export
#'
utils::globalVariables(c("Sequence"))

DE.analysis <- function(M, IN, S, OUT = getOption("OUT", "DE_analysis"), LF = getOption("LF", "NONE"), FD = getOption("FD", "NONE")){
args = commandArgs(trailingOnly=TRUE)
  if (M == "DESeq2"){

    Read_Matrix <- read.delim(IN, header = TRUE)
    Read_Matrix <- as.data.frame(Read_Matrix)

    S <- S

    Sample = as.integer(strsplit(S, ":")[[1]])
    S_1 <- Sample[1]
    S_2 <- Sample[2]

    metaData <- data.table(id = colnames(Read_Matrix[,2:ncol(Read_Matrix)]),
                           dex = c(rep("Control", S_1), rep("Treated", S_2)))

    Read_col <- length(colnames(Read_Matrix[,2:ncol(Read_Matrix)]))
    Sample_col <-length(c(rep("Control", S_1), rep("Treated", S_2)))

    if (Read_col != Sample_col){
      stop("Please Correct the sample numbers! Their sum must be ", Read_col, ", corresponding to the read columns in provided input file.\n",
           call.=FALSE)
    }

    dds <- DESeqDataSetFromMatrix(countData = Read_Matrix,
                                  colData = metaData,
                                  design = ~dex, tidy = TRUE)

    dds <- DESeq(dds)

    res <- results(dds)
    res <- as.data.frame(res)

    res1 <- cbind(Sequence=row.names(res), res)
    row.names(res1) <- NULL

    write.table(res1, file= paste0(OUT, "_DESeq2.txt"), row.names = FALSE, sep="\t", quote=FALSE)

    if(LF == "NONE" && FD == "NONE"){
      print ("Attention: As no options are provided for filter,  the DESeq2 result will be filtered for LFC value of 1/-1 and adj-pvalue of 0.05.")
      up_res1 <- subset(res1, res1$log2FoldChange >= 1 & res1$padj <= 0.05)
      down_res1 <- subset(res1, res1$log2FoldChange <= -1 & res1$padj <= 0.05)
      DF_res1 <- rbind(up_res1, down_res1)
    } else if(LF == "NONE"){
      PV = as.numeric(FD)
      DF_res1 <- subset(res1, res1$padj <= PV)
    }else if(FD == "NONE"){
      LFC = as.numeric(LF)
      Negative_LFC = LFC * (-1)
      up_res1 <- subset(res1, res1$log2FoldChange >= LFC)
      down_res1 <- subset(res1, res1$log2FoldChange <= Negative_LFC)
      DF_res1 <- rbind(up_res1, down_res1)
    } else{
      LFC = as.numeric(LF)
      PV = as.numeric(FD)
      Negative_LFC = LFC * (-1)
      up_res1 <- subset(res1, res1$log2FoldChange >= LFC & res1$padj <= PV)
      down_res1 <- subset(res1, res1$log2FoldChange <= Negative_LFC & res1$padj <= PV)
      DF_res1 <- rbind(up_res1, down_res1)
    }
    write.table(DF_res1, file= paste0(OUT, "_DESeq2_LFC_", LF, "_and_padj_", FD, ".txt"), row.names = FALSE, sep='\t', quote=FALSE)

    ##Generate plots

    pdf(paste0(OUT, "_DESeq2_plots.pdf"))

    if (nrow(DF_res1) >= 1){
      hist(DF_res1$pvalue, breaks=20, col="grey", xlab = "p-value", main = "Histogram of p-vlaue of DE sequences")
      hist(DF_res1$log2FoldChange, breaks=10, col="grey", xlab = "log2FoldChange",  main = "Histogram of log2FoldChange of DE sequences")
    }

    plotDispEsts(dds, ylim = c(1e-6, 1e1), main = "Dispersion Plot")

    res1_max_LFC <- max(res1$log2FoldChange, na.rm = TRUE)
    res1_min_LFC <- min(res1$log2FoldChange, na.rm = TRUE)
    plotMA(dds, ylim = c(res1_min_LFC, res1_max_LFC), main = "MA Plot" )

    rld <- rlog(dds)
    # Heatmap with gene clustering
    topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ), 20 )
    heatmap.2( assay(rld)[ topVarGenes, ], scale="row", srtCol = 45, trace="none", dendrogram="column", margins=c(5, 20), col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100), main = "                            Heatmap showing clustering of \n                     top-20 most variable genes", keysize = 1.1)

    # produce rlog-transformed data
    rld <- rlogTransformation(dds, blind=TRUE)
    # create a distance matrix between the samples : hierarchical clustering of the samples.
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(mat, trace="none", srtCol = 45, col = rev(hmcol), margins=c(13, 13), main = "                            Heatmap showing clustering of the samples", keysize = 1.1)

    rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="topright", main="PCA plot", textcx=1, ...) {
      #require(genefilter)
      #require(calibrate)
      #require(RColorBrewer)
      rv = rowVars(assay(rld))
      select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
      pca = prcomp(t(assay(rld)[select, ]))
      fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
      if (is.null(colors)) {
        if (nlevels(fac) >= 3) {
          colors = brewer.pal(nlevels(fac), "Paired")
        }   else {
          colors = c("black", "red")
        }
      }
      pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
      pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
      pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
      pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
      plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
      with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
      legend(legendpos, legend=levels(fac), col=colors, pch=20)
    }

    rld_pca(rld, intgroup="dex", xlim=c(-25, 25))

    dev.off()

  }
  else if (M == "edgeR"){

    Read_Matrix <- read.delim(IN, header = TRUE)

    #Remove duplicate entries, if any
    Read_Matrix <- Read_Matrix %>% group_by(Sequence) %>%
    summarise_at(c(colnames(Read_Matrix[,2:ncol(Read_Matrix)])),sum, na.rm = TRUE)

    seq_num <- nrow(Read_Matrix)

    Sample = as.integer(strsplit(S, ":")[[1]])
    S_1 <- Sample[1]
    S_2 <- Sample[2]

    Read_col <- length(colnames(Read_Matrix[,2:ncol(Read_Matrix)]))
    Sample_col <-length(c(rep("Control", S_1), rep("Treated", S_2)))

    if (Read_col != Sample_col){
      stop("Please Correct the sample numbers!\n",
           call.=FALSE)
    }

    dex = c(rep("Control", S_1), rep("Treated", S_2))

    data_edgeR <- DGEList(counts = Read_Matrix, group =factor(dex))

    normalized_data_edgeR <- calcNormFactors(data_edgeR)

    Disp <- estimateCommonDisp(normalized_data_edgeR, verbose=T)
    Disp <- estimateTagwiseDisp(Disp)

    DE_edgeR_et <- exactTest(Disp, pair=c(1,2))
    DE_edgeR_et_results <- topTags(DE_edgeR_et, n = seq_num)

    DE_edgeR_et_results <- as.data.frame(DE_edgeR_et_results)

    write.table(DE_edgeR_et_results, file= paste0(OUT, "_edgeR.txt"), row.names = FALSE, sep="\t", quote=FALSE)

    DE_edgeR_et_Filter <- decideTests(DE_edgeR_et, adjust.method="BH", p.value=0.05)

    if(LF == "NONE" && FD == "NONE"){
      print ("Attention: As no options are provided for filter,  the edgeR result will be filtered for LFC value of 1/-1 and FDR value of 0.05.")
      up_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$logFC >= 1 & DE_edgeR_et_results$FDR <= 0.05)
      down_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$logFC <= -1 & DE_edgeR_et_results$FDR <= 0.05)
      DF_res1 <- rbind(up_res1, down_res1)
    } else if(LF=="NONE"){
      FDR_value = as.numeric(FD)
      DF_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$FDR <= FDR_value)
    } else if(FD=="NONE"){
      LFC = as.numeric(LF)
      Negative_LFC = LFC * (-1)
      up_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$logFC >= LFC)
      down_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$logFC <= Negative_LFC)
      DF_res1 <- rbind(up_res1, down_res1)
    } else{
      LFC = as.numeric(LF)
      FDR_value = as.numeric(FD)
      Negative_LFC = LFC * (-1)
      up_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$logFC >= LFC & DE_edgeR_et_results$FDR <= FDR_value)
      down_res1 <- subset(DE_edgeR_et_results, DE_edgeR_et_results$logFC <= Negative_LFC & DE_edgeR_et_results$FDR <= FDR_value)
      DF_res1 <- rbind(up_res1, down_res1)
    }
    write.table(DF_res1, file= paste0(OUT, "_edgeR_LFC_", LF, "_and_FDR_", FD, ".txt"), row.names = FALSE, sep='\t', quote=FALSE)


    pdf(paste0(OUT, "_edgeR_plots.pdf"))

    if (nrow(DF_res1) >= 1){
      hist(DF_res1$PValue, breaks=20, col="grey", xlab = "p-value", main = "Histogram of p-vlaue of DE sequences")
      hist(DF_res1$logFC, breaks=10, col="grey", xlab = "log2FoldChange",  main = "Histogram of log2FoldChange of DE sequences")
    }

    plot <- plotMDS(normalized_data_edgeR, method ="logFC", col=as.numeric(normalized_data_edgeR$samples$group))
    legend("topright", as.character(unique(normalized_data_edgeR$samples$group)), col=1:2, pch=20)

    plotBCV(Disp)

    de1tags12 <- rownames(Disp)[as.logical(DE_edgeR_et_Filter)]
    plotSmear(DE_edgeR_et, de.tags=de1tags12)
    abline(h = c(-2, 2), col = "blue")

    dev.off()

  }
  else {
    stop("Specifiy correct method either 'DESeq2' or 'edgeR'!\n",
         call.=FALSE)
  }
  message("Analysis finished successfully! Check your working directory for text and PDF outputs.")
}
