require("EpiDISH")
require("CellMix")
require("estimate")
require("MCPcounter")
require("xCell")
require(pheatmap)
require(magrittr)
library(wesanderson)

# data in the form of a matrix, gene names in row names and sample names in column names
combinedEstimates <- function(data, LM22 = NULL, AbbasGenes = NULL, datadir = NULL, titleprefix = "", 
                              outdir = getwd(), fileprefix = outdir, outfn = FALSE,
                              annotfactor = NULL, annotname = "", sampnames = colnames(data)) {
  if(missing(datadir) & (missing(AbbasGenes) | missing(LM22))) { 
    stop("Must specify both LM22 & AbbasGenes, or datadir that contains the RDS files") 
  }
  if(class(annotfactor) != "factor") {
    stop("annotfactor is not a factor")
  }
  if(fileprefix != "") {
    fileprefix %<>% paste0("_")
  }
  if(outfn == "" | outfn == TRUE) {
    outfn <- paste0(fileprefix, "combined.txt")
  }
  # get LM22 and AbbasGenes
  if(is.null(datadir)) {
    if(class(LM22) == "character") {
      LM22=readRDS(LM22)
    }
    if(class(AbbasGenes) == "character") {
      AbbasGenes=readRDS(AbbasGenes)
    }
  }
  else {
    LM22=readRDS(file.path(datadir, "LM22.RDS"))
    AbbasGenes=readRDS(file.path(datadir, "AbbasGenes.RDS"))
  }
  
  # function to pick regular or annotated pheatmap function based on whether annotation is specified
  pickpheatmap <- function(annotcheck, toplot, titlesuffix) {
    plottitle <- paste0(titleprefix, titlesuffix)
    if(is.null(annotcheck)) {
      pheatmap(toplot, cluster_cols = F, main = plottitle, treeheight_row = 0)
    }
    else {
      annotatedpheatmap(toplot, annotfactor, annotname, sampnames, main = plottitle)
    }
  }
  
  ##### NNLS from CellMix #####
  message("NNLS")
  cm=intersect(rownames(data), rownames(AbbasGenes))
  nnlsLM22=ged(data[cm,], AbbasGenes[cm,],'lsfit', fit = 'nnls')
  # nnlsLM22@fit@H %>% dim
  pickpheatmap(annotfactor, nnlsLM22@fit@H, "Immune infiltration: NNLS")
  
  
  ##### CIBERSORT from EpiDISH #####
  message("CIBERSORT")
  cm=intersect(rownames(data), rownames(LM22))
  CS22=epidish(data[cm,], LM22[cm,], method="CBS")
  # t(CS22$estF) %>% dim
  pickpheatmap(annotfactor, t(CS22$estF), "Immune infiltration: CIBERSORT")
  
  
  ##### ESTIMATE #####
  message("ESTIMATE")
  data[unique(rownames(data)), ] %>%
    write.table(file=file.path(outdir,paste0(fileprefix, "data.txt")), quote = F, sep = '\t')
  # unify different number of genes per microarray platforms to 10,412 common genes
  filterCommonGenes(
    input.f=file.path(outdir,paste0(fileprefix, "data.txt")), 
    output.f=file.path(outdir,paste0(fileprefix, "10412genes.gct")), 
    id="GeneSymbol")
  # get stromal, immune, and ESTIMATE scores
  estimateScore(
    file.path(outdir,paste0(fileprefix, "10412genes.gct")), 
    file.path(outdir,paste0(fileprefix, "estimate_score.gct")), 
    platform="illumina")
  # read in results
  data_estimate <- read.table(file.path(outdir,paste0(fileprefix, "estimate_score.gct")), header = T, sep = '\t', quote = "", skip = 2)
  rownames(data_estimate) <- data_estimate$NAME
  data_estimate %<>% extract( , c(-1,-2))
  colnames(data_estimate) %<>% gsub("^X", "", .)
  # plot
  pickpheatmap(annotfactor, data_estimate, "Immune infiltration: ESTIMATE")
  
  
  ##### MCP-counter #####
  message("MCP-counter")
  data_mcp <- MCPcounter.estimate(data,featuresType="HUGO_symbols")
  # visualize
  pickpheatmap(annotfactor, data_mcp, "Immune infiltration: MCP-counter")
  
  
  ##### xCell #####
  message("xCell")
  if(nrow(data) < 5000) {
    data_xcell <- NULL
    message("Not enough genes for xCell analysis")
  }
  else {
    data_xcell <- xCellAnalysis(data)
    # visualize
    pickpheatmap(annotfactor, data_xcell, "Immune infiltration: xCell")
  }
  
  
  ##### Combine the datasets #####
  print(all(sapply(list(colnames(t(CS22$estF)), colnames(data_estimate), colnames(data_mcp), colnames(data_xcell)), FUN = identical, colnames(nnlsLM22@fit@H))))
  combined <- rbind(NNLS=nnlsLM22@fit@H, CIBERSORT=t(CS22$estF), ESTIMATE=data_estimate, MCPcounter=data_mcp, xCell=data_xcell)
  if(outfn != FALSE) {
    write.table(combined, file = file.path(outdir,outfn), sep = "\t")
  }
  combined
}

annotatedpheatmap <- function(toplot, annotateVar, varName, sampNames, ...) {
  if(length(annotateVar) != length(sampNames)) {
    stop("Length of annotation variable does not match the number of samples.")
  }
  
  # annotation data frame
  hm.annotate <- annotateVar %>%
    as.factor %>%
    list %>% setNames(varName) %>%
    do.call(data.frame, args = .)
  rownames(hm.annotate) <- sampNames
  # annotation colors
  # mypalette <- brewer.pal(length(levels(annotateVar)), "Set3")
  mypalette <- wes_palette(n=length(levels(annotateVar)), name="Moonrise3")
  mycolors <- mypalette %>%
    setNames(levels(annotateVar)) %>%
    list %>% setNames(varName) %>%
    do.call(list, args = .)
  
  pheatmap(
    toplot, 
    # cluster_cols = F, 
    treeheight_row = 0,
    annotation_col = hm.annotate,
    annotation_colors = mycolors,
    ...
  )
}
