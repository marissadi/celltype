# data: gene (rows) x sample (cols) matrix

combinedEstimates(
      data = data, 
      LM22 = readRDS("LM22.RDS"),
      AbbasGenes = readRDS("AbbasGenes.RDS"),
      outdir = "immune",
      annot = annotation$Factor,
      annotname = "Factor name",
      sampnames = annotation$sampleID,
      titleprefix = "Plot title\n",
      fileprefix = "immune-output"
)
