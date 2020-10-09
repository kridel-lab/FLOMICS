doSeuratProc <- function(exp, samp){

    print(samp)

    #1. create seurat object
    object_fl <- CreateSeuratObject(counts = exp,
      project = "FL",
      min.cells = 3, min.features = 200)

    object_fl$sample <- samp

    #2. calculate mitochondrial percentage
    object_fl[["percent.mt"]] <- PercentageFeatureSet(object_fl, pattern = "^MT-")

    #3. Filter cells that have unique feature counts over 2,500 or less than 200
    object_fl <- subset(object_fl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

    #4. Normalize the Seurat objects
    object_fl <- NormalizeData(object_fl)

    #5. FindVariableFeatures
    object_fl <- FindVariableFeatures(object_fl, selection.method = "vst", nfeatures = 2000)

    return(object_fl)

}








































  # Readout proportion of mitochondrial genes
  mito.genes <- grep(pattern = "^MT-",
                     x = rownames(x = Sobj@data), value = T)
  pc.mito <- Matrix::colSums(Sobj@raw.data[mito.genes, ])/
   Matrix::colSums(Sobj@raw.data)

  # Readout proportion of ribosomal genes
  ribo.genes <- grep(pattern = "^RP",
                    x = rownames(x = Sobj@data), value = T)
  pc.ribo <- Matrix::colSums(Sobj@raw.data[ribo.genes, ])/
   Matrix::colSums(Sobj@raw.data)

  Sobj <- AddMetaData(Sobj, metadata = pc.mito, col.name = "pc.mito")
  Sobj <- AddMetaData(Sobj, metadata = pc.ribo, col.name = "pc.ribo")

  par(mfrow = c(1, 2))

  GenePlot(object = Sobj, gene1 = "nUMI", gene2 = "pc.mito", cex.use = 0.1)
  GenePlot(object = Sobj, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.1)

  par(mfrow = c(1, 1))

  # Filter cells with exceedingly high or low number of genes and high percentage of mitochondrial genes
  Sobj <- FilterCells(object = Sobj, subset.names = c("nGene", "pc.mito"),
                      low.thresholds = c(lowthresh, 0),
                      high.thresholds = c(highthresh, highmito))


  # Normalize function with default scale factor of 10000
  Sobj <- NormalizeData(object = Sobj, normalization.method = "LogNormalize",
                        display.progress = F)

  # Identify variable genes
  cat("Plotting dispersion versus average expression")

  Sobj <- FindVariableGenes(Sobj, mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.5)

  if(regress.cc %in% c(1, 2)){

    # Assign and regress cell cycle
    Sobj <- CellCycleScoring(Sobj, s.genes = s.genes,
                             g2m.genes = g2m.genes, set.ident = T)
    ## Plot cell cycle of all cells
    ggplot(Sobj@meta.data, aes(x=G2M.Score, y=S.Score, color=Phase))+
      geom_point()+theme_bw()

    if(regress.cc==1) {

      Sobj <- ScaleData(Sobj, display.progress = F, do.par = T, num.cores = 3,
                        vars.to.regress = c("S.Score", "G2M.Score", "nUMI","pc.mito"))

     }else{

      Sobj@meta.data$CC.Difference <- Sobj@meta.data$S.Score - Sobj@meta.data$G2M.Score
      Sobj <- ScaleData(Sobj, display.progress = T, do.par = T, num.cores = 3,
                        vars.to.regress = c("CC.Difference", "nUMI", "pc.mito"))
    }

    ## cell cycle effects strongly mitigated in PCA
    Sobj <- RunPCA(Sobj, pc.genes = c(s.genes, g2m.genes),
                   do.print = FALSE)

    PCAPlot(Sobj)

  }

  # Run principal component analysis
  Sobj <- RunPCA(object = Sobj,
               pc.genes = Sobj@var.genes, do.print = TRUE,
               pcs.print = 1:5, genes.print = 10)


  # Looking at SD of PCs
  PCElbowPlot(object = Sobj)

  # Perform shared nearest neighbor (SNN) clustering
  Sobj <- FindClusters(Sobj, reduction.type = "pca",
                     dims.use = ClustDims,  resolution = r, print.output = 0,
                     save.SNN = TRUE, force.recalc = TRUE)

   # Perform T-distributed stochastic neighborhood embedding (t-SNE)
   Sobj <- RunTSNE(Sobj, theta=0.5,
                dims.use = tDims, do.fast = TRUE, perplexity=perplexity)

  return(Sobj)

}
