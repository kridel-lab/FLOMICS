doSeuratProc <- function(exp, samp){

    print(samp)

    #1. create seurat object
    object_fl <- CreateSeuratObject(counts = exp,
      project = "FL",
      min.cells = 3, min.features = 200)

    object_fl$sample <- samp

    #2. calculate mitochondrial percentage
    object_fl[["percent.mt"]] <- PercentageFeatureSet(object_fl, pattern = "^MT-")

    #3. Visualize QC metrics as a violin plot
    vlnplot = VlnPlot(object_fl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+
     ggtitle(object_fl@meta.data$sample[1])
    print(vlnplot)

    feature_plot = FeatureScatter(object_fl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    ggtitle(object_fl@meta.data$sample[1])
    print(feature_plot)

    #4. Filter cells that have unique feature counts over 2,500 or less than 200
    object_fl <- subset(object_fl, subset = nFeature_RNA > 200 & percent.mt < 5 & nFeature_RNA < 5000)
    #change nFeature_RNA < 5000 instead of 2500

    #5. Normalize the Seurat objects
    object_fl <- NormalizeData(object_fl)

    #6. FindVariableFeatures
    object_fl <- FindVariableFeatures(object_fl, selection.method = "vst", nfeatures = 5000) #changed nfeatures from 2000 to 5000

    return(object_fl)

}
