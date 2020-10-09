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
