#library(data.table)

doSeuratProc <- function(exp, samp, mito_rm, nc_rm, norm_type){

    print(samp)
    print(mito_rm)
    print(dim(exp))
    print(paste("genes in matrix =", dim(exp)[1]))

    if(nc_rm == "yes"){ #remove non coding genes
      z = which(rownames(exp) %in% pc_genes)
      exp = exp[z,]
    }

    print(paste("genes in matrix =", dim(exp)[1]))

    #1. create seurat object
    object_fl <- CreateSeuratObject(counts = exp,
      project = "FL",
      min.cells = 3, min.features = 200)

    object_fl$sample <- samp
    print(object_fl)

    #2. calculate mitochondrial percentage
    object_fl[["percent.mt"]] <- PercentageFeatureSet(object_fl, pattern = "^MT-")
    print(summary(object_fl[["percent.mt"]]))

    #3. Visualize QC metrics as a violin plot
    vlnplot = VlnPlot(object_fl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+
     ggtitle(object_fl@meta.data$sample[1])
    print(vlnplot)

    feature_plot = FeatureScatter(object_fl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    ggtitle(object_fl@meta.data$sample[1])
    print(feature_plot)

    #4. Filter cells that have unique feature counts over 5,000 or less than 200
    object_fl <- subset(object_fl, subset = nFeature_RNA > 200 & percent.mt < 5 & nFeature_RNA < 5000)
    print(object_fl)

    if(mito_rm == "yes"){
      counts <- GetAssayData(object_fl, assay = "RNA")
      mito=which(str_detect(rownames(counts), "MT-"))
      counts <- counts[-mito,]
      object_fl <- subset(object_fl, features = rownames(counts))
    }

    #5. Normalize the Seurat objects
    object_fl <- NormalizeData(object_fl)

    if(norm_type == "SC"){
      object_fl <- SCTransform(object_fl)
    }

    #6. FindVariableFeatures
    if(!(norm_type == "SC")){
      object_fl <- FindVariableFeatures(object_fl, selection.method = "vst", nfeatures = 2000) #changed nfeatures from 5000 to 2000
    }

    return(object_fl)

}
