library(data.table)

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


#cell types
#gene marker versus cell type data
cells = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/PanglaoDB_markers_27_Mar_2020.tsv")
cells = as.data.table(filter(cells, species %in% c("Mm Hs", "Hs"), organ=="Immune system"))
colnames(cells)[c(2,3,5, 6, 9)] = c("gene", "cell", "ubiquitousness_index", "product", "germlayer")
cells = as.data.table(cells %>% select(gene, cell, ubiquitousness_index, product, sensitivity_human, specificity_human) %>%
	filter(sensitivity_human > 0.1, specificity_human > 0.1))

imsig_cells = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/imsig_gene_cells.csv")

cell_markers = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Human_cell_markers_HRBMUS_EDU.txt")
cell_markers = as.data.table(filter(cell_markers, tissueType=="Blood"))
cell_markers = cSplit(cell_markers, "cellMarker", ",", direction = "long")
cell_markers = cell_markers %>% select(tissueType, cellType, cellName, cellMarker)
colnames(cell_markers)[4] = "gene"

cell_markers_single_cell = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Single_cell_markers_HRBMU_EDU.txt")
cell_markers_single_cell = as.data.table(filter(cell_markers_single_cell, tissueType=="Blood"))
cell_markers_single_cell = cSplit(cell_markers_single_cell, "proteinName", ",", direction = "long")
cell_markers_single_cell = cell_markers_single_cell %>% select(tissueType, cellType, cellName, proteinName)
colnames(cell_markers_single_cell)[4] = "gene"

#data from nature methods paper for CellAssign
fl_combined_wov = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/41592_2019_529_MOESM4_ESM_HGSC_FL_combined_edited.csv")
fl_cells = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/41592_2019_529_MOESM4_ESM_FL_light_chain.csv")
fl_lightchain = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/41592_2019_529_MOESM4_ESM_FL_celltype.csv")

fl_combined_wov = melt(fl_combined_wov)
fl_combined_wov$data = "fl_combined_wov"

fl_cells = melt(fl_cells)
fl_cells$data = "fl_cells"

fl_lightchain = melt(fl_lightchain)
fl_lightchain$data = "fl_lightchain"

all_cellassign = rbind(fl_lightchain, fl_combined_wov, fl_cells)
colnames(all_cellassign)[1:2] = c("gene", "cell")
all_cellassign = as.data.table(filter(all_cellassign, value==1))

#data from christian steidl paper (manually put into spreadsheet)
steidl = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/FL_cells_types_aoki_cancerdiscovery.csv")
