

if(FALSE){
  devtools::install_github(repo = "hhoeflin/hdf5r")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  devtools::install_github(repo = "satijalab/seurat", ref = "loom")
  
  
  
  # devtools::install_github(repo = "hhoeflin/hdf5r")
  # devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
}

library(loomR)
library(Seurat)
library(hdf5r)

cellxgene <- H5File$new("/media/mahogny/TOSHIBA EXT/collab/mattias/spleen.cellxgene.h5ad", mode = "r")
#cellxgene$close_all()

cellxgene[["X"]][1:5,1:5]

cellxgene$ls()
cellxgene[["uns"]][["patient_categories"]]
cellxgene[["uns"]][["Celltypes_categories"]]

cellxgene[["uns"]][["Celltypes_categories"]][]



h5attr(cellxgene,"uns")
  

mca <- connect(filename = "/media/mahogny/TOSHIBA EXT/collab/mattias/spleen.cellxgene.h5ad", mode = "r+")

mca$attr_get_number()



mca[["col_attrs"]]

unique(mca[["col_attrs/cell_suspension.provenance.document_id"]][])
unique(mca[["col_attrs/derived_organ_label"]][])

mca[["row_attrs"]]

mca[["layers"]]

NormalizeData(object = mca, chunk.size = 1000, scale.factor = 10000, display.progress = FALSE)

FindVariableGenes(object = mca)
hv.genes <- head(x = GetVariableGenes(object = mca)$index, n = 1000)

# Pull the indices of the mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = mca[["row_attrs/gene_names"]][], value = FALSE)
# Calculate percent.mito and store directly to the loom object Similar to
# creating a vector with the percentage of expression from mitochondrial
# genes, then using AddMetaData to put in to an object
mca$apply(name = "col_attrs/percent_mito", FUN = function(mat) {
  return(rowSums(x = mat[, mito.genes])/rowSums(x = mat))
}, MARGIN = 2, dataset.use = "matrix")
ScaleData(object = mca, genes.use = hv.genes, chunk.size = 20, display.progress = FALSE,
          vars.to.regress = "percent_mito")



RunPCA(object = mca, pc.genes = hv.genes, online.pca = FALSE, pcs.compute = 100,
       do.print = TRUE, pcs.print = 1:5, genes.print = 5, display.progress = FALSE)


PCElbowPlot(object = mca, num.pc = 100)

PCHeatmap(mca, pc.use = c(1:3, 70:75), cells.use = 500, do.balanced = TRUE)

FindClusters(object = mca, reduction.type = "pca", dims.use = 1:75, resolution = 3,
             save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)



RunTSNE(object = mca, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE",
        max_iter = 2000, nthreads = 4, overwrite = TRUE)
p1 <- DimPlot(object = mca, reduction.use = "tsne", no.legend = TRUE, do.return = TRUE,
              pt.size = 0.1, ident.use = "col_attrs/res.3", vector.friendly = TRUE) +
  ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = mca, reduction.use = "tsne", no.legend = TRUE, do.return = TRUE,
              pt.size = 0.1, ident.use = "col_attrs/tissue", vector.friendly = TRUE) +
  ggtitle("Tissue") + theme(plot.title = element_text(hjust = 0.5))
plot_grid(p1, p2)



FeaturePlot(mca, c("S100a9", "Sftpc"), reduction.use = "tsne", dark.theme = TRUE,
            pt.size = 0.1, vector.friendly = TRUE)


