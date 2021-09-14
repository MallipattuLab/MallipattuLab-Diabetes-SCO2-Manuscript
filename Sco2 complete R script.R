###comparing all four groups (WT,Sco2KO/KI,db/db,Sco2KO/KI;db/db)
library(dplyr)
library(Seurat)
library(Matrix)
library(cowplot)
library(rmarkdown)


#Create Seurat objects for each sample first
data_dir <- 'filtered_feature_bc_matrix_SC7_2/'
WT_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_SC7_2/')
WT_Data <- CreateSeuratObject(counts = WT_Data,
                              min.cells = 3,
                              min.features = 200,
                              project = "WT")

data_dir <- 'filtered_feature_bc_matrix_SC7_3/'
DBDB_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_SC7_3/')
DBDB_Data <- CreateSeuratObject(counts = DBDB_Data,
                                min.cells = 3,
                                min.features = 200,
                                project = "db/db")

data_dir <- 'filtered_feature_bc_matrix_SC7_1/'
KOKIDBDB_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_SC7_1/')
KOKIDBDB_Data <- CreateSeuratObject(counts = KOKIDBDB_Data,
                                    min.cells = 3,
                                    min.features = 200,
                                    project = "Sco2KO/KI;db/db")

data_dir <- 'filtered_feature_bc_matrix_SC7_4/'
KOKI_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_SC7_4/')
KOKI_Data <- CreateSeuratObject(counts = KOKI_Data,
                                min.cells = 3,
                                min.features = 200,
                                project = "Sco2KO/KI")




#Creates new merged seurat object
SC7.combined <- merge(KOKIDBDB_Data, y = c(WT_Data, DBDB_Data, KOKI_Data),   add.cell.ids = c("WT", "DBDB", "KOKIDBDB", "KOKI"), project = "Combined")
SC7 <- SC7.combined

# notice the cell names now have an added identifier
head(colnames(SC7))

#See original cell numbers
table(SC7$orig.ident)



# Get number of cells per cluster and per sample of origin
table(data.combined@meta.data[["seurat_clusters"]], data.combined@meta.data$orig.ident)

#Assess Mito%

SC7[["percent.mt"]] <- PercentageFeatureSet(SC7, pattern = "^mt-")

# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(SC7, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SC7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#mitochondrial content less than 5 and other QC features (will change based on the graphs above (Data))
SC7<- subset(SC7, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

#Identify Variable Features 
obj.list <- SplitObject(SC7, split.by = "orig.ident")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Perform Integration and Identify Anchors
data.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:30)

# Perform an integrated analysis
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)

#UMAP and Clustering
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(data.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE)
plot_grid(p2)

DimPlot(data.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Get number of cells per cluster and per sample of origin
table(data.combined@meta.data[["seurat_clusters"]], data.combined@meta.data$orig.ident)

#Get average expression of genes in each cluster
cluster.averages <- AverageExpression(data.combined, assays = "RNA", add.ident = "orig.ident")
write.csv(cluster.averages[["RNA"]], file = "CombinedClusterAveragespremRNA.csv")

#save file(please run this)
saveRDS(SC7, file = "SC7.rds")

#load saved file
data.combined<- readRDS("SC7.rds")

#set default assay to RNA (as opposed to integrated)
DefaultAssay(data.combined) <- "RNA"


#markers for cluster identification
new.cluster.idss <- c("0", "1","2", "4", "5", "6","7", "8", "9", "10", "11", "13",  "15", "16", "17", "18", "19", "20")
for(i in new.cluster.idss) {
  clus.markers <- FindMarkers(data.combined, ident.1 = i, ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
  head(clus.markers)
  str1 = i
  str2 = "deg.csv"
  print(i)
  result = paste(str1,str2, sep = ".")
  write.csv(clus.markers, file = result)
}
C0 <- read.csv("0.deg.csv")
C1 <- read.csv("1.deg.csv")
C2 <- read.csv("2.deg.csv")
C3 <- read.csv("3.deg.csv")
C4 <- read.csv("4.deg.csv")
C5 <- read.csv("5.deg.csv")
C6 <- read.csv("6.deg.csv")
C7 <- read.csv("7.deg.csv")
C8 <- read.csv("8.deg.csv")
C9 <- read.csv("9.deg.csv")
C10 <- read.csv("10.deg.csv")
C11 <- read.csv("11.deg.csv")
C12 <- read.csv("12.deg.csv")
C13 <- read.csv("13.deg.csv")
C14 <- read.csv("14.deg.csv")
C15 <- read.csv("15.deg.csv")
C16 <- read.csv("16.deg.csv")
C17 <- read.csv("17.deg.csv")
C18 <- read.csv("18.deg.csv")
C19 <- read.csv("19.deg.csv")
C20 <- read.csv("20.deg.csv")



C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)


list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2,  "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,  "C13" = C13,  "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20 )
#list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2,  "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,  "C13" = C13,  "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20)
write.xlsx(list_of_datasets, file = "Summary.DEGforclusters.Genes.xlsx")


#to identify differential genes between any two genotype (choose ident.1 and ident.2)

for(i in new.cluster.ids) {
  a = paste(i, "WT", sep = "_")
  print(a)
  b = paste(i, "KOKI", sep = "_")  
  c = paste(i, "dbdb", sep = "_")
  d = paste(i, "KOKIdbdb", sep = "_")
  
  MC.response <- FindMarkers(data.combined, ident.1 = a, ident.2 = b, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}


C0 <- read.csv("0comparison.csv")
C1 <- read.csv("1comparison.csv")
C2 <- read.csv("2comparison.csv")
C3 <- read.csv("3comparison.csv")
C4 <- read.csv("4comparison.csv")
C5 <- read.csv("5comparison.csv")
C6 <- read.csv("6comparison.csv")
C7 <- read.csv("7comparison.csv")
C8 <- read.csv("8comparison.csv")
C9 <- read.csv("9comparison.csv")
C10 <- read.csv("10comparison.csv")
C11 <- read.csv("11comparison.csv")
C12 <- read.csv("12comparison.csv")
C13 <- read.csv("13comparison.csv")
C14 <- read.csv("14comparison.csv")
C15 <- read.csv("15comparison.csv")
C16 <- read.csv("16comparison.csv")
C17 <- read.csv("17comparison.csv")
C18 <- read.csv("18comparison.csv")
C19 <- read.csv("19comparison.csv")
C20 <- read.csv("20comparison.csv")


C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)


list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3,"C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,"C12" = C12, "C13" = C13, "C14" = C14,"C15" = C15, "C16" = C16, "C17" = C17,"C18" = C18, "C19" = C19, "C20" = C20)
write.xlsx(list_of_datasets, file = "xxxx vs xxx.xlsx")


# label cluster
data.combined <- RenameIdents(data.combined, `0` = "PT(S3)-1", `1` = "PT(S1-S2)-1", `2` = "PT(S1-S2)2", 
                                `3` = "PT(S3)-2", `4` = "LH(AL)", `5` = "DCT", `6` = "Mes", `7` = "PT(S3)/LH(DL)",`8` = "CNT", `9` = "Endo", 
                                `10` = "DCT/CNT-1", `11` = "IC-B", `12` = "DCT/CNT-2", '13' = "CD-PC", '14'= "LH(AL)", '15' = "Pod", 
                                '16' = "IC-A", '17' = "LH(DL)", '18' = "Novel", '19' = "MØ-1", '20' = "MØ-2")

#To remove the cluster only present in Sco2KO/KIdbdb
levels(x = data.combined.withoutnovel) <- c(       "PT(S1-S2)-1" ,  "PT(S1-S2)-2" ,   "PT(S3)"  , "PT(S3)/LH(DL)",   "LH(DL)" , "LH(AL)"  ,  
                                                     "DCT"    ,    "DCT/CNT"   ,    "CNT"  ,      "CD-PC"  ,    "IC-A"     ,    "IC-B" ,        
                                                     "MØ-1"        ,  "MØ-2" ,   "Mes"  ,    "Pod" , "Endo"   )

#For dot plot combined
markers.to.plot <- c(  "Slc7a7","Slc5a12", "Slc13a3", "Cyp7b1", "Cp" , "Slc12a1", "Slc12a3","Slc8a1","Phactr1","Mgat4c","Clnk","Insrr", "Cd74","Diaph3","Cfh","Nphs1", "Flt1")
DotPlot(data.combined.withoutnovel, features = rev(markers.to.plot), cols = c("Red", "white"), 
        dot.scale = 8) + RotatedAxis()


#####Comparing Sco2KO/KI;db/db vs db/db
data_dir <- 'filtered_feature_bc_matrix_SC7_1/'
KOKIdbdb_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_SC7_1/')
KOKIdbdb_Data <- CreateSeuratObject(counts = KOKIdbdb_Data,
                                    min.cells = 3,
                                    min.features = 200,
                                    project = "KOKIdbdb")

data_dir <- 'filtered_feature_bc_matrix_SC7_3/'
DBDB_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_SC7_3/')
DBDB_Data <- CreateSeuratObject(counts = DBDB_Data,
                                min.cells = 3,
                                min.features = 200,
                                project = "dbdb")


#Creates new merged seurat object for Sco2KO/KI;dbdb and dbdb
SC713.combined <- merge(KOKIdbdb_Data, y = DBDB_Data,   add.cell.ids = c("KOKIdbdb", "DBDB"), project = "Combined")
SC713 <- SC713.combined


#Assess mt content
SC713[["percent.mt"]] <- PercentageFeatureSet(SC713, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(SC713, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
plot1 <- FeatureScatter(SC713, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SC713, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Filter cells, only keep the cells with mitochondrial content less than 5, and more than 200 and less than 5500 features.
SC713<- subset(SC713, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 5)

#Identify Variable Features
obj.list <- SplitObject(SC713, split.by = "orig.ident")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Perform Integration and identify anchors
data.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:30)

#Perform an intergrated analysis
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)

#UMAP and Clustering
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5)

# Visualization of UMAP
p1 <- DimPlot(data.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

#show Sco2KO/KI;db/db and db/db side by side
DimPlot(data.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)


# Get number of cells per cluster and per sample of origin
table(data.combined@meta.data[["seurat_clusters"]], data.combined@meta.data$orig.ident)

#To output average expression of genes in each cluster by group
cluster.averages <- AverageExpression(data.combined, assays = "RNA", add.ident = "orig.ident")
write.csv(cluster.averages[["RNA"]], file = "CombinedClusterAveragespremRNA.csv")


#For differential expressed genes between a single cluster and all the other cells

new.cluster.idss <- c("0", "1","2", "4", "5", "6","7", "8", "9", "10", "11", "13",  "15", "16", "17")
for(i in new.cluster.idss) {
  clus.markers <- FindMarkers(data.combined, ident.1 = i, ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
  head(clus.markers)
  str1 = i
  str2 = "deg.csv"
  print(i)
  result = paste(str1,str2, sep = ".")
  write.csv(clus.markers, file = result)
}
C0 <- read.csv("0.deg.csv")
C1 <- read.csv("1.deg.csv")
C2 <- read.csv("2.deg.csv")
C3 <- read.csv("3.deg.csv")
C4 <- read.csv("4.deg.csv")
C5 <- read.csv("5.deg.csv")
C6 <- read.csv("6.deg.csv")
C7 <- read.csv("7.deg.csv")
C8 <- read.csv("8.deg.csv")
C9 <- read.csv("9.deg.csv")
C10 <- read.csv("10.deg.csv")
C11 <- read.csv("11.deg.csv")
C12 <- read.csv("12.deg.csv")
C13 <- read.csv("13.deg.csv")
C14 <- read.csv("14.deg.csv")
C15 <- read.csv("15.deg.csv")
C16 <- read.csv("16.deg.csv")
C17 <- read.csv("17.deg.csv")



C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)


list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2,  "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,  "C13" = C13,  "C15" = C15, "C16" = C16, "C17" = C17 )
write.xlsx(list_of_datasets, file = "Summary.DEGforclusters.Genes.xlsx")


#to identify differential genes between db/db vs Sco2KO/KI;db/db

for(i in new.cluster.ids) {
  a = paste(i, "dbdb", sep = "_")
  print(a)
  b = paste(i, "KOKIdbdb", sep = "_")

  MC.response <- FindMarkers(data.combined, ident.1 = a, ident.2 = b, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}




C0 <- read.csv("0comparison.csv")
C1 <- read.csv("1comparison.csv")
C2 <- read.csv("2comparison.csv")
C3 <- read.csv("3comparison.csv")
C4 <- read.csv("4comparison.csv")
C5 <- read.csv("5comparison.csv")
C6 <- read.csv("6comparison.csv")
C7 <- read.csv("7comparison.csv")
C8 <- read.csv("8comparison.csv")
C9 <- read.csv("9comparison.csv")
C10 <- read.csv("10comparison.csv")
C11 <- read.csv("11comparison.csv")
C12 <- read.csv("12comparison.csv")
C13 <- read.csv("13comparison.csv")
C14 <- read.csv("14comparison.csv")
C15 <- read.csv("15comparison.csv")
C16 <- read.csv("16comparison.csv")
C17 <- read.csv("17comparison.csv")
C18 <- read.csv("18comparison.csv")


C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)


list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3,"C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,"C12" = C12, "C13" = C13, "C14" = C14,"C15" = C15, "C16" = C16, "C17" = C17)


# label cluster

data.combined <- RenameIdents(data.combined, `0` = "PT-S1-S2-1", `1` = "PT-S3", `2` = "PT-S1-S2-2", 
                                `3` = "CNT", `4` = "DCT", `5` = "Endo", `6` = "LH(AL)", `7` = "LH(AL)",`8` = "MC", `9` = "IC-B", 
                                `10` = "PT-S3/LH(DL)", `11` = "IC-A", `12` = "CD/PC", '13' = "Pod", '14'= "Endo",
                                '15' = "LH(DL)", '16' = "Mix", '17'= "MØ")


data.combined$celltype.orig.ident<- paste(Idents(data.combined), data.combined$orig.ident, sep = "_")

data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "orig.ident"
head(data.combined@meta.data)

#For Dot plot (combined)
Idents(data.combined)<-"celltype.orig.ident"
data.combined.endo<-subset(data.combined, idents = c("Endo_dbdb", "Endo_KOKIdbdb") )

#For violin plot for up-regulated genes in the endo cells
#Violinplot of genes in individual clusters in specific genotype (make sure you split the object)
plots <- VlnPlot(data.combined, idents = c("Endo_dbdb", "Endo_KOKIdbdb"),features = c( "Kdr", "Gpx3","Zbtb16","Hspg2", "Plvap","Igfbp5" ), split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

#For Violin plot for down-regulated genes in the endo cells
VlnPlot(data.combined, idents = c("Endo_dbdb", "Endo_KOKIdbdb"),features = c(  "Malat1","Mecom", "Eln","Thsd4","St3gal4", "Ptprj"), split.by = "orig.ident", 
        cols = c("white", "white") ,  pt.size = 0) + stat_summary(fun.y = median.stat, geom='point', size = 10, colour = "blue") 
wrap_plots(plots = plots, ncol = 3)

#Idents(data.combined) <- factor(Idents(data.combined), levels = c("PT (S1-1)", "PT (S1-2)", "CNT", "DCT", "EC/IMCT", "LOH-Thick", "PT/Intermedullary CT",  "IC/IMCT",  "Intercalated B", "PT",  "PT (S2-S3)",  "CT",  "Intercalated A",  "IntercalatedA",  "?",  "PT-S1",  "only in SC7_3"))
markers.to.plot <- c("Flt4", "Kdr","Flt1","Cat", "Txnrd1", "Gclc", "Gsr","Gstt2","Gpx1","Gpx3")
DotPlot(data.combined, features = rev(markers.to.plot), cols = c("yellow", "red"), 
        dot.scale = 8) + RotatedAxis()

###Subclustering of EC cluster 

data.combined.EC <- subset(data.combined, idents = c("Endo_dbdb", "Endo_KOKIdbdb"))
DefaultAssay(object = data.combined.EC) <- "integrated"

data.combined.EC <- FindNeighbors(data.combined.EC, dims = 1:10)
data.combined.EC <- FindClusters(data.combined.EC, resolution = 0.5)

# Generate a new column called sub_cluster in the metadata
data.combined$sub_cluster <- as.character(Idents(data.combined))

# Change the information of cells containing sub-cluster information
data.combined$sub_cluster[Cells(data.combined.EC)] <- paste("EC",Idents(data.combined.EC))
DimPlot(data.combined, group.by = "sub_cluster")

DimPlot(data.combined.EC, reduction = "umap", group.by = "sub_cluster")
DimPlot(data.combined.EC, reduction = "umap", label = FALSE)
DimPlot(data.combined.EC, reduction = "umap", label = TRUE)
DimPlot(data.combined.EC, reduction = "umap", label = TRUE, split.by = "orig.ident")
DimPlot(data.combined.EC, reduction = "umap", label = FALSE, split.by = "orig.ident")

#To compare one cluster with other cluster (in combined object)
C0 <- FindMarkers(data.combined.EC, ident.1 = "0", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(data.combined.EC, ident.1 = "1", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(data.combined.EC, ident.1 = "2", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(data.combined.EC, ident.1 = "3", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(data.combined.EC, ident.1 = "4", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(data.combined.EC, ident.1 = "5", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5)
write.xlsx(list_of_datasets, file = "EC_subcluster_markers.xlsx")

#To compare one genotype vs another (after you split the clusters  into their respective genotype)


data.combined.EC$celltype.orig.ident <- paste(Idents(data.combined.EC), data.combined.EC$orig.ident, sep = "_")
data.combined.EC$celltype <- Idents(data.combined.EC)
Idents(data.combined.EC) <- "celltype"

#to look at the metadata headings
head(data.combined.EC@meta.data)


C0 <- FindMarkers(data.combined.EC, ident.1 = "0_KOKIdbdb", ident.2 = "0_dbdb", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(data.combined.EC, ident.1 = "1_KOKIdbdb", ident.2 = "1_dbdb", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(data.combined.EC, ident.1 = "2_KOKIdbdb", ident.2 = "2_dbdb", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(data.combined.EC, ident.1 = "3_KOKIdbdb", ident.2 = "3_dbdb", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(data.combined.EC, ident.1 = "4_KOKIdbdb", ident.2 = "4_dbdb", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(data.combined.EC, ident.1 = "5_KOKIdbdb", ident.2 = "5_dbdb", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)



saveRDS(data.combined, file = "EC.rds")

DefaultAssay(object = data.combined.EC) <- "RNA"

# Change the information of cells containing sub-cluster information
data.combined.EC$subcluster <- paste(Idents(data.combined.EC), data.combined.EC$orig.ident, sep = "_")
Idents(data.combined.EC) <- "subcluster"
levels(x = ) <- c("0_dbdb", "0_KOKIdbdb", "1_dbdb", "1_KOKIdbdb", "2_dbdb", "2_KOKIdbdb", "3_dbdb", "3_KOKIdbdb", "4_dbdb", "4_KOKIdbdb", "5_dbdb", "5_KOKIdbdb")

markers.to.plot <- c("Plat","Tgfbr2","Ehd3", "Cdh5","Plvap", "Tek", "Kdr" )
DotPlot(data.combined.EC, features = rev(markers.to.plot), cols = c("yellow", "red"), 
        dot.scale = 8) + RotatedAxis()

plots <-VlnPlot(object = data.combined, features = c("Kdr", "Gpx3", "Zbtb16", "Hspg2","Plvap", "Igfbp5"), cols = c("White", "Grey"), ncol = 1, split.by = "celltype.orig.ident", group.by = "seurat_clusters", split.plot=TRUE, pt.size = 0, combine = FALSE, log=T)
CombinePlots(plots, ncol = 1)

