setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/2_Neutrophil_scRNAseq_ECM_Enrichment/')
options(stringsAsFactors = F)

library(Seurat)
library(beeswarm)
library(readxl)
library(ggplot2)
library(UCell)

# 2022-10-14
# Analyze ECM gene expression in scRNAseq data

# 2022-12-13
# Add analysis by subset

####################################################################################################
# Load Kim et al, 2022 dataset
load('2022-04-12_10x_BM_Ntph_USC_Xie_SingleCellNet_predictions_Seurat_object.RData')
ntph.singlet.cl
# An object of class Seurat 
# 22402 features across 6025 samples within 3 assays 
# Active assay: SCT (11199 features, 3000 variable features)
# 2 other assays present: RNA, HTO
# 2 dimensional reductions calculated: pca, umap

####################################################################################################
# Load and prep ECM genesets for UCell

# read lists
my.genesets <- read_xlsx('SupplementalTable1.xlsx',sheet = 1)
colnames(my.genesets) <- c( "Division","Category", "Gene_Symbol")

# sub chaaraacters for R
my.genesets$Division  <- gsub(" |-","_",my.genesets$Division)
my.genesets$Category  <- gsub(" |-","_",my.genesets$Category)

# get categories and remove meaningless ones
divisions  <- setdiff(unique(my.genesets$Division),"Retired")
categories <- setdiff(unique(my.genesets$Category),"n/a")

# initialize and populate
gs.list.divisions         <- vector(mode = "list", length = length(divisions))
names(gs.list.divisions)  <- divisions
gs.list.categories        <- vector(mode = "list", length = length(categories))
names(gs.list.categories) <- categories


for (i in 1:length(divisions)) {
  gs.list.divisions[[i]] <- my.genesets$Gene_Symbol[my.genesets$Division == divisions[i]]
}

for (i in 1:length(categories)) {
  gs.list.categories[[i]] <- my.genesets$Gene_Symbol[my.genesets$Category == categories[i]]
}

gs.ecm.list <- c(gs.list.divisions,gs.list.categories)

# Score ECM gene lists with UCell
DefaultAssay(ntph.singlet.cl) <- "SCT"
ntph.singlet.cl <- AddModuleScore_UCell(ntph.singlet.cl, features=gs.ecm.list, name=NULL)

pdf(paste0(Sys.Date(),"_DotPLot_ECM_GeneLists_UCell_Scores.pdf"),width = 7, height = 5)
DotPlot(ntph.singlet.cl, 
        features = rev(c("Collagens"              ,
                         "Core_matrisome"         ,
                         "Matrisome_associated"   ,
                         "Secreted_Factors"       ,
                         "ECM_affiliated_Proteins",
                         "ECM_Glycoproteins"      ,      
                         "ECM_Regulators"         ,
                         "Proteoglycans" )), 
        group.by = "Sample_ID", cols = c("grey", "red"),col.min = -0.5,col.max = 0.5) + coord_flip()
dev.off()

wilcox.test(Collagens               ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 0.3222
wilcox.test(Core_matrisome          ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 9.93e-15
wilcox.test(Matrisome_associated    ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 1.312e-06
wilcox.test(Secreted_Factors        ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 0.6313
wilcox.test(ECM_affiliated_Proteins ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 0.004643
wilcox.test(ECM_Glycoproteins       ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 2.717e-14
wilcox.test(ECM_Regulators          ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 3.799e-06
wilcox.test(Proteoglycans           ~ orig.ident, data = ntph.singlet.cl@meta.data) # p-value = 2.257e-10


# Analysis as a function of neutrophil subset
#### by subset
pdf(paste0(Sys.Date(),"_DotPLot_ECM_GeneLists_UCell_Scores_by_XieSubset.pdf"),width = 7, height = 5)
DotPlot(ntph.singlet.cl, 
        features = rev(c("Collagens"              ,
                         "Core_matrisome"         ,
                         "Matrisome_associated"   ,
                         "Secreted_Factors"       ,
                         "ECM_affiliated_Proteins",
                         "ECM_Glycoproteins"      ,      
                         "ECM_Regulators"         ,
                         "Proteoglycans" )), 
        group.by = "SingleCellNet_Xie", cols = c("grey", "red"),col.min = -0.5,col.max = 0.5) + coord_flip()
dev.off()


#############################################################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()


