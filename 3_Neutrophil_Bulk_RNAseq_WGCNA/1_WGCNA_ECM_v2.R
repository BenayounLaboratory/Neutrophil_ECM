setwd("/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/3_Neutrophil_Bulk_RNAseq_WGCNA/")
options(stringsAsFactors = FALSE)

library("WGCNA")
library("pheatmap")

#===============================================================================
#  WGCNA for Neutrophil ECM Collaboration
#===============================================================================

################################################################################################
#  1. Read the gene counts table and plot the sample tree

# Read data and update column names
expressionList           <- read.csv('2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt', header = TRUE, sep = "\t", row.names = NULL);
colnames(expressionList) <- c("gene_id","Female_4m_F1a",	"Female_4m_F1b",	"Female_4m_F2a",	"Female_4m_F2b",	"Female_20m_F3a",	"Female_20m_F3b",	"Female_20m_F4a",	"Female_20m_F4b",	"Male_4m_M1a"	,"Male_4m_M1b"	,"Male_4m_M2a"	,"Male_4m_M2b"	,"Male_20m_M3a"	,"Male_20m_M3b",	"Male_20m_M4a",	"Male_20m_M4b"	)

## Prepare and clean data
datExpr0              <- as.data.frame(t(expressionList[,-c(1)]));
colnames(datExpr0)    <- expressionList$gene_id;
rownames(datExpr0)    <- names(expressionList)[-c(1)];

# Check that all genes and samples have sufficiently low numbers of missing values. #13858 genes left
gsg = goodSamplesGenes(datExpr0, verbose = 3);
# Flagging genes and samples with too many missing values...
# ..step 1
# ..Excluding 926 genes from the calculation due to too many missing samples or zero variance.
# ..step 2

# Removing genes that did not pass filter
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]

## Clustering
sampleTree = hclust(dist(datExpr0), method = "average");

pdf(file = "1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

# Clustering is decent, all samples are kept
datExpr  = datExpr0
nGenes   = ncol(datExpr)
nSamples = nrow(datExpr)

collectGarbage();

save(datExpr, file = "WGCNA_data_input.RData")

## Modules construction
# enableWGCNAThreads()

################################################################################################
#  2. Choose soft threshold parameter

# set a saerch space
powers = c(c(1:10), seq(from = 12, to=20, by=2)) #to find the best power

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

## Set Power
power     = sft$powerEstimate #estimates the power if you're not sure from above #10
softPower = 10; # from above
adjacency = adjacency(datExpr, power = softPower);

################################################################################################
#  3. Turn data expression into topological overlap matrix

# run similarity analysis
TOM     = TOMsimilarity(adjacency); # may take a long time but it's ok
dissTOM = 1-TOM
save(TOM, file = "TOM.RData")

# load("TOM.RData")

# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");

pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()


################################################################################################
#  4. Construct modules

# Setting the minimum module size relatively high, but modules can be made smaller:
minModuleSize = 100;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

pdf(file = "4-module_tree_100.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()


# Call an automatic merging function
MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()

write.table(merge$oldMEs,file="oldMEs.txt");
write.table(merge$newMEs,file="newMEs.txt");

################################################################################################
#  5. Export of networks to external software

# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = names(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("original_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("original_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}

# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = names(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("merged_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("merged_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}


################################################################################################
#  6. Plot the heatmap of module eigen-genes and samples

# Heatmap of old module eigen-genes and samples
pdf(file="6-oldMEs.pdf")
row.names(merge$oldMEs)=row.names(datExpr0)
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=F,show_rownames=T,show_colnames=T,fontsize=6)
dev.off()

# Heatmap of new module eigen-genes and samples
pdf(file="7-newMEs.pdf")
row.names(merge$newMEs)=row.names(datExpr0)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=F,show_rownames=T,show_colnames=T,fontsize=6)
dev.off()

clust <- pheatmap(merge$newMEs,cluster_col=T,cluster_row=F,show_rownames=T,show_colnames=T,fontsize=6)
module.order <- clust$tree_col$labels[clust$tree_col$order]

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, dynamicColors, file = "Data-networkConstruction.RData")


################################################################################################
# 7.  compare modules with each other 

# Do not plot full TOM (not useful, takes too long)
# # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
# plotTOM = dissTOM^7;
# # Set diagonal to NA for a nicer plot
# diag(plotTOM) = NA;
# # Call the plot function
# # sizeGrWindow(9,9)
# # TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#
nSelect = 500

# For reproducibility, we set the random seed
set.seed(10);
select    = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];

# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree   = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];

# Open a graphical window
sizeGrWindow(9,9)

# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^10;
diag(plotDiss) = NA;

pdf(file="8-TOMPlot_select.pdf")
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# # Plot the relationships among the eigengenes and the trait
# sizeGrWindow(5,7.5);
# par(cex = 0.9)
# plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)

################################################################################################
#  8.  Module-trait analysis 

# Read meta data
traitData = read.csv("ClinicalTraits.txt", header = T, sep = "\t");
dim(traitData)
names(traitData)


# Form a multi-set structure that will hold the traits.
Traits = vector(mode="list", length = 1);
setSamples = rownames(datExpr);
traitRows  = match(setSamples, traitData$Mice);
Traits[[1]] = list(data = traitData[traitRows, -1]);
rownames(Traits[[1]]$data) = traitData[traitRows, 1];

collectGarbage();

# Define numbers of genes and samples
nGenes   = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs  = orderMEs(MEs0, order = module.order)
colnames(MEs) <- module.order
moduleTraitCor = cor(MEs, traitData[,-1], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));



# Display the correlation values within a heatmap plot (with pheatmap order)
pdf(file="9-Module_Trait_analysis.pdf")
labeledHeatmap(Matrix        = moduleTraitCor,
               xLabels       = names(traitData)[-1],
               yLabels       = names(MEs),
               ySymbols      = names(MEs),
               colorLabels   = FALSE,
               colors        = blueWhiteRed(50),
               textMatrix    = textMatrix,
               setStdMargins = FALSE,
               cex.text      = 0.5,
               zlim          = c(-1,1),
               main          = paste("Module-trait relationships"))
dev.off()

################################################################################################


#######################
sink(file = paste(Sys.Date(),"_WGCNA_session_Info.txt", sep =""))
sessionInfo()
sink()


