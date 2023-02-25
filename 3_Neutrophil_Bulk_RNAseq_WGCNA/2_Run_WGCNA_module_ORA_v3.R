setwd("/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/3_Neutrophil_Bulk_RNAseq_WGCNA/")
options(stringsAsFactors = FALSE)

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

library(ReactomePA)
# ReactomePA v1.38.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
#   If you use ReactomePA in published research, please cite:
#   Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479

# analyze Functional enrichments for thee top female biased cluster
# salmon is Sex-biased (female high)
# magenta is Sex-biased (male high)


################################################################################################
# Read module gene lists
my.universe   <- read.table('2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt', sep = "\t", header = T)
my.salmon     <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-salmon.txt', sep = "\t", header = T)
my.magenta    <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-magenta.txt', sep = "\t", header = T)


# genes lists for enrichment analysis
genes.universe   <- rownames(my.universe)
genes.salmon     <- my.salmon$nodeName   
genes.magenta    <- my.magenta$nodeName   

# Convert to ENTREZ ID
entrezID.genes.salmon       <- bitr(genes.salmon  , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.genes.magenta      <- bitr(genes.magenta  , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.genes.universe     <- bitr(genes.universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


########################################################################
#### A. REACTOME over-representation test
rct.salmon <- enrichPathway(gene          = entrezID.genes.salmon$ENTREZID,
                            universe      = entrezID.genes.universe$ENTREZID,
                            readable      = TRUE,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.05,
                            organism      = 'mouse',
                            minGSSize     = 10  ,
                            maxGSSize     = 5000)

# write results to file
write.table(rct.salmon@result[rct.salmon@result$p.adjust < 0.05,]    ,  file = paste(Sys.Date(),"Neutrophil_WGCNA_REACTOME_Salmon_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")

# make some plots
pdf(paste(Sys.Date(),"Neutrophil_WGCNA_REACTOME_Salmon_Analysis_dotplot_FDR5.pdf", sep = "_"))
dotplot(rct.salmon,  x = "Count", title = "Salmon REACTOME gene sets")
dev.off()

pdf(paste(Sys.Date(),"Neutrophil_WGCNA_REACTOME_Salmon_Analysis_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 10, height = 10)
cnetplot(rct.salmon,  categorySize="pvalue")
dev.off()


rct.magenta <- enrichPathway(gene          = entrezID.genes.magenta$ENTREZID,
                             universe      = entrezID.genes.universe$ENTREZID,
                             readable      = TRUE,
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             organism      = 'mouse',
                             minGSSize     = 10  ,
                             maxGSSize     = 5000)

# write results to file
write.table(rct.magenta@result[rct.magenta@result$p.adjust < 0.05,]    ,  file = paste(Sys.Date(),"Neutrophil_WGCNA_REACTOME_magenta_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")

# make some plots
pdf(paste(Sys.Date(),"Neutrophil_WGCNA_REACTOME_magenta_Analysis_dotplot_FDR5.pdf", sep = "_"))
dotplot(rct.magenta,  x = "Count", title = "magenta REACTOME gene sets")
dev.off()

pdf(paste(Sys.Date(),"Neutrophil_WGCNA_REACTOME_magenta_Analysis_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 10, height = 10)
cnetplot(rct.magenta,  categorySize="pvalue")
dev.off()


########################################################################
#### B. GO over-representation test
GO.salmon <- enrichGO(gene          = entrezID.genes.salmon$ENTREZID,
                      universe      = entrezID.genes.universe$ENTREZID,
                      ont           = "ALL",
                      OrgDb         = "org.Mm.eg.db",
                      readable      = TRUE,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      minGSSize     = 10  ,
                      maxGSSize     = 5000)

write.table(GO.salmon@result[GO.salmon@result$p.adjust < 0.05,]    ,  file = paste(Sys.Date(),"Neutrophil_WGCNA_GO_Salmon_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


pdf(paste(Sys.Date(),"Neutrophil_WGCNA_GO_Salmon_Analysis_dotplot_FDR5.pdf", sep = "_"), width = 4.5, height = 3.5)
dotplot(GO.salmon,  x = "Count", title = "Salmon GO gene sets")+scale_x_continuous(limits = c(0,8))+scale_size(breaks = c(0.01,0.02), name="Gene Ratio")
dev.off()


pdf(paste(Sys.Date(),"Neutrophil_WGCNA_GO_Salmon_Analysis_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 5, height = 5)
cnetplot(GO.salmon,  categorySize="pvalue", layout = "circle")
dev.off()

##
GO.magenta <- enrichGO(gene          = entrezID.genes.magenta$ENTREZID,
                       universe      = entrezID.genes.universe$ENTREZID,
                       ont           = "ALL",
                       OrgDb         = "org.Mm.eg.db",
                       readable      = TRUE,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

write.table(GO.magenta@result[GO.magenta@result$p.adjust < 0.05,]    ,  file = paste(Sys.Date(),"Neutrophil_WGCNA_GO_magenta_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


pdf(paste(Sys.Date(),"Neutrophil_WGCNA_GO_magenta_Analysis_dotplot_FDR5.pdf", sep = "_"), width = 4.5, height = 3.5)
dotplot(GO.magenta,  x = "Count", title = "magenta GO gene sets")+scale_x_continuous(limits = c(0,8))+scale_size(breaks = c(0.01,0.02), name="Gene Ratio")
dev.off()


pdf(paste(Sys.Date(),"Neutrophil_WGCNA_GO_magenta_Analysis_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 5, height = 5)
cnetplot(GO.magenta,  categorySize="pvalue", layout = "circle")
dev.off()
################################################################################################


#######################
sink(paste0(Sys.Date(),"_Enrichment_analysis_sessionInfo.txt"))
sessionInfo()
sink()
