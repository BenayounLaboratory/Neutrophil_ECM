setwd("/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/3_Neutrophil_Bulk_RNAseq_WGCNA/")
options(stringsAsFactors = FALSE)
library('grDevices')
library('bitops')
library('pheatmap')



################################################################################################
#### load expression data #######
tissue.cts   <- read.table('2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt')
my.salmon     <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-salmon.txt', sep = "\t", header = T)

pdf(paste0(Sys.Date(),"_Slamon_module_expression_heatmap.pdf"),width = 5, height = 5, onefile = F)
pheatmap(tissue.cts[my.salmon$nodeName,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F,
         scale="row",
         main = paste0("Salmon module ", length(my.salmon$nodeName)) )
dev.off()

################################################################################################


#######################
sink(file = paste(Sys.Date(),"_WGCNA_Heatmap_session_Info.txt", sep =""))
sessionInfo()
sink()

