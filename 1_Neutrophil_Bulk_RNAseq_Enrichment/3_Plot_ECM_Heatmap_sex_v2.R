setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/1_Neutrophil_Bulk_RNAseq_Enrichment')
options(stringsAsFactors = F)

library('grDevices')
library('bitops')
library('pheatmap')

# 2021-09-06
# get ECM associated genes heatmap

################################################################################################
#### load expression data #######
tissue.cts <- read.table('2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt')

# load gene lists
my.1.ECM     <- read.table('./GeneSets/CoreMatrisome.txt', sep = "\t")
my.2.ECM     <- read.table('./GeneSets/Collagens.txt', sep = "\t")
my.3.ECM     <- read.table('./GeneSets/ECMaffiliated.txt', sep = "\t")
my.4.ECM     <- read.table('./GeneSets/ECMglycoproteins.txt', sep = "\t")
my.5.ECM     <- read.table('./GeneSets/ECMregulators.txt', sep = "\t")
my.6.ECM     <- read.table('./GeneSets/Proteoglycans.txt', sep = "\t")
my.7.ECM     <- read.table('./GeneSets/Matrisomeassociated.txt', sep = "\t")
my.8.ECM     <- read.table('./GeneSets/Secretedfactors.txt', sep = "\t")

# Cannot plot genes with null variance
my.var.null <- apply(tissue.cts,1,var)

my.1.ECM     <- intersect(unique(my.1.ECM$V1) , rownames(tissue.cts)[my.var.null != 0])
my.2.ECM     <- intersect(unique(my.2.ECM$V1) , rownames(tissue.cts)[my.var.null != 0])
my.3.ECM     <- intersect(unique(my.3.ECM$V1) , rownames(tissue.cts)[my.var.null != 0])
my.4.ECM       <- intersect(unique(my.4.ECM$V1)   , rownames(tissue.cts)[my.var.null != 0])
my.5.ECM       <- intersect(unique(my.5.ECM$V1)   , rownames(tissue.cts)[my.var.null != 0])
my.6.ECM       <- intersect(unique(my.6.ECM$V1)   , rownames(tissue.cts)[my.var.null != 0])
my.7.ECM       <- intersect(unique(my.7.ECM$V1)   , rownames(tissue.cts)[my.var.null != 0])
my.8.ECM       <- intersect(unique(my.8.ECM$V1)   , rownames(tissue.cts)[my.var.null != 0])

pdf(paste0(Sys.Date(),"_CoreMatrisome.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.1.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="CoreMatrisome")
dev.off()

pdf(paste0(Sys.Date(),"_Collagens.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.2.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Collagens")
dev.off()

pdf(paste0(Sys.Date(),"_ECMaffiliated.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.3.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="ECM affiliated")
dev.off()

pdf(paste0(Sys.Date(),"_ECMglycoproteins.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.4.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="ECM glycoproteins")
dev.off()

pdf(paste0(Sys.Date(),"_ECMregulators.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.5.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="ECM regulators")
dev.off()

pdf(paste0(Sys.Date(),"_Proteoglycans.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.6.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Proteoglycans")
dev.off()

pdf(paste0(Sys.Date(),"_Matrisomeassociated.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.7.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Matrisome associated")
dev.off()

pdf(paste0(Sys.Date(),"_Secretedfactors.pdf"),width = 10, height = 15, onefile = F)
pheatmap(tissue.cts[my.8.ECM,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Secreted Factors")
dev.off()


################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Heatmap_session_Info.txt", sep =""))
sessionInfo()
sink()

