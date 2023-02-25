setwd("/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2022_Neutrophil_ECM_Collaboration/Code/4_Xlinked_gene_enrichment_analysis")
options(stringsAsFactors = FALSE)

library(Vennerable)        # Vennerable_3.1.0.9000

# Analyze enrichment of X-linked genes for female-biased Salmon module

####################################################################################################
# 1. Load in DESeq2 results and WGCNA Salmon module

# load DEseq2 results from Lu et al, 2021
load('2020-05-21_Neutrophils_SEX.RData')

res.sex          <- data.frame(res.sex)
res.sex$GeneName <- rownames(res.sex)

# load Salmon module genes
my.salmon     <- read.table('../3_Neutrophil_Bulk_RNAseq_WGCNA/Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-salmon.txt', sep = "\t", header = T)$nodeName

# load in X-linked gene list downloaded from ENSMBL Biomart
Xgenes <- read.csv("2023-02-02_Ens108_mouse_xlinked_genes.txt", header = T, sep ="\t")

# Prepare TERM2GENE file for enrichment
colnames(Xgenes) <- c("GenestableID","gene",	"term")
Xlinkedgenes     <- unique(Xgenes[, c(3,2)]) # 2623

# Only genes detected in RNAseq for fair comparison
Xlinkedgenes <- Xlinkedgenes[Xlinkedgenes$gene %in% res.sex$GeneName, ] #380



##### Venn Diagrams using Vennerable
my.salmon.X <- list("X-linked"  = Xlinkedgenes$gene   ,
                        "Salmon"    = my.salmon  )
my.Venn <- Venn(my.salmon.X)

# without scaling
pdf(paste0(Sys.Date(),"_Salmon_Venn.pdf"))
plot(my.Venn, doWeights=T, show = list(FaceText = "weight", SetLabels = TRUE, Faces = FALSE))
dev.off()

# test enrichment (overlap greater than by chance)
fisher.test(matrix(c(6, 
                     374,
                     312,
                     nrow(res.sex)-6-374-312),
                   2,2), 
            alternative = 'greater') # test if overlap is larger than expected
# Fisher's Exact Test for Count Data
# data:  matrix(c(6, 374, 312, nrow(res.sex) - 6 - 374 - 312), 2, 2)
# p-value = 0.9153
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.2671956       Inf
# sample estimates:
#   odds ratio 
# 0.6215153 

####################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()
