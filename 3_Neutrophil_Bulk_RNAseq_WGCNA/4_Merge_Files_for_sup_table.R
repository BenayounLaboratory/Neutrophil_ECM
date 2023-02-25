setwd("/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/3_Neutrophil_Bulk_RNAseq_WGCNA/")
options(stringsAsFactors = FALSE)

my.salmon     <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-salmon.txt', sep = "\t", header = T)
my.brown      <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-brown.txt', sep = "\t", header = T)
my.red        <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-red.txt', sep = "\t", header = T)
my.purple     <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-purple.txt', sep = "\t", header = T)
my.magenta    <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-magenta.txt', sep = "\t", header = T)
my.black      <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-black.txt', sep = "\t", header = T)
my.yellow     <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-yellow.txt', sep = "\t", header = T)
my.tan        <- read.table('./Output/1_WGCNA_run/Module_GeneLists/merged_CytoscapeInput-nodes-tan.txt', sep = "\t", header = T)

my.merged <- rbind(my.salmon   , 
                   my.brown    , 
                   my.red      , 
                   my.purple   , 
                   my.magenta  , 
                   my.black    , 
                   my.yellow   , 
                   my.tan       )

my.merged <- my.merged[,-2]
colnames(my.merged) <- c("GeneSymbol","Module")

write.table(my.merged   ,  file = paste(Sys.Date(),"Neutrophil_WGCNA_Modules.txt", sep = "_"), quote = F, sep = "\t", row.names = F)
