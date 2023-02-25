# setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2022/2022_Neutrophil_ECM_Collaboration/Code/1_Neutrophil_Bulk_RNAseq_Enrichment')
options(stringsAsFactors = F)

library(DESeq2)
library(phenoTest)

# 2022-09-12
# Analysis with ECM related gene lists looking at differences in Aging
# update GSEA parameters for reproducibility (seed, number of permutations)


################### A. Load DEseq2 results for analysis ########################
# load DEseq2 results
load('2020-05-21_Neutrophils_Aging_BOTH.RData')
my.ECM.age           <- data.frame(res.age)
my.ECM.age$GeneName  <- rownames(my.ECM.age )

########## B. Prepare GeneLists using DEseq2 t-statistic to rank genes #########
ECM.age.geneList         <- my.ECM.age$stat
names(ECM.age.geneList)  <- my.ECM.age$GeneName
ECM.age.geneList         <- sort(ECM.sex.geneList , decreasing = TRUE)

######################## C. Prep ECM gene Sets #################################
my.1.ECM     <- read.table('./GeneSets/CoreMatrisome.txt', sep = "\t", header = T)
my.2.ECM     <- read.table('./GeneSets/Collagens.txt', sep = "\t", header = T)
my.3.ECM     <- read.table('./GeneSets/ECMaffiliated.txt', sep = "\t", header = T)
my.4.ECM     <- read.table('./GeneSets/ECMglycoproteins.txt', sep = "\t", header = T)
my.5.ECM     <- read.table('./GeneSets/ECMregulators.txt', sep = "\t", header = T)
my.6.ECM     <- read.table('./GeneSets/Proteoglycans.txt', sep = "\t", header = T)
my.7.ECM     <- read.table('./GeneSets/Matrisomeassociated.txt', sep = "\t", header = T)
my.8.ECM     <- read.table('./GeneSets/Secretedfactors.txt', sep = "\t", header = T)

my.ECM.curated.gs <- list("CoreMatrisome"       = my.1.ECM$GeneSymbol,
                          "Collagens"           = my.2.ECM$GeneSymbol,
                          "ECMaffiliated"       = my.3.ECM$GeneSymbol,
                          "ECMglycoproteins"    = my.4.ECM$GeneSymbol,
                          "ECMregulators"       = my.5.ECM$GeneSymbol,
                          "Proteoglycans"       = my.6.ECM$GeneSymbol,
                          "Matrisomeassociated" = my.7.ECM$GeneSymbol,
                          "Secretedfactors"     = my.8.ECM$GeneSymbol)


######################## C. Gene Set Enrichment Analysis #######################
# set seed to stabilize output
set.seed(123456789)
# run phenotest GSEA
gsea.data <- gsea( x         =  ECM.age.geneList ,
                   gsets     =  my.ECM.curated.gs,
                   mc.cores  =  1                 ,
                   logScale  =  FALSE             ,
                   B         =  10000            ,
                   minGenes  =  5                 ,
                   maxGenes  =  5000              ,
                   center=TRUE)
my.summary <- data.frame(gsea.data$significance$summary)


# write results to file
my.outfile <- paste(Sys.Date(), "Neutrophil_age_ECM_Gene_Lists_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary, file = my.outfile, quote = F, sep = "\t")

save(gsea.data , file = paste(Sys.Date(),"Neutrophil_age_ECM_Gene_Lists_GSEA.RData", sep = "_") )
# no significant pathways


#######################
sink(file = paste(Sys.Date(),"_GSEA_AGING_session_Info.txt", sep =""))
sessionInfo()
sink()
