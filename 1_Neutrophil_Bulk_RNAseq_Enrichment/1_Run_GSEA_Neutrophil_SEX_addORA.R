# setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2022_Neutrophil_ECM_Collaboration/Code/1_Neutrophil_Bulk_RNAseq_Enrichment/')
options(stringsAsFactors = F)

library(DESeq2)
library(phenoTest)

library(ggplot2) 
library(scales) 
theme_set(theme_bw())

library(bitops)

# 2022-09-12
# Analysis with ECM related gene lists based on Sex
# update GSEA parameters for reproducibility (seed, number of permutations)

# 2023-04-01
# add ORA for female biased genes

################################################################################################
######################## 1. Run GSEA analysis

######################## A. Load DEseq2 results for analysis ########################
# load DEseq2 results from Lu et al, 2021
load('2020-05-21_Neutrophils_SEX.RData')
my.ECM.sex           <- data.frame(res.sex)
my.ECM.sex$GeneName  <- rownames(my.ECM.sex )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ########################
ECM.sex.geneList         <- my.ECM.sex$stat
names(ECM.sex.geneList)  <- my.ECM.sex$GeneName
ECM.sex.geneList         <- sort(ECM.sex.geneList , decreasing = TRUE)


######################## C. Prep ECM gene Sets ########################
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


######################## C. Gene Set Enrichment Analysis ########################
# set random seed to stabilize output
set.seed(123456789)

# run phenotest GSEA
gsea.data <- gsea( x         =  ECM.sex.geneList ,
                   gsets     =  my.ECM.curated.gs,
                   logScale  =  FALSE             ,
                   B         =  8000              ,
                   minGenes  =  5                 ,
                   maxGenes  =  5000              ,
                   center = TRUE)
my.summary <- data.frame(gsea.data$significance$summary)


# write results to file
my.outfile <- paste(Sys.Date(), "Neutrophil_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary, file = my.outfile, quote = F, sep = "\t")

save(gsea.data , file = paste(Sys.Date(),"Neutrophil_sex_ECM_Gene_Lists_GSEA.RData", sep = "_") )

pdf(paste(Sys.Date(), "Neutrophil_ECM_CoreMatrisome_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='CoreMatrisome', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_Collagens_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='Collagens', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_ECMaffiliated_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='ECMaffiliated', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_ECMglycoproteins_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='ECMglycoproteins', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_ECMregulators_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='ECMregulators', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_Proteoglycans_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='Proteoglycans', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_Matrisomeassociated_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='Matrisomeassociated', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Neutrophil_ECM_Secretedfactors_GSEA_plot.pdf", sep = "_"))
plot.gseaData(gsea.data, es.nes='nes',selGsets='Secretedfactors', color = "purple")
dev.off()
################################################################################################

################################################################################################
######################## 2. Plot GSEA results

# prep for bubble plot plotting
my.gsea.sex <- my.summary

# create -log10 FDR for plotting
my.gsea.sex$minlog10fdr  <- -log10(my.gsea.sex$fdr + 1e-30)
my.gsea.sex$Description  <- rownames(my.gsea.sex)

# Sort results
my.sorting         <- sort(my.gsea.sex$nes, index.return = T, decreasing = T)
my.gsea.sex.sorted <- my.gsea.sex[my.sorting$ix,]

# create and preserve wanted display order
my.gsea.sex.sorted$sex <- ifelse(my.gsea.sex.sorted$nes < 0, "Male", "Female")  # male/female avg flag
my.gsea.sex.sorted <- my.gsea.sex.sorted[order(my.gsea.sex.sorted$sex), ]

my.max.char <- max(nchar(my.gsea.sex.sorted$Description))

my.gsea.sex.sorted$Description <- factor(my.gsea.sex.sorted$Description, levels = rev(unique(my.gsea.sex.sorted$Description)))
my.gsea.sex.sorted$CellType <- factor(rep("Neutrophils",length(  my.gsea.sex.sorted$Description)))

# Female/Male color scale
my.max <- max(my.gsea.sex.sorted$nes)
my.min <- min(my.gsea.sex.sorted$nes)

my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("white","lavenderblush","plum1","orchid1","deeppink")

my.plot <- ggplot(my.gsea.sex.sorted, aes(x=CellType,y=Description,colour=nes,size=minlog10fdr)) + theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("gsea Analysis") + labs(x = "-log10(pvalue)", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, na.value = "grey50", guide = "colourbar", values = my.scaled, limits=c(0,2.5))
my.plot


my.pdfname <- paste(Sys.Date(),"GSEA_sex_BALLOON_plot_ECM_genesets.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 6, width=max(4,my.max.char/8) )
print(my.plot)
dev.off()  
################################################################################################



################################################################################################
######################## 3. Run ORA analysis

### A. Load DEseq2 results for analysis ########################
# load DEseq2 results from Lu et al, 2021
load('2020-05-21_Neutrophils_SEX.RData')

my.F.up <- rownames(res.sex)[bitAnd(res.sex$padj < 0.05,res.sex$log2FoldChange > 0) > 0] # 1073
my.M.up <- rownames(res.sex)[bitAnd(res.sex$padj < 0.05,res.sex$log2FoldChange < 0) > 0] # 661

### B. Prep ECM gene Sets 
my.1.ECM     <- read.table('./GeneSets/CoreMatrisome.txt', sep = "\t", header = T)
my.2.ECM     <- read.table('./GeneSets/Collagens.txt', sep = "\t", header = T)
my.3.ECM     <- read.table('./GeneSets/ECMaffiliated.txt', sep = "\t", header = T)
my.4.ECM     <- read.table('./GeneSets/ECMglycoproteins.txt', sep = "\t", header = T)
my.5.ECM     <- read.table('./GeneSets/ECMregulators.txt', sep = "\t", header = T)
my.6.ECM     <- read.table('./GeneSets/Proteoglycans.txt', sep = "\t", header = T)
my.7.ECM     <- read.table('./GeneSets/Matrisomeassociated.txt', sep = "\t", header = T)
my.8.ECM     <- read.table('./GeneSets/Secretedfactors.txt', sep = "\t", header = T)

# only detected genesfor Fisher
my.ECM.curated.gs <- list("CoreMatrisome"       = intersect(rownames(res.sex),my.1.ECM$GeneSymbol) ,
                          "Collagens"           = intersect(rownames(res.sex),my.2.ECM$GeneSymbol) ,
                          "ECMaffiliated"       = intersect(rownames(res.sex),my.3.ECM$GeneSymbol) ,
                          "ECMglycoproteins"    = intersect(rownames(res.sex),my.4.ECM$GeneSymbol) ,
                          "ECMregulators"       = intersect(rownames(res.sex),my.5.ECM$GeneSymbol) ,
                          "Proteoglycans"       = intersect(rownames(res.sex),my.6.ECM$GeneSymbol) ,
                          "Matrisomeassociated" = intersect(rownames(res.sex),my.7.ECM$GeneSymbol) ,
                          "Secretedfactors"     = intersect(rownames(res.sex),my.8.ECM$GeneSymbol) )

### C. overepresentation Analysis 
my.fisher.res <- as.data.frame(matrix(0,8,6))
rownames(my.fisher.res) <- c("CoreMatrisome"      ,
                             "Collagens"          ,
                             "ECMaffiliated"      ,
                             "ECMglycoproteins"   ,
                             "ECMregulators"      ,
                             "Proteoglycans"      ,
                             "Matrisomeassociated",
                             "Secretedfactors"    )
colnames(my.fisher.res) <- c("n","N","FoldEnrich","enrichPval", "enrichFDR","genes")

# Female
for (i in 1:8) {
  
  # calculate contigency table
  a <- length(intersect(my.F.up,my.ECM.curated.gs[[i]]) ) 
  b <- length(setdiff(my.F.up,my.ECM.curated.gs[[i]])   )
  c <- length(setdiff(my.ECM.curated.gs[[i]],my.F.up)   )
  d <- length(rownames(res.sex)) - a - b -c
  
  # get fisher test
  tmp.test <- fisher.test(matrix(c(a,b,c,d),2,2), alternative = 'greater') # test if overlap is larger than expected
  my.fisher.res$enrichPval[i] <- tmp.test$p.value
  
  # calculate fold enrichment
  my.fisher.res$FoldEnrich[i] <- (a/(a+b))  / (c/(c+d))
  my.fisher.res$n[i]            <- a # number of overlapping
  my.fisher.res$N[i]            <- length(my.ECM.curated.gs[[i]])
  my.fisher.res$genes[i]        <- paste(intersect(my.F.up,my.ECM.curated.gs[[i]]),collapse = ";")
  
}

my.fisher.res$enrichFDR <- p.adjust(my.fisher.res$enrichPval, method = "BH")

# write results to file
my.outfile <- paste(Sys.Date(), "Serum_FemaleBias_ECM_Gene_Lists_Fisher_Analysis_table.txt", sep = "_")
write.table(my.fisher.res, file = my.outfile, quote = F, sep = "\t")
################################################################################################

################################################################################################
##### 4. Plot ORA results

# prep for bubble plot plotting
my.ora.sex <- my.fisher.res

# create -log10 FDR for plotting
my.ora.sex$minlog10fdr  <- -log10(my.ora.sex$enrichFDR)
my.ora.sex$Description  <- rownames(my.ora.sex)

# Sort results
my.sorting         <- sort(my.ora.sex$FoldEnrich, index.return = T, decreasing = T)
my.ora.sex.sorted <- my.ora.sex[my.sorting$ix,]

# create and preserve wanted display order
my.max.char <- max(nchar(my.ora.sex.sorted$Description))

my.ora.sex.sorted$Description <- factor(my.ora.sex.sorted$Description, levels = rev(unique(my.ora.sex.sorted$Description)))
my.ora.sex.sorted$Sample <- factor(rep("RNA",length(  my.ora.sex.sorted$Description)))

# Female/Male color scale
my.max <- max(my.ora.sex.sorted$FoldEnrich)
my.min <- min(my.ora.sex.sorted$FoldEnrich)

my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("white","lavenderblush","plum1","orchid1","deeppink")

my.plot <- ggplot(my.ora.sex.sorted, aes(x=Sample,y=Description,colour=FoldEnrich,size=minlog10fdr)) + theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("ora Analysis") + labs(x = "-log10(pvalue)", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, na.value = "grey50", guide = "colourbar", values = my.scaled, limits=c(0,8))
my.plot


my.pdfname <- paste(Sys.Date(),"ORA_FemaleNTPH_BALLOON_plot_ECM_genesets.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 6, width=max(4,my.max.char/8) )
print(my.plot)
dev.off()  
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_GSEA_ORA_Sex_session_Info.txt", sep =""))
sessionInfo()
sink()

