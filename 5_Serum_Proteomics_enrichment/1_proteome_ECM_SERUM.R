setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2022_Neutrophil_ECM_Collaboration/Code/5_Serum_Proteomics_enrichment')
options(stringsAsFactors = F)

library('readxl')
library(limma)
library('pheatmap')
library(phenoTest)
library(bitops)

# 2023-03-27
# Analysis with ECM related gene lists
# Proteomics datasets

### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  ###
###               Mouse female vs. male serum proteomics /// Aumailley et al, 2021                   ###
### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  ###

########################################################################################################
# 1. Preprocess data

# read proteomic quantification
serum.prot <- read_xlsx('pr1c00542_si_006.xlsx', sheet = 2)

## LFQ /// LFQ (label-free quantitation): protein intensities are normalized to exclude some outliers to best represent the ratio changes of different samples.
WT.LFQ.int.col <- intersect(grep("FWT|MWT",colnames(serum.prot)),grep("LFQ",colnames(serum.prot)))

# get data from WT Male and Female (don't care about VitC)
serum.prot.cl <- serum.prot[,c(1:4,WT.LFQ.int.col)]

# rename columns
colnames(serum.prot.cl) <- gsub(" ", "_",colnames(serum.prot.cl))
colnames(serum.prot.cl) <- gsub("LFQ_intensity_", "",colnames(serum.prot.cl))

# collapse based on gene name
serum.prot.cl.2 <- aggregate(serum.prot.cl[,-c(1:4)], by = list(serum.prot.cl$Gene_names), FUN = "mean")
serum.prot.cl.2 <- serum.prot.cl.2[!is.na(serum.prot.cl.2$Group.1),]

# do MDS analysis
mds.result <- cmdscale(1-cor(serum.prot.cl.2[,-1],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.colors.sex <- rep("",length(colnames(serum.prot.cl.2))-1)
my.colors.sex[grep("FGL|FWT",colnames(serum.prot.cl.2)[-1])] <- "deeppink"
my.colors.sex[grep("MGL|MWT",colnames(serum.prot.cl.2)[-1])] <- "deepskyblue"

my.pch.geno <- rep(NA,length(colnames(serum.prot.cl.2)))
my.pch.geno[grep("FGL|MGL",colnames(serum.prot.cl.2)[-1])] <- 15
my.pch.geno[grep("FWT|MWT",colnames(serum.prot.cl.2)[-1])] <- 16

pdf(paste(Sys.Date(),"Aumailley_SerumProteome_2021_MDS_plot_ALL.pdf",sep="_"))
plot(x, y,
     pch = my.pch.geno, col = my.colors.sex,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",cex=3)
legend("topright", c("Female","Male"), col = c("deeppink","deepskyblue"), pch = 16, bty = 'n', pt.cex = 2)
dev.off()


########################################################################################################
# 2. Limma differential expression analysis

# use sex as covariate
rownames(serum.prot.cl.2) <- serum.prot.cl.2$Group.1

my.sex <- rep("",length(colnames(serum.prot.cl.2))-1)
my.sex[grep("FGL|FWT",colnames(serum.prot.cl.2)[-1])] <- "F"
my.sex[grep("MGL|MWT",colnames(serum.prot.cl.2)[-1])] <- "M"
my.sex <- factor(my.sex, levels=c("M","F")) # male bias is negative logFC

my.treat <- as.numeric(sapply(strsplit(gsub("FWT|MWT","",colnames(serum.prot.cl.2)[-1]),"_"),"[[",1))

# build model matrix
my.model <- model.matrix( ~ my.sex + my.treat)

# fit limma model
my.fit <- lmFit(serum.prot.cl.2[,-c(1)], my.model)
my.fit <- eBayes(my.fit)

# get output table
my.sig      <- topTable(my.fit, coef = "my.sexF", p.value = 1 , number = Inf)
sig.prots   <- rownames(my.sig)[my.sig$adj.P.Val < 0.05]
n.sig.prots <- length(sig.prots) # 115

# output result table
my.output <- paste(Sys.Date(),"Aumailley_SerumProteome_2021.limma.results_ALL.txt",sep="_")
write.table(my.sig, file = my.output, sep = "\t", quote = F)

# get heatmap
my.heatmap.out <- paste(Sys.Date(),"Aumailley_SerumProteome_2021_Heatmap_significant_proteins_ALL.pdf",sep="_")

pdf(my.heatmap.out, height = 5, width = 9, onefile = F)
my.heatmap.title <- paste("Sex serum significant (FDR < 5%), ",n.sig.prots, " proteins",sep="")
pheatmap(serum.prot.cl.2[rownames(my.sig)[my.sig$adj.P.Val < 0.05],-1],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 10,
         border = NA)
dev.off()

# get significant proteins
my.F.up <- rownames(my.sig)[bitAnd(my.sig$adj.P.Val < 0.05,my.sig$logFC >0)>0]
my.M.up <- rownames(my.sig)[bitAnd(my.sig$adj.P.Val < 0.05,my.sig$logFC <0)>0]


########################################################################################################
# 3. ORA analysis with Fisher

### A. prep gene list
ECM.sex.geneList         <- my.sig$t
names(ECM.sex.geneList)  <- rownames(my.sig)
ECM.sex.geneList         <- sort(ECM.sex.geneList , decreasing = TRUE)

### B. Prep ECM gene Sets 
my.1.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/CoreMatrisome.txt', sep = "\t", header = T)
my.2.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/Collagens.txt', sep = "\t", header = T)
my.3.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/ECMaffiliated.txt', sep = "\t", header = T)
my.4.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/ECMglycoproteins.txt', sep = "\t", header = T)
my.5.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/ECMregulators.txt', sep = "\t", header = T)
my.6.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/Proteoglycans.txt', sep = "\t", header = T)
my.7.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/Matrisomeassociated.txt', sep = "\t", header = T)
my.8.ECM     <- read.table('../1_Neutrophil_Bulk_RNAseq_Enrichment/GeneSets/Secretedfactors.txt', sep = "\t", header = T)

# only detected
my.ECM.curated.gs <- list("CoreMatrisome"       = intersect(rownames(my.sig),my.1.ECM$GeneSymbol) ,
                          "Collagens"           = intersect(rownames(my.sig),my.2.ECM$GeneSymbol) ,
                          "ECMaffiliated"       = intersect(rownames(my.sig),my.3.ECM$GeneSymbol) ,
                          "ECMglycoproteins"    = intersect(rownames(my.sig),my.4.ECM$GeneSymbol) ,
                          "ECMregulators"       = intersect(rownames(my.sig),my.5.ECM$GeneSymbol) ,
                          "Proteoglycans"       = intersect(rownames(my.sig),my.6.ECM$GeneSymbol) ,
                          "Matrisomeassociated" = intersect(rownames(my.sig),my.7.ECM$GeneSymbol) ,
                          "Secretedfactors"     = intersect(rownames(my.sig),my.8.ECM$GeneSymbol) )


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
colnames(my.fisher.res) <- c("n","N","FoldEnrich","enrichPval", "enrichFDR","proteins")

# Female
for (i in 1:8) {
  
  # calculate contigency table
  a <- length(intersect(my.F.up,my.ECM.curated.gs[[i]]) ) 
  b <- length(setdiff(my.F.up,my.ECM.curated.gs[[i]])   )
  c <- length(setdiff(my.ECM.curated.gs[[i]],my.F.up)   )
  d <- length(rownames(my.sig)) - a - b -c
  
  # get fisher test
  tmp.test <- fisher.test(matrix(c(a,b,c,d),2,2), alternative = 'greater') # test if overlap is larger than expected
  my.fisher.res$enrichPval[i] <- tmp.test$p.value
  
  # calculate fold enrichment
  my.fisher.res$FoldEnrich[i] <- (a/(a+b))  / (c/(c+d))
  my.fisher.res$n[i]            <- a # number of overlapping
  my.fisher.res$N[i]            <- length(my.ECM.curated.gs[[i]])
  my.fisher.res$proteins[i]        <- paste(intersect(my.F.up,my.ECM.curated.gs[[i]]),collapse = ";")
  
}

my.fisher.res$enrichFDR <- p.adjust(my.fisher.res$enrichPval, method = "BH")

# write results to file
my.outfile <- paste(Sys.Date(), "Serum_FemaleBias_ECM_Gene_Lists_Fisher_Analysis_table.txt", sep = "_")
write.table(my.fisher.res, file = my.outfile, quote = F, sep = "\t")
################################################################################################

################################################################################################
##### 4. Plot ORA results
library(scales) 

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
my.ora.sex.sorted$Sample <- factor(rep("Serum",length(  my.ora.sex.sorted$Description)))

# Female/Male color scale
my.max <- max(my.ora.sex.sorted$FoldEnrich)
my.min <- min(my.ora.sex.sorted$FoldEnrich)

my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("white","lavenderblush","plum1","orchid1","deeppink")

my.plot <- ggplot(my.ora.sex.sorted, aes(x=Sample,y=Description,colour=FoldEnrich,size=minlog10fdr)) + theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("ora Analysis") + labs(x = "-log10(pvalue)", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, na.value = "grey50", guide = "colourbar", values = my.scaled, limits=c(0,6))
my.plot


my.pdfname <- paste(Sys.Date(),"ORA_FemaleSerum_BALLOON_plot_ECM_genesets.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 6, width=max(4,my.max.char/8) )
print(my.plot)
dev.off()  
################################################################################################

#######################
sink(file = paste(Sys.Date(),"_SERUM_proteomics_Sex_session_Info.txt", sep =""))
sessionInfo()
sink()
