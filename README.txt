## README ## 

Code and data for the manuscript:

Cassandra McGill, Collin Y. Ewald and Bérénice A. Benayoun
"Sex-dimorphic expression of extracellular matrix genes in mouse bone marrow neutrophils"
         (biorXiv preprint at: https://www.biorxiv.org/content/10.1101/2023.02.25.530027v1)

The code is arranged by analysis:

	# 1_Neutrophil_Bulk_RNAseq_Enrichment:
			- 1_Run_GSEA_Neutrophil_SEX_addORA.R     : GSEA of curated ECM related genesets as a function of sex + Fisher enrichement
			- 2_Run_GSEA_Neutrophil_ECM_Aging_BOTH.R : GSEA of curated ECM related genesets as a function of age
			- 3_Plot_ECM_Heatmap_sex_v2.R            : Expression Heatmap plotting for ECM related genesets 
			
			# Previous data files from Lu et al, Nat Aging, 2021: (provided)
				%%% 2020-05-21_Neutrophils_SEX.RData
				%%% 2020-05-21_Neutrophils_Aging_BOTH.RData
				%%% 2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt

         	# Gene sets curated by Ewald lab:
         		%%% GeneSets
         
         
	# 2_Neutrophil_scRNAseq_ECM_Enrichment:
			- 2022-10-14_analyze_ECM_genes_scRNAseq_v2.R : UCell analysis of curated ECM related genesets in neutrophil scRNAseq
						
			# Previous Seurat file from Kim et al, Sci Data, 2022:
				%%% 2022-04-12_10x_BM_Ntph_USC_Xie_SingleCellNet_predictions_Seurat_object.RData,
					 available at https://figshare.com/articles/dataset/Annotated_Seurat/19623978
			
			
	# 3_Neutrophil_Bulk_RNAseq_WGCNA: 
			- 1_WGCNA_ECM_v2.R           : Run WGCNA on bulk RNAseq
			- 2_Run_WGCNA_module_ORA.R   : Run Overrepresentation analysis for key WGCNA modules
			
			
	# 4_Xlinked_gene_enrichment_analysis:
			- 1_X_linked_enrich.R           : Calculate fisher exact's test for X-linked genes in Salmon cluster

			# Previous data files from Lu et al, Nat Aging, 2021: (provided)
				%%% 2020-05-21_Neutrophils_SEX.RData

			# Downloaded annotation from ENSEMBL Biomart: (provided)
				%%% 2023-02-02_Ens108_mouse_xlinked_genes.txt

	# 5_Serum_Proteomics_enrichment:
			- 1_proteome_ECM_SERUM.R    : Analysis of serum proteomics for sex-biased proteins and ECM enrichments
			
			# Previous data files from Aumailley et al, J. Proteome Res., 2021: 
				available at https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.1c00542/suppl_file/pr1c00542_si_006.xlsx
