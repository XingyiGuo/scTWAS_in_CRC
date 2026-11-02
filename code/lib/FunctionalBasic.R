#conda activate twasguo
library("stringr")
library("data.table")
library("dplyr")
library("readxl")

### CRISPRGeneEffect.csv [https://depmap.org/portal/; https://www.nature.com/articles/s41588-019-0395-x]
# Pipeline: Achilles
# _Post-Chronos_
# Gene effect estimates for all models, integrated using Chronos.  Copy number corrected, scaled, and screen quality corrected.
# - Columns: Gene - Rows: ModelID

w_k="/data/l2_bioinfo1/liq17/LQ_references/GeneFunctioalSupports/FunctionalDatabases/"

get_depmap_diseases <- function(){
  depmap_CRISPRCAS9_meta=paste0(w_k, "depmap/Model_26Q1.csv")
  depmap_CRISPRCAS9_filename=paste0(w_k, "depmap/CRISPRGeneEffect_26Q1.csv")
  
  #Get a list contains cell lines for the disease
  depmap_meta=fread(depmap_CRISPRCAS9_meta)
  return(unique(depmap_meta$OncotreePrimaryDisease))
}

#Colorectal Adenocarcinoma
get_depmap_ceres_scores <- function(genelist, disease){
  depmap_CRISPRCAS9_meta=paste0(w_k, "depmap/Model_26Q1.csv")
  depmap_CRISPRCAS9_filename=paste0(w_k, "depmap/CRISPRGeneEffect_26Q1.csv")
  
  #Get a list contains cell lines for the disease
  depmap_meta=fread(depmap_CRISPRCAS9_meta)
  depmap_disease_modelID=unlist(depmap_meta[depmap_meta$OncotreePrimaryDisease == disease, "ModelID"])
  cat("The numebr of cell lines related to ", disease, " is ", as.character(length(depmap_disease_modelID)), "\n")
  
  #Get CERES socres for genes in the genelist
  deepmap_effect=fread(depmap_CRISPRCAS9_filename)
  deepmap_effect_modelIDs = filter(deepmap_effect, V1 %in% depmap_disease_modelID) 
  cat("The numebr of cell lines DepMap scores related to ", disease, " is ", as.character(dim(deepmap_effect_modelIDs)[1]), "\n")
  
  new_colnames <- str_extract(colnames(deepmap_effect_modelIDs), "\\w+(?=\\s\\(\\d+\\))")
  # Replace NAs with original names
  new_colnames[is.na(new_colnames)] <- colnames(deepmap_effect_modelIDs)[is.na(new_colnames)]
  colnames(deepmap_effect_modelIDs) <- new_colnames
  
  common_genes=intersect(as.character(new_colnames), as.character(genelist))
  deepmap_effect_modelIDs_genes=select(deepmap_effect_modelIDs, all_of(common_genes))
  
   #Mean values
  #deepmap_effect_modelIDs_genes_mean = as.data.frame(sapply(deepmap_effect_modelIDs_genes, mean, na.rm = TRUE))
  #colnames(deepmap_effect_modelIDs_genes_mean) = c("depmap_mean_CERES") 
  #genes_df = data.frame("Genes"=genelist)
  #deepmap_effect_modelIDs_genes_mean["Genes"]=rownames(deepmap_effect_modelIDs_genes_mean)
  #genes_scores_df = merge(genes_df, deepmap_effect_modelIDs_genes_mean, by="Genes", all.x=TRUE)

   #Raw values
   genes_df = data.frame("Genes"=genelist)
   gene_depmap = as.data.frame(t(deepmap_effect_modelIDs_genes))
   colnames(gene_depmap) <- paste0(disease, "_", seq_along(colnames(gene_depmap)))
   gene_depmap['Genes']=rownames(gene_depmap)
   genes_scores_df = merge(genes_df, gene_depmap, by="Genes", all.x=TRUE)
   return(genes_scores_df)
}

get_sc_de_genes <- function(){

}

get_cis_gene_annotations <- function(){

}

get_distal_chromatin_chromatin_interaction <- function(){

}