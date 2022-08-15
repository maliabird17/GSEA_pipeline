# Loading Libraries #
library(tidyverse)
library(janitor)
library(AnnotationDbi)
library(readxl)
library("org.Dr.eg.db")
library("org.Hs.eg.db")
library(dplyr)

# OPTIMIZATION OF ZEBRAFISH TO HUMAN MAPPING #
column_name = "gene_id"

# LOADING DATA FROM DROPBOX 
dropbox_to_dataframe <- function(link, h = TRUE) {
  raw_link <- gsub("dl=0", "raw=1", link) # exracting raw file from dropbox link
  df <- data.frame(read.csv(raw_link, sep = "\t", header = h,  stringsAsFactors = FALSE)) # load as dataframe
  df # return dataframe
}

# ENSEMBL ID METHOD
mapping_efficiency_fish_to_human_ENSEMBL <- function(fish_data, column_name, export_name) {
  
  # Load from dropbox 
  df <- dropbox_to_dataframe(fish_data)
  
  samples <- c(colnames(df))
  preprocessed <- nrow(df) # number of genes input
  
  # Load ensembl ortholog map (zebrafish and human ensembl ids)
  ensembl_data <- data.frame(read.csv("https://www.dropbox.com/s/93kjxvac9jjnyxq/ensembl_human_zebra_map.txt?raw=1", sep = ",", header = TRUE)) %>% 
    janitor::clean_names() %>% 
    dplyr::select(gene_stable_id, human_gene_stable_id) %>% 
    distinct()
  
  # Merging ensembl database onto input transcriptomic data
  data_mapped <- df %>%
    merge(ensembl_data, by.x = column_name , by.y = "gene_stable_id") %>% 
    na_if("") %>%
    filter(!is.na(human_gene_stable_id))
  
  distinct(data_mapped, human_gene_stable_id, .keep_all = TRUE)
  
  # Calculating efficiency of mapping 
  post_processed <- nrow(data_mapped) # number of genes mapped
  efficiency = (post_processed/preprocessed)*100 # efficiency of mapping
  print(c('Percent of genes mapped to human genes:', efficiency))

}

# ZFIN ID METHOD
mapping_efficiency_fish_to_human_ZFIN <- function(fish_data, column_name, export_name) {
  
  # Load from dropbox 
  df <- dropbox_to_dataframe(fish_data) %>% 
    rename('ensembl_id' = 'gene_id')
  
  samples <- c(colnames(df))
  preprocessed <- nrow(df) # number of genes input
  
  # LOAD ZFIN ORTHOLOG/HOMOLOG DATA (HUMAN:ZEBRAFISH VIA ZFIN)
  zfin_ortholog_data <- read.csv("https://raw.githubusercontent.com/maliabird17/Telomere-Lab/main/human_orthos_2022.06.07.txt", skip = 1, sep = "\t", header = TRUE)
  DataFrame(zfin_ortholog_data)
  zfin_data <- zfin_ortholog_data %>% 
    clean_names() %>% 
    dplyr::select(zfin_id, gene_id) %>% 
    distinct()
  
  ### RETREIVE HUMAN ENSEMBL IDS ASSOCIATED WITH ZFIN IDS 
  gene_keys <- as.character(zfin_data$gene_id)
  zfin_data['ensembl_id'] <- mapIds(org.Hs.eg.db, keys = gene_keys, keytype = "ENTREZID", column = "ENSEMBL")
  
  ensembl <- df$ensembl_id
  
  df['zfin_id'] <- mapIds(org.Dr.eg.db, keys = ensembl, keytype = 'ENSEMBL', column = 'ZFIN') 
  df <- merge(df, zfin_data, by = 'zfin_id') 
  
  distinct(df, ensembl_id.y, .keep_all = T)
  
  # Calculating efficiency of mapping 
  post_processed <- nrow(df) # number of genes mapped
  efficiency = (post_processed/preprocessed) *100 # efficiency of mapping
  print(c('Percent of genes mapped to human genes:', efficiency))
  
}

# TESTIS 
fish_data = "https://www.dropbox.com/s/26v2vmdy2zn3cmk/Batch_noCre_WT.normalized_countsTESTIS.txt?dl=0"
mapping_efficiency_fish_to_human_ENSEMBL(fish_data, column_name, export_name) # 78.6%
mapping_efficiency_fish_to_human_ZFIN(fish_data, column_name, export_name) # 62.3%

#KIDNEY
fish_data = "https://www.dropbox.com/s/qjy69mqpljzizba/Batch_noCre_WT.normalized_counts.txt?dl=0"
mapping_efficiency_fish_to_human_ENSEMBL(fish_data, column_name, export_name) # 78.7%
mapping_efficiency_fish_to_human_ZFIN(fish_data, column_name, export_name) # 62.3%

#GUT
fish_data = "https://www.dropbox.com/s/jvvlg867uwaufhw/Batch_Zebrafish_Gut.normalized_counts_nocre_wt.txt?dl=0"
mapping_efficiency_fish_to_human_ENSEMBL(fish_data, column_name, export_name) # 86.2%
mapping_efficiency_fish_to_human_ZFIN(fish_data, column_name, export_name) # 69.7%

