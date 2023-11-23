# This file describes the dependencies of the whole workflow. For now I copy/pasted Emily's setup from the 1_Dowload_seqs.R script

# SETUP
package.setup <- function() {
  if (!require("rentrez", quietly = TRUE)){
    install.packages("rentrez")
  } 
  if (!require("seqinr", quietly = TRUE)){
    install.packages("seqinr")
  }
  if (!require("stringr", quietly = TRUE)){
    install.packages("stringr")
  }
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  if (!require("Biostrings", quietly = TRUE)){
    BiocManager::install("Biostrings")
  }
  library(rentrez)
  library(seqinr)
  library(stringr)
  library("Biostrings")
  api_key <- "4bb20e27b9f0e52e14014832f00d2f139a08"
  set_entrez_key(api_key)
}
package.setup() # one function to install and load all required packages & libraries & set API key (N.B. don't share your key)
setwd("/Users/emilyhodgson/Documents/Autophylo/")  # change the path for your own use