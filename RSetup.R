package.setup <- function(workingdir = NULL) {
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
  if(!is.null(workingdir)){
    setwd(workingdir)
    message(paste0("The working directory was set to ", workingdir))
  }
}



