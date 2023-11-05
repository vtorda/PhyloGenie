## Emily Hodgson
## Last modified: 4.11.2023

# A script to download extra information uploaded on GenBank associated with the sequences we want to download. 

# Eventually merge with the script that downloads seqs?

package.setup <- function() {
  if (!require("readr", quietly = TRUE)){
    install.packages("readr")
  } 
  if (!require("seqinr", quietly = TRUE)){
    install.packages("seqinr")
  }
  if (!require("stringr", quietly = TRUE)){
    install.packages("stringr")
  }
  if (!require("rentrez", quietly = TRUE)){
    install.packages("rentrez")
  }
  if (!require("restez", quietly = TRUE)){
    install.packages("restez") # # this is the package containing the gb_extract() function, the main tool in this script for extracting the seq. info!
  }
  library(readr)
  library(seqinr)
  library(stringr)
  library(rentrez)
  library(restez)
}
package.setup()

# Remove ">" glyph function
remove_glyph <- function(x) {
  gsub(pattern = ">", replacement = "", x)
}

seq_dir <- "/Users/emilyhodgson/Documents/Autophylo/Loci_sorted_sequences/"
setwd(seq_dir)
seq_fasta_files <- list.files(path = seq_dir)

# Get basic info from annotations first
# Fast

metadata_df <- vector()

for(x in 1:length(seq_fasta_files)) {
  locus <- str_remove(string = seq_fasta_files[x], pattern = "_Res.fasta")
  seq_fastas <- read.fasta(file = seq_fasta_files[x])
  seq_annots <- getAnnot(seq_fastas)
  seq_annots_2 <- unlist(lapply(seq_annots, remove_glyph))
  annot_split <- str_split(string = seq_annots_2, pattern = " ")
  seq_accessions <- sapply(annot_split, "[[", 1)
  locus_df <- as.data.frame(cbind(locus, seq_accessions, seq_annots_2))
  metadata_df <- rbind(metadata_df, locus_df)
}

# Download the NCBI data and fetch the features
# Slow!

ncbi_data <- vector()
features <- c()

for(x in 1:nrow(metadata_df)) {
  ncbi_data[x] <- rentrez::entrez_fetch(db = "nuccore",
                                             id = metadata_df$seq_accessions[x],
                                             rettype = "gb",
                                             retmode = "text")
  features[x] <- gb_extract(record = ncbi_data[x], what = "features")
  cat("Accession ", x, " of ", nrow(metadata_df), "found")
}

# Warnings: In features[x] <- gb_extract(record = ncbi_data[x], what = "features"): number of items to replace is not a multiple of replacement length

# setwd("/Users/emilyhodgson/Documents/Autophylo/")
# save(ncbi_data, features, file = "NCBI_data.RData")

# work out what all the feature options are
unique_features <- vector()
for(x in 1:length(features)) {
  unique_features <- unique(c(names(features[[x]]), unique_features))
}

metadata_df[unique_features] <- NA

# check 
length(features) == nrow(metadata_df)

for(x in 1:length(unique_features)) {
  for(z in 1:length(features)) {
    if(is.null(features[[z]][[unique_features[x]]]) == TRUE) {
      metadata_df[z, unique_features[x]] <- "NULL"
    } else {
      metadata_df[z, unique_features[x]] <- features[[z]][[unique_features[x]]]
    }
  }
}
