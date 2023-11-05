## Emily Hodgson
## Last modified: 3.11.2023

# This script downloads all available sequence data for a given list of genera. (not including metadata, that's another script, but these scripts could maybe be combined to save time!)

# I found searching by genus easier & more effective than searching by species.

# Inputs required:
## Personal API key (increases your e-utils limit to 10 requests/second)
## List of genera (I used csv file)

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
setwd("/Users/emilyhodgson/Documents/Autophylo/")
path_to_output_dir <- "/Users/emilyhodgson/Documents/Autophylo/Available_Seqs/"

# ADD FUNCTION TO ENVIRONMENT
# This function outputs a FASTA file with all the available sequences in NCBI for your search term.
# Set 'path_to_output_dir' (above) as the place you want your file to go
single_attempt_NCBI_seq_search <- function(search_term) {
  master_fasta_list <- list()
  r_search <- entrez_search(db = "taxonomy", 
                            term = paste0(search_term, "[subtree]"),
                            retmax = 999, 
                            use_history = TRUE) # gets all IDs for genus
  cat("\nNumber of ", search_term, " IDs in NCBI taxonomy database: ", length(r_search$ids), "\n")
  if(length(r_search$ids) != 0) {
      # FETCHING SEQUENCE DATA FOR EVERY ID
      loop <- 1
      loop_v <- vector()
      all_recs_list <- list()
      unavailable_list <- matrix(ncol = 2)
      colnames(unavailable_list) <- c("Tax_ID", "Skipped_taxa")
      for (i in r_search$ids) {
        cat("\nID ", loop, ":\t", i, "\t")
        loop_v[loop] <- loop
        upload <- entrez_post(db = "taxonomy", 
                              id = i) # getting the taxonomy with IDs
        fetch_id <- entrez_link(dbfrom = "taxonomy", 
                                db = "nuccore", 
                                web_history = upload) # linking between the nucleotide & taxonomy databases - this is when the inconsistencies come in
        fetch_id2 <- fetch_id$links$taxonomy_nuccore
        if(!is.null(fetch_id2)) {
          all_recs_list[[loop]] <- entrez_fetch(db = "nuccore",
                                                id = fetch_id2, 
                                                rettype = "fasta")
          loop <- loop + 1
        }
        if(is.null(fetch_id2)) { 
          cat("no sequence data")
          unavailable_list <- rbind(unavailable_list, 
                                    c(i, 
                                      entrez_summary(db = "taxonomy", 
                                                     id = i)$scientificname))}
      }
      fasta <- str_split(all_recs_list, pattern = "\n")
      if(isEmpty(fasta) == FALSE) {
        names_v <- vector()
        fasta_list <- list()
        count <- 1
        for(i in 1:length(fasta)){
          fasta_temp <- fasta[[i]]
          idx <- str_detect(fasta_temp, ">")
          if(sum(idx) == 1){
            names_v[count] <- fasta_temp[idx]
            fasta_list[[count]] <- str_c(fasta_temp[!idx], collapse = "")
            count <- count + 1
          }
          if(sum(idx) > 1){
            idx2 <- which(idx)
            idx2 <- c(idx2, length(idx) + 1)
            j <- 2
            for(j in 1:(length(idx2) - 1)){
              names_v[count] <- fasta_temp[idx2[j]]
              fasta_list[[count]] <- str_c(fasta_temp[(idx2[j] + 1):(idx2[j+1] - 1)], collapse = "")
              count <- count + 1
            }
          }
        }
        names(fasta_list) <- str_sub(names_v, start = 2)
      }
  }
  if(length(r_search$ids) == 0) {cat("\nThis genus has 0 accessions in NCBI taxonomy database\n")}
  cat("\n\n", length(fasta_list), " sequences found for ", paste0(search_term), "\n", sep = "")
  setwd(path_to_output_dir)
  write.fasta(sequences = fasta_list, names = names(fasta_list), file.out = paste0(search_term, "_available_seqs.fasta"))
}

# TEST THE FUNCTION WITH A SMALL GENUS
single_attempt_NCBI_seq_search("Glaziella")

# GET LIST OF GENERA
setwd("/Users/emilyhodgson/Documents/Autophylo/")
genera_df <- read.csv(file = 'Otideaceae_genera_list.csv')

# RUN FOR YOUR LIST OF GENERA 
for(g in genera_df$Genus) {
  ten_NCBI_seq_search(g)
}