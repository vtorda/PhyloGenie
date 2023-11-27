## Emily Hodgson
## Last modified: 3.11.2023

# This script downloads all available sequence data for a given list of genera. (not including metadata, that's another script, but these scripts could maybe be combined to save time!)

# I found searching by genus easier & more effective than searching by species.

# Inputs required:
## Personal API key (increases your e-utils limit to 10 requests/second)
## List of genera (I used csv file)
source("RSetup.R")
package.setup(workingdir = "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/")
search_term <- "Glaziella"
api_key <- NULL
path_to_output_dir <- "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/"
repeats <- 5

# ADD FUNCTION TO ENVIRONMENT
# This function outputs a FASTA file with all the available sequences in NCBI for your search term.
# Set 'path_to_output_dir' (above) as the place you want your file to go
NCBI_seq_search <- function(search_term, api_key = NULL, path_to_output_dir, repeats = 10 ) {
  # set api key within the function
  if(!is.null(api_key)){
    set_entrez_key(api_key)
  }else{
    warning("Request an API key by registering NCBI to get a faster download!")
  }
  r_search <- entrez_search(db = "taxonomy", 
                            term = paste0(search_term, "[subtree]"),
                            retmax = 999, 
                            use_history = TRUE) # gets all IDs for genus
  ID_n <- length(r_search$ids) # reducing repetitions by reusing this variable
  cat("\nNumber of ", search_term, " IDs in NCBI taxonomy database: ", ID_n, "\n")
  if(ID_n != 0) {
    write(">0\n", file = paste0(path_to_output_dir, "/", search_term, "_available_seqs.fasta"))
    k <- 2
    for(k in 1:repeats){
    # FETCHING SEQUENCE DATA FOR EVERY ID
    loop <- 1
    all_recs_list <- list()
    unavailable_list <- matrix(ncol = 2)
    colnames(unavailable_list) <- c("Tax_ID", "Skipped_taxa")
    for (i in 1:ID_n) {
      ID <- r_search$ids[i]
      cat("\nID ", i, ":\t", ID, "\t")
      upload <- entrez_post(db = "taxonomy", 
                            id = ID) # getting the taxonomy with IDs
      fetch_id <- entrez_link(dbfrom = "taxonomy", 
                              db = "nuccore", 
                              web_history = upload) # linking between the nucleotide & taxonomy databases - this is when the inconsistencies come in
      fetch_id2 <- fetch_id$links$taxonomy_nuccore
      if(!is.null(fetch_id2)) {
        all_recs_list[[loop]] <- entrez_fetch(db = "nuccore",
                                              id = fetch_id2, 
                                              rettype = "fasta")
        loop <- loop + 1
      }else{
        cat("no sequence data")
        unavailable_list <- rbind(unavailable_list, 
                                  c(ID, 
                                    entrez_summary(db = "taxonomy", 
                                                   id = ID)$scientificname))
        }
    }
    write(unlist(all_recs_list), file = paste0(path_to_output_dir, "/", search_term, "_available_seqs_temp.fasta"))
    temp <- read.fasta(paste0(path_to_output_dir, "/", search_term, "_available_seqs_temp.fasta"))
    all_fasta <- read.fasta(paste0(path_to_output_dir, "/", search_term, "_available_seqs.fasta"))
    temp_names <- unlist(getAnnot(temp))
    all_names <- unlist(getAnnot(all_fasta))
    select_idx <- which(!temp_names %in% all_names)
    all_fasta2 <- c(all_fasta, temp[select_idx])
    write.fasta(all_fasta2, names = str_remove(unlist(getAnnot(all_fasta2)), ">"), file.out = paste0(path_to_output_dir, "/", search_term, "_available_seqs.fasta"))
    }
  }
  #delete the empty first sequences
  all_fasta <- read.fasta(paste0(path_to_output_dir, "/", search_term, "_available_seqs.fasta"))
  all_fasta <- all_fasta[-1]
  n_seq <- length(all_fasta)
  if(length(r_search$ids) == 0) {cat("\nThis genus has 0 accessions in NCBI taxonomy database\n")}
  cat("\n\n", n_seq, " sequences found for ", paste0(search_term), "\n", sep = "")
  write.fasta(all_fasta2, names = getAnnot(all_fasta2), file.out = paste0(path_to_output_dir, "/", search_term, "_available_seqs.fasta"))
}

# TEST THE FUNCTION WITH A SMALL GENUS
ten_NCBI_seq_search("Diehliomyces")

# GET LIST OF GENERA
setwd("/Users/emilyhodgson/Documents/Autophylo/")
genera_df <- read.csv(file = 'Otideaceae_genera_list.csv')

# RUN FOR YOUR LIST OF GENERA 
for(g in genera_df$Genus) {
  ten_NCBI_seq_search(g)
}
# This can be pretty slow! And can hit 'HTTP failure' a lot :/ If it does, just retry.
# Rather than redoing all genera you can just go from where it stopped