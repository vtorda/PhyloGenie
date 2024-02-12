## Emily Hodgson
## Last modified: 3.11.2023

# This script downloads all available sequence data for a given list of genera. (not including metadata, that's another script, but these scripts could maybe be combined to save time!)

# I found searching by genus easier & more effective than searching by species.

###
## ToDO
## Write a report file out

# Inputs required:
## Personal API key (increases your e-utils limit to 10 requests/second)
## List of genera (I used csv file)
# source("RSetup.R")
# package.setup(workingdir = "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/")
# search_term <- "Otidea"
# api_key <- "NULL"
# path_to_output_dir <- "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/"
# # minlength <- 150
# # maxlength <- 4000

#entrez_db_searchable("taxonomy")
# ADD FUNCTION TO ENVIRONMENT
# This function outputs a FASTA file with all the available sequences in NCBI for your search term.
# Set 'path_to_output_dir' (above) as the place you want your file to go
NCBI_seq_fetch <- function(search_term, api_key = NULL, path_to_output_dir) {
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
    # linking all taxid with nuccore data 
    link_history <- entrez_link(dbfrom = "taxonomy", 
                                db = "nuccore",
                                web_history = r_search$web_history,
                                cmd = "neighbor_history")
    # get the info of all linked noccure records
    summary_result <- entrez_summary("nuccore", 
                                     web_history = link_history$web_histories$taxonomy_nuccore, retmode = "xml")
    meta <- lapply(summary_result, function(x) x[1:27]) # it seems that the first 27 element always the same info
    #changing NULL to NA https://stackoverflow.com/questions/22870198/is-there-a-more-efficient-way-to-replace-null-with-na-in-a-list
    nullToNA <- function(x) {
      x[sapply(x, is.null)] <- NA
      return(x)
    }
    meta2 <- lapply(meta, nullToNA)
    meta3 <- lapply(meta2, unlist) 
    meta_df <- as.data.frame(do.call(rbind, meta3))
    # download sequences in batches using web history data
    max_seq <- nrow(meta_df)
    seq_start <- seq(1,max_seq,50)
    batch_n <- length(seq_start)
    for(j in 1:batch_n){
      recs <- rentrez::entrez_fetch(db="nuccore", web_history=link_history$web_histories$taxonomy_nuccore,
                           rettype="fasta", retmax=50, retstart=seq_start[j]-1) # if restez package is loaded need to define rentrez package here
      cat(recs, file=paste0(path_to_output_dir, search_term, ".fasta"), append=TRUE)
      if(j < batch_n){
        cat("\n",seq_start[j]+49, "sequences downloaded\n")
      }else{
        cat("\n", max_seq, "sequences downloaded\n")
      }
    }
    # write out metadata
    write.table(x = meta_df, file = paste0(path_to_output_dir, search_term, ".metadata.tsv"),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  }else{
    warning("No taxon ID was found")
  }
}
    

