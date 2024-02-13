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
# package.setup(workingdir = "")
# search_term <- "Otidea"
# api_key <- ""
# path_to_output_dir <- ""
# minlength <- 150
# maxlength <- 5000
# TechFilter <- c("wgs", "targeted", "tsa", "est")
#entrez_db_searchable("taxonomy")
# ADD FUNCTION TO ENVIRONMENT
# This function outputs a FASTA file with all the available sequences in NCBI for your search term.
# Set 'path_to_output_dir' (above) as the place you want your file to go
NCBI_seq_fetch <- function(search_term, api_key = NULL, path_to_output_dir, minlength = 100, maxlength = 4000, TechFilter = NULL, chunk = 200) {
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
    cat("\nGathering sequence information and prefiltering the dataset have been started\n")
    link_history <- entrez_link(dbfrom = "taxonomy", 
                                db = "nuccore",
                                web_history = r_search$web_history,
                                cmd = "neighbor_history")
    # get the info of all linked nuccore records
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
    # filter out sequences
    idx <- (as.numeric(meta_df$Slen) <= maxlength & as.numeric(meta_df$Slen) >= minlength)
    meta_keep <- meta_df[idx, ]
    if(!is.null(TechFilter)){
      idx <- meta_keep$Tech %in% TechFilter
      meta_keep <- meta_keep[!idx,]
    }
    meta_dropped <- meta_df[!meta_df$Caption %in% meta_keep$Caption,]
    if(nrow(meta_dropped) > 0){
      write.table(x = meta_dropped, file = paste0(path_to_output_dir, search_term, ".metadata.skipped.tsv"),
                  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      cat(nrow(meta_dropped), " number of sequences were filtered out. The metadata of these sequences can be found in the \n",
          paste0(path_to_output_dir, search_term, ".metadata.skipped.tsv"),"\n File\n")
    }else{
      cat("No sequences were filtered out")
    }
    # download sequences in batches and fetch extra meta info
    
    Caption <- meta_keep$Caption
    max_seq <- length(Caption)
    seq_start <- seq(1,max_seq,chunk)
    batch_n <- length(seq_start)
    gb_info_list <- list()
    cat("\nDownloading ", max_seq, " number of sequences and their metadata has been started\n")
    for(j in 1:batch_n){
      if(j != batch_n){
        upload <- entrez_post(db = "nuccore",
                              id = Caption[seq_start[j]:(seq_start[j]+(chunk-1))])
        gb_info_temp <- rentrez::entrez_fetch(db="nuccore",
                                              web_history = upload,
                                              rettype="gb",
                                              retmode = "text")
        gb_info_list[[j]] <- gb_info_temp
        recs <- rentrez::entrez_fetch(db="nuccore",
                                      web_history=upload,
                                      rettype="fasta") # if restez package is loaded need to define rentrez package here
        cat(recs, file=paste0(path_to_output_dir, search_term, ".fasta"), append=TRUE)
        cat("\n",seq_start[j]+chunk-1, "specimens info has been fetched\n")
      }
      if(j == batch_n){
        final_add <- max_seq - seq_start[j]
        upload <- entrez_post(db = "nuccore",
                              id = Caption[seq_start[j]:(seq_start[j]+final_add)])
        gb_info_temp <- rentrez::entrez_fetch(db="nuccore",
                                              web_history = upload,
                                              rettype="gb",
                                              retmode = "text")
        gb_info_list[[j]] <- gb_info_temp
        recs <- rentrez::entrez_fetch(db="nuccore",
                                      web_history=upload,
                                      rettype="fasta") # if restez package is loaded need to define rentrez package here
        cat(recs, file=paste0(path_to_output_dir, search_term, ".fasta"), append=TRUE)
        cat("\n", max_seq, "specimens info has been fetched\n")
      }
    }
    
    gb_info <- paste0(unlist(gb_info_list), collapse = "")
    gb_info_list2 <- strsplit(gb_info, split = "(?<=LOCUS)", perl = TRUE) # split along LOCUS, but keep LOCUS in result too
    
    # I need to refer to the 1st element of the list where I have all the sequences. I need to skip the first among the sequences which is just one word: LOCUS
    feature <- lapply(2:length(gb_info_list2[[1]]), function(x) gb_extract(record = gb_info_list2[[1]][x], what = "features")[[1]])
    acc <- lapply(2:length(gb_info_list2[[1]]), function(x) gb_extract(record = gb_info_list2[[1]][x], what = "accession")[[1]])
    extra_info <- unique(unlist(lapply(feature, names)))
    all_info <- matrix(ncol = length(extra_info))
    for(i in 1:length(feature)){
      n <- unlist(feature[[i]])
      names(n) <- names(feature[[i]])
      n2 <- n[match(extra_info, names(n))] # with this I expand the vector to the sice of extra_info, NA will be for no matches
      # names(n2) <- extra_info
      all_info <- rbind(all_info, n2)
    }
    all_info <- all_info[-1,]
    colnames(all_info) <- extra_info
    all_info2 <- as.data.frame(all_info)
    all_info2$accession_no <- unlist(acc)
    if(nrow(meta_keep) == nrow(all_info2)){
      df2 <- cbind(meta_keep, all_info2)
    }else{
      cat("You couldn't download all the sequences, try to increase the chunk size")
    }

    write.table(x = df2, file = paste0(path_to_output_dir, search_term, ".metadata.tsv"),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

  }else{
    warning("No taxon ID was found")
  }
}