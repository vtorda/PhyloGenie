## Torda based on Emily Hodgson's script
## Last modified: 12.02.2024

# A script to download extra information uploaded on GenBank associated with the sequences we want to download. 

# Eventually merge with the script that downloads seqs?

# path_to_output_dir <- "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/"
# meta_df <- "Glaziella.metadata.tsv"
extr_info_fetch <- function(path_to_output_dir, meta_df){
df <- readr::read_tsv(paste0(path_to_output_dir, meta_df))
# create a web history object in batches
max_seq <- nrow(df)
chunk <- 300
seq_start <- seq(1,max_seq,chunk)
batch_n <- length(seq_start)
gb_info_list <- list()
for(j in 1:batch_n){
  if(j != batch_n){
  upload <- entrez_post(db = "nuccore",
                        id = df$Caption[seq_start[j]:(seq_start[j]+(chunk-1))])
  gb_info_temp <- rentrez::entrez_fetch(db="nuccore",
                                        web_history = upload,
                                        rettype="gb",
                                        retmode = "text")
  gb_info_list[[j]] <- gb_info_temp
  cat("\n",seq_start[j]+chunk-1, "specimens info has been fatched\n")
  }
  if(j == batch_n){
  final_add <- max_seq - seq_start[j]
  upload <- entrez_post(db = "nuccore",
                          id = df$Caption[seq_start[j]:(seq_start[j]+final_add)])
  gb_info_temp <- rentrez::entrez_fetch(db="nuccore",
                                        web_history = upload,
                                        rettype="gb",
                                        retmode = "text")
  gb_info_list[[j]] <- gb_info_temp
  cat("\n", max_seq, "specimens info has been fatched\n")
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
df2 <- cbind(df, all_info2)
out <- str_remove(meta_df, ".metadata.tsv")
write.table(x = df2, file = paste0(path_to_output_dir, out, ".metadata_long.tsv"),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}
