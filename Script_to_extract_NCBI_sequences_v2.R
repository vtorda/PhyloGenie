# Install packages i.e. stringr under Biostrings
install.packages("rentrez")
install.packages("seqinr")
install.packages("Biostrings")

# Check installed packages 
library(rentrez)
library(seqinr)
library("Biostrings")
library(stringr)

# Obtain ids from taxonomy database of NCBI and their DNA sequences -------
# retrieve all ids under Pezizaceae family
# term (search term), use_history = T to return a web_history object for use in later calls to the NCBI
r_search <- entrez_search(db="taxonomy", term="Pezizaceae[subtree]",retmax = 999,use_history = T)  #or [SBTR]

# loop through every id in Pezizaceae to retrieve info, may need to truncate into batches as it may crash
loop <- 1
loop_v <- vector()
all_recs_list <- list()
unavailable_list <- matrix(ncol = 2)
colnames(unavailable_list) <- c("Tax_ID", "Skipped_taxa")

for (i in r_search$ids[1:length(r_search$ids)]){  
  cat(loop, "\t", i, "\n")
  loop_v[loop] <- loop 
  upload <- entrez_post(db="taxonomy", id=i)
  fetch_id <- entrez_link(dbfrom = "taxonomy", db = "nuccore", web_history = upload)
  fetch_id2 <- fetch_id$links$taxonomy_nuccore
  if(!is.null(fetch_id2)){
    all_recs_list[[loop]] <- entrez_fetch(db = "nuccore",id = fetch_id2, rettype = "fasta")
    loop <- loop + 1
  }
  if (is.null(fetch_id2)){
    unavailable_list <- rbind(unavailable_list, c(i, entrez_summary(db = "taxonomy", id = i)$scientificname))}
}

# write out unavailable list in csv (those id without info) and all_recs_list in fasta file (ids and their DNA sequences)
write.csv(unavailable_list, file = "./Unavailable_Rlist.csv")
write.fasta(sequences = all_recs_list, names = names(all_recs_list), file.out = "./Available_Rlist.fasta")

# filter out long whole genome shotgun sequences and short reads - remove less than 200bp and more than 5000bp
# also filter out those with the word "whole genome shotgun sequence"
all_recs_list_non_filtered <- readDNAStringSet(file = "./Available_Rlist_all.fasta")
all_recs_list_filtered <- all_recs_list_non_filtered[nchar(all_recs_list_non_filtered) >= 200 & nchar(all_recs_list_non_filtered) <= 5000]
all_recs_list_filtered <- all_recs_list_filtered[!grepl("whole genome shotgun sequence", names(all_recs_list_filtered), ignore.case = TRUE)]
writeXStringSet(all_recs_list_filtered, "./Available_Rlist_filtered.fasta")
