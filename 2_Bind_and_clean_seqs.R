## Emily Hodgson
## Last modified: 3.11.2023

# Setup
package.setup <- function() {
  if (!require("stringr", quietly = TRUE)){
    install.packages("stringr")
  } 
  if (!require("seqinr", quietly = TRUE)){
    install.packages("seqinr")
  }
  library(seqinr)
  library(stringr)
}
package.setup()
path_input <- "/Users/emilyhodgson/Documents/Autophylo/"
files <- list.files(path = paste0(path_input, "Available_Seqs"), pattern = ".fasta")
genera_df <- read.csv(file = paste0(path_input, 'Otideaceae_genera_list.csv'))

# Data check:
# Is there a FASTA file for every taxa in your list?
if(all(genera_df$Genus %in% str_remove(files, "_available_seqs.fasta")) == TRUE) {
  cat("FASTA files are ready")
} else {
  cat("STOP & CHECK, FASTA files missing:\n")
  idx <- genera_df$Genus %in% str_remove(files, "_available_seqs.fasta")
  genera_df$Genus[!idx]
}

# Are there any extra FASTA files for taxa not in your list?
if(all(str_remove(files, "_available_seqs.fasta") %in% genera_df$Genus) == TRUE) {
  cat("FASTA file taxa are all in your list :)")
} else {
  cat("N.B. There are extra FASTA files here not in your list & will not be included downstream:\n")
  idx <- str_remove(files, "_available_seqs.fasta") %in% genera_df$Genus
  cat(paste0("\n", files[!idx]), sep = "\n")
  files <- files[idx] # remove extra FASTA files
}

# Binding files
all_fasta <- list()
for(i in 1:length(files)){
  temp_fasta <- read.fasta(paste0(path_input, "Available_Seqs/", files[i]))
  all_fasta <- c(all_fasta, temp_fasta)
}

length(names(all_fasta)) == length(unique(names(all_fasta))) 
# checking for duplicates
length(names(all_fasta))
length(unique(names(all_fasta)))
# yes there are some duplicated accessions. Remove duplicated ones to give fasta2
all_fasta <- all_fasta[!duplicated(names(all_fasta))]
sum(str_detect(sapply(all_fasta, function(x) attr(x, "Annot")), "genome"))
idx <- str_detect(sapply(all_fasta, function(x) attr(x, "Annot")), "genome")
all_fasta2 <- all_fasta[!idx]

# filter out sequences shorter than 200bp and more than 5000bp to give fasta3
all_fasta3 <- all_fasta2[nchar(all_fasta2) >= 200 & nchar(all_fasta2) <= 5000] 
setwd(path_input)
names(all_fasta3) <- sapply(all_fasta3, function(x) attr(x, "Annot"))
names(all_fasta3) <- str_remove(names(all_fasta3), ">")

# save filtered sequences, specify your own path
write.fasta(sequences = all_fasta3, 
            names = names(all_fasta3), 
            file.out = "./Available_Seqs_filtered.fasta") 
