# Tests of scripts can run in a separate file
###
# setting up R
source("../RSetup.R")
package.setup(workingdir = "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/")
#######
#### 1_Download_seqs.R aka single_attempt_NCBI_seq_search
#######
source("../1_Download_seqs.R")
# TEST THE FUNCTION WITH A SMALL GENUS
single_attempt_NCBI_seq_search(search_term = "Glaziella", api_key = NULL, path_to_output_dir = "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/")
list.files(".")
Glaziella_fasta <- read.fasta("./Glaziella_available_seqs.fasta")
Glaziella_fasta_origin <- read.fasta("./Glaziella_available_seqs_original.fasta")
identical(Glaziella_fasta, Glaziella_fasta_origin)

##### I haven't run the lines below (Torda)
# GET LIST OF GENERA
setwd("/Users/emilyhodgson/Documents/Autophylo/")
genera_df <- read.csv(file = 'Otideaceae_genera_list.csv')

# RUN FOR YOUR LIST OF GENERA 
for(g in genera_df$Genus) {
  ten_NCBI_seq_search(g)
}