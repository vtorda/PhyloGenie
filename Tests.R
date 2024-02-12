# Tests of scripts can run in a separate file
###
# setting up R
source("./RSetup.R")
package.setup(workingdir = "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/")
#######
#### 1_Download_seqs.R aka single_attempt_NCBI_seq_search
#######
source("../1_Download_seqs.R")
# TEST THE FUNCTION WITH A SMALL GENUS
NCBI_seq_fetch(search_term = "Otidea",
               api_key = "ce946cc5385e86927f230c4aea4a5f68ac08",
               path_to_output_dir = "/Users/varga/OneDrive/Documents/GitHub/PhyloGenie/TestFolder/")
list.files(".")
Glaziella_fasta <- read.fasta("./Glaziella_available_seqs.fasta")


##### I haven't run the lines below (Torda)
# GET LIST OF GENERA
setwd("/Users/emilyhodgson/Documents/Autophylo/")
genera_df <- read.csv(file = 'Otideaceae_genera_list.csv')

# RUN FOR YOUR LIST OF GENERA 
for(g in genera_df$Genus) {
  ten_NCBI_seq_search(g)
}