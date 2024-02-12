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
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  if (!require("Biostrings", quietly = TRUE)){
    BiocManager::install("Biostrings")
  }
  library(seqinr)
  library(stringr)
  library("Biostrings")
}
package.setup()
path_input <- "/Users/emilyhodgson/Documents/Autophylo/"
path_baits <- "/Users/emilyhodgson/Documents/Autophylo/Emily_baits/"
path_outputs <- "/Users/emilyhodgson/Documents/Autophylo/Loci_sorted_sequences/"

# Bait-sorting and pattern-searching script

setwd(path_baits)

# SETTING THE BAIT FILES
ITS <- read.fasta("./ITS_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq
LSU <- read.fasta("./LSU_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq
SSU <- read.fasta("./SSU_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq
RPB2 <- read.fasta("./RPB2_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq
EF1a <- read.fasta("./EF1a_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq
Btub <- read.fasta("./Btub_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq
RPB1 <- read.fasta("./RPB1_bait.fasta", seqtype = "DNA", as.string = TRUE) # read the bait seq

# can combine those of poorer resolution, like SSU and ITS as one, or RPB1 and RPB2 
# ??? ^
ITS <- BString(ITS[[1]][[1]])
LSU <- BString(LSU[[1]][[1]])
SSU <- BString(SSU[[1]][[1]])
RPB2 <- BString(RPB2[[1]][[1]])
EF1a <- BString(EF1a[[1]][[1]])
Btub <- BString(Btub[[1]][[1]])
RPB1 <- BString(RPB1[[1]][[1]])
Bait_list <- list(ITS = ITS, LSU = LSU, SSU = SSU, RPB2 = RPB2, EF1a = EF1a, Btub = Btub, RPB1 = RPB1)
# storing as BString objects. Useful for sequence manipulation.

setwd(path_input)

all_recs_list_uncategorized <- read.fasta(file = "./Available_Seqs_filtered.fasta", as.string = T)
ITS_master <- list()
LSU_master <- list()
SSU_master <- list()
RPB2_master <- list()
EF1a_master <- list()
Btub_master <- list()
RPB1_master <- list()
uncategorised_master <- list()

loop <- 1
loop_v <- vector()

for(i in seq_along(all_recs_list_uncategorized)) {
  cat("1st for loop\n")
  sink("./temp.fasta")
  cat(paste0(getAnnot(all_recs_list_uncategorized[[i]]), 
             "\n", 
             as.character(all_recs_list_uncategorized[[i]]), 
             "\n")) # print sequence annotations for our FASTA file and create a character vector? does sink() send it to the new temp FASTA file?
  sink() 
  seqs <- read.fasta("./temp.fasta", seqtype = "DNA", as.string = TRUE)
  loci_v <- vector()
  names_v <- vector()
  cat(loop, "\t", i, "\n")
  loop_v[loop] <- loop
  for(j in seq_along(seqs)) {
    pair_wise <- sapply(Bait_list, function(x) { # applies the following function over the Bait_list
      align <- pairwiseAlignment(x, 
                                 BString(seqs[[j]][[1]]), 
                                 gapOpening=0, 
                                 gapExtension=-5,
                                 type = "local")
      score <- score(align)
      pid <- pid(align, type="PID1")
      cbind(pid, score)
    }, simplify = FALSE, USE.NAMES = TRUE)
    pair_wise_matrix <- t(matrix(unlist(pair_wise), 
                                 ncol = 2, 
                                 byrow = TRUE))
    rownames(pair_wise_matrix) <- c("pid", "score")
    colnames(pair_wise_matrix) <- names(Bait_list)
    # Find the index of the alignment with the highest score 
    best_index <- which.max(pair_wise_matrix[2,])
    # Use the best alignment to assign the sequence to a locus
    
    loci_v[j] <- colnames(pair_wise_matrix)[best_index]
    # Use the best alignment to assign the full name of the sequence
    #names_v[j] <- gsub(".*\\[|\\].*", "", names(seqs[j]))
  }
  names(seqs) <- loci_v
  for(k in seq_along(seqs)) {
    if(str_detect(names(seqs)[k], pattern = "LSU")) {
      if(str_detect(getAnnot(seqs[[k]]), pattern = c("large subunit ribosomal RNA gene|28S|LSU|25S"))){ 
        #next
        #} else {
        LSU_master <- c(LSU_master, seqs[k])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs)[k], pattern = "ITS")){
      if(str_detect(getAnnot(seqs[[k]]), pattern = c("internal transcribed spacer|ITS|5.8S|internal"))){
        #next
        #} else {
        ITS_master <- c(ITS_master, seqs[k])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs)[k], pattern = "SSU")){
      if(str_detect(getAnnot(seqs[[k]]), pattern = c("small subunit ribosomal RNA gene|18S"))){
        #next
        #} else {
        SSU_master <- c(SSU_master, seqs[k])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs)[k], pattern = "RPB2")){
      if(str_detect(getAnnot(seqs[[k]]), pattern = c("RNA polymerase II|RPB2|rpb2"))){
        #next
        #} else {
        RPB2_master <- c(RPB2_master, seqs[k])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs)[k], pattern = "EF1a")){
      if(str_detect(getAnnot(seqs[[k]]), pattern = c("elongation|EF1a|TEF1a|TEF|tef"))){
        #next
        #} else {
        EF1a_master <- c(EF1a_master, seqs[k])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs)[k], pattern = "Btub")){
      if(str_detect(getAnnot(seqs[[k]]), pattern = c("tubulin|beta|Tub"))){
        #next
        #} else {
        Btub_master <- c(Btub_master, seqs[k])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs)[k], pattern = "RPB1")){
      if(str_detect(getAnnot(seqs[[k]]), pattern =c("RNA polymerase II largest subunit|RPB1|rpb1"))){
        #next
        #} else {
        RPB1_master <- c(RPB1_master, seqs[k])
        loop <- loop + 1
      }}
    if(!str_detect(names(seqs)[k],pattern = "LSU|ITS|SSU|RPB2|EF1a|Btub|RPB1")){
      uncategorised_master <- c(uncategorised_master, seqs[k])
      loop <- loop + 1
    }
  }
}

setwd(path_outputs)

# ORGANISED SEQUENCES GO INTO OUTPUT FOLDER
write.fasta(sequences = lapply(LSU_master, toupper), names = gsub(">", "", getAnnot(LSU_master)), file.out = "./LSU_Res.fasta")
write.fasta(sequences = lapply(ITS_master, toupper), names = gsub(">", "", getAnnot(ITS_master)), file.out = "./ITS_Res.fasta")
write.fasta(sequences = lapply(SSU_master, toupper), names = gsub(">", "", getAnnot(SSU_master)), file.out = "./SSU_Res.fasta")
write.fasta(sequences = lapply(RPB2_master, toupper), names = gsub(">", "", getAnnot(RPB2_master)), file.out = "./RPB2_Res.fasta")
write.fasta(sequences = lapply(EF1a_master, toupper), names = gsub(">", "", getAnnot(EF1a_master)), file.out = "./EF1a_Res.fasta")
write.fasta(sequences = lapply(Btub_master, toupper), names = gsub(">", "", getAnnot(Btub_master)), file.out = "./Btub_Res.fasta")
write.fasta(sequences = lapply(RPB1_master, toupper), names = gsub(">", "", getAnnot(RPB1_master)), file.out = "./RPB1_Res.fasta")

# SUMMARY
ITS_length <- length(ITS_master)
LSU_length <- length(LSU_master)
SSU_length <- length(SSU_master)
RPB2_length <- length(RPB2_master)
EF1a_length <- length(EF1a_master)
Btub_length <- length(Btub_master)
RPB1_length <- length(RPB1_master)
uncategorised_length <- length(uncategorised_master)
cat("Total categorised sequences = ", sum(ITS_length + LSU_length + SSU_length + RPB2_length + EF1a_length + Btub_length + RPB1_length), "\nTotal uncategorised sequences = ", uncategorised_length)
