# Categorise DNA sequences according to loci ------------------------------
# read in the selected bait files, opprtunity to modify them
ITS <- read.fasta("./ITS_bait.fasta", seqtype = "DNA", as.string = TRUE) 
LSU <- read.fasta("./LSU_bait.fasta", seqtype = "DNA", as.string = TRUE) 
SSU <- read.fasta("./SSU_bait.fasta", seqtype = "DNA", as.string = TRUE) 
RPB2 <- read.fasta("./RPB2_bait.fasta", seqtype = "DNA", as.string = TRUE) 
EF1a <- read.fasta("./EF1a_bait.fasta", seqtype = "DNA", as.string = TRUE) 
Btub <- read.fasta("./Btub_bait.fasta", seqtype = "DNA", as.string = TRUE) 
RPB1 <- read.fasta("./RPB1_bait.fasta", seqtype = "DNA", as.string = TRUE) 

# create general container for each loci 
ITS <- BString(ITS[[1]][[1]])
LSU <- BString(LSU[[1]][[1]])
SSU <- BString(SSU[[1]][[1]])
RPB2 <- BString(RPB2[[1]][[1]])
EF1a <- BString(EF1a[[1]][[1]])
Btub <- BString(Btub[[1]][[1]])
RPB1 <- BString(RPB1[[1]][[1]])
Bait_list <- list(ITS = ITS, LSU = LSU, SSU = SSU, RPB2 = RPB2, EF1a = EF1a, Btub = Btub, RPB1 = RPB1)

# read in saved fasta file containing all DNA sequences
all_recs_list_uncategorized <- read.fasta(file = "./Available_Rlist_filtered.fasta", as.string = T)

# Prepare list() for the loops later
ITS_master <- list()
LSU_master <- list()
SSU_master <- list()
RPB2_master <- list()
EF1a_master <- list()
Btub_master <- list()
RPB1_master <- list()

#  a loop to categorise sequences including a search for their exact names 
loop <- 1
loop_v <- vector()
for(i in seq_along(all_recs_list_uncategorized)){
  sink("./temp.fasta")
  cat(paste0(getAnnot(all_recs_list_uncategorized[[i]]), "\n", as.character(all_recs_list_uncategorized[[i]]), "\n"))   # renew per loop
  sink()
  seqs <- read.fasta("./temp.fasta", seqtype = "DNA", as.string = TRUE) # create temporary fasta file to reduce error
  loci_v <- vector()
  # names_v <- vector()
  cat(loop, "\t", i, "\n")
  loop_v[loop] <- loop
  
  for(j in seq_along(seqs)){
    pair_wise <- sapply(Bait_list, function(x){
      align <- pairwiseAlignment(x, BString(seqs[[j]][[1]]), gapOpening=0, gapExtension=-5,type = "local")
      score <- score(align)
      pid <- pid(align, type="PID1")
      cbind(pid, score)
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    pair_wise_matrix <- t(matrix(unlist(pair_wise), ncol = 2, byrow = TRUE))
    rownames(pair_wise_matrix) <- c("pid", "score")
    colnames(pair_wise_matrix) <- names(Bait_list)
    
    # Find the index of the alignment with the highest score and highest pid
    best_index <- which.max(pair_wise_matrix[1,] + pair_wise_matrix[2,])
    # Use the best alignment to assign the sequence to a locus
    loci_v[j] <- names(pair_wise)[best_index]
    # Use the best alignment to assign the full name of the sequence
    # names_v[j] <- gsub(".*\\[|\\].*", "", names(seqs[j]))
  }
  
  names(seqs[[j]]) <- str_c(loci_v[j] ,sep = "_")
  for(j in seq_along(seqs)){
    if(str_detect(names(seqs[[j]]), pattern = "LSU")){
      if(str_detect(getAnnot(seqs[[j]]), pattern = c("large subunit ribosomal RNA gene|28S|LSU|25S"))){
        LSU_master <- c(LSU_master, seqs[j])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs[[j]]), pattern = "ITS")){
      if(str_detect(getAnnot(seqs[[j]]), pattern = c("internal transcribed spacer|ITS|5.8S|internal"))){
        ITS_master <- c(ITS_master, seqs[j])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs[[j]]), pattern = "SSU")){
      if(str_detect(getAnnot(seqs[[j]]), pattern = c("small subunit ribosomal RNA gene|18S"))){
        SSU_master <- c(SSU_master, seqs[j])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs[[j]]), pattern = "RPB2")){
      if(str_detect(getAnnot(seqs[[j]]), pattern = c("RNA polymerase II|RPB2|rpb2"))){
        RPB2_master <- c(RPB2_master, seqs[j])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs[[j]]), pattern = "EF1a")){
      if(str_detect(getAnnot(seqs[[j]]), pattern = c("elongation|EF1a|TEF1a|TEF|tef"))){
        EF1a_master <- c(EF1a_master, seqs[j])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs[[j]]), pattern = "Btub")){
      if(str_detect(getAnnot(seqs[[j]]), pattern = c("tubulin|beta|Tub"))){
        Btub_master <- c(Btub_master, seqs[j])
        loop <- loop + 1
      }}
    if(str_detect(names(seqs[[j]]), pattern = "RPB1")){
      if(str_detect(getAnnot(seqs[[j]]), pattern =c("RNA polymerase II largest subunit|RPB1|rpb1"))){
        RPB1_master <- c(RPB1_master, seqs[j])
        loop <- loop + 1
      }}
    if(!str_detect(names(seqs[[j]]),pattern = "LSU|ITS|SSU|RPB2|EF1a|Btub|RPB1")){
      uncategorised_master <- c(uncategorised_master, seqs[j])
      loop <- loop + 1
    }
  }}

# write out fasta files of each loci
write.fasta(sequences = lapply(LSU_master, toupper), names = gsub(">", "", getAnnot(LSU_master)), file.out = "./LSU_Res.fasta")
write.fasta(sequences = lapply(ITS_master, toupper), names = gsub(">", "", getAnnot(ITS_master)), file.out = "./ITS_Res.fasta")
write.fasta(sequences = lapply(SSU_master, toupper), names = gsub(">", "", getAnnot(SSU_master)), file.out = "./SSU_Res.fasta")
write.fasta(sequences = lapply(RPB2_master, toupper), names = gsub(">", "", getAnnot(RPB2_master)), file.out = "./RPB2_Res.fasta")
write.fasta(sequences = lapply(EF1a_master, toupper), names = gsub(">", "", getAnnot(EF1a_master)), file.out = "./EF1a_Res.fasta")
write.fasta(sequences = lapply(Btub_master, toupper), names = gsub(">", "", getAnnot(Btub_master)), file.out = "./Btub_Res.fasta")
write.fasta(sequences = lapply(RPB1_master, toupper), names = gsub(">", "", getAnnot(RPB1_master)), file.out = "./RPB1_Res.fasta")

# remove duplicates
all_seq_list <- list(read.fasta("LSU_Res.fasta", seqtype = "DNA", as.string = TRUE),
                     read.fasta("ITS_Res.fasta", seqtype = "DNA", as.string = TRUE), 
                     read.fasta("SSU_Res.fasta", seqtype = "DNA", as.string = TRUE), 
                     read.fasta("RPB1_Res.fasta", seqtype = "DNA", as.string = TRUE), 
                     read.fasta("RPB2_Res.fasta", seqtype = "DNA", as.string = TRUE), 
                     read.fasta("EF1a_Res.fasta", seqtype = "DNA", as.string = TRUE), 
                     read.fasta("Btub_Res.fasta", seqtype = "DNA", as.string = TRUE))
new_names <- c("LSU", "ITS", "SSU", "RPB1", "RPB2", "EF1a", "Btub")
for (i in seq_along(all_seq_list)){
  all_seq_list[[i]] <- unique(all_seq_list[[i]])
  # Write out to new fasta file
  write.fasta(sequences = lapply(all_seq_list[[i]], toupper),
              names = gsub(">", "", getAnnot(all_seq_list[[i]])),
              file.out = paste0("new_sequences_", new_names[i], ".fasta"))
}
