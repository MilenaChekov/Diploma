if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('Rfast', force = TRUE)
library('Rfast')
BiocManager::install("protr")
library("protr")
BiocManager::install("doParallel")
library("doParallel")
BiocManager::install("fastmatch")
library("fastmatch")
BiocManager::install("ShortRead")
library("ShortRead")

aLign = function(seq_1, seq_2, mtx = protr::AABLOSUM62, aas = colnames(protr::AABLOSUM62)) {
  sum(mtx[cbind(fmatch(unlist(strsplit(seq_1, ""), use.names = F), aas), fmatch(unlist(strsplit(seq_2, ""), use.names = F), aas))])
}

peptides = read.delim("EXAMPLE.csv", sep = ",", stringsAsFactors = F) #the file was downloaded from https://adaptivepublic.blob.core.windows.net/publishedproject-supplements/covid-2020/ImmuneCODE-MIRA-Release002.zip. The file can be downloaded from Mendeley Data (dx.doi.org/10.17632/63xr979845.1).
peptide_seqs = unique(unlist(strsplit(peptides$index, ",")))
peptide_seqs = peptide_seqs[nchar(peptide_seqs) == 9]
proteome = unlist(readFASTA("uniprot-proteome_UP000005640+reviewed_yes.fasta")) #The file can be downloaded from Mendeley Data (dx.doi.org/10.17632/63xr979845.1).
proteome = proteome[nchar(proteome) >= 9]
all_nonamers = sapply(proteome, FUN = function(x) substring(x, 1:(nchar(x) - 8), 9 : nchar(x)))
all_nonamers = unique(unlist(all_nonamers, use.names = F))
c1 = makeCluster(7) 
clusterExport(c1, c("aLign", "AABLOSUM62", "fmatch"))
self_all_nonamers = parSapply(c1, all_nonamers, FUN = function(x) aLign(seq_1 = x, seq_2 = x))
self_peptide_seqs = parSapply(c1, peptide_seqs, FUN = function(x) aLign(seq_1 = x, seq_2 = x))
aas = colnames(AABLOSUM62)
mtx = AABLOSUM62
clusterExport(c1, c("mtx", "aas"))
similarity = lapply(peptide_seqs, FUN = function(x) {
  indices_x = fmatch(substring(x,1:9,1:9), aas)
  clusterExport(c1, c("indices_x"), envir = environment())
  seq_to_all = parSapply(c1, all_nonamers, FUN = function(y) sum(mtx[cbind(indices_x, fmatch(substring(y,1:9,1:9), aas))]))
  similarities = seq_to_all/sqrt(self_peptide_seqs[x]*self_all_nonamers)
  similarities = sort(similarities, decreasing = T)
  return(similarities[1:10])
}) #It is a long-time running


write.csv(similarity, "res_proteome_EXAMPLE.csv")

