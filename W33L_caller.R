#--------------------------------------------------------------------#
#
## R script to call W33L from Sanger sequencing of VH186.2 locus S####
#
# CC0 1.0 Universal
#
# Authors:
# Ed Carr, edward.carr at babraham.ac.uk
# Danika Hill, danika.hill at babraham.ac.uk
# Alyssa Silva-Cayetano, Alyssa.Silva-Cayetano at babraham.ac.uk
#
# Linterman lab, Babraham Institute, Cambridge, UK
# 
# Released June 2022
#
#--------------------------------------------------------------------#
## Load libraries ####
#--------------------------------------------------------------------#
library(sangerseqR)
library(Biostrings)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)

#--------------------------------------------------------------------#
## Set parameters #####
#--------------------------------------------------------------------#
indir <- "/path/to/dir/with/sequencing/files"
outdir <- "/path/for/reports/"
dir.create(outdir)
ratio <- 0.5
pwd <- getwd()
## W33L parameters ##
max.mismatch <- 2
with.indels <- T
#--------------------------------------------------------------------#
## Import the sequence traces ####
#--------------------------------------------------------------------#
fls <- list.files(path = indir, recursive = T, full.names = T, pattern = ".ab1$")
big.data <- lapply(fls, function (x) { readsangerseq(filename = x) })
# #View first one:
chromatogram(big.data[[1]], width = 200, height = 2, trim5 = 50, trim3 = 50, showcalls = "both")

#--------------------------------------------------------------------#
## Call base calls (any secondary calls < ratio are ignored.) ####
#--------------------------------------------------------------------#
big.data.calls <- lapply(big.data, function(x) {makeBaseCalls(x, ratio = ratio)})

# Now we can extract the primary sequence data:
big.DNAstr <- lapply(big.data.calls, function(x) {x@primarySeq})
big.DNAstrset <- lapply(big.DNAstr, function(x) {DNAStringSet(x)})

save(big.DNAstrset, file = "big.DNAstrset.RData")
#--------------------------------------------------------------------#
## Summary statistics ####
#--------------------------------------------------------------------#
## determine length of each string
seq.length <- function(x){
  nchar(x) #determine the length of the 'string'
  }
how.long <- unlist(lapply(big.DNAstr, seq.length))

#determine summary statistics
summary(how.long)

#some PCR fragments are long, some are very short. 
hist(how.long, breaks = 100)
which(how.long <300) #two are less than 300 - should be removed? 

#remove short sequences
big.DNAstrset <- big.DNAstrset[-which(how.long <300)]
seq.names <- basename(fls[-which(how.long <300)]) #update file names

# what about sequence quality? 
N.length <- function(x){
  str_count(x, pattern = "N") #determine the length of the 'string'
}
how.many.N <- unlist(lapply(big.DNAstr, N.length))

summary(how.many.N)
which(how.many.N >1) ## consider removing 
big.DNAstrset[which(how.many.N >1)] # are the N's at the beginning? this code outputs the start and end of the sequence

#remove any sequences with an "N" called
big.DNAstrset <- big.DNAstrset[-(which(how.many.N >1))]
seq.names <- seq.names[-(which(how.many.N >1))] #update file names

#--------------------------------------------------------------------#
## W33L mutation calling ####
#--------------------------------------------------------------------#
################################
# W33L mutation identification #
# This is in CDR1              #
# TGG is W (germline)          #
# TTG is L                     #
################################
#  Preceding sequence:
# ACCAGCTACT
#  Postceding sequence:
# GATGCACTGG
# Raw sequence needs to be reverse complemented.
####

W33L.caller <- function(DNAstrset,
                        dna_sequence = "ACCAGCTACTNNATGCACTGG",
                        max.mismatch = max.mismatch,
                        with.indels = with.indels){
  m1 <- matchPattern(dna_sequence, 
                     reverseComplement(DNAstrset[[1]]), 
                     max.mismatch=max.mismatch, 
                     min.mismatch = 0, fixed = FALSE, with.indels = with.indels)
  
  W33L <- ifelse(length(m1) == 0, "Locus_not_found",
                 ifelse(as.character(substr(m1[[1]], 10,12)) == "TGG", "W", 
                        ifelse(as.character(substr(m1[[1]], 10,12)) == c("TTG", "TTA"), "L", "other")))
  
  if(length(m1) == 1) {a <- cbind("nt_sequence" = as.character(m1[[1]]), "W33L_call" = W33L)}
  if(length(m1) == 0) {a <- cbind("nt_sequence" = "locus_not_found", "W32_call" = "locus_not_found")}

  return(a)
}

W33L <- lapply(big.DNAstrset, W33L.caller)

W33L <- do.call(rbind, W33L) # makes into a dataframe/matrix

rownames(W33L) <- seq.names # apply original sequencing file information



#--------------------------------------------------------------------#
## Tabulate the results ####
#--------------------------------------------------------------------#
######################################
###                                 ##
### adding data to W33L table       ##
###                                 ##
######################################

W33L <- as.data.frame(W33L, stringsAsFactors = F)
## this substr selection will depend on the filenames used in project ##
W33L$short.name<-substr(rownames(W33L),11,23)
W33L$mouse.ID <- substr(W33L$short.name,11,13)
W33L$genotype <- substr(W33L$mouse.ID,1,2)

#Split into WT and KO
W33L.KO <- as.data.frame(W33L[grepl(rownames(W33L), pattern = "KO"),], stringsAsFactors = F)
W33L.WT <- as.data.frame(W33L[grepl(rownames(W33L), pattern = "WT"),], stringsAsFactors = F)

# outputs summary statistics (means etc)
summary(factor(W33L.KO[,"W33L_call"]))
summary(factor(W33L.WT[,"W33L_call"]))

#show proportion of "L" out of all rows other than "locus not found" 
summary(factor(W33L.KO[,"W33L_call"]))[1]/(nrow(W33L.KO)-summary(factor(W33L.KO[,"W33L_call"]))[2])
summary(factor(W33L.WT[,"W33L_call"]))[1]/(nrow(W33L.WT)-summary(factor(W33L.WT[,"W33L_call"]))[2])

KO.tally <- W33L.KO %>%
  group_by(mouse.ID, W33L_call) %>% 
  count()
  
WT.tally <- W33L.WT %>%
  group_by(mouse.ID, W33L_call) %>% 
  count() 

write.csv(W33L, file = paste0(outdir,"W33L_calls.csv"))
write.csv(KO.tally, file = paste0(outdir,"W33L.KO.tally.csv"))
write.csv(WT.tally, file = paste0(outdir,"W33L.WT.tally.csv"))

save.image(file = "W33L.mutation.call.RData")



