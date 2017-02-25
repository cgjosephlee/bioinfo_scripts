# extract sequences with given IDs
# Hsinhan Lee, 2017

# Usuage
# 0. copy this script to where the fasta file at
# 1. open this script with RStudio (must, test in v1.0.136)
# 2. edit the following two filenames (supposing they are in the same directory)
#    in_fa      : your fasta file
#    wanted_list: a text file with IDs, one ID a line
in_fa <- "test.fa"
wanted_list <- "test.list.txt"
# 3. hit Ctrl+Alt+R (run all!), or Cmd+Opt+R for macOS
# 4. find the output "filter.fasta"

### main ###
if (!require("Biostrings", character.only = TRUE))
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
  if(!require("Biostrings", character.only = TRUE)) stop("Package not found")
}

# edit by your favor if you are smart
out_fa <- "filtered.fasta"
path <- dirname(rstudioapi::getActiveDocumentContext()$path) # where the script at

fa <- readBStringSet(filepath = file.path(path, in_fa))
readlist <- file(file.path(path, wanted_list), open = "r")
list <- readLines(readlist)
close(readlist)
fa <- fa[names(fa) %in% list]
writeXStringSet(fa, filepath = file.path(path, out_fa))
