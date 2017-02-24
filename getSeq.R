# extract sequences from a fasta file
# Hsinhan Lee, 2017

# Usuage
# 0. copy this script to where the fasta file at
# 1. open this script with RStudio (recommand)
# 2. edit the following two filenames (supposing they are in the same directory)
#    in_fa: your fasta file
#    wanted_list: a text file with IDs, one ID a line
in_fa <- "test.aa.fa"
wanted_list <- "test.list.txt"
# 3. hit Ctrl+Alt+R (run all!), or Cmd+Opt+R for macOS
# 4. find the output "filter.fasta"

# install the reguired library, only run at the fisrt time
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")

### main ###
library(Biostrings)

# edit by your favor if your are smart
out_fa <- "filtered.fasta"
path <- getwd() # where the script at

fa <- readBStringSet(filepath = file.path(path, in_fa))
readlist <- file(file.path(path, wanted_list), open = "r")
list <- readLines(readlist)
close(readlist)
fa <- fa[names(fa) %in% list]
writeXStringSet(fa, filepath = file.path(path, out_fa))
