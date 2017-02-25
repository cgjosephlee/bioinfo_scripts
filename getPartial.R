# extract sequence with given coordinates
# Hsinhan Lee, 2017

# Usuage
# 0. copy this script to where the fasta file at
# 1. open this script with RStudio (must, test in v1.0.136)
# 2. edit the following parameters
#    in_fa: your fasta file (supposing they are in the same directory)
#    ID   : sequence title
#           if (ID <- ""), default is the first sequence in fasta
#    coord: sequence coordinates c(start, end)
#           if start > end, output will be reverse complentment (DNA only)
in_fa <- "test.fa"
ID <- "PNOA_0000900.1"
coord <- c(50, 21)
# 3. hit Ctrl+Alt+R (run all!), or Cmd+Opt+R for macOS
# 4. find the output "sub_seq.fasta"

### main ###
if (!require("Biostrings", character.only = TRUE))
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
  if(!require("Biostrings", character.only = TRUE)) stop("Package not found")
}

# edit by your favor if you are smart
out_fa <- "sub_seq.fasta"
path <- dirname(rstudioapi::getActiveDocumentContext()$path) # where the script at

fa <- readBStringSet(filepath = file.path(path, in_fa))
start <- coord[1]; end <- coord[2]

if(ID == ""){
  ID <- names(fa[1])
}

if(start < end){ # normal situation
  out_seq <- subseq(fa[ID], start, end)
  names(out_seq) <- paste0(ID, ":", start, "-", end, ";len=", width(out_seq))
} else if(start > end){ # only support DNA strings
  fa <- DNAStringSet(fa)
  out_seq <- reverseComplement(subseq(fa[ID], end, start))
  names(out_seq) <- paste0(ID, "_rev:", start, "-", end, ";len=", width(out_seq))
}

writeXStringSet(out_seq, filepath = file.path(path, out_fa))
