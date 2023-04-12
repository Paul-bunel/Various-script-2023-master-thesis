library(tidyverse)
library(ggtext)
library(normentR)

# Load genomic positions of genes of interest

genes_of_interest <- data.frame(
  MGAM = c(7, 141607613, 141806547),
  SI = c(3, 164696686, 164796284),
  TREH = c(11, 118528026, 118550359),
  SLC2A2 = c(3, 170714137, 170744539),
  SLC2A5 = c(1, 9095166, 9148537),
  SLC5A1 = c(22, 32439248, 32509016),
  TAS1R3 = c(1, 1266660, 1270694),
  LCT = c(2, 136545420, 136594754),
  AMY1A = c(1, 104198141, 104207176),
  AMY1B = c(1, 104230037, 104239075),
  AMY1C = c(1, 104292276, 104301314),
  AMY2A = c(1, 104159999, 104168402),
  AMY2B = c(1, 104096437, 104122156)
)
window_size <- 50000

# Load genomic position of SNPs of interest

SNPs_of_interest <- read.delim(
  "SNPs/SNPs_of_interest/SNPs_of_interest.map",
  header = FALSE,
  col.names = c("CHR", "ID", "#", "POS")
)

# Load every .fst.var file in SNPs/pops_of_interest + use files names to get
# populations names

FSTs <- list()
pops <- list()
dir_path <- "SNPs/pops_of_interest"
fst_files <- list.files(path = dir_path, pattern = "*.fst.var")
for (file in fst_files) {
  pops <- append(pops, strsplit(file, ".", fixed = TRUE)[[1]][2:3])
  path <- paste0("SNPs/pops_of_interest/", file)
  FSTs[[file]] <- read.delim(
    path,
    col.names = c("CHROM", "POS", "ID", "NOBS", "FST")
  )
  
  FSTs[[file]]$FST[FSTs[[file]]$FST < 0] <- 0
}

pops <- unique(pops)
