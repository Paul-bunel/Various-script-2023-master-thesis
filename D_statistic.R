# This script is written by Paul Bunel (paulbunel34@gmail.com)
# GitHub: https://github.com/Paul-bunel
#
# Output: The script compute D statistic for each SNP between every group of
# population present in a directory containing PLINK .fst.var files, then
# generate a plot at the following path :
# plots/D_statistic/D_statistic_<pop>.png
#
# How tu use: put the path of the directory containing your PLINK .fst.var
# in the dir_path variable, then run the script
# If you don't want to highlight SNPs or genes of interest, just make the
# corresponding variables empty

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

# Loop over populations names to compute D statistic and genereate plot for each

for (pop in pops) {
  pop_indexes <- grepl(pop, names(FSTs))
  na <- c()
  na <- sapply(FSTs[pop_indexes], \(x) {
    return(x$ID[is.na(x$FST)])
  })
  na <- unlist(na, use.names = FALSE)
  
  tmp_pops <- sapply(FSTs[pop_indexes], \(x) {
    return(x[!is.element(x$ID, na),])
  }, simplify = FALSE)
  
  print(names(tmp_pops))
  
  desired_length <- nrow(FSTs[[1]]) - length(unique(na))
  if (sum(sapply(tmp_pops, nrow) == rep(desired_length, 3)) != 3) {
    print("ATTENTION MAUVAISE TAILLE DE TMP_POPS")
  }
  
  tmp_SNPs_of_interest <- SNPs_of_interest[!is.element(SNPs_of_interest$ID, na),]

  pop_D <- rowSums(sapply(tmp_pops, \(x) { (x$FST - mean(x$FST)) / sd(x$FST) }))
  
  pop_df <- data.frame(
    CHR = tmp_pops[[1]]$CHROM,
    POS = tmp_pops[[1]]$POS,
    ID = tmp_pops[[1]]$ID,
    D = pop_D
  )
  
  # Following code for manhattan plot come from
  # https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

  pop_cumul <- pop_df %>%
    group_by(CHR) %>%
    summarise(max_bp = max(POS)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHR, bp_add)
  
  pop_df <- pop_df %>%
    inner_join(pop_cumul, by = "CHR") %>%
    mutate(bp_cumul = POS + bp_add)
  
  pop_axis_set <- pop_df %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cumul))
  
  pop_genes_of_interest_window <- pop_df[as.logical(rowSums(
    sapply(
      1:ncol(genes_of_interest),
      \(x) {
        pop_df$CHR == genes_of_interest[1, x] & (
          between(pop_df$POS,
            genes_of_interest[2, x] - window_size,
            genes_of_interest[2, x]
          ) |
          between(pop_df$POS,
            genes_of_interest[3, x],
            genes_of_interest[3, x] + window_size
          )
        )
      }
    )
  )),]
  
  manplot <- ggplot(pop_df, aes(x = bp_cumul, y = D)) +
    geom_point(alpha = 0.75, aes(colour = ifelse(CHR %% 2 == 0, "1", "2"))) +
    geom_point(
      data = pop_genes_of_interest_window,
      aes(x = bp_cumul, y = D, color = "4")
    ) +
    geom_point(
      data = pop_df[pop_df$ID %in% tmp_SNPs_of_interest$ID,],
      aes(x = bp_cumul, y = D, color = "0")
    ) +
    geom_hline(
      yintercept = quantile(pop_df$D, 0.995),
      color = "red"
      linetype = "dashed"
    ) +
    scale_x_continuous(label = pop_axis_set$CHR, breaks = pop_axis_set$center) +
    scale_y_continuous(limits = c(min(pop_df$D), max(pop_df$D))) +
    scale_color_manual(values = setNames(
      c("orange", "yellow", "#183059", "#276FBF"),
      c("4", "0", "1", "2")
    )) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )

  plot_name <- paste0("plots/D_statistic/D_statistic_", pop, ".png")
  ggsave(plot_name, width = 12, height = 8, bg = "white")
}
