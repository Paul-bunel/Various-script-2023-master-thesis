# This script is written by Paul Bunel (paulbunel34@gmail.com) 2O23
# GitHub: https://github.com/Paul-bunel
#
# Output: The script computes D statistic from Akey et al. (2010) for all SNP in
# each populations present in a directory containing PLINK .fst.var files, then
# generates a plot at the following path :
# results/D_statistic/D_statistic_<pop>.png
# Also saves a tsv file at the same path containing the computed D statistic
#
# How to use: put the path of the directory containing your PLINK .fst.var
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

# SNPs_of_interest <- read.delim(
#   "Africaneo_dataset/SNPs_of_interest/SNPs_of_interest.map",
#   header = FALSE,
#   col.names = c("CHR", "ID", "#", "POS")
# )

SNPs_of_interest <- read.delim(
  "dataset2/SNPs_of_interest.bim",
  header = FALSE,
  col.names = c("CHR", "ID", "#", "POS", "A", "B")
)

# Load file containing names of genes associated with SNP if needed

# Load every .fst.var file in dir_path + use files names to get
# populations names

FSTs <- list()
pops <- list()
# dir_path <- "Africaneo_dataset/pops_of_interest"
dir_path <- "dataset2/Pops/FST"
fst_files <- list.files(path = dir_path, pattern = "*.fst.var")
for (file in fst_files) {
  pops <- append(pops, strsplit(file, ".", fixed = TRUE)[[1]][2:3])
  path <- paste(dir_path, file, sep="/")
  FSTs[[file]] <- read.delim(
    path,
    col.names = c("CHROM", "POS", "ID", "NOBS", "FST")
  )

  FSTs[[file]]$FST[FSTs[[file]]$FST < 0] <- 0
}

pops <- unique(pops)

# Loop over populations names to compute D statistic and genereate plot for each

# pops <- c("Senegal")
comp <- c(
  "Namibia_TsumkweKung",
  "SouthAfrica_Khomani",
  "DRC",
  "East_Africa",
  "SA",
  "Senegal"
)

for (pop in pops) {
  pop_indexes <- grepl(
    paste0('\\.', pop, '\\.'), names(FSTs)) &
    !grepl("\\.Namibia\\.|Sudan", names(FSTs)
  )

  # Remove NA values

  na <- c()
  na <- sapply(FSTs[pop_indexes], \(x) {
    return(x$ID[is.na(x$FST)])
  })
  na <- unlist(na, use.names = FALSE)
  
  tmp_SNPs_of_interest <- SNPs_of_interest[!is.element(SNPs_of_interest$ID, na),]

  tmp_pops <- sapply(FSTs[pop_indexes], \(x) {
    return(x[!is.element(x$ID, na),])
  }, simplify = FALSE)

  # Compute D statistic for each SNP and store it in a list

  pop_D <- rowSums(sapply(tmp_pops, \(x) { (x$FST - mean(x$FST)) / sd(x$FST) }))
  
  pop_df <- data.frame(
    CHR = tmp_pops[[1]]$CHROM,
    POS = tmp_pops[[1]]$POS,
    ID = tmp_pops[[1]]$ID,
    D = pop_D
  )
  

  # Following code for plot was inspired by the code from:
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
      seq_len(ncol(genes_of_interest)),
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
  )), ]

  res_file_name <- paste0("results/tmp/D_statistic_", pop, ".tsv")
  res_SOI <- pop_df[
    pop_df$ID %in% tmp_SNPs_of_interest$ID &
      pop_df$D > quantile(pop_df$D, 0.995),
    c("CHR", "POS", "D")
  ]
  write.table(res_SOI,
              file = res_file_name,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE
  )

  manplot <- ggplot(pop_df, aes(x = bp_cumul, y = D)) +
    geom_point(alpha = 0.75, aes(colour = ifelse(CHR %% 2 == 0, "even", "odd"))) +
    geom_point(
      data = pop_genes_of_interest_window,
      aes(x = bp_cumul, y = D, color = "genes_window")
    ) +
    geom_point(
      data = pop_df[pop_df$ID %in% tmp_SNPs_of_interest$ID,],
      aes(x = bp_cumul, y = D, color = "SNPs_interest")
    ) +
    geom_hline(
      linetype = "dashed",
      aes(
        yintercept = quantile(D, 0.995),
        color = "quantile"
      )
    ) +
    scale_x_continuous(label = pop_axis_set$CHR, breaks = pop_axis_set$center) +
    scale_y_continuous(limits = c(-10, 100)) +
    scale_color_manual(
      breaks = c(
        "even",
        "odd",
        "genes_window",
        "SNPs_interest",
        "quantile"
      ),
      values = c(
        "even" = "#183059",
        "odd" = "#276FBF",
        "genes_window" = "orange",
        "SNPs_interest" = "yellow",
        "quantile" = "red"
      ),
      labels = c(
        "Even chromosome",
        "Odd Chromosome",
        str_wrap(paste0("Window around genes of interest (", window_size, " bp)"), 25),
        "SNPs of interest",
        "Genome-wide 0.995 quantile"
      ),
      guide = guide_legend(title = "Legend", override.aes = list(
        linetype = c(rep("blank", 4), "dashed"),
        shape = c(rep(19, 4), NA)
      ))
    ) +
    geom_text(
      data = pop_df[pop_df$ID %in% tmp_SNPs_of_interest$ID,],
      aes(label = ifelse(
        D > quantile(pop_df$D, 0.995),
        ID,
        ''
      )),
      hjust = 0, vjust = 0
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.justification = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 9),
      panel.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
      panel.grid.major = element_line(color = "white", linewidth = 1),
      panel.grid.minor = element_line(color = "white", linewidth = 0.75),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
    ) +
    labs(
      title = paste(pop, "D statistic for each SNP"),
      x = "Position on genome"
    ) +
    annotate(
      "rect",
      xmin=32627244+1061782806,
      xmax=32636160+1061782806,
      ymin=-Inf,
      ymax=+Inf,
      color="transparent",
      fill="green",
      alpha=0.2
    )# +
    # annotate(
    #   "rect",
    #   xmin=50712358+492250183,
    #   xmax=51421629+492250183,
    #   ymin=-Inf,
    #   ymax=+Inf,
    #   # color="transparent",
    #   fill="green",
    #   alpha=0.2
    # ) +
    # coord_cartesian(xlim = c(1093419927, 1095419927))# +
    # geom_text(
    #   # data = trio_df[trio_df$ID %in% tmp_SNPs_of_interest$ID,],
    #   aes(label = ifelse(
    #     D > quantile(D, 0.995) & D > 60,
    #     ID,
    #     ''
    #   )),
    #   hjust = 0, vjust = 0
    # )

    # Commented parts after "labs" instruction are for higlighting the CD36 gene
    # region and zoom on it, respectively

  plot_name <- paste0("results/tmp/D_statistic_", pop, ".png")
  ggsave(plot_name, width = 14, height = 8, bg = "white")
}
