# This script is written by Paul Bunel (paulbunel34@gmail.com) 2O23
# GitHub: https://github.com/Paul-bunel
#
# Output: The script compute D statistic for each SNP between every group of
# population present in a directory containing PLINK .fst.var files, then
# generate a plot at the following path :
# results/D_statistic/D_statistic_<pop>.png
# Also save a tsv file at the same path containing the D statistic computed
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

SNPs_of_interest <- read.delim(
  "Africaneo_dataset/SNPs_of_interest/SNPs_of_interest.map",
  header = FALSE,
  col.names = c("CHR", "ID", "#", "POS")
)

# dbSNP <- read.delim("Africaneo_dataset/h3achip-annotated-with-gene.bed",
#   header = F,
#   col.names = c("CHR", "N", "POS", "GENE")
#   # colClasses = rep("character", 5)
# )

# Load every .fst.var file in dir_path + use files names to get
# populations names

FSTs <- list()
pops <- list()
dir_path <- "Africaneo_dataset/pops_of_interest"
fst_files <- list.files(path = dir_path, pattern = "*.fst.var")
for (file in fst_files) {
  pops <- append(pops, strsplit(file, ".", fixed = TRUE)[[1]][2:3])
  path <- paste0("Africaneo_dataset/pops_of_interest/", file)
  FSTs[[file]] <- read.delim(
    path,
    col.names = c("CHROM", "POS", "ID", "NOBS", "FST")
  )
  
  FSTs[[file]]$FST[FSTs[[file]]$FST < 0] <- 0
}

pops <- unique(pops)

# Loop over populations names to compute D statistic and genereate plot for each

# pops <- c("Yoruba")
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
  
  tmp_SNPs_of_interest <- SNPs_of_interest[!is.element(SNPs_of_interest$ID, na),]

  pop_D <- rowSums(sapply(tmp_pops, \(x) { (x$FST - mean(x$FST)) / sd(x$FST) }))
  
  pop_df <- data.frame(
    CHR = tmp_pops[[1]]$CHROM,
    POS = tmp_pops[[1]]$POS,
    ID = tmp_pops[[1]]$ID,
    D = pop_D
  )
  
  res_file_name <- paste0("results/tmp/D_statistic_", pop, ".tsv")
  write.table(pop_df,
    file = res_file_name,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
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
  )),]
  
  # sof <- pop_df[
  #   pop_df$CHR == 7 &
  #     pop_df$bp_cumul > 1310000000 &
  #     pop_df$bp_cumul < 1320000000 &
  #     pop_df$D > quantile(pop_df$D, 0.995)
  # , ]
  
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
    # geom_vline(linetype="dashed", xintercept = 1312567881, color = "yellow") +
    # geom_vline(linetype="dashed", xintercept = 1312733432, color = "green") +
    scale_x_continuous(label = pop_axis_set$CHR, breaks = pop_axis_set$center) +
    scale_y_continuous(limits = c(min(pop_df$D), max(pop_df$D))) +
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
    geom_rect(
      # data=matrix_amy1,
      # inherit.aes=FALSE,
      aes(xmin=79998891+1232462990, xmax=80308593+1232462990, ymin=-Inf, ymax=+Inf),
      color="transparent",
      fill="green",
      alpha=0.3
    )
    
    # coord_cartesian(xlim = c(1310000000, 1320000000))
    # geom_text(aes(label = ifelse(
    #       POS %in% sof$POS & CHR == 7,
    #       dbSNP$GENE[dbSNP$CHR == 7 & dbSNP$POS %in% sof$POS],
    #       ''
    #     )),
    #   hjust = 0, vjust = 0
    # )

  plot_name <- paste0("results/tmp/D_statistic_", pop, ".png")
  ggsave(plot_name, width = 14, height = 8, bg = "white")
}
