# This script is written by Paul Bunel (paulbunel34@gmail.com) 2023
# GitHub: https://github.com/Paul-bunel
#
# Output: The script computes PBS making trios of population from the reference
# and outgroup variable, along with the targets list, then generates a plot at
# the following path : results/D_statistic/PBS_<target>_<reference>_<outgroup>.png
# Also saves a tsv file at the same path containing the PBS computed
#
# How to use: put the path of the directory containing your PLINK .fst.var
# in the dir_path variable, then put the name of the reference and an outgroup
# population you want in the corresponding variables, then run the script.
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

# Load file containing names of genes associated with SNP if needed

# dbSNP <- read.delim("Africaneo_dataset/h3achip-annotated-with-gene.bed",
#                     header = F,
#                     col.names = c("CHR", "N", "POS", "GENE")
#                     # colClasses = rep("character", 5)
# )

# Load .fst.var files in dir_path that contain either the reference or outgroup
# population in their name + use files names to get targets populations names

reference <- "Yoruba"
outgroup <- "Europe"

FSTs <- list()
targets <- list()
dir_path <- "Africaneo_dataset/pops_of_interest"
fst_files <- list.files(path = dir_path, pattern = "*.fst.var")
for (file in fst_files) {
  pops <- strsplit(file, ".", fixed = TRUE)[[1]][2:3]
  if (reference %in% pops || outgroup %in% pops) {
    targets <- append(targets, pops)
    path <- paste0("Africaneo_dataset/pops_of_interest/", file)
    FSTs[[file]] <- read.delim(
      path,
      col.names = c("CHROM", "POS", "ID", "NOBS", "FST")
    )
    
    FSTs[[file]]$FST[FSTs[[file]]$FST < 0] <- 0
  }
}

targets <- unique(targets)
targets <- targets[targets != reference]
targets <- targets[targets != outgroup]

targets <- "Baka"

BC_index <- grepl(reference, names(FSTs)) & grepl(outgroup, names(FSTs))

# Loop over targets to compute PBS and generate a plot for each

for (target in targets) {
  # Here, target, reference and outgroup population will be designed with
  # letters A, B and C, respectively.
  
  AB_index <- grepl(target, names(FSTs)) & grepl(reference, names(FSTs))
  AC_index <- grepl(target, names(FSTs)) & grepl(outgroup, names(FSTs))

  # Remove NA values

  na <- c()
  na <- sapply(FSTs[AB_index | AC_index | BC_index], \(x) {
    return(x$ID[is.na(x$FST)])
  })
  na <- unlist(na, use.names = FALSE)
  
  tmp_SNPs_of_interest <- SNPs_of_interest[!is.element(SNPs_of_interest$ID, na),]
  
  AB_FST <- FSTs[AB_index][[1]][!is.element(FSTs[AB_index][[1]]$ID, na),]
  AC_FST <- FSTs[AC_index][[1]][!is.element(FSTs[AC_index][[1]]$ID, na),]
  BC_FST <- FSTs[BC_index][[1]][!is.element(FSTs[BC_index][[1]]$ID, na),]

  # Compute PBS for each SNP and store it in a list
  
  T_AB <- -log(1 - AB_FST$FST)
  T_AC <- -log(1 - AC_FST$FST)
  T_BC <- -log(1 - BC_FST$FST)

  PBS <- (T_AB + T_AC + T_BC) / 2
  
  trio_df <- data.frame(
    CHR = AB_FST$CHROM,
    POS = AB_FST$POS,
    ID = AB_FST$ID,
    PBS = PBS
  )
  
  res_file_name <- paste0(
    "results/tmp/PBS_",
    paste(target, reference, outgroup, sep = "_"),
    ".tsv"
  )
  write.table(trio_df,
    file = res_file_name,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  # Following code for plot was inspired by the code from:
  # https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
  
  trio_cumul <- trio_df %>%
    group_by(CHR) %>%
    summarise(max_bp = max(POS)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHR, bp_add)
  
  trio_df <- trio_df %>%
    inner_join(trio_cumul, by = "CHR") %>%
    mutate(bp_cumul = POS + bp_add)
  
  trio_axis_set <- trio_df %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cumul))
  
  trio_genes_of_interest_window <- trio_df[as.logical(rowSums(
    sapply(
      seq_len(ncol(genes_of_interest)),
      \(x) {
        trio_df$CHR == genes_of_interest[1, x] & (
          between(trio_df$POS,
                  genes_of_interest[2, x] - window_size,
                  genes_of_interest[2, x]
          ) |
            between(trio_df$POS,
                    genes_of_interest[3, x],
                    genes_of_interest[3, x] + window_size
            )
        )
      }
    )
  )), ]
    
  manplot <- ggplot(trio_df, aes(x = bp_cumul, y = PBS)) +
    geom_point(alpha = 0.75, aes(colour = ifelse(CHR %% 2 == 0, "even", "odd"))) +
    geom_point(
      data = trio_genes_of_interest_window,
      aes(x = bp_cumul, y = PBS, color = "genes_window")
    ) +
    geom_point(
      data = trio_df[trio_df$ID %in% tmp_SNPs_of_interest$ID,],
      aes(x = bp_cumul, y = PBS, color = "SNPs_interest")
    ) +
    geom_hline(
      linetype = "dashed",
      aes(
        yintercept = quantile(PBS, 0.995),
        color = "quantile"
      )
    ) +
    scale_x_continuous(label = trio_axis_set$CHR, breaks = trio_axis_set$center) +
    scale_y_continuous(limits = c(min(trio_df$PBS), max(trio_df$PBS))) +
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
      title = paste(target, "PBS for each SNP, with", reference,
                    "as reference and", outgroup, "as outgroup."),
      x = "Position on genome"
    ) +
    annotate(
      "rect",
      xmin=31699382+492250183,
      xmax=32119072+492250183,
      ymin=-Inf,
      ymax=+Inf,
      # color="transparent",
      fill="green",
      alpha=0.2
    ) +
    annotate(
      "rect",
      xmin=50712358+492250183,
      xmax=51421629+492250183,
      ymin=-Inf,
      ymax=+Inf,
      # color="transparent",
      fill="green",
      alpha=0.2
    ) +
    coord_cartesian(xlim = c(492350902, 592350902)) +
    geom_text(
      # data = trio_df[trio_df$ID %in% tmp_SNPs_of_interest$ID,],
      aes(label = ifelse(
        PBS > quantile(trio_df$PBS, 0.995) & PBS > 0.75,
        ID,
        ''
      )),
      hjust = 0, vjust = 0
    )

    # Parts after "labs" instruction are highlighting a gene region, zooming on
    # it and adding text, respectively. They are optional and can be commented

  plot_name <- paste0(
    "results/tmp/PBS_",
    paste(target, reference, outgroup, sep = "_"),
    ".png"
  )
  ggsave(plot_name, width = 14, height = 8, bg = "white")
}

