library(tidyverse)
library(ggtext)
library(normentR)

# Compute D statistic for each SNP between two groups of population

# Load every .fst.var file in SNPs/pops_of_interest
# Then for each population, compute D statistic

FSTs <- list()
pops <- list()
fst_files <- list.files(path = "SNPs/pops_of_interest", pattern = "*.var")
for (file in fst_files) {
  pops <- append(pops, strsplit(file, ".", fixed = TRUE)[[1]][2:3])
  path <- paste0("SNPs/pops_of_interest/", file)
  FSTs[[file]] <- read.delim(
    path,
    col.names = c("CHROM", "POS", "ID", "NOBS", "FST")
  )
}

pops <- unique(pops)
na <- vector("list", length(pops))
names(na) <- pops

read_plink2_fst_file <- function(pop1, pop2) {
  path <- paste(
    "SNPs/pops_of_interest/pops_of_interest_fst",pop1,pop2,"fst.var",
    sep = "."
  )
  if (!file.exists(path)) {
    path <- paste(
      "SNPs/pops_of_interest/pops_of_interest_fst",pop2,pop1,"fst.var",
      sep = "."
    )
    if (!file.exists(path)) {
      return(data.frame())
    }
  }
  p1p2 <- read.delim(path, col.names = c("CHROM", "POS", "ID", "NOBS", "FST"))
  na <<- append(na, p1p2$ID[is.na(p1p2$FST)])
  p1p2 <- p1p2[!is.na(p1p2$FST),]
  p1p2$FST[p1p2$FST < 0] <- 0
  
  return(p1p2)
}

jap_yor <- read_plink2_fst_file("Japanese", "Yoruba")
jap_drc <- read_plink2_fst_file("Japanese", "DRC")
jap_bak <- read_plink2_fst_file("Japanese", "Baka")

jap_yor <- jap_yor[!is.element(jap_yor$ID, na),]
jap_drc <- jap_drc[!is.element(jap_drc$ID, na),]
jap_bak <- jap_bak[!is.element(jap_bak$ID, na),]

mean_jap_yor <- mean(jap_yor$FST)
sd_jap_yor <- sd(jap_yor$FST)

mean_jap_drc <- mean(jap_drc$FST)
sd_jap_drc <- sd(jap_drc$FST)

mean_jap_bak <- mean(jap_bak$FST)
sd_jap_bak <- sd(jap_bak$FST)


jap_D <- ((jap_yor$FST - mean_jap_yor) / sd_jap_yor) +
  ((jap_drc$FST - mean_jap_drc) / sd_jap_drc) +
  ((jap_bak$FST - mean_jap_bak) / sd_jap_bak)

jap <- data.frame(
  CHR = jap_bak$CHROM, POS = jap_bak$POS, ID = jap_bak$ID,
  jap_bak$NOBS, jap_bak$FST,
  jap_drc$NOBS, jap_drc$FST,
  jap_yor$NOBS, jap_yor$FST,
  D = jap_D
)

SNPs_of_interest <- read.delim(
  "SNPs/SNPs_of_interest/SNPs_of_interest.map",
  header = FALSE,
  col.names = c("CHR", "ID", "#", "POS")
)
SNPs_of_interest <- SNPs_of_interest[!is.element(SNPs_of_interest$ID, na),]

# Following code for manhattan plot come from
# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

jap_cumul <- jap %>%
  group_by(CHR) %>%
  summarise(max_bp = max(POS)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

jap <- jap %>%
  inner_join(jap_cumul, by = "CHR") %>%
  mutate(bp_cumul = POS + bp_add)

jap_axis_set <- jap %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cumul))

manplot <- ggplot(jap, aes(x = bp_cumul, y = D)) +
  geom_point(alpha = 0.75, aes(colour = ifelse(CHR %% 2 == 0, "1", "2"))) +
  geom_point(
    data = jap[jap$ID %in% SNPs_of_interest$ID,],
    aes(x = bp_cumul, y = D, color = "0")) +
  geom_hline(
    yintercept = quantile(jap$D, 0.995),
    color = "red",
    linetype = "dashed"
  ) +
  scale_x_continuous(label = jap_axis_set$CHR, breaks = jap_axis_set$center) +
  scale_y_continuous(limits = c(min(jap$D), max(jap$D))) +
  scale_color_manual(values = setNames(
    c("yellow", "#183059", "#276FBF"),
    c("0", "1", "2")
  )) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
# ggsave("plots/D_statistic/D_statistic_JAP.png")

