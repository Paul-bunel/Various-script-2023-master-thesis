library(reshape2)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(htmltools)
library(qqman)

pop_datas <- read.csv("Africaneo_dataset_simplified_c.csv")
pop_datas <- pop_datas[c("Group", "Population_label", "N", "main_lifestyle")]
SNPs <- read.table("SNPs_of_interest_analysis_clst.frq.strat", header = T,
                   sep = "\t", colClasses = c(rep("character", 4), rep("numeric", 3)))

SNP_by_population <- dcast(SNPs, CLST ~ CODE, value.var = "MAF")

SNP_by_population <- merge(
  x = pop_datas,
  y = SNP_by_population,
  by.x = "Population_label",
  by.y = "CLST"
)

SNPs_of_interest <- read.csv("SNPs_of_interest.csv", colClasses = rep("character", 5))

MGAM <- SNPs_of_interest[grep("MGAM", SNPs_of_interest$GENE), c("CHR", "POS")]
MGAM <- paste0(MGAM[,1], "_", MGAM[,2])

SI <- SNPs_of_interest[grep("SI", SNPs_of_interest$GENE), c("CHR", "POS")]
SI <- paste0(SI[,1], "_", SI[,2])

SLC2A5 <- SNPs_of_interest[grep("SLC2A5", SNPs_of_interest$GENE), c("CHR", "POS")]
SLC2A5 <- paste0(SLC2A5[,1], "_", SLC2A5[,2])

AMY2B <- SNPs_of_interest[grep("AMY2B", SNPs_of_interest$GENE), c("CHR", "POS")]
AMY2B <- paste0(AMY2B[,1], "_", AMY2B[,2])

LCT <- SNPs_of_interest[grep("LCT", SNPs_of_interest$GENE), c("CHR", "POS")]
LCT <- paste0(LCT[,1], "_", LCT[,2])

SLC5A1 <- SNPs_of_interest[grep("SLC5A1", SNPs_of_interest$GENE), c("CHR", "POS")]
SLC5A1 <- paste0(SLC5A1[,1], "_", SLC5A1[,2])

SLC2A2 <- SNPs_of_interest[grep("SLC2A2", SNPs_of_interest$GENE), c("CHR", "POS")]
SLC2A2 <- paste0(SLC2A2[,1], "_", SLC2A2[,2])

farmers_means <- cbind(
  colnames(SNP_by_population[5:85]),
  colMeans(SNP_by_population[SNP_by_population$main_lifestyle == "farmer",5:85])
)
pastoralists_means <- cbind(
  colnames(SNP_by_population[5:85]),
  colMeans(SNP_by_population[SNP_by_population$main_lifestyle == "pastoralist",5:85])
)
foragers_means <- cbind(
  colnames(SNP_by_population[5:85]),
  colMeans(SNP_by_population[SNP_by_population$main_lifestyle == "foragers",5:85])
)
freq_mean_by_lifestyle <- rbind(farmers_means, pastoralists_means, foragers_means)
freq_mean_by_lifestyle <- cbind(freq_mean_by_lifestyle, c(
    rep("farmer", 81),
    rep("pastoralist", 81),
    rep("foragers", 81)
  )
)
colnames(freq_mean_by_lifestyle) <- c("CODE", "MAF", "LIFESTYLE")
freq_mean_by_lifestyle <- data.frame(freq_mean_by_lifestyle)
freq_mean_by_lifestyle$MAF <- as.numeric(freq_mean_by_lifestyle$MAF)
freq_mean_by_lifestyle$LIFESTYLE <- as.factor(freq_mean_by_lifestyle$LIFESTYLE)

plot(
  rowMeans(SNP_by_population[5:85]),
  col=ifelse(
    SNP_by_population$main_lifestyle == "farmer", "green",
    ifelse(
      SNP_by_population$main_lifestyle == "pastoralist", "red",
      ifelse(
        SNP_by_population$main_lifestyle == "foragers", "blue", "black"
      )
    )
  )
)

plot(
  rowMeans(SNP_by_population[MGAM]),
  col=ifelse(
    SNP_by_population$main_lifestyle == "farmer", "green",
    ifelse(
      SNP_by_population$main_lifestyle == "pastoralist", "red",
      ifelse(
        SNP_by_population$main_lifestyle == "foragers", "blue", "black"
      )
    )
  )
)

plot(
  rowMeans(SNP_by_population["7_141738362"]),
  col=ifelse(
    SNP_by_population$main_lifestyle == "farmer", "green",
    ifelse(
      SNP_by_population$main_lifestyle == "pastoralist", "red",
      ifelse(
        SNP_by_population$main_lifestyle == "foragers", "blue", "black"
      )
    )
  )
)

# plot(MAF ~ LIFESTYLE, freq_mean_by_lifestyle)
# plot.default(
#   freq_mean_by_lifestyle$LIFESTYLE,
#   freq_mean_by_lifestyle$MAF
# )

p <- ggplot(freq_mean_by_lifestyle, aes(
    x = LIFESTYLE,
    y = MAF,
    text = paste("CODE :", CODE)
  )) +
  geom_line(aes(group = CODE), linewidth=0.1) +
  geom_point(aes(x = LIFESTYLE), size=0.5)

  # geom_label(aes(label = CODE))
  # geom_segment(aes(x = LIFESTYLE))

gg <- ggplotly(p, tooltip = c("text", "MAF"))
saveWidget(gg, file="plot.html")


freq_mean_by_lifestyle %>%
  select(MAF, LIFESTYLE, CODE) %>%
  pivot_longer(-MAF&LIFESTYLE, names_to = "stat", values_to = "val") %>%
  ggplot(aes(x = val, y = MAF)) +
  geom_point() +
  geom_segment(aes(
    x = freq_mean_by_lifestyle$LIFESTYLE == "farmer",
    xend = freq_mean_by_lifestyle$LIFESTYLE == "foragers",
    yend = MAF
  ))




