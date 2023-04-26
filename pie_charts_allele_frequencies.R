library(dplyr)
library(leaflet)
library(htmlwidgets)
library(leaflet.minicharts)
library(htmltools)

coordinates <- read.csv("Africaneo_dataset/Africaneo_dataset_simplified_c.csv")
SNPs <- read.table(
  "Africaneo_dataset/SNPs_of_interest/SNPs_of_interest_analysis_clst.frq.strat",
  header = T,
  sep = "\t",
  colClasses = c(rep("character", 4), rep("numeric", 3))
)

SNP_of_interest <- "7_141761438"

SNP <- subset(SNPs, CODE == SNP_of_interest)
colnames(SNP)[colnames(SNP) == "MAC"] = "MA1"
SNP$MA2 <- SNP$NCHROBS - SNP$MA1
SNP_coordinates <- merge(
  x = SNP,
  y = coordinates,
  by.x = "CLST",
  by.y = "Population_label"
)

dbSNP <- read.csv(
  "Africaneo_dataset/SNPs_of_interest/dbSNP_retrieved.csv",
  header = F,
  col.names = c("GENE", "RSID", "CHR", "POS", "ALLELES"),
  colClasses = rep("character", 5)
)
CHR_of_interest <- strsplit(SNP_of_interest, "_")[[1]][1]
POS_of_interest <- strsplit(SNP_of_interest, "_")[[1]][2]
GENE_of_interest <- dbSNP$GENE[dbSNP$CHR == CHR_of_interest & dbSNP$POS == POS_of_interest]

chart_datas <- SNP_coordinates[, c("MA1", "MA2")]

tag.map.title <- tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 15px;
  }
"))

title <- tags$div(
  tag.map.title, HTML(paste0(
    "GENE OF INTEREST : ", GENE_of_interest, "<br>",
    "CHROMOSOME : ", CHR_of_interest, "<br>",
    "POSITION : ", POS_of_interest
  ))
)

lifestyles <- c("farmer", "pastoralist", "foragers")

map <- leaflet() %>%
  addTiles(
    "https://tile.thunderforest.com/cycle/{z}/{x}/{y}.png?apikey=688b32b79cbf44edb274836064f80b41"
  ) %>%
  addCircleMarkers(
    lat = SNP_coordinates$Latitude, lng = SNP_coordinates$Longitude,
    opacity = ifelse(SNP_coordinates$main_lifestyle %in% lifestyles, 1, 0),
    weight = 3,
    color = ifelse(
      SNP_coordinates$main_lifestyle == "farmer", 'green',
      ifelse(
        SNP_coordinates$main_lifestyle == "pastoralist", 'red', 'blue'
      )
    ),
    fillOpacity = 0,
    radius = sqrt(SNP_coordinates$N) / sqrt(max(SNP_coordinates$N)) * 20,
  ) %>%
  addMinicharts(
    lat = SNP_coordinates$Latitude, lng = SNP_coordinates$Longitude,
    type = "pie",
    chartdata = chart_datas,
    width = sqrt(SNP_coordinates$N) / sqrt(max(SNP_coordinates$N)) * 40,
    legend = FALSE
  ) %>%
  addCircleMarkers(
    lat = SNP_coordinates$Latitude, lng = SNP_coordinates$Longitude,
    label = SNP_coordinates$CLST,
    opacity = 0,
    fillOpacity = 0,
    radius = sqrt(SNP_coordinates$N) / sqrt(max(SNP_coordinates$N)) * 20,
    popup = paste0(
              "Echantillon : ", SNP_coordinates$NCHROBS / 2, "<br>",
              SNP_coordinates$A1[1], " : ", SNP_coordinates$MA1, "<br>",
              SNP_coordinates$A2[1], " : ", SNP_coordinates$MA2
            )
  ) %>%
  addLegend(
    position = "topright",
    labels = SNP_coordinates[1,c("A1", "A2")],
    colors = d3.schemeCategory10[1:length(chart_datas)],
    opacity = 1,
    layerId = "minichartsLegend",
  ) %>%
  addControl(title, position = "topleft", className="map-title")
map

saveWidget(map, file=paste0("results/", SNP_of_interest, ".html"))
