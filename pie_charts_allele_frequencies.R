library(dplyr)
library(leaflet)
library(htmlwidgets)
library(leaflet.minicharts)
library(htmltools)

coordinates <- read.csv2("Africaneo_dataset_geographical_coordinates.csv")
SNPs <- read.table("SNPs_of_interest_analysis_clst.frq.strat", header = T,
                   sep = "\t", colClasses = c(rep("character", 4), rep("numeric", 3)))
SNP_of_interest <- "1_9115700"
SNP <- subset(SNPs, CODE == SNP_of_interest)
# SNP <- read.table("retrieved_MGAM_first_SNP.frq.strat", header = T, sep="\t",
#                   colClasses = c(rep("character", 4), rep("numeric", 3)))
colnames(SNP)[colnames(SNP) == "MAC"] = "MA1"
SNP$MA2 <- SNP$NCHROBS - SNP$MA1
SNP_coordinates <- merge(
  x = SNP,
  y = coordinates,
  by.x = "CLST",
  by.y = "Population.label"
)

chart_datas <- SNP_coordinates[, c("MA1", "MA2")]

map <- leaflet() %>%
  addTiles(
    "https://tile.thunderforest.com/cycle/{z}/{x}/{y}.png?apikey=688b32b79cbf44edb274836064f80b41"
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
              "Echantillon : ", SNP_coordinates$N,
              "<br>", SNP_coordinates$A1[1], " : ", SNP_coordinates$MA1,
              "<br>", SNP_coordinates$A2[1], " : ", SNP_coordinates$MA2
            )
  ) %>%
  addLegend(
    position = "topright",
    labels = SNP_coordinates[1,c("A1", "A2")],
    colors = d3.schemeCategory10[1:length(chart_datas)],
    opacity = 1,
    layerId = "minichartsLegend",
  )
map

saveWidget(map, file="map.html")
