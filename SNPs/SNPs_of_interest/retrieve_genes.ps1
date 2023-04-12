grep -E `
    "\sMGAM$|\sSI$|\sTREH$|\sSLC2A2$|\sSLC2A5$|\sSLC5A1$|\sLCT$|\sAMY[0-9][AB]$" `
    h3achip-annotated-with-gene.bed > retrieved_genes.bed
