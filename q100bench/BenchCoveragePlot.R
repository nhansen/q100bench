#setwd("/Users/nhansen/HG002_diploid_benchmark/plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "test")
outputdir <- ifelse(!( is.na(args[2])), args[2], "./")
resourcedir <- ifelse(!( is.na(args[3])), args[3], "./")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("Benchmark aligned coverage vs. ", genomename), sep="", collapse=""))
genomefile <- paste(c(resourcedir, "/v1.0.1.genome.bed"), sep="", collapse="")
benchgenome <- toGRanges(genomefile)

maternalchroms <- paste0("chr", c(1:22, "X"), "_MATERNAL")
paternalchroms <- paste0("chr", c(1:22, "Y"), "_PATERNAL")
chroms <- c(maternalchroms, paternalchroms)

benchcoveredfile <- paste(c(outputdir, "/", genomename, ".benchcovered.truth.merged.bed"), sep="", collapse="")
benchcoveredranges <- toGRanges(benchcoveredfile)
benchcovereddf <- read.table(benchcoveredfile, header=FALSE, sep="\t")

nlocfile <- paste(c(resourcedir, "/v1.0.1.ncoords.bed"), sep="", collapse="")
nlocranges <- toGRanges(nlocfile)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 100
pp$topmargin <- 300

plotname <- paste(c(outputdir, "/", genomename, ".benchcovered.truth.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 11.0)
aligncolor <- ifelse(str_detect(benchcovereddf$V1, "MAT"), "green", "blue")
kp <- plotKaryotype(genome=benchgenome, plot.type=1, chromosomes=chroms, main=plottitle, cex=0.5, plot.params=pp)
kpAddBaseNumbers(kp, tick.dist = 20000000, tick.len = 10, cex=0.3)
kpRect(kp, data=benchcoveredranges, col=aligncolor, data.panel="ideogram", y0=rep(0, length(benchcoveredranges)), y1=rep(1, length(benchcoveredranges)))

kpRect(kp, data=nlocranges, col="red", data.panel="ideogram", y0=rep(0, length(nlocranges)), y1=rep(1, length(nlocranges)))

legend("bottomright", c("Maternal", "Paternal", "Ns"), fill=c("green", "blue", "red"))

dev.off()

