#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readaligns")

args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "illumina2x250mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.0.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("Homopolymer length concordance for ", readsetname, " vs ", genomename), sep="", collapse=""))

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="")
  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  mndf <- mndf[mndf$readlength != -1, ]
  mndf$reflength <- as.integer(mndf$reflength)
  mndf$readlength <- as.integer(mndf$readlength)
  mndf$numcorrect <- as.integer(mndf$numcorrect)
  mndf$numerror <- as.integer(mndf$numerror)
  
  return(mndf)
}

readmononuccounts <- function(filename) {
  mndf <- read.table(filename, header=FALSE, sep="")
  names(mndf) <- c("count", "base", "reflength", "readlength", "type")
  mndf <- mndf[mndf$readlength != "COMPLEX", ]
  
  return(mndf)
}

mndiffs <- function(filename, i) {
  mndiff <- readmononuchist(filename);
  
  if (i==0) {
    diffcount <- sum(mndiff[mndiff$readlength==mndiff$reflength, "numcorrect"])
  }
  else {
    diffcount <- sum(mndiff[mndiff$readlength-mndiff$reflength==i, "numerror"])
  }
  return(diffcount)
}

onebaseindelratio <- function(filename) {
  onebaseratio <- mndiffs(filename, 1)*1.0/mndiffs(filename, -1)
  
  return(onebaseratio)
}

diffarray <- function(filename, maxdiff=4) {
  platformdiffs <- sapply(seq(-1*maxdiff, maxdiff), function(i) {mndiffs(filename, i)})
  platformdiffs <- platformdiffs/sum(platformdiffs)
  
  return(platformdiffs)
}

makemultiplot <- function(maxdiff=4) {
  plotcolors <- c("blue", "lightblue", "lightgreen", "darkgreen", "violet", "darkviolet")
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  barplot(diffarray("sequelII_15k.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[1], cex.names=1.0, main="SequelII_15k", xlab="Difference from benchmark", ylab="Fraction of reads")
  barplot(diffarray("sequelII_20kb_20210812.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[2], cex.names=1.0, main="SequelII_20kb_20210812", xlab="Difference from benchmark")
  barplot(diffarray("hifidcv1.1.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[3], cex.names=1.0, main="HiFi/DCv1.1", xlab="Difference from benchmark")
  barplot(diffarray("hifirevio3cells.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[4], cex.names=1.0, main="HiFi Revio", xlab="Difference from benchmark")
  barplot(diffarray("illumina2x250mat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[5], cex.names=0.7, main="Illumina 2x250 (maternal)", xlab="Difference from benchmark")
  barplot(diffarray("element_avitilongmat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[6], cex.names=0.7, main="Element Aviti (maternal)", xlab="Difference from benchmark")
  mtext("Homopolymer length concordance", side = 3, line = 1, outer = TRUE)
}

plotmononuchist <- function(file, plottitle="", maxdiff=4){
  differences <- diffarray(file, maxdiff=maxdiff)
  accrate <- as.integer(differences[maxdiff+1]*1000)/10
  barplot(differences, names=seq(-1*maxdiff,maxdiff), col="blue", cex.names=1.0, main=plottitle, xlab="Difference from benchmark", ylab="Fraction of reads")
  mtext(paste(c("Accuracy: ", as.character(accrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-4)
  indelratio <- as.integer(onebaseindelratio(file)*1000)/1000
  mtext(paste(c("1bp Ins/Dels: ", as.character(indelratio)), sep="", collapse=""), side=3, adj=0.85, line=-5 )
}

mononucfile <- paste(c(outputdir, "/", readsetname, ".mononuchist.txt"), sep="", collapse="")
plotname <- paste(c(outputdir, "/", readsetname, ".mononucsizehist.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
plotmononuchist(mononucfile, plottitle, 4)
dev.off()

