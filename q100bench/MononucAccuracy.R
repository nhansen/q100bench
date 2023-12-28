#setwd("/Users/nhansen/HG002_diploid_benchmark/plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
#BiocManager::install("karyoploteR")

args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "ash1")
outputdir <- ifelse(!( is.na(args[2])), args[2], ".")
resourcedir <- ifelse(! ( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c(genomename, " accuracy of mononucleotide runs"), sep="", collapse=""))

mononucsitefile <- paste(c(outputdir, "/", genomename, ".mononucstats.txt"), sep="", collapse="")
mononucsites <- read.table(mononucsitefile, header=FALSE, sep="\t")
names(mononucsites) <- c("name", "base", "reflength", "assemblylength", "type")

consensuserrors <- mononucsites[(mononucsites$assemblylength != -1) & (mononucsites$type == "CONSENSUS"), ]
noncomplexcovered <- mononucsites[mononucsites$assemblylength != -1, ]

mononucerrorrate <- length(consensuserrors$reflength)/length(noncomplexcovered$reflength)
mononucerrorperc <- round(mononucerrorrate*100, 2)

consensuserrorcounts <- hist(consensuserrors$reflength, plot=FALSE, breaks=seq(10, 100, 1))
noncomplexcovcounts <- hist(noncomplexcovered$reflength, plot=FALSE, breaks=seq(10, 100, 1))

accrate <- 1.0 - consensuserrorcounts$counts/noncomplexcovcounts$counts

plotname <- paste(c(outputdir, "/", genomename, ".mononuc_accuracy.v1.0.1.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)

plot(consensuserrorcounts$mids, accrate, xlim=c(10,40), pch=16, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=paste(c("Accuracy of mononucleotide runs in ", genomename), sep="", collapse=""))
text(20, 0.5, labels= paste(c("Overall error rate: ", mononucerrorperc, "%"), sep="", collapse=""))

dev.off()

#consensuserrorlengths <- consensuserrors$assemblylength - consensuserrors$reflength
#consensuslengthdiffs <- consensuserrors$altlength - consensuserrors$reflength

#v0.8length_hist_notail <- v0.8length_hist[v0.8length_hist$V1>=-4 & v0.8length_hist$V1<=4,]
#hprcyear1length_hist_notail <- hprcyear1length_hist[hprcyear1length_hist$V1>=-4 & hprcyear1length_hist$V1<=4,]
#hprccuratedlength_hist_notail <- hprccuratedlength_hist[hprccuratedlength_hist$V1>=-4 & hprccuratedlength_hist$V1<=4,]

#barplot(lengthdiffs, names.arg=v0.8length_hist_notail$V1, main="Assembly error lengths", xlab="length(v0.9)-length(assembly)", ylab="Number of errors")

#barplot(rbind(v0.8length_hist_notail$V2, hprcyear1length_hist_notail$V2, hprccuratedlength_hist_notail$V2), main="Assembly error lengths", beside=TRUE, names.arg=v0.8length_hist_notail$V1, xlab="Error length (bases)", ylab="# Errors", col=c("brown", "blue", "darkgreen"), cex.lab=1.1, cex.axis=1.1)
#legend("topright", c("v0.8 Verkko", "HPRC Year1", "Jarvis et al"), col=c("brown", "blue", "darkgreen"), pch=15, cex=1.2)

