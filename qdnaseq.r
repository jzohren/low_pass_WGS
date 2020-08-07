library(QDNAseq)
library(glue)

args <- commandArgs(trailingOnly = T)
sample <- args[1]
bin_size <- as.numeric(args[2])
out_dir <- args[3]

setwd(out_dir)
bam <- glue("{sample}_merged.bam")

bins <- getBinAnnotations(binSize = bin_size)

readCounts <- binReadCounts(bins, bamfiles = bam)
p3 <- readCounts

readCountsFiltered <- applyFilters(readCounts)
p4 <- readCountsFiltered

readCountsCorrected <- estimateCorrection(readCountsFiltered)
p5 <- readCountsCorrected

readCountsFiltered <- applyFilters(readCountsCorrected, chromosomes = NA)

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
p1 <- copyNumbersSegmented

copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
p2 <- copyNumbersSegmented


pdf(glue("{out_dir}/{sample}_{bin_size}.pdf")
	plot(p1, main = "SQRT")
	plot(p2, main = "Log2") 
	plot(p3, logTransform = F)
	isobarPlot(p4)
	noisePlot(p5)
dev.off()
	
