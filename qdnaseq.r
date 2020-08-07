library(QDNAseq)

args <- commandArgs(trailingOnly = T)
slurm_id <- as.numeric(args[1])
binSize <- as.numeric(args[2])
bams <- args[3]
dir <- as.character(args[4])
chr <- args[5]
out <- as.character(args[6]) 

setwd(dir)
bamfiles <- scan(bams, what = character())
file <- bamfiles[slurm_id]

load(paste0("/camp/lab/turnerj/working/resources/QDNAseq/mm10_qdnaseq_bins_", binSize, ".RData"))

tmp <- rev(strsplit(file, "/")[[1]])[1]
lib <- strsplit(tmp, "\\.")[[1]][1]

readCounts <- binReadCounts(bins, bamfiles = file)
p3 <- readCounts

readCountsFiltered <- applyFilters(readCounts)
p4 <- readCountsFiltered

readCountsCorrected <- estimateCorrection(readCountsFiltered)
p5 <- readCountsCorrected

if (chr == "T") {
	readCountsFiltered <- applyFilters(readCountsCorrected, chromosomes = NA)
} else {
	readCountsFiltered <- applyFilters(readCountsCorrected, chromosomes = c(2:19))
}


copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
p1 <- copyNumbersSegmented

copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
p2 <- copyNumbersSegmented

#idx <- chromosomes(copyNumbersSegmented) == 6
#p1.2 <- plot(copyNumbersSegmented[idx])

pdf(paste0(out, lib, "_", binSize, ".pdf"))
	plot(p1, main = "SQRT")
	plot(p2, main = "Log2") 
	plot(p3, logTransform = F, ylim = c(-100, 1500))
	if (slurm_id > 9) {
		highlightFilters(p3, logTransform = F, residual = T, blacklist = T) 
	}
	isobarPlot(p4)
	noisePlot(p5)
dev.off()

print(paste0("Written file ", lib, "_", binSize, ".pdf", " to directory ", out))
