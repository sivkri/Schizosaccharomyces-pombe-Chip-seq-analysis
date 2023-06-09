# Import libraries
library(GenomicFeatures)
library(biomaRt)
library(ChIPseeker)
library(GenomicRanges)
library(AnnotationDbi)
library(clusterProfiler)
library(rtracklayer)
library(ggplot2)
library(ChIPpeakAnno)
library(ggupset)
library(ggimage)
library(ReactomePA)
library(DOSE)
library(meshes)

# Import the annotation from biomart
txdb <- makeTxDbFromBiomart(biomart = "fungi_mart", dataset = "spombe_eg_gene", host = "fungi.ensembl.org")

# Get information about the TxDb object
columns(txdb)
keytypes(txdb)

# Import the files from the directory
samplefiles <- list.files(".", pattern = ".bed", full.names = TRUE)

# Read peak file
peak <- readPeakFile(samplefiles[[2]])
peak

# Plot coverage
covplot(peak, weightCol = "V5")

# Get promoters and create tag matrix
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peak, windows = promoter)

# Plot tag heatmap
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = "red")

# Plot peak heatmap
peakHeatmap(samplefiles[[2]], TxDb = txdb, upstream = 3000, downstream = 3000, color = "red")

# Plot average profile
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")

# Plot average profile with confidence interval and resampling
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95, resample = 1000)

# Get tag matrix with binning
tagMatrix_binning <- getTagMatrix(peak = peak, TxDb = txdb, upstream = 3000, downstream = 3000, type = "start_site", by = "gene", weightCol = "V5", nbin = 800)

# Plot peak profile by gene body with binning
plotPeakProf2(peak = peak, upstream = rel(0.2), downstream = rel(0.2), conf = 0.95, by = "gene", type = "body", nbin = 800, TxDb = txdb, weightCol = "V5", ignore_strand = FALSE)

# Get tag matrix for gene body without flank extension
genebody <- getBioRegion(TxDb = txdb, by = "gene", type = "body")
matrix_no_flankextension <- getTagMatrix(peak, windows = genebody, nbin = 800)

# Plot peak profile without flank extension
plotPeakProf(matrix_no_flankextension, conf = 0.95)

# Get tag matrix for gene body with actual extension
matrix_actual_extension <- getTagMatrix(peak, windows = genebody, nbin = 800, upstream = 1000, downstream = 1000)
plotPeakProf(matrix_actual_extension, conf = 0.95)

# Get tag matrix for 5' UTR body
five_UTR_body <- getTagMatrix(peak = peak, TxDb = txdb, upstream = rel(0.2), downstream = rel(0.2), type = "body", by = "5UTR", weightCol = "V5", nbin = 50)

# Plot peak profile for 5' UTR body
plotPeakProf(tagMatrix = five_UTR_body, conf = 0.95)

# Get tag matrix for transcription termination site (TTS)
TTS_matrix <- getTagMatrix(peak = peak, TxDb = txdb, upstream = 3000, downstream = 3000, type = "end_site", by = "gene", weightCol = "V5")

# Plot peak profile for TTS
plotPeakProf(tagMatrix = TTS_matrix, conf = 0.95)

# Annotate peak
peakAnno <- annotatePeak(samplefiles[[2]], tssRegion = c(-3000, 3000), TxDb = txdb)

# Plot annotation pie chart and bar plot
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)

# Plot venn diagram pie chart
vennpie(peakAnno)

# Plot upset plot
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie = TRUE)

# Plot distribution of transcription factor-binding loci relative to TSS
plotDistToTSS(peakAnno, title = "Distribution of transcription factor-binding loci\nrelative to TSS")

# Plot peak profiles by gene start site (facet row)
plotPeakProf2(samplefiles, upstream = 3000, downstream = 3000, conf = 0.95, by = "gene", type = "start_site", TxDb = txdb, facet = "row")
plotPeakProf2(samplefiles, upstream = 3000, downstream = 3000, conf = 0.95, by = "gene", type = "start_site", TxDb = txdb, facet = "row", nbin = 800)

# Plot peak profiles by gene body (facet row)
plotPeakProf2(samplefiles, upstream = rel(0.2), downstream = rel(0.2), conf = 0.95, by = "gene", type = "body", TxDb = txdb, facet = "row", nbin = 800)

# Annotate peaks for each sample file
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)

# Plot annotation bar plot for multiple samples
plotAnnoBar(peakAnnoList)

# Plot distribution of transcription factor-binding loci for multiple samples
plotDistToTSS(peakAnnoList)

# Extract gene IDs from peak annotation list
genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

# Plot venn diagram for gene IDs
vennplot(genes)
