setwd("F:/mehdi-ada")
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


txdb<-makeTxDbFromBiomart(biomart ="fungi_mart" ,dataset="spombe_eg_gene" ,host="fungi.ensembl.org")
columns(txdb)
keytypes(txdb)

samplefiles <- list.files(".", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
print(samplefiles)

peak <- readPeakFile(samplefiles[[2]])
peak

covplot(peak, weightCol="V5")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
#data("tagMatrixList")
#tagMatrix <- tagMatrixList[[4]]
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
peakHeatmap(samplefiles[[2]], TxDb=txdb, upstream=3000, downstream=3000, color="red")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

tagMatrix_binning <- getTagMatrix(peak = peak, TxDb = txdb, 
                                  upstream = 3000, downstream = 3000, 
                                  type = "start_site", by = "gene", 
                                  weightCol = "V5", nbin = 800)


plotPeakProf2(peak = peak, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)


genebody <- getBioRegion(TxDb = txdb, by = "gene", type = "body")
matrix_no_flankextension <- getTagMatrix(peak,windows = genebody, nbin = 800)

plotPeakProf(matrix_no_flankextension,conf = 0.95)

matrix_actual_extension <- getTagMatrix(peak,windows = genebody, nbin = 800,
                                        upstream = 1000,downstream = 1000)
plotPeakProf(matrix_actual_extension,conf = 0.95)

five_UTR_body <- getTagMatrix(peak = peak, TxDb = txdb, upstream = rel(0.2),
                              downstream = rel(0.2), type = "body",
                              by = "5UTR", weightCol = "V5",
                              nbin = 50)

plotPeakProf(tagMatrix = five_UTR_body, conf = 0.95)


TTS_matrix <- getTagMatrix(peak = peak, 
                           TxDb = txdb,
                           upstream = 3000,
                           downstream = 3000, 
                           type = "end_site",
                           by = "gene",
                           weightCol = "V5")
plotPeakProf(tagMatrix = TTS_matrix, conf = 0.95)

peakAnno <- annotatePeak(samplefiles[[2]], tssRegion=c(-3000, 3000), TxDb=txdb)

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)

plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotPeakProf2(samplefiles, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row")

plotPeakProf2(samplefiles, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row", nbin = 800)

plotPeakProf2(samplefiles, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

getGEOspecies()
getGEOgenomeVersion()

#txdb1 <- makeTranscriptDbFromGTF("Schizosaccharomyces_pombe.ASM294v2.22.gtf", format="gtf")
#??makeTranscriptDbFromGFF
#link <- "ftp://ftp.ensemblgenomes.org/pub/release-34/fungi/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.34.gff3.gz" 
#txdb <- makeTxDbFromGFF(link, format = "gff3", organism = "Schizosaccharomyces pombe", taxonomyId = "4896")
#txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", dataSource = "TAIR", organism = "Arabidopsis thaliana")
#files <- list.files(path="./", pattern="*!!.narrowPeak$")
#files <- files[1:13]
#files
#peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000), verbose=FALSE)
#peak <- readPeakFile("C2_peaks.narrowPeak", header=FALSE)
#peakAnno <- annotatePeak(C2_peaks.sorted.narrowPeak, tssRegion = c(-500, 2000), TxDb=txdb)
#print(peakAnno)
#write.table(peakAnno, paste("annot_", i, sep="") , sep="\t", col.names=T, row.names = F, quote=F)
#peakAnnoList <- lapply(C1.df, annotatePeak, TxDb=txdb, verbose=FALSE)
#keytypes(txdb)
#?annotatePeak
