## create Arch object /data/ArchR-Proj
library(ArchR)
library(parallel)

inputFiles <- c("~/data/500K_HT_pseudo_fragments.tsv.gz")
names(inputFiles) <- c("HT_pseudo")
addArchRGenome("hg38")
addArchRThreads(threads = 4)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 1, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "~/Save_HT_pseudo",
  copyArrows = TRUE
)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    n.start = 5),
  varFeatures = 25000, 
  dimsToUse = 1:5
)

proj <- addClusters(proj, reducedDims = "IterativeLSI",
                    method = "Seurat", name = "Clusters")
n <- nCells(proj)
proj$Clusters <- factor(c(rep("C1", ceiling(n/2)), rep("C2", floor(n/2))))
coverageFiles <- getGroupBW(ArchRProj = proj, groupBy = "Clusters")

proj <- addUMAP(proj, reducedDims = "IterativeLSI",  
                nNeighbors = 2,
                force = TRUE)
#######################################################
valid_cells <- intersect(proj$cellNames, rownames(proj@embeddings$UMAP$df))
proj@cellColData <- proj@cellColData[valid_cells, , drop = FALSE]
#######################################################

proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = "~/Save-HT_pseudo", load = TRUE)

library(Seurat)
library(ShinyCell2)
#archr_obj <- loadArchRProject("~/Save-HT_pseudo")
scConf2 <- createConfig(proj)
# Columns you want to delete
to_del <- c("ReadsInTSS", "ReadsInPromoter", "ReadsInBlacklist",
            "NucleosomeRatio", "nMultiFrags", "nMonoFrags", "nFrags",
            "nDiFrags", "ReadsInPeaks")
# Keep only those that exist in your config
to_del_existing <- to_del[to_del %in% colnames(scConf2$meta)]
scConf2 <- delMeta(scConf2, to_del_existing)

### this failed
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Clusters",
  peakMethod = "Tiles",
  method = "p",
  minCells = 2,
  cutOff = 0.01,
)

## Go with
### psesudo peaks
featureDF <- getFeatures(proj, useMatrix = "TileMatrix")
myPeaks <- GRanges(
  seqnames = featureDF@seqnames,
  ranges   = IRanges(start = featureDF@ranges@start, 
                     end = featureDF@ranges@start + 499),
  name     = featureDF$idx
)

proj <- addPeakSet(proj, peakSet = myPeaks)
proj <- addPeakMatrix(proj)

makeShinyFiles(proj, scConf = scConf2, dimred.to.use = "UMAP", 
               bigWigGroup = c("Clusters"), shiny.prefix = "scATAC",
               shiny.dir = "~/Shiny_HT_pseudo/", default.gene1 = "IRF1", 
               default.multigene = NA, default.dimred = "UMAP", chunkSize =  500)

makeShinyCodes(shiny.title="DEMO of Heart Data: scATAC shinycell2",shiny.prefix=c("scATAC"),
               shiny.headers = c("scATAC-seq Data"),
               shiny.dir="~/Shiny_HT_pseudo/")



#proj$cellColData <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
### meaningless for HT_pseudo
#p1 <- plotFragmentSizes(ArchRProj = projHeme1)
#p2 <- plotTSSEnrichment(ArchRProj = projHeme1)

projHeme1 <- proj
paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

p1 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges",
  baseSize = 10
)
p1
p2 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)
p3 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  baseSize = 10
)
p4 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)






write it as a standalone Rscript, that took input parameters as input file
for example: inputFiles <- c("~/HT_pseudos.tsv.gz")
and save output as 
~/ShinyCell2_HT_pseudo
