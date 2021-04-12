library(MotifDb)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(VariantAnnotation)
library(ghdb)
source("~/github/fimoService/batchMode/fimoBatchTools.R")
ghdb <- GeneHancerDB()

snps <- c("rs4575098", "rs10194375", "rs73223431")
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   igv <- start.igv("PTK2B", "hg19")
   showGenomicRegion(igv, "CCDC25")     # chr8:27,730,074-27,768,330
   showGenomicRegion(igv, "PTK2B")      # chr8:27,310,478-27,460,386
   snp.loc <- 27219987 # hg19
   tbl.snp <- data.frame(chrom="chr8", start=snp.loc, end=snp.loc, stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("rs73223431", tbl.snp, color="red")
   displayTrack(igv, track)

   tbl.gh.ptk2b <- retrieveEnhancersFromDatabase(ghdb, "PTK2B", tissues="all")
   track <- DataFrameQuantitativeTrack("gh.ptk2b", tbl.gh.ptk2b[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)

   tbl.gh.ccdc25 <- retrieveEnhancersFromDatabase(ghdb, "CCDC25", tissues="all")
   dim(tbl.gh.ccdc25)
   track <- DataFrameQuantitativeTrack("gh.ccdc25", tbl.gh.ccdc25[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="darkgreen")
   displayTrack(igv, track)

   for(i in seq_len(nrow(tbl.fimo))){
      tf <- tbl.fimo$tf[i]
      track <- DataFrameQuantitativeTrack(tf, tbl.fimo[i, c("chrom", "start", "end", "score")],
                                          autoscale=FALSE, min=0, max=15, color="random")
      displayTrack(igv, track)
      } # for i

   vcf <- readVcf("ptk2b-region-hg19.vcf.gz", "hg19")
   track <- VariantTrack("rosmap", vcf)
   displayTrack(igv, track)


} # viz
#----------------------------------------------------------------------------------------------------
findBrokenMotifs <- function()
{

    snpsById(SNPlocs.Hsapiens.dbSNP150.GRCh38, snps)

#       seqnames       pos strand |   RefSNP_id alleles_as_ambig
#              1 161185602      * |   rs4575098                R
#              2 127082205      * |  rs10194375                M
#              8  27362470      * |  rs73223431                Y

    snps.gr <- snps.from.rsid(rsid = snps[3],
                              dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                              search.genome=BSgenome.Hsapiens.UCSC.hg38)

    motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core"))
    length(motifs) # 1305
    meme.file <- "tmp.meme"
    export(motifs, con=meme.file, format="meme")

    tbl.loc <- data.frame(chrom="chr8", start=27362470-15, end=27362470+15, stringsAsFactors=FALSE)
    tbl.fimo <- fimoBatch(tbl.loc, matchThreshold=1e-4, genomeName="hg38", pwmFile=meme.file)
    dim(tbl.fimo)

    motif.ids <- unique(tbl.fimo$motif_id)
    length(motif.ids)
    motifs.selected <- MotifDb[motif.ids]

    results <- motifbreakR(snpList = snps.gr,
                           filterp = TRUE,
                           pwmList = motifs.selected,
                           show.neutral=FALSE,
                           method = c("ic", "log", "notrans")[2],
                           bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                           BPPARAM = BiocParallel::bpparam(),
                           verbose=TRUE)
    length(results)
    length(results[results$effect=="strong"])   # 13
    names(results) <- NULL
    tbl.broken <- as.data.frame(results)
    colnames(tbl.broken)[1] <- "chrom"
    tbl.broken$chrom <- as.character(tbl.broken$chrom)

}
#----------------------------------------------------------------------------------------------------
vcf <- function()
{
     # on khaleesi, slice off a bit of the full vcf
   dir <- "./"
   file.exists(dir)
   file <- "ptk2b-region.vcf.gz.tbi"
   full.path <- file.path(dir, file)
   file.exists(full.path)

   vcf <- readVcf(full.path, "hg19")

   #rosmap.file <-
   regions <- GRanges(seqnames="Y", IRanges(start=13902747, end=16161938))
   vcf.file <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_Y.recalibrated_variants.annotated.vcf.gz"
   vcf <- readVcf(vcf.file, "hg19", regions)
   dim(geno(vcf)$GT)   # NULL
   head(samples(scanVcfHeader(vcf.file))) # character(0)

} # vcf
#----------------------------------------------------------------------------------------------------

