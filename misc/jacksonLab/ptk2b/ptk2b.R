library(MotifDb)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(VariantAnnotation)
library(ghdb)
library(trena)

library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/fimoService/batchMode/fimoBatchTools.R")

#----------------------------------------------------------------------------------------------------
khaleesi.check <- function()
{
  suppressWarnings(
    db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE))
    )

 if(length(db.access.test) == 0)
   stop("khaleesi database server unavailable")

} # khaleesi.check
#----------------------------------------------------------------------------------------------------
if(!exists("ghdb")){
    khaleesi.check()
    ghdb <- GeneHancerDB()
    }

snps <- c("rs4575098", "rs10194375", "rs73223431")

if(!exists("hg38.to.hg19.chain")){
   chain.filename <- "hg38ToHg19.over.chain"
   if(!file.exists(chain.filename)){
     system("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz")
     system(sprintf("gunzip %s.gz", chain.filename))
     }
   hg38.to.hg19.chain <- import.chain(chain.filename)
   }

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_to.hg19()

} # runTests
#----------------------------------------------------------------------------------------------------
# liftover hg38 coordinates to hg19
to.hg19 <- function(tbl.hg38)
{
   if(!grepl("^chr", tbl.hg38$chrom[1])){
      #message(sprintf("prepending 'chr' to chrom column"))
      tbl.hg38$chrom <- paste0("chr", tbl.hg38$chrom)
      }

   gr.hg38 <- GRanges(tbl.hg38)
   seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
   x <- liftOver(gr.hg38, hg38.to.hg19.chain)
   gr.hg19 <- unlist(x)
   tbl.hg19 <- as.data.frame(gr.hg19)
   colnames(tbl.hg19)[1] <- "chrom"
   tbl.hg19$chrom <- as.character(tbl.hg19$chrom)

   strand.column <- grep("strand", colnames(tbl.hg19))
   if(length(strand.column) == 1)
       tbl.hg19 <- tbl.hg19[, -strand.column]

   width.column <- grep("width", colnames(tbl.hg19))
   if(length(width.column) == 1)
       tbl.hg19 <- tbl.hg19[, -width.column]

   tbl.hg19

} # to.hg19
#----------------------------------------------------------------------------------------------------
test_to.hg19 <- function()
{
   khaleesi.check()
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, "PTK2B", tissues="all")
   dim(tbl.gh)
   tbl.gh.hg19 <- to.hg19(tbl.gh)
   dim(tbl.gh.hg19)
   checkEquals(dim(tbl.gh), dim(tbl.gh.hg19))
   checkEquals(colnames(tbl.gh), colnames(tbl.gh.hg19))
   # checkEquals(tbl.gh$chrom, tbl.gh.hg19$chrom)
   checkTrue(length(intersect(tbl.gh$start, tbl.gh.hg19$start)) == 0)
   checkTrue(length(intersect(tbl.gh$end, tbl.gh.hg19$end)) == 0)

} # test_to.hg19
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   igv <- start.igv("PTK2B", "hg19")
   showGenomicRegion(igv, "CCDC25")     # chr8:27,730,074-27,768,330
   showGenomicRegion(igv, "PTK2B")      # chr8:27,310,478-27,460,386
   snp.loc <- 27219987 # hg19
   tbl.snp <- data.frame(chrom="chr8", start=snp.loc-1, end=snp.loc, stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("rs73223431", tbl.snp, color="red")
   displayTrack(igv, track)

   tbl.gh.ptk2b <- retrieveEnhancersFromDatabase(ghdb, "PTK2B", tissues="all")
   tbl.gh.ptk2b.hg19 <- to.hg19(tbl.gh.ptk2b)

   track <- DataFrameQuantitativeTrack("gh.ptk2b", tbl.gh.ptk2b.hg19[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)

   tbl.gh.ccdc25 <- retrieveEnhancersFromDatabase(ghdb, "CCDC25", tissues="all")
   tbl.gh.ccdc25.hg19 <- to.hg19(tbl.gh.ccdc25)
   track <- DataFrameQuantitativeTrack("gh.ccdc25", tbl.gh.ccdc25.hg19[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="darkgreen")
   displayTrack(igv, track)

   tbl.stamTF <- get(load("~/github/geneReg-2021/stam-dhs/tbl.stamTFs-hg19.RData"))
   tbl.stam <- subset(tbl.stamTF, chrom=="chr8" & start > 26043685 & end < 28915481)
   dim(tbl.stam)
   track <- DataFrameAnnotationTrack("stam", tbl.stam[, c("chrom", "start", "end", "tf")],
                                     color="darkgreen")
   displayTrack(igv, track)

   for(i in seq_len(nrow(tbl.fimo))){
      tf <- tbl.fimo$tf[i]
      track <- DataFrameQuantitativeTrack(tf, tbl.fimo[i, c("chrom", "start", "end", "score")],
                                          autoscale=FALSE, min=0, max=15, color="random")
      displayTrack(igv, track)
      } # for i

   vcf <- readVcf("ptk2b-236k-region-hg19.vcf.gz", "hg19")
   track <- VariantTrack("rosmap", vcf)
   displayTrack(igv, track)
   length(vcf)

} # viz
#----------------------------------------------------------------------------------------------------
buildTrenaModel <- function()
{
  meme.file <- "jaspar2018-hocomocoCore.meme"
  motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core"))
  export(motifs, con=meme.file, format="meme")
  mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))

  dim(tbl.gh.hg19)

  tbl.fimo <- fimoBatch(tbl.gh.hg19[, c("chrom", "start", "end")],
                        matchThreshold=1e-3, genomeName="hg19", pwmFile=meme.file)
  dim(tbl.fimo)
  head(tbl.fimo)
  sort(table(tbl.fimo$tf))
  tfs <- unique(tbl.fimo$tf)
  length(tfs)  # 290
  tfs.oi <- intersect(tfs, rownames(mtx))
  length(tfs.oi) # 160
  targetGene <- "PTK2B"
  targetGene %in% rownames(mtx)
  cors <- lapply(tfs.oi, function(tf) cor(mtx[targetGene,], mtx[tf,]))
  names(cors) <- tfs.oi
  fivenum(as.numeric(cors))
  tfs.final <- names(cors)[which(abs(as.numeric(cors)) > 0.25)]
  length(tfs.final) # 56

  solver <- EnsembleSolver(mtx,
                           targetGene=targetGene,
                           candidateRegulators=tfs.final,
                           solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"))
  tbl.out <- run(solver)
  dim(tbl.out)
  tbl.out <- tbl.out[order(abs(tbl.out$spearmanCoeff), decreasing=TRUE),]
    #   gene    betaLasso   betaRidge spearmanCoeff pearsonCoeff   rfScore     xgboost
    #  MEF2C  0.095804097  0.05861931     0.5164458    0.5386553 58.367362 0.166214223
    #   RARA  0.271161463  0.15000484     0.5151084    0.5168680 83.738846 0.192812221
    # PKNOX2  0.131712422  0.08745204     0.4971322    0.5265822 40.804790 0.053275951
    #  VEZF1 -0.085680775 -0.06477160    -0.4900935   -0.4909379 37.565401 0.056665307
    #   KLF3  0.000000000 -0.02421796    -0.4845785   -0.5093787 20.458958 0.006273603
    #    SP3 -0.057353515 -0.05351812    -0.4596344   -0.5055334 26.437693 0.047660863
    #   KLF5  0.052558011  0.07046639     0.3733406    0.3931583 27.686257 0.047774154
    #    MAF -0.008169798 -0.06641919    -0.2927013   -0.3026912  6.982959 0.015473276
    # ZBTB7A  0.044219125  0.07955861     0.2655768    0.3223686 13.348027 0.023715614
    #   HSF4  0.115512098  0.08937000     0.2265415    0.2621282 16.290378 0.036750569

  head(tbl.fimo)

  for(gene in tbl.out$gene){
      tbl <- subset(tbl.fimo, tf==gene)[, c("chrom", "start", "end")]
      track <- DataFrameAnnotationTrack(gene, tbl, color="random", trackHeight=25)
      displayTrack(igv, track)
      }


} # buildTrenaModel
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

    snps.gr <- snps.from.rsid(rsid = c(snps[3],"rs17057043"),
                              dbSNP=SNPlocs.Hsapiens.dbSNP144.GRCh37,
                              search.genome=BSgenome.Hsapiens.UCSC.hg19)


    motifs.selected <- query(MotifDb, c("sapiens", "RARA"), c("jaspar2018", "HOCOMOCOv11-core"))
    length(motifs) # 1305
    meme.file <- "tmp.meme"
    export(motifs, con=meme.file, format="meme")

    #tbl.loc <- data.frame(chrom="chr8", start=27362470-15, end=27362470+15, stringsAsFactors=FALSE)
    #tbl.fimo <- fimoBatch(tbl.loc, matchThreshold=1e-4, genomeName="hg38", pwmFile=meme.file)
    #dim(tbl.fimo)

    #motif.ids <- unique(tbl.fimo$motif_id)
    #length(motif.ids)
    #motifs.selected <- MotifDb[motif.ids]

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

