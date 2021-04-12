library(SKAT)
library(RUnit)
library(VariantAnnotation)
#----------------------------------------------------------------------------------------------------
tbl.cov <- get(load("~/github/kats/podkat/rosmap/tbl.covAll-rosmap.RData"))
rownames(tbl.cov) <- tbl.cov$individualID
cogdx.simple <- rep(0, nrow(tbl.cov))
pure.ad <- which(tbl.cov$cogdx == 4)
cogdx.simple[pure.ad] <- 10
tbl.cov$cogdx.simple <- cogdx.simple

vcfFile <- "ptk2b-region-hg19.vcf.gz"
chrom.loc <-  "8"
   # min and max from vcf
start.loc <- 27215027
end.loc   <- 27224807
snps.rs <- c("rs73223431", #  chr8 27219987
             "rs17057043")  #: chr8 27220310
ad.snps <- c("8:27219987_C/T", "8:27220310_G/A")

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createGeno()
   test_skatifyGenotypeMatrix()

} # runTests
#----------------------------------------------------------------------------------------------------
createGeno <- function(vcf.region=GRanges(seqnames=chrom.loc, IRanges(start=start.loc, end=end.loc)))
{
   vcf <- readVcf(vcfFile, "hg19") #, vcf.region)
   mtx.geno <- t(geno(vcf)$GT)


  #      readGenotypeMatrix(file, regions, subset,
  #                         noIndels=TRUE, onlyPass=TRUE,
  #                         na.limit=1, MAF.limit=1,
  #                         na.action=c("impute.major", "omit", "fail"),
  #                         MAF.action=c("invert", "omit", "ignore", "fail"),
  #                         sex=NULL)

   keepers <- grep("SM-", rownames(mtx.geno))
   mtx.geno <- mtx.geno[keepers,]

   name.match.indices <-  match(rownames(mtx.geno), tbl.cov$specimenID)
   na.matches <- which(is.na(name.match.indices))
   length(na.matches)
   mtx.geno <- mtx.geno[-na.matches,]
   name.match.indices <-  match(rownames(mtx.geno), tbl.cov$specimenID)
   stopifnot(length(which(is.na(name.match.indices))) == 0)

   rownames(mtx.geno) <- tbl.cov$individualID[name.match.indices]
   mtx.geno

} # createGeno
#----------------------------------------------------------------------------------------------------
test_createGeno <- function()
{
    message(sprintf("--- test_createGeno"))

    mtx.geno <- createGeno()
    checkTrue(is.matrix(mtx.geno))
    # khaleesi, 2kb checkEquals(dim(mtx.geno), c(1144, 32))
    # hagfish, 100kb, 2288
    colnames(mtx.geno)
    checkTrue(all(grepl("8:", colnames(mtx.geno))))
    checkTrue(all(grepl("^R", rownames(mtx.geno))))
      # get just the AD-related positions, hand check the variants

      # AD vcf's in hg19:
      #   rs429358   19:45411941 (GRCh37)   ref: T
      #   rs7412     19:45412079 (GRCh37)   ref: C
    all(ad.snps %in% colnames(mtx.geno))
    mtx.adSnps <- mtx.geno[, ad.snps]
                                                          #  0/0  0/1   1/1

    checkEquals(as.numeric(table(mtx.adSnps[, ad.snps[1]])), c(487, 514, 143))
    checkEquals(as.numeric(table(mtx.adSnps[, ad.snps[2]])), c(485, 516, 143))

       # pick a double snp, make sure the patient is in tbl.cov
    pts.double <- sort(intersect(names(which(mtx.adSnps[,1] == "0/1")),
                                 names(which(mtx.adSnps[,2] == "0/1"))))
    checkEquals(length(pts.double), 514)
    pt.double.first <- head(pts.double, n=1)
    checkEquals(pt.double.first, "R1015854")
    sampleID <- subset(tbl.covAll, individualID==pt.double.first)$specimenID
    checkEquals(sampleID, "SM-CTDQP")

       # go back to the tradtional VariantAnnotation
    # checkTrue(sampleID %in% vcf.sampleIDs)
    # mtx.orig <- geno(vcf)$GT
    # khaleesi, 2kb: checkEquals(dim(mtx.orig), c(34, 1894))
    # hagfish, 100kb: checkEquals(dim(mtx.orig), c(2754, 1894))
    # roi <- unlist(lapply(ad.snps, function(loc) grep(loc, rownames(mtx.orig))))
    # checkEquals(as.character(mtx.orig[roi, sampleID]), c("0/1", "0/1"))


} # test_createGeno
#----------------------------------------------------------------------------------------------------
# create a numeric genotype matrix with each row as a different individual and each column
# as a separate gene/snp.  Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA,
# Aa, aa, and missing, where A is a major allele and a is a minor allele.  Missing
# genotypes will be imputed by the simple Hardy-Weinberg equilibrium (HWE) based
# imputation.
skatifyGenotypeMatrix <- function(mtx.raw)
{
   #stopifnot(all(rownames(mtx.raw) %in% tbl.cov$individualID))

   noReads <- which(mtx.raw == "./.")
   AA <- which(mtx.raw == "0/0")
   aa <- which(mtx.raw == "1/1")
   Aa <- grep("0/[1-9]", mtx.raw)
   head(setdiff(seq_len(length(mtx.raw)), c(AA, noReads, aa, Aa)))
      # see this discussion of unusual genotypes in rosmap vcf file
      # ~/github/notes/log
      # --- skat vs actual genotype from rosmap: 98% fit, 2% do not, as
      #  discovered in ~/github/SKAT/inst/unitTests/test_rosmapSkat.R

   #stopifnot(sum(length(AA),
   #              length(noReads),
   #              length(aa),
   #              length(Aa)) > 0.98 * length(mtx.raw))
   # wt <- c(hets, noReads)
   # variants <- setdiff(seq_len(length(mtx.raw)), wt)

   mtx.out <- matrix(0, nrow=nrow(mtx.raw), ncol=ncol(mtx.raw))
   mtx.out[noReads] <- 9
   mtx.out[AA] <- 0
   mtx.out[Aa] <- 1
   mtx.out[aa] <- 2
   rownames(mtx.out) <- rownames(mtx.raw)
   colnames(mtx.out) <- colnames(mtx.raw)

   #patientIDs <- tbl.cov$patientID[match(colnames(mtx.out), tbl.cov$specimen)]
   #browser()
   #stopifnot(length(patientIDs) == ncol(mtx.out))
   #colnames(mtx.out) <- patientIDs

   invisible(mtx.out)

} # skatifyGenotypeMatrix
#----------------------------------------------------------------------------------------------------
test_skatifyGenotypeMatrix <- function()
{
   message(sprintf("--- test_skatifyGenotypeMatrix"))

   mtx.geno <- createGeno()
   checkEquals(dim(mtx.geno), c(1144, 203))
   sample.rows <- c("R9880904", "R9677385", "R1236313", "R9794121", "R6692433", "R8261694")
   sample.cols <- c("8:27219865_G/A",  "8:27219987_C/T", "8:27220128_C/A")
   mtx.test <- mtx.geno[sample.rows, sample.cols]

   mtx.test [2,3] <- "./."   # inject a no-reads value


                                             # ./. 0/0 0/1 1/1
   checkEquals(as.numeric(table(mtx.test)), c(1, 12,  3, 2))

   mtx <- skatifyGenotypeMatrix(mtx.test)
                                          #  0  1  2  9
   checkEquals(as.numeric(table(mtx)), c(12, 3, 2, 1))

   checkEquals(dim(mtx), dim(mtx.test))
   checkEquals(rownames(mtx), sample.rows)
   checkEquals(colnames(mtx), sample.cols)

   mtx.test <- mtx.geno[1000:1020, 1:34]
   checkEquals(dim(mtx.test), c(21, 34))
   mtx <- skatifyGenotypeMatrix(mtx.test)
   checkEquals(dim(mtx), dim(mtx.test))
   checkEquals(as.numeric(table(mtx)), c(712, 2)) # AA, Aa

      # now the entire mtx.geno
   mtx <- skatifyGenotypeMatrix(mtx.geno)
   checkEquals(dim(mtx), dim(mtx.geno))
                                      #       0     1     2  9
   checkEquals(as.numeric(table(mtx)), c(226450, 3847, 1933, 2))

   checkEquals(head(colnames(mtx), n=3), c("8:27215027_C/T", "8:27215044_C/A", "8:27215048_A/G"))
   checkEquals(head(rownames(mtx), n=3), c("R1977848", "R6478102", "R4567280"))

} # test_skatifyGenotypeMatrix
#----------------------------------------------------------------------------------------------------
runSingleSkato <- function()
{
  mtx.geno <- createGeno()
  dim(mtx.geno)   # 1144 203
  head(rownames(mtx.geno))
  checkTrue(all(rownames(mtx.geno) %in% tbl.cov$individualID))

  mtx.skat <- skatifyGenotypeMatrix(mtx.geno)

  null.model <- SKAT_Null_Model(cogdx ~ age_death + msex, tbl.cov)
  x <- SKAT(mtx.skat, null.model, method="SKATO")
  x$p.value # [1] 0.02725835


} # runSingleSkato
#----------------------------------------------------------------------------------------------------
to.chrom.loc <- function(col.title)
{
   tokens <- strsplit(col.title, ":|_")[[1]]
   list(chrom=tokens[1], pos=as.integer(tokens[2]))
}
#----------------------------------------------------------------------------------------------------
runBlockedSkato <- function()
{
   mtx.geno <- createGeno()

   count <- 10

   starts <- as.integer(seq(1, ncol(mtx.geno), length.out=count+1))
   ends <- starts-1
   starts <- starts[-length(starts)]
   ends <- ends[-1]
   ends[length(ends)] <- ncol(mtx.geno)
   stopifnot(sum(ends - starts + 1) == ncol(mtx.geno))

   tbls.skat <- list()

   for(i in seq_len(length(starts))){
      start <- starts[i]
      end   <- ends[i]
      coi <- colnames(mtx.geno)[start:end]
      mtx.skat <- skatifyGenotypeMatrix(mtx.geno[, coi])
      mtx.skat <- mtx.skat[rownames(tbl.cov),]
      null.model <- SKAT_Null_Model(cogdx.simple ~ age_death + msex, tbl.cov)
      x <- SKAT(mtx.skat, null.model, method="SKATO")
      chrom.loc.start <- to.chrom.loc(colnames(mtx.geno)[start])
      chrom.loc.end   <- to.chrom.loc(colnames(mtx.geno)[end])
      score <- -log10(x$p.value)
      tbls.skat[[i]] <- data.frame(chrom=chrom.loc.start$chrom,
                                   start=chrom.loc.start$pos,
                                   end=chrom.loc.end$pos,
                                   score = score,
                                   pval = x$p.value,
                                   stringsAsFactors=FALSE)
      if(x$p.value < 0.05) printf("%d - %d: %f", start, end, x$p.value)
      }

   tbl.skat <- do.call(rbind, tbls.skat)
   dim(tbl.skat)

} # runBlockedSkato
#----------------------------------------------------------------------------------------------------
#   cogdx: final clinical consensus diagnosis, blinded to all postmortem data
#          value coding
#          1     NCI, No cognitive impairment (No impaired domains)
#          2     MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI
#          3     MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI
#          4     AD, Alzheimer's disease and NO other cause of CI (NINCDS PROB AD)
#          5     AD, Alzheimer's disease AND another cause of CI (NINCDS POSS AD)
#          6     Other dementia. Other primary cause of dementia
runHandCraftedRegionsSkato <- function()
{
   cogdx.fake <- rep(0, nrow(tbl.cov))
   cogdx.fake[which(mtx.skat != 0)] <- 1
   tbl.cov$cogdx.fake <- cogdx.fake

   cogdx.lumped <- rep(0, nrow(tbl.cov))
   ad <- which(tbl.cov$cogdx > 1)
   length(ad)
   cogdx.lumped[ad] <- 1
   tbl.cov$cogdx.lumped <- cogdx.lumped

   tbl.cov$cogdx.fake
   null.model <- SKAT_Null_Model(cogdx.lumped ~ age_death + msex, tbl.cov, out_type="C")

   snp.01 <- "chr8:27,219,964-27,220,023"
   snp.02 <- "chr8:27,220,308-27,220,311"
   snps.rs <- c("rs73223431", #  chr8 27219987
                "rs17057043")  #: chr8 27220310
   ad.snps <- c("8:27219987_C/T",   # col 104
                "8:27220310_G/A")   # col 111
   coi <- 104:111
   mtx.skat <- skatifyGenotypeMatrix(mtx.geno[, c(104), drop=FALSE]) # , 111)])
   mtx.skat <- skatifyGenotypeMatrix(mtx.geno[, c(99), drop=FALSE]) # , 111)])
   mtx.skat <- skatifyGenotypeMatrix(mtx.geno[, c(95), drop=FALSE]) # , 111)])
   x <- SKAT(mtx.skat, null.model, method="SKATO")
   x$p.value


} # runHandCraftedRegionsSkato
#----------------------------------------------------------------------------------------------------
viz <- function()
{
    igv <- start.igv("PTK2B", "hg19")
    track <- DataFrameAnnotationTrack("rs73223431",
                                      data.frame(chrom=c("chr8", "chr8"),
                                                 start=c(27219987-1, 27220310-1),
                                                 end=c(27219987,27220310),
                                                 name=c("rs73223431", "rs17057043"),
                                                 stringsAsFactors=FALSE),
                                      color="red")
    displayTrack(igv, track)

    track <- DataFrameQuantitativeTrack("SKATO", tbl.skat, autoscale=TRUE, color="brown")
    displayTrack(igv, track)

    track <- VariantTrack("rosmap", vcf)
    displayTrack(igv, track)


} # viz
#----------------------------------------------------------------------------------------------------
