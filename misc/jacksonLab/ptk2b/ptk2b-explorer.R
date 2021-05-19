source("../RegExplorer.R")
library(RUnit)

if(!exists("rex")){
    rex <- RegExplorer$new(targetGene="PTK2B", genome="hg19")
    rex$igv.init()
    }

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    checkTrue(all(c("RegExplorer", "R6") %in% class(rex)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_geneHancer <- function()
{
    tbl.gh <- rex$runGeneHancer()
    tbl.transformed <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
    tbl.transformed$combinedscore <- asinh(tbl.transformed$combinedscore)
    rex$displayQuantitativeTrack("gh", tbl.transformed, color="red")
    rex$displayAnnotationTrack("HiC", subset(tbl.gh, !is.na(hic)), color="brown")
    shoulder <- 5000
    big.loc <- sprintf("%s:%d-%d", tbl.gh$chrom[1], min(tbl.gh$start) - shoulder,
                       max(tbl.gh$end) + shoulder)
    showGenomicRegion(rex$get.igv(), big.loc)

} # test_geneHancer
#----------------------------------------------------------------------------------------------------
test_displayVariants <- function()
{
   snp.1 <- "rs73223431" # 8:27219987 (GRCh37) C>T
   snp.2 <- "rs17057043" # 8:27220310 (GRCh37) G>A   (300 bp downstream)
   snp.3 <- "rs566760173"
   snp.loc.1 <- 27219987
   snp.loc.2 <- 27220310
   snp.loc.3 <- 27220128

   tbl.snps <- data.frame(chrom=rep("chr8", 3),
                          start= c(snp.loc.1-1, snp.loc.2-1, snp.loc.3-1),
                          end=c(snp.loc.1, snp.loc.2, snp.loc.3),
                          name=c(snp.1, snp.2, snp.3),
                          stringsAsFactors=FALSE)
   rex$displayAnnotationTrack("snps", tbl.snps, color="darkgreen")

} # test_displayVariants
#----------------------------------------------------------------------------------------------------
test_score.TFBS.snp <- function()
{
   message(sprintf("--- test_score.TFBS.snp"))

   motifs <- query(MotifDb, c("sapiens", "RARA"), c("jaspar2018"))[2]
   snp.loc <- 27219222
   snp <- "rs73223431"   # 8:27219987 (GRCh37) C>T
   snp.2 <- "rs17057043" # 8:27220310 (GRCh37) G>A   (300 bp downstream)

   results.0 <- rex$score.TFBS.snp(snp, motifs, method="log") #, "chr8", snp.loc-10, snp.loc+10)
   results.1 <- rex$score.TFBS.snp(snp, motifs, method="default") #, "chr8", snp.loc-10, snp.loc+10)
   plotMB(results=results.0, rsid = snp, effect = "strong")
   plotMB(results=results.1, rsid = snp, effect = "strong")

   motifs.2 <- query(MotifDb, c("sapiens", "TBR1"), c("jaspar2018", "HOCOMOCOv11-core"))
   results.2 <- rex$score.TFBS.snp(snp.2, motifs, method="log") #, "chr8", snp.loc-10, snp.loc+10)
   results.3 <- rex$score.TFBS.snp(snp.2, motifs.2, method="log") #, "chr8", snp.loc-10, snp.loc+10)


} # test_score.TFBS.snp
#----------------------------------------------------------------------------------------------------
test_trena <- function()
{
    mtx.expr <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))

    chrom.loc <- "chr8"
    start.loc <- 27166944
    end.loc   <- 27274015

    tbl.regions <- data.frame(chrom="chr8", start=27219862, end=27220442, stringsAsFactors=FALSE)
    tbl.regions <- data.frame(chrom=chrom.loc, start=start.loc, end=end.loc, stringsAsFactors=FALSE)

    tbl.regions <- rex$getGeneHancerTable()
    tbl.trena <- rex$build.trena.model(tbl.regions,
                                       mtx.expr,
                                       fimo.threshold=1e-3,
                                       orderBy="spearmanCoeff"
                                       )
    rownames(tbl.trena) <- NULL
    head(tbl.trena, n=10)
    tbl.fimo <- rex$getFimoTable()
    ntfs <- tbl.trena$gene[1:10]
    for(tf.name in tfs){
      tbl.sub <- subset(tbl.fimo, tf==tf.name)
      rex$displayAnnotationTrack(tf.name, tbl.sub[, 1:4], color="random")
      }

    tbl.fimo.rara <- subset(tbl.fimo, tf=="RARA")
    nrow(tbl.fimo)


    rex$viz()

    x <- rex$score.TFBS.snp("rs73223431", query(MotifDb, c("Hsapiens", "jaspar2018"), "MA0729.1"), method="log")
    x <- calculatePvalue(x)

    library(RPostgreSQL)
    db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")

    query <- sprintf("select * from chipseq where chrom='%s' and start >= %d and endpos <= %d",
                     chrom.loc, start.loc, end.loc)
    printf(query)
    tbl.hits <- dbGetQuery(db, query)
    printf("ChIP-seq hits: %d", nrow(tbl.hits))

} # test_trena
#----------------------------------------------------------------------------------------------------
test_raraEverywhere <- function()
{
   motifs <- query(MotifDb, c("sapiens", "RARA"), c("jaspar2018", "HOCOMOCOv11-core"))
   tbl.gh <- rex$getGeneHancerTable()
   threshold = 1e-2
   tbl.fimo <- rex$runFimo(tbl.gh, motifs, threshold)

   rex$displayAnnotationTrack("RARA*", tbl.fimo[, 1:4], color="darkblue")
   rex$displayAnnotationTrack("RARA*new", tbl.fimo[, 1:4], color="darkblue")
   tfs <- unique(tbl.fimo$tf)
   for(tf.name in tfs){
      tbl.sub <- subset(tbl.fimo, tf==tf.name)
      rex$displayAnnotationTrack(tf.name, tbl.sub[, 1:4], color="random")
      }

   vcf <- readVcf("ptk2b-236k-region-hg19.vcf.gz", "hg19")
   track <- VariantTrack("rosmap", vcf)
   displayTrack(rex$get.igv(), track)

} # test_raraEverywhere
#----------------------------------------------------------------------------------------------------
# look downstream
test_otherVariants <- function()
{
   snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, "8:27220127-27220130") #

   snp <- "rs566760173"
   results <- rex$score.TFBS.snp(snp, query(MotifDb, c("Hsapiens", "jaspar2018"), "RARA"), method="log")
   results <- calculatePvalue(results)

    # seqnames                                    chr8                                chr8                                chr8                                chr8
    # start                                   27220128                            27220128                            27220128                            27220128
    # end                                     27220128                            27220128                            27220128                            27220128
    # width                                          1                                   1                                   1                                   1
    # strand                                         -                                   -                                   -                                   -
    # SNP_id                               rs566760173                         rs566760173                         rs566760173                         rs566760173
    # REF                                            C                                   C                                   C                                   C
    # ALT                                            A                                   A                                   A                                   A
    # varType                                      SNV                                 SNV                                 SNV                                 SNV
    # motifPos                                  -15, 1                              -15, 2                              -15, 2                              -15, 1
    # geneSymbol                           RARA(var.2)                          RARA::RXRG                                RARA                          RARA::RXRA
    # dataSource                            jaspar2018                          jaspar2018                          jaspar2018                          jaspar2018
    # providerName                            MA0730.1                            MA1149.1                            MA0729.1                            MA0159.1
    # providerId                              MA0730.1                            MA1149.1                            MA0729.1                            MA0159.1
    # seqMatch     ccaccttcaccctgaccccAgcagtgttttctt   ccaccttcaccctgaccccAgcagtgttttcttt  ccaccttcaccctgaccccAgcagtgttttcttt  ccaccttcaccctgaccccAgcagtgttttctt
    # pctRef                                   0.72063                           0.8404957                           0.6767596                           0.7966474
    # pctAlt                                 0.6522503                           0.7778234                           0.6591926                           0.7165348
    # scoreRef                               -7.126781                            6.299785                           -22.24321                            4.491062
    # scoreAlt                               -13.65075                            2.848347                           -24.65005                         -0.04153728
    # Refpvalue                            0.002533439                         0.000163018                         0.007249348                        0.0004864033
    # Altpvalue                             0.01566023                          0.00161657                          0.01175797                          0.00634619
    # altPos                                         1                                   1                                   1                                   1
    # alleleDiff                             -6.523968                           -3.451438                           -2.406837                           -4.532599
    # effect                                    strong                              strong                              strong                              strong


} # test_otherVariants
#----------------------------------------------------------------------------------------------------
# create two mtx.expr, one for patients with the variant, one for those without
test_twoTrenaModels <- function()
{
   vcf <- readVcf("ptk2b-236k-region-hg19.vcf.gz", "hg19")
   # 8:27219987  rs73223431
   mtx.geno <- geno(vcf)$GT
   dim(mtx.geno)
   mtx.geno[1:10, 1:10]
   grep("8:27219987", rownames(mtx.geno)) # 2260
   snp.vec <- mtx.geno[2260,]
   table(snp.vec)  #  0/0 0/1 1/1
                   #  796 860 238
   mtx.expr <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))
   tbl.cov <- get(load("~/github/kats/podkat/rosmap/tbl.covAll-rosmap.RData"))

   mtx.new <- read.table("ROSMAP_RNAseq_FPKM_gene.tsv", sep="\t", as.is=TRUE, header=TRUE)
   dim(mtx.new)
   mtx.new[1:10, 1:10]

   pt.00 <- names(which(snp.vec=="0/0"))
   length(pt.00)
   pt.01 <- names(which(snp.vec=="0/1"))
   pt.11 <- intersect(names(which(snp.vec=="1/1")), tbl.cov$specimenID)
   length(pt.11)  # 143
   length(intersect(pt.11, colnames(mtx.expr)))

   tbl.md <- read.table("ROSMAP_assay_rnaSeq_metadata.csv", sep=",", header=TRUE, as.is=TRUE)

} # test_twoTrenaModels
#----------------------------------------------------------------------------------------------------
