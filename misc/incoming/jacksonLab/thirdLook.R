library(TrenaProjectAD)
library(TrenaMultiScore)
library(GenomicRanges)
library(rtracklayer)
library(ghdb)
library(GenomicScores)
library(phastCons7way.UCSC.hg38)
library(phastCons100way.UCSC.hg38)
library(httr)
library(EnsDb.Hsapiens.v79)


#--------------------------------------------------
# make sure the databases on khaleesi are available
#--------------------------------------------------
suppressWarnings(
  db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
if(length(db.access.test) == 0)
   stop("khaleesi database server unavailable")

if(!exists("tpad"))
   tpad <- TrenaProjectAD()

if(!exists("ghdb"))
   ghdb <- GeneHancerDB()

table.file <- "Top3_Candidate_Regions_Non_Coding_Variants.txt"
table.file <- "Top3_Candidate_Regions_Non_Coding_Variants-2.24.txt"
if(!exists("tbl.nc3")){
  tbl.nc3 <- read.table(table.file, sep="\t", as.is=TRUE, header=TRUE, quote="", fill=TRUE, nrow=-1)
  colnames(tbl.nc3)[1:4] <- c("chromLoc", "chrom", "start", "end")
  tbl.nc3$chrom <- sprintf("chr%s", tbl.nc3$chrom)
  }

# 3 non-coding variants with AD associations
# rs73223431: C>T chr8:27362470 (GRCh38)
#   regulates increased expression and splicing of PTK2B in brain
#   ptk2b intron variant
#   Connecting the dots: potential of data integration to identify regulatory SNPs in late-onset
#     Alzheimer's disease GWAS findings. 2014. shows RegulomeDB in action
#   Genetic Association of FERMT2, HLA-DRB1, CD2AP, and PTK2B Polymorphisms With Alzheimer's
#     Disease Risk in the Southern Chinese Population. associated with early onset AD (EOAD)
# rs10194375 is linked to increased BIN1 expression in the brain and other tissues,
#              in a region that is active in oligodendrocytes.
#   rs73223431 was published by Phil deJager.


#------------------------------------------------------------------------------------------------------------------------
# first intro on PTK2B
# 27355483
# 27376578
# 27361392
# 27365825
rs73223431 <- function()
{
   uri <- "http://www.ebi.ac.uk/eqtl/api/associations/rs73223431"
   response <- GET(uri)
   tbls.raw <- fromJSON(httr::content(response, as="text"))
   tbls <- lapply(tbls.raw$`_embedded`$associations, function(item) as.data.frame(item))
   tbl.effects <- do.call(rbind, tbls)
   tbl.effects <- tbl.effects[order(tbl.effects$neg_log10_pvalue, decreasing=TRUE),]
   ensgs <- tbl.effects$gene_id
   symbols <- lapply(ensgs, function(ensg) unique(select(EnsDb.Hsapiens.v79, key=ensg, keytype="GENEID")$SYMBOL))
   stopifnot(length(ensgs) == length(symbols))
   tbl.effects$gene = unlist(symbols)
   coi <- c("chromosome", "position", "rsid", "ref", "alt", "pvalue", "gene", "maf", "r2", "neg_log10_pvalue", "qtl_group", "median_tpm", "se", "type", "an", "beta", "study_id", "ac", "variant", "molecular_trait_id", "gene_id", "tissue")
   tbl.effects <- tbl.effects[, coi]
   rownames(tbl.effects) <- NULL
   length(ensg) # 20
   length(symbols) # 12
   igv <- start.igv("PTK2B", "hg38")
   chrom <- "chr8"
   snp.loc <- 27362470
   tbl.snp <- data.frame(chrom=chrom, start=snp.loc, end=snp.loc, stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("rs73223431", tbl.snp, color="red")
   displayTrack(igv, track)

   shoulder <- 10000
   tbl.gh <- queryByRegion(ghdb, chrom, start=snp.loc-shoulder, end=snp.loc+shoulder)
   dim(tbl.gh)
   tbl.gh <- subset(tbl.gh, combinedscore > 5)
   track <- DataFrameAnnotationTrack("gh", tbl.gh[, c("chrom", "start", "end")], color="red")
   displayTrack(igv, track)

   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, "PTK2B", tissues="all")
   track <- DataFrameQuantitativeTrack("gh.ptk2b", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)
   loc <- getGenomicRegion(igv)
   tbl.intersecting.eqtl <- subset(tbl.gh, start > loc$start & end < loc$end)

   tbl.gh <- getEnhancers(tpad, maxSize=30000)
   with(tbl.gh, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom[1], min(start)-1000, max(end)+ 1000)))

   tbl.gh$score <- asinh(tbl.gh$combinedscore)
   track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c("chrom", "start", "end", "score")], color="brown",
                                       autoscale=TRUE)
   displayTrack(igv, track)


   loc <- getGenomicRegion(igv)
   starts <- with(loc, seq(start, end, by=5))
   ends <- starts + 5
   count <- length(starts)
   tbl.blocks <- data.frame(chrom=rep(loc$chrom, count), start=starts, end=ends, stringsAsFactors=FALSE)
   gr.blocks <- GRanges(tbl.blocks)

   tbl.cons7 <- as.data.frame(gscores(phastCons7way.UCSC.hg38, gr.blocks), stringsAsFactors=FALSE)
   tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
   tbl.cons7 <- tbl.cons7[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast7", tbl.cons7, autoscale=TRUE, color="red")
   displayTrack(igv, track)

   tbl.cons100 <- as.data.frame(gscores(phastCons100way.UCSC.hg38, gr.blocks), stringsAsFactors=FALSE)
   tbl.cons100$chrom <- as.character(tbl.cons100$seqnames)
   tbl.cons100 <- tbl.cons100[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast100", tbl.cons100, autoscale=TRUE, color="blue")
   displayTrack(igv, track)

} # region.three
#------------------------------------------------------------------------------------------------------------------------
