source("../RegExplorer.R")
library(RUnit)

if(!exists("rex")){
    rex <- RegExplorer$new(targetGene="PTK2B", genome="hg19")
    rex$igv.init()
    }

#------------------------------------------------------------------------------------------------------------------------
gh <- function()
{
    tbl.gh <- rex$runGeneHancer()
    #tbl.gh <- subset(tbl.gh, elite)
    tbl.transformed <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
    tbl.transformed$combinedscore <- asinh(tbl.transformed$combinedscore)
    rex$displayQuantitativeTrack("gh.all", tbl.transformed, color="red")
    #rex$displayAnnotationTrack("HiC", subset(tbl.gh, !is.na(hic)), color="green")
    shoulder <- 10000
    big.loc <- sprintf("%s:%d-%d", tbl.gh$chrom[1], min(tbl.gh$start) - shoulder,
                       max(tbl.gh$end) + shoulder)
    showGenomicRegion(rex$get.igv(), big.loc)

} # gh
#------------------------------------------------------------------------------------------------------------------------
snps.in.LD <- function()
{
   # ptk2b variants:
   #   rs73223431:      1.29e-11
   #   rs17057043:      8.34e-10


   # 8	27337604	0.97	0.99	rs28834970	T	C	0.22	0.30	0.23	0.35		BLD, GI	FAT, BLD, THYM	6 tissues		6 altered motifs	1 hit		3 hits	PTK2B	intronic
   # 8	27347529	0.95	0.98	rs57735330	AG	A	0.46	0.31	0.23	0.35			BLD			Foxa,Pou5f1			3 hits	PTK2B	intronic
   # 8	27350609	0.95	0.98	rs6987305	G	A	0.46	0.31	0.23	0.35		BLD	12 tissues	4 tissues		Hoxc6,Pou2f2			3 hits	PTK2B	intronic
   # 8	27354393	0.97	0.99	rs2322599	G	A	0.45	0.31	0.23	0.35			9 tissues	BLD,BLD,BLD					6 hits	PTK2B	intronic
 ### 8	27362470	1	1	rs73223431	C	T	0.22	0.29	0.23	0.35		4 tissues	19 tissues	15 tissues	7 bound proteins	Brachyury,CEBPD,Eomes			3 hits	PTK2B	intronic
   # 8	27362793	0.99	1	rs17057043	G	A	0.47	0.31	0.23	0.35		BLD, GI	12 tissues	6 tissues	IRF1	4 altered motifs		2 hits	3 hits	PTK2B	intronic
   # 8	27369273	0.8	0.97	rs755951	A	C	0.52	0.35	0.24	0.39

  tbl.snps <- data.frame(rsid=c("rs28834970", "rs57735330", "rs6987305", "rs2322599", "rs73223431",
                                "rs17057043", "rs755951"),
                         chrom=rep("chr8", 7),
                         hg38=c(27337604, 27347529, 27350609, 27354393, 27362470, 27362793, 27369273),
                         hg19=c(27195121, 27205046, 27208126, 27211910, 27219987, 27220310, 27226790),
                         rSquared=c(0.97, 0.95, 0.95, 0.97, 1, 0.99, 0.8),
                         dPrime=c(0.99, 0.98, 0.98, 0.99, 1, 1, 0.97),
                         ref=c("T", "AG", "G", "G", "C", "G", "A"),
                         alt=c("C", "A", "A", "A", "T", "A", "C"),
                         stringsAsFactors=FALSE)

  save(tbl.snps, file="tbl.snps.RData")

  tbl.track <- tbl.snps[, c("chrom", "hg19", "hg19", "rSquared")]
  colnames(tbl.track) <- c("chrom", "start", "end", "score")
  tbl.track$start <- tbl.track$start - 1;
  rex$displayQuantitativeTrack("LD snps", tbl.track, color="brown")
  tbl.track <- tbl.snps[, c("chrom", "hg19", "hg19", "rsid")]
  colnames(tbl.track) <- c("chrom", "start", "end", "rsid")
  tbl.track$start <- tbl.track$start - 1
  rex$displayAnnotationTrack("rsid", tbl.track, color="blue")

} # snps.in.LD
#------------------------------------------------------------------------------------------------------------------------
buildTrenaModel <- function()
{
   mtx.expr <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/rosmap.14235x632.RData"))

   #chrom.loc <- "chr8"
   #start.loc <- 27166944
   #end.loc   <- 27274015
   #tbl.regions <- data.frame(chrom="chr8", start=27219862, end=27220442, stringsAsFactors=FALSE)
   #tbl.regions <- data.frame(chrom=chrom.loc, start=start.loc, end=end.loc, stringsAsFactors=FALSE)
      # intersect genehancer and hint20 footprints
   tbl.gh <- rex$runGeneHancer()

   tbl.regions <- data.frame(chrom="chr8", start=27160441, end=27284347, stringsAsFactors=FALSE)

   library(RPostgreSQL)
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_20", host="khaleesi")

   roi <- getGenomicRegion(rex$get.igv())

   query <- sprintf("select * from regions where chrom='chr8'  and start >= %d and endpos <= %d",
                    roi$start, roi$end)

   tbl.hint <- dbGetQuery(db, query)
   tbl.hint <- tbl.hint[, c("chrom", "start", "endpos")]
   colnames(tbl.hint)[3] <- "end"
   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.gh[,c("chrom", "start", "end")]),
                                        GRanges(tbl.hint)))
   dim(tbl.ov)
   tbl.hint.gh <- tbl.hint[unique(tbl.ov$subjectHits),]
   dim(tbl.hint.gh)   # 8611 3
   track <- DataFrameAnnotationTrack("hint+gh", tbl.hint.gh, color="black")
   displayTrack(rex$get.igv(), track)

   tbl.trena <- rex$build.trena.model(tbl.hint.gh,
                                      mtx.expr,
                                      fimo.threshold=1e-2,
                                      orderBy="spearmanCoeff"
                                      )
   rownames(tbl.trena) <- NULL
   tbl.20 <- head(tbl.trena, n=20)
   tbl.20$betaLasso <- round(tbl.20$betaLasso, digits=3)
   tbl.20$betaRidge <- round(tbl.20$betaRidge, digits=3)
   tbl.20$spearmanCoeff <- round(tbl.20$spearmanCoeff, digits=3)
   tbl.20$pearsonCoeff <- round(tbl.20$pearsonCoeff, digits=3)
   tbl.20$rfScore <- round(tbl.20$rfScore, digits=3)
   tbl.20$xgboost <- round(tbl.20$xgboost, digits=3)
   tbl.20$tfbs.e3 <- unlist(lapply(tbl.20$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value < 1e-3))))
   tbl.20$tfbs.e4 <- unlist(lapply(tbl.20$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value < 1e-4))))
   tbl.20$tfbs.e5 <- unlist(lapply(tbl.20$gene, function(gene) nrow(subset(tbl.fimo, tf==gene & p.value < 1e-5))))

   save(tbl.20, file="tbl.20.RData")

   gr.snps <- GRanges(seqnames=tbl.snps$chrom, IRanges(start=tbl.snps$hg19))
   tbl.fimo.select <- subset(tbl.fimo, tf %in% tbl.20$gene)
   gr.fimo <- with(tbl.fimo.select, GRanges(seqnames=chrom, IRanges(start=start-5, end=end+5)))
   tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.fimo))
   dim(tbl.ov)
   tbl.ov
   tbl.snps[tbl.ov$queryHits,]
   tbl.fimo.select[tbl.ov$subjectHits,]

} # buildTrenaModel
#------------------------------------------------------------------------------------------------------------------------
score.variants <- function()
{
   tfs <- tbl.trena$gene[1:10]
   motifs <- query(MotifDb, c("sapiens", "jaspar2018"), tfs)
   length(motifs)  # 10

   tbls.break <- list()
   list.break <- list()

   for(rsid in tbl.snps$rsid){
     tryCatch({
         x <- rex$score.TFBS.snp(rsid, motifs, method="ic")
         x <- calculatePvalue(x)
         tbl.x <- as.data.frame(x[x$effect=="strong"], row.names=NULL)
         tbls.break[[rsid]] <- tbl.x
         list.break[[rsid]] <- x
         print(tbl.x[, c("motifPos", "geneSymbol", "alleleDiff")])
         },
       error=function(e){printf("error with %s", rsid)})
     } # for rsid

   tbl.breaks <- do.call(rbind, tbls.break)

} # score.variants
#------------------------------------------------------------------------------------------------------------------------
showBindingSitesForTopTFs <- function()
{

  tbl.fimo <- rex$getFimoTable()
  fivenum(subset(tbl.fimo, tf=="CUX2")$p.value)
  roi <- getGenomicRegion(rex$get.igv())
    # $chrom: [1] "chr8"
    # $start: [1] 27190094
    # $end: [1] 27230658
    # $string: [1] "chr8:27,190,094-27,230,658"

  showGenomicRegion(rex$get.igv(), "chr8:27,190,094-27,230,658")

  top.tfs <- head(tbl.trena, n=10)$gene
  for(TF in top.tfs){
    tbl.tf <- subset(tbl.fimo, tf==TF) #  & chrom=="chr8" & start >= roi$start & end <= roi$end & p.value <= 1e-5)
    printf("--- %s: %d", TF, nrow(tbl.tf))
    if(nrow(tbl.tf) > 0)
       rex$displayAnnotationTrack(TF, tbl.tf[, c("chrom", "start", "end")], color="random")
    } # for TF

} # showBindingSitesForTopTFs
#------------------------------------------------------------------------------------------------------------------------
#   chrom  start endpos     tf                  name strand peakStart peakEnd
# 1  chr1 941011 941305 ARID1B GSE69566.ARID1B.HEPG2      .    941200  941201
# 2  chr1 941767 942084 ARID1B GSE69566.ARID1B.HEPG2      .    941890  941891
# 3  chr1 966645 967083 ARID1B GSE69566.ARID1B.HEPG2      .    966904  966905

chip <- function()
{
  db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
  roi <- getGenomicRegion(rex$get.igv())

  tfs <- tbl.10$gene
  for(tf in tfs){
     query <- sprintf("select * from chipseq where tf='%s' and chrom='%s' and start >= %d and endpos <= %d",
                       tf, roi$chrom, roi$start, roi$end)
     tbl <- dbGetQuery(db, query)
     if(nrow(tbl) > 0){
        colnames(tbl)[3] <- "end"
        track <- DataFrameAnnotationTrack(tf, tbl[, 1:3], color="random")
        displayTrack(rex$get.igv(), track)
        }
     } # for tf


} # chip
#----------------------------------------------------------------------------------------------------
