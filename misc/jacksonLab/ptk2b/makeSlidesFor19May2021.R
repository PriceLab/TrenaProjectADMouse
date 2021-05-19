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
    tbl.transformed <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
    tbl.transformed$combinedscore <- asinh(tbl.transformed$combinedscore)
    rex$displayQuantitativeTrack("gh", tbl.transformed, color="red")
    #rex$displayAnnotationTrack("HiC", subset(tbl.gh, !is.na(hic)), color="green")
    shoulder <- 10000
    big.loc <- sprintf("%s:%d-%d", tbl.gh$chrom[1], min(tbl.gh$start) - shoulder,
                       max(tbl.gh$end) + shoulder)
    showGenomicRegion(rex$get.igv(), big.loc)

} # gh
#------------------------------------------------------------------------------------------------------------------------
snps.in.LD <- function()
{
   # 8	27337604	0.97	0.99	rs28834970	T	C	0.22	0.30	0.23	0.35		BLD, GI	FAT, BLD, THYM	6 tissues		6 altered motifs	1 hit		3 hits	PTK2B	intronic
   # 8	27347529	0.95	0.98	rs57735330	AG	A	0.46	0.31	0.23	0.35			BLD			Foxa,Pou5f1			3 hits	PTK2B	intronic
   # 8	27350609	0.95	0.98	rs6987305	G	A	0.46	0.31	0.23	0.35		BLD	12 tissues	4 tissues		Hoxc6,Pou2f2			3 hits	PTK2B	intronic
   # 8	27354393	0.97	0.99	rs2322599	G	A	0.45	0.31	0.23	0.35			9 tissues	BLD,BLD,BLD					6 hits	PTK2B	intronic
   # 8	27362470	1	1	rs73223431	C	T	0.22	0.29	0.23	0.35		4 tissues	19 tissues	15 tissues	7 bound proteins	Brachyury,CEBPD,Eomes			3 hits	PTK2B	intronic
   # 8	27362793	0.99	1	rs17057043	G	A	0.47	0.31	0.23	0.35		BLD, GI	12 tissues	6 tissues	IRF1	4 altered motifs		2 hits	3 hits	PTK2B	intronic
   # 8	27369273	0.8	0.97	rs755951	A	C	0.52	0.35	0.24	0.39

  tbl.snps <- data.frame(chrom=rep("chr8", 7),
                         hg38=c(27337604, 27347529, 27350609, 27354393, 27362470, 27362793, 27369273),
                         hg19=c(27195121, 27205046, 27208126, 27211910, 27219987, 27220310, 27226790),
                         rSquared=c(0.97, 0.95, 0.95, 0.97, 1, 0.99, 0.8),
                         dPrime=c(0.99, 0.98, 0.98, 0.99, 1, 1, 0.97),
                         rsid=c("rs28834970", "rs57735330", "rs6987305", "rs2322599", "rs73223431",
                                "rs17057043", "rs755951"),
                         ref=c("T", "AG", "G", "G", "C", "G", "A"),
                         alt=c("C", "A", "A", "A", "T", "A", "C"),
                         stringsAsFactors=FALSE)

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
   tbl.regions <- data.frame(chrom="chr8", start=27160441, end=27284347, stringsAsFactors=FALSE)

   library(RPostgreSQL)
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_20", host="khaleesi")
   query <- "select count(*) from regions where chrom='chr14' and start >= 74955470 and endpos <= 74955570"
   query <- "select * from regions where chrom='chr8'  and start >= 27164032 and endpos <= 27284691"
   tbl.regions <- dbGetQuery(db, query)
   tbl.regions <- tbl.regions[, c("chrom", "start", "endpos")]
   colnames(tbl.regions) <- c("chrom", "start", "end")
   #tbl.regions <- rex$getGeneHancerTable()
   tbl.trena <- rex$build.trena.model(tbl.regions,
                                      mtx.expr,
                                      fimo.threshold=1e-3,
                                      orderBy="spearmanCoeff"
                                      )
   rownames(tbl.trena) <- NULL
   tbl.10 <- head(tbl.trena, n=10)
   tbl.10$betaLasso <- round(tbl.10$betaLasso, digits=3)
   tbl.10$betaRidge <- round(tbl.10$betaRidge, digits=3)
   tbl.10$spearmanCoeff <- round(tbl.10$spearmanCoeff, digits=3)
   tbl.10$pearsonCoeff <- round(tbl.10$pearsonCoeff, digits=3)
   tbl.10$rfScore <- round(tbl.10$rfScore, digits=3)
   tbl.10$xgboost <- round(tbl.10$xgboost, digits=3)


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
    tbl.tf <- subset(tbl.fimo, tf==TF & chrom=="chr8" & start >= roi$start & end <= roi$end & p.value <= 1e-5)
    printf("--- %s: %d", TF, nrow(tbl.tf))
    if(nrow(tbl.tf) > 0)
       rex$displayAnnotationTrack(TF, tbl.tf[, c("chrom", "start", "end")], color="random")
    } # for TF

} # showBindingSitesForTopTFs
#------------------------------------------------------------------------------------------------------------------------
