library(GenomicRanges)
library(rtracklayer)
library(ghdb)
if(!exists("ghdb"))
   ghdb <- GeneHancerDB()

table.file <- "Top3_Candidate_Regions_Non_Coding_Variants.txt"
table.file <- "Top3_Candidate_Regions_Non_Coding_Variants-2.24.txt"
if(!exists("tbl.nc3")){
  tbl.nc3 <- read.table(table.file, sep="\t", as.is=TRUE, header=TRUE, quote="", fill=TRUE, nrow=-1)
  colnames(tbl.nc3)[1:4] <- c("chromLoc", "chrom", "start", "end")
  tbl.nc3$chrom <- sprintf("chr%s", tbl.nc3$chrom)
  }


if(!exists("hg19.to.hg38.chain")){
   chain.filename <- "hg19ToHg38.over.chain"
   if(!file.exists(chain.filename)){
     system("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
     system(sprintf("gunzip %s.gz", chain.filename))
     }
   hg19.to.hg38.chain <- import.chain(chain.filename)
   }


gr.hg19 <- with(tbl.nc3, GRanges(seqnames=chrom, IRanges(start, end)))
x <- liftOver(gr.hg19, hg19.to.hg38.chain)
length(x)   # 3
gr.hg38 <- unlist(x)
length(gr.hg38)
seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
tbl.hg38 <- as.data.frame(gr.hg38)
colnames(tbl.hg38)[1] <- "chrom"
tbl.hg38$chrom <- as.character(tbl.hg38$chrom)
dim(tbl.hg38)

igv.hg38 <- start.igv("chr1:161185278-161186278", "hg38")


tbl.nc3[1, "human.target.genes"]

strsplit(tbl.nc3[1,"human.target.genes"], ",")[[1]]

track <- DataFrameAnnotationTrack("nc1", tbl.hg38[, 1:3], color="blue")
displayTrack(igv.hg38, track)

shoulder <- 10000

i <- 1
i <- 2
i <- 3
showGenomicRegion(igv.hg38, with(tbl.hg38, sprintf("%s:%d-%d", chrom[i], start[i]-shoulder, end[i]+shoulder)))

tbl.gh <- with(tbl.hg38, queryByRegion(ghdb, chrom[1], start[1]-shoulder, end[1]+shoulder))

gh.hits <- as.data.frame(findOverlaps(GRanges(tbl.hg38[1, 1:3]), GRanges(tbl.gh[, 1:3]), type="any"))$subjectHits

tbl.gh2 <- tbl.gh[gh.hits, c(1,2,3,8,6)]
tbl.gh2$combinedscore <- asinh(tbl.gh2$combinedscore)
for(r in 1:nrow(tbl.gh2)){
    track <- DataFrameAnnotationTrack(sprintf("GH-%d", r), tbl.gh2[r,], color="random")
    displayTrack(igv.hg38, track)
    } # for r


dim(tbl.gh)
   tbl.elite <- subset(tbl.gh, elite)


if(!exists("igv.mm10"))
    igv.mm10 <- start.igv("Ankrd45", "mm10")

if(!exists("igv.hg38"))
    igv.hg38 <- start.igv("ANKRD45", "hg38")

if(!exists("mouse.to.human.chain")){
   chain.filename <- "mm10ToHg38.over.chain"
   if(!file.exists(chain.filename)){
      system("curl -O http://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/mm10ToHg38.over.chain.gz")
      system(sprintf("gunzip %s.gz", chain.filename))
      }
   mouse.to.human.chain <- import.chain(chain.filename)
   }

if(!exists("human.to.mouse.chain")){
   chain.filename <- "hg38ToMm10.over.chain"
   if(!file.exists(chain.filename)){
      system("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz")
      system(sprintf("gunzip %s.gz", chain.filename))
      }
   human.to.mouse.chain <- import.chain(chain.filename)
   }

#------------------------------------------------------------------------------------------------------------------------
displayMouseLocation <- function(gene, chrom, start, end)
{
   loc.string <- sprintf("%s:%d-%d", chrom, start-5000, end+5000)
   showGenomicRegion(igv.mm10, loc.string)
   tbl <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack(gene, tbl, color="random")
   displayTrack(igv.mm10, track)

} # displayMouseLocation
#------------------------------------------------------------------------------------------------------------------------
test_dispayMouseLocation <- function()
{
    message(sprintf("--- test_dispayMouseLocation"))
    with(tbl.nc3[1,], displayMouseLocation("Ankrd45", chrom, start, end))

} # test_dispayMouseLocation
#------------------------------------------------------------------------------------------------------------------------
getLiftedGeneHancer <- function(chrom, start, end, shoulder=10000)
{
   browser()
   gr.mouse <- GRanges(seqnames=chrom, IRanges(start-shoulder, end+shoulder))
   x <- liftOver(gr.mouse, mouse.to.human.chain)
   length(x)   # 3
   gr.hg38 <- unlist(x)
   length(gr.hg38)
   seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
   tbl.hg38 <- as.data.frame(gr.hg38)
   colnames(tbl.hg38)[1] <- "chrom"
   tbl.hg38$chrom <- as.character(tbl.hg38$chrom)

   shoulder <- 5000
   start <- min(tbl.hg38$start - shoulder)
   end   <- max(tbl.hg38$end   + shoulder)
   chrom <- sub("chr", "", tbl.hg38$chrom[1])

   tbl.gh <- queryByRegion(ghdb, chrom, start, end)
   tbl.elite <- subset(tbl.gh, elite)
   dim(tbl.elite)
   gr.hg38 <- GRanges(tbl.gh)
   seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]

   x <- liftOver(gr.hg38, human.to.mouse.chain)
   length(x)   # 3
   length(unlist(x))
   gr.hg38 <- unlist(x)
   length(gr.hg38)
   seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
   tbl.hg38 <- as.data.frame(gr.hg38)
   colnames(tbl.hg38)[1] <- "chrom"
   tbl.hg38$chrom <- as.character(tbl.hg38$chrom)


} # getLiftedGeneHancer
#------------------------------------------------------------------------------------------------------------------------
test_getLiftedGeneHancer <- function()
{
   message(sprintf("--- test_getLiftedGeneHancer"))
   tbl.gh <- getLiftedGeneHancer("chr1", 161150068, 161161068, shoulder=1000)
   tbl.gh <- getLiftedGeneHancer("chr1", 161155068, 161156068)


} # test_getLiftedGeneHancer
#------------------------------------------------------------------------------------------------------------------------
# with(tbl.nc3[1,], displayMouseLocation("ANKRD45", chrom, start, end))
#
# showGenomicRegion(igv, rownames(tbl.nc3)[2])    # Acoxl      human ACOX1   chr17:75,940,517-75,980,199   first enzyme of the fatty acid beta-oxidation pathway
# showGenomicRegion(igv, rownames(tbl.nc3)[3])    # Got1l1     human GOT1L1  chr8:37,932,548-37,941,503    Putative aspartate aminotransferase, cytoplasmic 2
#
# showGenomicRegion(igv, "Ankrd45")               # chr1:161,141,546-161,171,510
# showGenomicRegion(igv, "Acox1")                 # chr11:116,170,883-116,200,045
# showGenomicRegion(igv, "Got1l1")                # chr8:27,196,459-27,224,844
#
#
#
#
#
# gr.human <- GRanges(seqnames=c("chr1", "chr17", "chr8"),
#                     IRanges(start=c(173596894, 75940517, 37932548),
#                               end=c(173694763, 75980199, 37941503)))
# x2 <- liftOver(gr.human[1], human.chain)
# tbl.mm10.ankrd45 <- as.data.frame(unlist(x2))
# colnames(tbl.mm10.ankrd45)[1] <- "chrom"
# tbl.mm10.ankrd45$chrom <- as.character(tbl.mm10.ankrd45$chrom)
#
# track <- DataFrameAnnotationTrack("from human", tbl.mm10.ankrd45[, 1:3], color="random")
# displayTrack(igv, track)
# showGenomicRegion(igv, sprintf("%s:%d-%d", tbl.mm10.ankrd45$chrom[1], min(tbl.mm10.ankrd45$start) - 2000, max(tbl.mm10.ankrd45$end) + 2000))
#
#   #-----------------------------------------------------------------------------
#   # lay down a track showing the original area jackson lab reported for ankrd45
#   #-----------------------------------------------------------------------------
# track <- DataFrameAnnotationTrack("jackson lab", data.frame(chrom="chr1", start=173596894, end=173694763, stringsAsFactors=FALSE, color="random"))
# displayTrack(igv, track)
#
# # chr1:173,596,894-173,694,763
# # chr17:75,940,517-75,980,199
# # chr8:37,932,548-37,941,503
#
#
#
#
# x <- liftover(gr.mouse, chain)
# length(x)   # 3
# gr.hg38 <- unlist(x)
# gr.hg38
# length(gr.hg38)  # 142
# seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
#
# gr.mouse
#
# igv2 <- start.igv("ANKRD45", "hg38")
# showGenomicRegion(igv2, "ACOX1")
# showGenomicRegion(igv2, "Got1l1")
#
#
