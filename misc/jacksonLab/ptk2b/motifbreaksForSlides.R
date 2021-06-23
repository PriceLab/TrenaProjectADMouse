library(MotifDb)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
if(!exists("tbl.snps"))
   load("tbl.snps.RData")


snps.gr <- snps.from.rsid(rsid = tbl.snps$rsid,
                          dbSNP=SNPlocs.Hsapiens.dbSNP144.GRCh37,
                          search.genome=BSgenome.Hsapiens.UCSC.hg19)

snps.gr <- snps.from.rsid(rsid = c(tbl.snps$rsid, "rs796228316"),
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
                       pwmList = motifs,
                       #pwmList = motifs.selected,
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
