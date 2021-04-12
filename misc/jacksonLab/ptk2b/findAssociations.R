# findAssociations.R: look for association between variants and phenotype in ROSMAP data
# rs73223431: chr8 27219987
# rs17057043: chr8 27220310
#----------------------------------------------------------------------------------------------------
library(VariantAnnotation)
vcf.file <- "ptk2b-region-hg19.vcf.gz"
tbl.cov <- get(load("~/github/kats/podkat/rosmap/tbl.covAll-rosmap.RData"))
#----------------------------------------------------------------------------------------------------
getGenotype <- function()
{
  vcf <- readVcf(vcf.file)
  mtx.geno <- geno(vcf)$GT
  dim(mtx.geno)   # 203 1894
  snp.1 <- mtx.geno[grep("8:27219987", rownames(mtx.geno)),]
  snp.2 <- mtx.geno[grep("8:27220310", rownames(mtx.geno)),]
  length(which(snp.1 == snp.2)) # 1861
  length(which(snp.1 != snp.2)) # 33

  as.data.frame(sort(table(snp.1)))
     #   snp.1 Freq
     # 1   1/1  238
     # 2   0/0  796
     # 3   0/1  860

  as.data.frame(sort(table(snp.2)))
     #   snp.2 Freq
     # 1   1/1  249
     # 2   0/0  773
     # 3   0/1  872

} # getGenotype
#----------------------------------------------------------------------------------------------------
exploreCovariates <- function()
{
   dim(mtx.geno) # 203 1894
   sample.ids <- colnames(mtx.geno)
   length(grep("SM-", sample.ids))  # 1151
   dim(tbl.cov)    # 1144 19
   sample.ids <- intersect(sample.ids, tbl.cov$specimenID)
   length(sample.ids) # 1144

} # exploreCovariates
#----------------------------------------------------------------------------------------------------
simpleAsociations <- function()
{

  snp.1 <- mtx.geno[grep("8:27219987", rownames(mtx.geno)), sample.ids]
  snp.2 <- mtx.geno[grep("8:27220310", rownames(mtx.geno)), sample.ids]
  length(which(snp.1 == snp.2)) # 1142
  length(which(snp.1 != snp.2)) # 2

  snp.1 <- snp.1[sample.ids]
  snp.2 <- snp.2[sample.ids]

  snp.01 <- rep(0, length(sample.ids))
  hets <- which(snp.1 == "0/1")
  homs <- which(snp.1 == "1/1")
  snp.01[hets] <- 1
  snp.01[homs] <- 2
  table(snp.01)  #   487, 514, 143

  snp.02 <- rep(0, length(sample.ids))
  hets <- which(snp.2 == "0/1")
  homs <- which(snp.2 == "1/1")
  snp.02[hets] <- 1
  snp.02[homs] <- 2
  table(snp.02)  #   485, 516, 143

  sample.index.matches <- match(sample.ids, tbl.cov$specimenID)
  length(sample.index.matches)  # 1144
  cogdx <- tbl.cov$cogdx[sample.index.matches]
  table(cogdx)
      #   1   2   3   4   5   6
      # 358 269  20 420  55  21

  tbl.01 <- data.frame(snp.01=snp.01, snp.02=snp.02, cogdx=cogdx)
  deleters <- which(is.na(cogdx))
  tbl.01 <- tbl.01[-deleters,]
  tbl.01$snps <- tbl.01$snp.01 + tbl.01$snp.02
     #    0   1   2   4
     #  484   2 514 143


  cogdx.simple <- rep(0, nrow(tbl.01))
  ad <- which(tbl.01$cogdx %in% c(4,5))
  length(ad) # 475
  cogdx.simple[ad] <- 4
  table(cogdx.simple)   # 668 0, 475 4
  tbl.01$cogdx.simple <- cogdx.simple

  with(tbl.01, cor(cogdx, snp.01))
  with(tbl.01, cor(cogdx.simple, snp.01))
  with(tbl.01, cor(cogdx.simple, snp.02))

  summary(lm(cogdx.simple ~ 0 + snp.01 * snp.02, data=tbl.01))
    #               Estimate Std. Error t value Pr(>|t|)
    # snp.01         -1.5535     1.6035  -0.969   0.3328
    # snp.02          4.0000     1.5884   2.518   0.0119 *
    # snp.01:snp.02  -0.7967     0.1365  -5.835 6.99e-09 ***
    # ---
    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    #
    # Residual standard error: 2.246 on 1140 degrees of freedom
    # Multiple R-squared:  0.2431,	Adjusted R-squared:  0.2411
    # F-statistic:   122 on 3 and 1140 DF,  p-value: < 2.2e-16

  summary(lm(cogdx.simple ~ 0 + snp.01 + snp.02, data=tbl.01))
    #        Estimate Std. Error t value Pr(>|t|)
    # snp.01   -2.770      1.613  -1.717   0.0862 .
    # snp.02    4.000      1.611   2.483   0.0132 *
    # ---
    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    #
    # Residual standard error: 2.279 on 1141 degrees of freedom
    # Multiple R-squared:  0.2205,	Adjusted R-squared:  0.2191
    # F-statistic: 161.3 on 2 and 1141 DF,  p-value: < 2.2e-16

  summary(lm(cogdx.simple ~ 0 + snps, data=tbl.01))
    #       Estimate Std. Error t value Pr(>|t|)
    #  snps  0.61666    0.03462   17.81   <2e-16 ***
    #  ---
    #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    #
    #  Residual standard error: 2.282 on 1142 degrees of freedom
    #  Multiple R-squared:  0.2175,	Adjusted R-squared:  0.2168
    #  F-statistic: 317.3 on 1 and 1142 DF,  p-value: < 2.2e-16

  summary(lm(cogdx.simple ~ 0 + snp.01, data=tbl.01))
  summary(lm(cogdx.simple ~ 0 + snp.02, data=tbl.01))
  summary(lm(cogdx.simple ~ 0 + snp.01 + snp.02, data=tbl.01))
  summary(lm(cogdx.simple ~ 0 + snp.01 * snp.02, data=tbl.01))
  summary(lm(cogdx ~ 0 + snp.01 * snp.02, data=tbl.01))

  tbl.01$fake <- sample(0:4, size=nrow(tbl.01), replace=TRUE)
  summary(lm(fake ~ 0 + snp.01 + snp.02, data=tbl.01))

} # simpleAssociations
#----------------------------------------------------------------------------------------------------

