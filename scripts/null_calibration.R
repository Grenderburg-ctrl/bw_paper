suppressPackageStartupMessages(library(DESeq2))
raw <- read.csv("data/salmon.merged.gene_counts_blast2_with_ensembl_ids.csv",
                check.names=FALSE, fileEncoding="UTF-8-BOM", stringsAsFactors=FALSE)
samples <- c("Bw1","Bw2","Bw3","Bw4","R004","R005","R006","R007")
sym <- as.character(raw$gene_id)
m <- round(as.matrix(raw[,samples])); storage.mode(m) <- "integer"; rownames(m) <- sym
mito <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6","COX1","COX2","COX3","CYTB","ATP6","ATP8")
drop <- toupper(sym) %in% mito | (grepl("^Rn",sym) & nchar(sym)<8)
m <- m[!drop,]; m <- m[rowSums(m)>=10,]

# Balanced mock contrasts: each arm holds 2 blast + 2 control, so the real condition
# effect cancels. Any "DEG" here is a false positive.
mocks <- list(
  mockA = c("A","A","B","B","A","A","B","B"),   # Bw1,Bw2,R004,R005 vs Bw3,Bw4,R006,R007
  mockB = c("A","B","A","B","A","B","A","B"),
  mockC = c("A","B","B","A","B","A","A","B")
)
sg <- function(r) sum(!is.na(r$padj) & r$padj < 0.05)

cat(sprintf("%-8s %12s %12s\n", "design", "DESeq2_EB", "geneWiseMLE"))
for (nm in names(mocks)) {
  cd <- data.frame(row.names=samples,
                   condition=factor(mocks[[nm]], levels=c("A","B")))
  d <- DESeq(DESeqDataSetFromMatrix(m, cd, ~condition), quiet=TRUE)
  r_eb <- results(d, contrast=c("condition","B","A"))
  dG <- d
  gw <- mcols(d)$dispGeneEst; gw[is.na(gw)] <- dispersions(d)[is.na(gw)]
  dispersions(dG) <- pmax(gw, 1e-8)
  dG <- nbinomWaldTest(dG, quiet=TRUE)
  r_gw <- results(dG, contrast=c("condition","B","A"))
  cat(sprintf("%-8s %12d %12d\n", nm, sg(r_eb), sg(r_gw)))
}
