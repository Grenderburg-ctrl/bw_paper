suppressPackageStartupMessages(library(DESeq2))
raw <- read.csv("data/salmon.merged.gene_counts_blast2_with_ensembl_ids.csv",
                check.names=FALSE, fileEncoding="UTF-8-BOM", stringsAsFactors=FALSE)
samples <- c("Bw1","Bw2","Bw3","Bw4","R004","R005","R006","R007")
sym <- as.character(raw$gene_id)
m <- round(as.matrix(raw[,samples])); storage.mode(m) <- "integer"; rownames(m) <- sym
mito <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6","COX1","COX2","COX3","CYTB","ATP6","ATP8")
drop <- toupper(sym) %in% mito | (grepl("^Rn",sym) & nchar(sym)<8)
m <- m[!drop,]; m <- m[rowSums(m)>=10,]
cd <- data.frame(row.names=samples, condition=factor(c(rep("blast",4),rep("control",4)),
                                                     levels=c("control","blast")))
dds <- DESeq(DESeqDataSetFromMatrix(m, cd, ~condition), quiet=TRUE)
ct <- c("condition","blast","control")
res <- results(dds, contrast=ct)
md <- mcols(dds)

panel <- c("Tmc6","Mcoln2","Kcnk2","Piezo2","Ano8","Tmem63a","Pkd1","Ano6",
           "Piezo1","Tmem63b","Tmc7","Trpm4","Kcnk3","Myh6","Myh7","Col3a1")
ref_padj <- c(2.9e-20,3.6e-08,9.3e-08,4.7e-07,6.6e-06,1.0e-05,3.4e-05,4.8e-05,
              7.7e-05,1.0e-03,3.2e-03,4.0e-03,0.53,4.8e-06,6.3e-06,9.5e-07)

# Hypothesis: reference stays at (or near) the gene-wise Cox-Reid MLE, little/no EB shrinkage.
ddsG <- dds
gw <- md$dispGeneEst
gw[is.na(gw)] <- dispersions(dds)[is.na(gw)]
gw <- pmax(gw, 1e-8)
dispersions(ddsG) <- gw
ddsG <- nbinomWaldTest(ddsG, quiet=TRUE)
resG <- results(ddsG, contrast=ct)

i  <- match(panel, rownames(dds))
cmp <- data.frame(gene=panel,
                  padj_ref     = ref_padj,
                  padj_DESeq2  = signif(res$padj[i],3),
                  padj_genewise= signif(resG$padj[i],3),
                  disp_fit     = signif(md$dispFit[i],3),
                  disp_gene    = signif(md$dispGeneEst[i],3),
                  disp_DESeq2  = signif(dispersions(dds)[i],3))
cmp$log10_err_DESeq2  <- round(log10(cmp$padj_DESeq2 / ref_padj),2)
cmp$log10_err_genewise<- round(log10(cmp$padj_genewise/ ref_padj),2)
print(cmp, row.names=FALSE)

sg <- function(r) sum(!is.na(r$padj) & r$padj<0.05)
cat(sprintf("\nDEG DESeq2 (EB-shrunk, default): %d   [reference claims 6371]\n", sg(res)))
cat(sprintf("DEG with gene-wise MLE dispersion: %d\n", sg(resG)))
cat(sprintf("  down %d / up %d\n",
    sum(!is.na(resG$padj)&resG$padj<0.05&resG$log2FoldChange<0),
    sum(!is.na(resG$padj)&resG$padj<0.05&resG$log2FoldChange>0)))
cat(sprintf("\nmedian |log10 padj error| DESeq2 vs ref:   %.2f\n",
            median(abs(cmp$log10_err_DESeq2))))
cat(sprintf("median |log10 padj error| genewise vs ref: %.2f\n",
            median(abs(cmp$log10_err_genewise))))
cat(sprintf("\nSpearman rho(padj) DESeq2 vs ref:   %.4f\n",
            cor(log10(cmp$padj_DESeq2), log10(ref_padj), method="spearman")))
cat(sprintf("Spearman rho(padj) genewise vs ref: %.4f\n",
            cor(log10(cmp$padj_genewise), log10(ref_padj), method="spearman")))
