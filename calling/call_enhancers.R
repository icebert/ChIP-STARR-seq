library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(SummarizedExperiment)
library(ChIPseeker)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(data.table)




lib_size <- read.table('library_size.tsv', header=T)
rownames(lib_size) <- lib_size$lib

data <- read.table('sliding_window.bed')

rranges <- GRanges(seqnames = Rle(data$V1),
                   ranges = IRanges(start=data$V2+1, end=data$V3),
                   strand = Rle('*', nrow(data)))

seqinfo(rranges) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(rranges)]


#
# AA
#
data <- read.table('count/AA_D_1.bed')
cnt <- matrix(data$V4, nrow=nrow(data))
data <- read.table('count/AA_D_2.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/AA_R_1.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/AA_R_2.bed')
cnt <- cbind(cnt, data$V4)

data <- DataFrame(row.names=c('AA_D_1', 'AA_D_2', 'AA_R_1', 'AA_R_2'))
data$library_size <- lib_size[rownames(data), ]$lib_size
data$type <- factor(c('DNA', 'DNA', 'RNA', 'RNA'), levels=c('DNA', 'RNA'))


AA_se <- SummarizedExperiment(assays=list(counts=cnt), rowRanges=rranges, colData=data)

assays(AA_se)$tpm <- t(t(assay(AA_se)) / (AA_se$library_size / 1000000))

AA_se <- AA_se[assays(AA_se)$tpm[, 1] >= 1 & assays(AA_se)$tpm[, 2] >= 1, ]




#
# AK
#
data <- read.table('count/AK_D_1.bed')
cnt <- matrix(data$V4, nrow=nrow(data))
data <- read.table('count/AK_D_2.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/AK_R_1.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/AK_R_2.bed')
cnt <- cbind(cnt, data$V4)

data <- DataFrame(row.names=c('AK_D_1', 'AK_D_2', 'AK_R_1', 'AK_R_2'))
data$library_size <- lib_size[rownames(data), ]$lib_size
data$type <- factor(c('DNA', 'DNA', 'RNA', 'RNA'), levels=c('DNA', 'RNA'))


AK_se <- SummarizedExperiment(assays=list(counts=cnt), rowRanges=rranges, colData=data)

assays(AK_se)$tpm <- t(t(assay(AK_se)) / (AK_se$library_size / 1000000))

AK_se <- AK_se[assays(AK_se)$tpm[, 1] >= 1 & assays(AK_se)$tpm[, 2] >= 1, ]




#
# KA
#
data <- read.table('count/KA_D_1.bed')
cnt <- matrix(data$V4, nrow=nrow(data))
data <- read.table('count/KA_D_2.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/KA_R_1.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/KA_R_2.bed')
cnt <- cbind(cnt, data$V4)

data <- DataFrame(row.names=c('KA_D_1', 'KA_D_2', 'KA_R_1', 'KA_R_2'))
data$library_size <- lib_size[rownames(data), ]$lib_size
data$type <- factor(c('DNA', 'DNA', 'RNA', 'RNA'), levels=c('DNA', 'RNA'))


KA_se <- SummarizedExperiment(assays=list(counts=cnt), rowRanges=rranges, colData=data)

assays(KA_se)$tpm <- t(t(assay(KA_se)) / (KA_se$library_size / 1000000))

KA_se <- KA_se[assays(KA_se)$tpm[, 1] >= 1 & assays(KA_se)$tpm[, 2] >= 1, ]





#
# KK
#
data <- read.table('count/KK_D_1.bed')
cnt <- matrix(data$V4, nrow=nrow(data))
data <- read.table('count/KK_D_2.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/KK_R_1.bed')
cnt <- cbind(cnt, data$V4)
data <- read.table('count/KK_R_2.bed')
cnt <- cbind(cnt, data$V4)

data <- DataFrame(row.names=c('KK_D_1', 'KK_D_2', 'KK_R_1', 'KK_R_2'))
data$library_size <- lib_size[rownames(data), ]$lib_size
data$type <- factor(c('DNA', 'DNA', 'RNA', 'RNA'), levels=c('DNA', 'RNA'))


KK_se <- SummarizedExperiment(assays=list(counts=cnt), rowRanges=rranges, colData=data)

assays(KK_se)$tpm <- t(t(assay(KK_se)) / (KK_se$library_size / 1000000))

KK_se <- KK_se[assays(KK_se)$tpm[, 1] >= 1 & assays(KK_se)$tpm[, 2] >= 1, ]






AA_res <- DESeqDataSet(AA_se, design = ~ type)
AA_res <- DESeq(AA_res, fitType='local')
AA_res <- results(AA_res)

rowRanges(AA_se)$log2FoldChange <- AA_res$log2FoldChange
rowRanges(AA_se)$pvalue <- AA_res$pvalue
rowRanges(AA_se)$score <- AA_res$log2FoldChange
export.bw(rowRanges(AA_se), 'AA.log2FC.bw')

AA_res_neg <- rowRanges(AA_se)[AA_res$pvalue < 0.05 & AA_res$log2FoldChange < 0]
AA_res <- rowRanges(AA_se)[AA_res$pvalue < 0.05 & AA_res$log2FoldChange > 0]
export.bed(AA_res, 'AA.enhancer.bed')
export.bed(AA_res_neg, 'AA.negative.bed')



AK_res <- DESeqDataSet(AK_se, design = ~ type)
AK_res <- DESeq(AK_res, fitType='local')
AK_res <- results(AK_res)

rowRanges(AK_se)$log2FoldChange <- AK_res$log2FoldChange
rowRanges(AK_se)$pvalue <- AK_res$pvalue
rowRanges(AK_se)$score <- AK_res$log2FoldChange
export.bw(rowRanges(AK_se), 'AK.log2FC.bw')

AK_res_neg <- rowRanges(AK_se)[AK_res$pvalue < 0.05 & AK_res$log2FoldChange < 0]
AK_res <- rowRanges(AK_se)[AK_res$pvalue < 0.05 & AK_res$log2FoldChange > 0]
export.bed(AK_res, 'AK.enhancer.bed')
export.bed(AK_res_neg, 'AK.negative.bed')



KA_res <- DESeqDataSet(KA_se, design = ~ type)
KA_res <- DESeq(KA_res, fitType='local')
KA_res <- results(KA_res)

rowRanges(KA_se)$log2FoldChange <- KA_res$log2FoldChange
rowRanges(KA_se)$pvalue <- KA_res$pvalue
rowRanges(KA_se)$score <- KA_res$log2FoldChange
export.bw(rowRanges(KA_se), 'KA.log2FC.bw')

KA_res_neg <- rowRanges(KA_se)[KA_res$pvalue < 0.05 & KA_res$log2FoldChange < 0]
KA_res <- rowRanges(KA_se)[KA_res$pvalue < 0.05 & KA_res$log2FoldChange > 0]
export.bed(KA_res, 'KA.enhancer.bed')
export.bed(KA_res_neg, 'KA.negative.bed')



KK_res <- DESeqDataSet(KK_se, design = ~ type)
KK_res <- DESeq(KK_res, fitType='local')
KK_res <- results(KK_res)

rowRanges(KK_se)$log2FoldChange <- KK_res$log2FoldChange
rowRanges(KK_se)$pvalue <- KK_res$pvalue
rowRanges(KK_se)$score <- KK_res$log2FoldChange
export.bw(rowRanges(KK_se), 'KK.log2FC.bw')

KK_res_neg <- rowRanges(KK_se)[KK_res$pvalue < 0.05 & KK_res$log2FoldChange < 0]
KK_res <- rowRanges(KK_se)[KK_res$pvalue < 0.05 & KK_res$log2FoldChange > 0]
export.bed(KK_res, 'KK.enhancer.bed')
export.bed(KK_res_neg, 'KK.negative.bed')




