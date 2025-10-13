rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','ggplot2','ggrepel','DESeq2','edgeR')
for (i in ps) {library(i, character.only = TRUE)}

fc = 2.0
pval = 0.05

## data
count = read.csv('count_matrix.csv', row.names = 1)
clinic = read.csv('clinical_62.csv', header = T, row.names = 1)
estimate = read.csv('esitimate.csv', header = T, row.names = 1)
clinic = clinic[clinic$sample_type == 'ln', ]
count = count[, colnames(count) %in% clinic$sample_id]
estimate = estimate %>%
  t() %>% as.data.frame() %>% select(4) %>%
  tibble::rownames_to_column(var = 'sample_id')
df = merge(clinic, estimate, by = 'sample_id', all.y = F)
df = df %>% mutate(group = case_when(MPR == "NMPR" ~ 0, MPR == "MPR" ~ 1))
df$group = factor(df$group)
table(colnames(count) == df$sample_id)

## edgeR
mpr = factor(df$group)
levels(mpr)
y = DGEList(counts = count, group = mpr)
# y = y[, group %in% c("2", "4")]
# y$samples$group = factor(y$samples$group)
keep = filterByExpr(y, group = mpr)
y = y[keep, , keep.lib.sizes = F]
y = calcNormFactors(y)
y = estimateDisp(y)
et = exactTest(y)
et = topTags(et, n = 10000)
et = as.data.frame(et)
et = tibble::rownames_to_column(et, var = 'gene')
colnames(et) = c('gene','log2fc','log2cpm','pvalue','padj')
etSig = na.omit(et)
etSig$up_down = ifelse(
  etSig$log2fc > log2(fc) & etSig$padj < pval, 'up',
  ifelse(etSig$log2fc < (-log2(fc)) & etSig$padj < pval, 'down', 'normal'))
etSig = etSig %>% na.omit() %>% arrange(desc(pvalue))
listdown = etSig[etSig$up_down == 'down', ] %>% arrange(log2fc)
listup = etSig[etSig$up_down == 'up', ] %>% arrange(desc(log2fc))
## plot
panel = c('SCGB1A1','SEC14L3','SPP1','MARCO','TREM2','CCL18','VSIG4','CD163',
          'AKR1C2','IGFL4','AKR1C1','AKR1C3','ANXA10',
          'SLC13A5','CERS3','FSTL5','CASC9','TM4SF4')
genelist = etSig[etSig$gene %in% panel,]
genelist = subset(genelist, select = -up_down)
ggplot(data = etSig, aes(x = log2fc, y = -log10(padj), color = up_down)) + 
  geom_point(alpha = 0.5, size = 4) + 
  scale_x_continuous(limits = c(-8, 6), breaks = seq(-8, 6, 2)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  xlab('log2FC') + ylab('-log10(Pvalue)') +
  geom_vline(xintercept = c(-1,1), lty = 3, col = 'black', lwd = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), lty = 3, col = 'black',lwd = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c('#0077b6', '#d5dfe5', '#e63946')) +
  geom_label_repel(
    data = genelist,
    aes(x = log2fc, y = -log10(padj), label = gene),
    size = 3, color = 'black',
    box.padding = unit(0.4, 'lines'),
    segment.color = 'black',
    segment.size = 0.4
  )


## DESeq2
y = DESeqDataSetFromMatrix(countData = count, colData = df, design= ~group)
dds = DESeq(y, fitType = 'mean', minReplicatesForReplace = 7, parallel = F)
et = results(dds, contrast = c('group', '1', '0'))
etSig = data.frame(et, stringsAsFactors = FALSE, check.names = FALSE)
etSig = na.omit(etSig)
etSig$up_down = ifelse(
  etSig$log2FoldChange > log2(fc) & etSig$padj < pval, 'up',
  ifelse(etSig$log2FoldChange < (-log2(fc)) & etSig$padj < pval, 'down', 'normal'))
etSig = etSig %>% tibble::rownames_to_column(var = 'gene') %>% arrange(desc(padj))
listdown = etSig[etSig$up_down == 'down', ] %>% arrange(log2FoldChange)
listdown = etSig[etSig$up_down == 'down', ] %>% arrange(padj)
listup = etSig[etSig$up_down == 'up', ] %>% arrange(desc(log2FoldChange))
listup = etSig[etSig$up_down == 'up', ] %>% arrange(padj)
## plot
panel = c('STK32C','KCNH7','IER5','RGS20','GOS2',
          'MAGEA1','ITPRID1','MUC5B','DMBT1','AKR1C3')
genelist = etSig[etSig$gene %in% panel,]
genelist = subset(genelist, select = -up_down)
ggplot(data = etSig, aes(x = log2FoldChange, y = -log10(padj), color = up_down)) + 
  geom_point(alpha = 0.5, size = 4) + 
  scale_x_continuous(limits = c(-27, 9), breaks = seq(-27, 9, 9)) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  xlab('log2FC') + ylab('-log10(Pvalue)') +
  geom_vline(xintercept = c(-1,1), lty = 3, col = 'black', lwd = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), lty = 3, col = 'black',lwd = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c('#0077b6', '#d5dfe5', '#e63946')) +
  geom_label_repel(
    data = genelist,
    aes(x = log2FoldChange, y = -log10(padj), label = gene),
    size = 3, color = 'black',
    box.padding = unit(0.4, 'lines'),
    segment.color = 'black',
    segment.size = 0.4
  )


## edgeR-adjust
mpr = factor(df$group)
levels(mpr)
y = DGEList(counts = count, group = mpr)
keep = filterByExpr(y, group = mpr)
y = y[keep, , keep.lib.sizes = F]
y = calcNormFactors(y)
y$samples$purity = df$TumorPurity
design = model.matrix(~ group + purity, data = y$samples)
y = estimateDisp(y, design)
fit = glmFit(y, design)
coef = grep('1', colnames(design))
lrt = glmLRT(fit, coef = coef)
et = topTags(lrt, n=Inf)
et = as.data.frame(et)
et = tibble::rownames_to_column(et, var = 'gene')
colnames(et) = c('gene','log2fc','log2cpm','lr','pvalue','padj')
etSig = na.omit(et)
etSig$up_down = ifelse(
  etSig$log2fc > log2(fc) & etSig$padj < pval, 'up',
  ifelse(etSig$log2fc < (-log2(fc)) & etSig$padj < pval, 'down', 'normal'))
etSig = etSig %>% arrange(desc(pvalue))
listdown = etSig[etSig$up_down == 'down', ] %>% arrange(log2fc)
listdown = etSig[etSig$up_down == 'down', ] %>% arrange(padj)
listup = etSig[etSig$up_down == 'up', ] %>% arrange(desc(log2fc))
listup = etSig[etSig$up_down == 'up', ] %>% arrange(padj)
## plot
panel = c('CERS3','GPT2','RRM2','ATF3','ARRDC4','TUBA1B')
genelist = etSig[etSig$gene %in% panel,]
genelist = subset(genelist, select = -up_down)
ggplot(data = etSig, aes(x = log2fc, y = -log10(pvalue), color = up_down)) + 
  geom_point(alpha = 0.5, size = 4) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  xlab('log2FC') + ylab('-log10(Pvalue)') +
  geom_vline(xintercept = c(-1,1), lty = 3, col = 'black', lwd = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), lty = 3, col = 'black',lwd = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c('blue', 'black', 'red')) +
  geom_label_repel(
    data = genelist,
    aes(x = log2fc, y = -log10(pvalue), label = gene),
    size = 6, color = 'black',
    box.padding = unit(0.4, 'lines'),
    segment.color = 'black',
    segment.size = 0.4
  )

