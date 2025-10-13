rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','clusterProfiler','GSVA','DESeq2','GSEABase','enrichplot','msigdbr')
for (i in ps) {library(i, character.only = TRUE)}

## data
count = read.csv('count_matrix.csv', row.names = 1)
patient = read.csv('clinical_48.csv', header = T, row.names = 1)
clinic = read.csv('clinical_83.csv', header = T, row.names = 1)
clinic = clinic %>% .[.$sample_id %in% patient$sample_id, ] %>% arrange(sample_id)
clinic = clinic %>% mutate(group = case_when(mpr == "No-MPR" ~ 0, mpr == "MPR" ~ 1))
count = count[, colnames(count) %in% clinic$sample_id]
clinic$group = factor(clinic$group)
table(colnames(count) == clinic$sample_id)

## DEG
y = DESeqDataSetFromMatrix(countData = count, colData = clinic, design= ~group)
dds = DESeq(y, fitType = 'mean', minReplicatesForReplace = 7, parallel = F)
et = results(dds, contrast = c('group', '1', '0'))
etSig = data.frame(et, stringsAsFactors = FALSE, check.names = FALSE)
etSig = na.omit(etSig)
etSig = etSig %>% dplyr::select(c('log2FoldChange')) %>% tibble::rownames_to_column(var = 'SYMBOL')

## Pathway
gene.id = bitr(etSig$SYMBOL, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
df = merge(etSig, gene.id, by = 'SYMBOL', all = F)
df = df[order(df$log2FoldChange, decreasing = T), ]
gene.expr = df$log2FoldChange
names(gene.expr) = df$ENTREZID
# gene list
metab.set = read.gmt('Metabolism.v2024.7.Hs.gmt')
metab.set = metab.set[metab.set$gene != '', ]
colnames(metab.set) = c('term','SYMBOL')
metab.id = bitr(metab.set$SYMBOL, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
metab.df = merge(metab.set, metab.id, by = 'SYMBOL', all = F)
metab.df = metab.df[,-1]
metab.df = arrange(metab.df, term)
metab.df = metab.df[metab.df$term %in% c(
  'Steroid biosynthesis','Steroid hormone biosynthesis',
  'Metabolism of steroid hormones','Metabolism of steroids',
  'Glycosphingolipid metabolism',
  'Glycosphingolipid biosynthesis - lacto and neolacto series',
  'Glycosphingolipid biosynthesis - globo and isoglobo series'), ]
metab.df$term = gsub(' - ', '_', metab.df$term)
metab.df$term = gsub(' ', '_', metab.df$term)

## plot
gene.expr = gene.expr[abs(gene.expr) > 0.5]
y = GSEA(gene.expr, TERM2GENE = metab.df)
gseaplot2(
  y, geneSetID = c(3,4),
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  base_size = 15,
  ES_geom = 'line',
  color = c('#483D8B'),
  pvalue_table = TRUE
)

