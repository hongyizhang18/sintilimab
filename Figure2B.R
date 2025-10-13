rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','openxlsx','ggplot2','clusterProfiler','pheatmap','GSVA')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
count = read.csv('tpm_matrix.csv', header = T, row.names = 1)
patient = read.csv('clinical_48.csv', header = T, row.names = 1)
clinic = read.csv('clinical_83.csv', header = T, row.names = 1)
metab = read.csv('signature_metab.csv', header = T, row.names = 1)
metab = metab %>% t() %>% as.data.frame()
clinic = clinic %>% .[.$sample_id %in% patient$sample_id, ] %>% arrange(sample_id)
clinic = clinic %>% mutate(group = case_when(
  mpr == 'MPR' & sample_time == 'baseline' ~ 'MPR_baseline',
  mpr == 'MPR' & sample_time == 'surgery' ~ 'MPR_surgery',
  mpr == 'No-MPR' & sample_time == 'baseline' ~ 'NMPR_baseline',
  mpr == 'No-MPR' & sample_time == 'surgery' ~ 'NMPR_surgery'))
group = select(clinic, c('sample_id','group'))
count = count[, colnames(count) %in% clinic$sample_id]
table(colnames(count) == clinic$sample_id)
metab = metab[, colnames(metab) %in% clinic$sample_id]
table(colnames(metab) == clinic$sample_id)

## metabolism
metab_set = read.gmt('Metabolism.v2023.2.Hs.gmt')
metab_lists = list()
for (pathway in unique(metab_set$term)) {metab_lists[[pathway]] = metab_set[metab_set$term == pathway, 2]}
metab_score = gsva(as.matrix(count), metab_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
metab_score = metab_score %>% as.data.frame() %>% .[, colnames(metab)]

## compare
test = metab %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
test = merge(test, group, by = 'sample_id')
pval = lapply(colnames(test)[2:(ncol(test)-1)], function(x) {
  formula = as.formula(paste(x, "~ group"))
  kruskal.test(formula, data = test)
})
pval1 = do.call(rbind, lapply(pval, function(x) {data.frame(pval = as.numeric(x[['p.value']]))}))
pval1$var = colnames(test)[2:(ncol(test)-1)]

rownames(metab_score) = gsub('/', ' ', rownames(metab_score))
rownames(metab_score) = gsub('-', ' ', rownames(metab_score))
rownames(metab_score) = gsub('"', '', rownames(metab_score))
rownames(metab_score) = gsub(',', '', rownames(metab_score))
rownames(metab_score) = gsub('\\(', '', rownames(metab_score))
rownames(metab_score) = gsub(')', '', rownames(metab_score))
rownames(metab_score) = gsub(' +', '_', rownames(metab_score))
test = metab_score %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
test = merge(test, group, by = 'sample_id')
pval = lapply(colnames(test)[2:(ncol(test)-1)], function(x) {
  formula = as.formula(paste(x, "~ group"))
  kruskal.test(formula, data = test)
})
pval2 = do.call(rbind, lapply(pval, function(x) {data.frame(pval = as.numeric(x[['p.value']]))}))
pval2$var = colnames(test)[2:(ncol(test)-1)]

## plot
# pattern = paste(keywords, collapse = "|")
# df = metab_score[-grep(pattern, rownames(metab_score)), ]
# write.csv(rownames(df), file = 'metabolism_pathway.csv')
# write.csv(rownames(metab_score), file = 'metabolism_pathway_all.csv')
metab_pathway = read.csv('result/metabolism_pathway.csv', header = T)
order = lapply(metab_pathway[,1:4], function(x) {x = x[x != '']; na.omit(x)})
order = unlist(order)

cell = read.csv('result/TIMER2_immune_cell_XCELL.csv', header = T, row.names = 1)
anno = clinic %>% select(c('sample_id','group','cTNM','smoking','pathology.x','pdl1'))
anno = tibble::column_to_rownames(anno, var = 'sample_id')
anno = select(anno, c('pdl1','smoking','cTNM','pathology.x','group'))
anno = arrange(anno, group)
anno = anno[colnames(cell), ]
df = metab_score[order, ]
plot = df[, rownames(anno)]

## plot
bk = c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(
  as.matrix(log2(plot+1)), annotation_col = anno,
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  scale = 'row',
  breaks = bk,
  legend_breaks = seq(-2,2,1),
  color = c(colorRampPalette(colors = c("#1e96fc","#f8f9fa"))(length(bk)/2),
            colorRampPalette(colors = c("#f8f9fa","#f77f00"))(length(bk)/2)),
  annotation_colors = list(
    group = c('MPR_baseline' = '#EDB16E', 'MPR_surgery' = '#D51F27', 'NMPR_baseline' = '#B5D6E4', 'NMPR_surgery' = '#4278B0'),
    pdl1 = c('0%' = '#A3A8F0', '1-49%' = '#757DE8', '>50%' = '#3F51B5'),
    cTNM = c('IIIA' = '#b7e4c7', 'IIIB' = '#52b788', 'IIIC' = '#248277'),
    pathology.x = c('LUAD' = '#e27d60', 'LUSC' = '#f2cc8f', 'NSCLC' = '#f4f1de'),
    smoking = c('never' = '#f8eceb', 'smoking' = '#9786ad', 'quit' = '#c69fb1')),
  gaps_col = c(5, 31, 35),
  gaps_row = c(29, 66, 94),
  fontsize = 4,
  fontsize_col = 5,
  fontsize_row = 5,
  border = F)

pheatmap(
  as.matrix(scale(plot)), annotation_col = anno,
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  scale = 'row',
  breaks = bk,
  legend_breaks = seq(-2,2,1),
  color = c(colorRampPalette(colors = c("#1e96fc","#f8f9fa"))(length(bk)/2),
            colorRampPalette(colors = c("#f8f9fa","#f77f00"))(length(bk)/2)),
  annotation_colors = list(
    group = c('MPR_baseline' = '#EDB16E', 'MPR_surgery' = '#D51F27', 'NMPR_baseline' = '#B5D6E4', 'NMPR_surgery' = '#4278B0'),
    pdl1 = c('0%' = '#A3A8F0', '1-49%' = '#757DE8', '>50%' = '#3F51B5'),
    cTNM = c('IIIA' = '#b7e4c7', 'IIIB' = '#52b788', 'IIIC' = '#248277'),
    pathology.x = c('LUAD' = '#e27d60', 'LUSC' = '#f2cc8f', 'NSCLC' = '#f4f1de'),
    smoking = c('never' = '#f8eceb', 'smoking' = '#9786ad', 'quit' = '#c69fb1')),
  gaps_col = c(5, 31, 35),
  gaps_row = c(29, 66, 94),
  fontsize = 4,
  fontsize_col = 5,
  fontsize_row = 5,
  border = F)


