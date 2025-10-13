rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','pheatmap','ComplexHeatmap','circlize')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
count = read.csv('tpm_matrix.csv', header = T, row.names = 1)
patient = read.csv('clinical_48.csv', header = T, row.names = 1)
clinic = read.csv('clinical_83.csv', header = T, row.names = 1)
estimate = read.csv('esitimate.csv', header = T, row.names = 1)
immune = read.csv('TIMER2_immune_cell.csv', header = T, row.names = 1)
clinic = clinic %>% .[.$sample_id %in% patient$sample_id, ] %>% arrange(sample_id)
clinic = clinic %>% mutate(group = case_when(
    mpr == 'MPR' & sample_time == 'baseline' ~ 'MPR_baseline',
    mpr == 'MPR' & sample_time == 'surgery' ~ 'MPR_surgery',
    mpr == 'No-MPR' & sample_time == 'baseline' ~ 'NMPR_baseline',
    mpr == 'No-MPR' & sample_time == 'surgery' ~ 'NMPR_surgery'))
group = select(clinic, c('sample_id','group'))
count = count[, colnames(count) %in% clinic$sample_id]
table(colnames(count) == clinic$sample_id)
immune = immune[, colnames(immune) %in% clinic$sample_id]
table(colnames(immune) == clinic$sample_id)
estimate = estimate[, colnames(estimate) %in% clinic$sample_id]
table(colnames(estimate) == clinic$sample_id)

genelist = c(
  'GZMA','GZMB','GZMH','GZMK','IFNG','PRF1',
  'IL2','IL4','IL6','IL10','IL17A','TGFB1',
  'CXCL8','CXCL9','CXCL10','CXCL11','CXCL12','CXCL13','CXCR3','CXCR4','CXCR5','CCR7','CCR8',
  'PDCD1','CD274','CTLA4','HAVCR2','LAG3','TIGIT')
test = count[rownames(count) %in% genelist, ]
test = test %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
test = merge(test, group, by = 'sample_id')
pval = lapply(colnames(test)[2:(ncol(test)-1)], function(x) {
  formula = as.formula(paste(x, "~ group"))
  kruskal.test(formula, data = test)
})
pval = do.call(rbind, lapply(pval, function(x) {data.frame(pval = as.numeric(x[['p.value']]))}))
pval$var = colnames(test)[2:(ncol(test)-1)]
genelist = c(
  'GZMA','GZMK','GZMH','IFNG','IL17A','CXCL9',
  'CD274','PDCD1','CTLA4','HAVCR2','LAG3','TIGIT')
gene = count[rownames(count) %in% genelist, ]
gene = gene[genelist, ]

cell = immune[grepl('XCELL', rownames(immune)), ]
rownames(cell) = gsub('_XCELL', '', rownames(cell))
order = c(
  'T_cell_CD8_naive','T_cell_CD8_central_memory','T_cell_CD8_effector_memory',
  'T_cell_CD4_naive','T_cell_CD4_memory','T_cell_CD4_Th1','T_cell_CD4_Th2',
  'NK_cell',
  'B_cell','B_cell_naive','B_cell_memory','B_cell_plasma',
  'Monocyte','Macrophage','Macrophage_M1','Macrophage_M2','Myeloid_dendritic_cell',
  'Neutrophil','Eosinophil','Mast_cell')
cell = cell[order, ]
test = immune %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
test = merge(test, group, by = 'sample_id')
pval = lapply(colnames(test)[2:(ncol(test)-1)], function(x) {
  formula = as.formula(paste(x, "~ group"))
  kruskal.test(formula, data = test)
})
pval = do.call(rbind, lapply(pval, function(x) {data.frame(pval = as.numeric(x[['p.value']]))}))
pval$var = colnames(test)[2:(ncol(test)-1)]
pval = pval[pval$pval < 0.1, ]

anno = clinic %>% select(c('sample_id','group','cTNM','smoking','pathology.x','pdl1'))
anno = tibble::column_to_rownames(anno, var = 'sample_id')
pval = lapply(anno[-1], function(x) {fisher.test(anno[[1]], x)})
pval = do.call(rbind, lapply(pval, function(x) {data.frame(pval = as.numeric(x[['p.value']]))}))
test = estimate %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'sample_id')
test = merge(test, group, by = 'sample_id')
pval = lapply(colnames(test)[2:(ncol(test)-1)], function(x) {
  formula = as.formula(paste(x, "~ group"))
  kruskal.test(formula, data = test)
})
pval = do.call(rbind, lapply(pval, function(x) {data.frame(pval = as.numeric(x[['p.value']]))}))
pval$var = colnames(test)[2:(ncol(test)-1)]
pval = pval[pval$pval < 0.1, ]

cell = read.csv('TIMER2_immune_cell_XCELL.csv', header = T, row.names = 1)
estimate = read.csv('esitimate_48_scale.csv', header = T, row.names = 1)
rownames(cell) = gsub('_XCELL', '', rownames(cell))
cell = cell[order, ]

anno = anno[colnames(cell), ]
estimate = estimate[colnames(cell), ]
gene = gene[, colnames(cell)]

df = rbind(cell, gene)
anno = cbind(estimate, anno)
colnames(anno)[1:4] = c('stroma','immune','estimate','purity')
anno = select(anno, c('stroma','immune','estimate','pdl1','smoking','cTNM','pathology.x','group'))

# cell = read.csv('TIMER2_immune_cell_scale_XCELL.csv', header = T, row.names = 1)
# estimate = read.csv('esitimate_48_scale.csv', header = T, row.names = 1)
# table(colnames(cell) == colnames(estimate))
# order = colnames(cell)
# estimate = estimate[order, ]
# gene = gene[, order]
# anno = anno[order, ]
# gene = scale(gene)

## plot
bk = c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(
  as.matrix(log2(df+1)), annotation_col = anno,
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
    smoking = c('never' = '#f8eceb', 'smoking' = '#9786ad', 'quit' = '#c69fb1'),
    purity = c('white', '#2E5278'),
    estimate = c('white', '#2E5278'),
    immune = c('white', '#2E5278'),
    stroma = c('white', '#2E5278')
    ),
  gaps_row = c(12, 20, 26),
  gaps_col = c(5, 31, 35),
  fontsize = 4,
  fontsize_col = 5,
  fontsize_row = 5,
  border = F)

cell = cell %>% t() %>% as.data.frame() %>% scale() %>% t() %>% as.data.frame()
col_fun = colorRamp2(c(-1,1), c("white", "#2E5278"))
col_fun = colorRamp2(c(-2, 0,2), c("#1e96fc", "#f8f9fa", "#f77f00"))
col_fun = colorRamp2(c(0,1), c("white", "#BCAB79"))
col_fun = colorRamp2(c(0,1), c("white", "#BCAB79"))
Heatmap(
  cell,
  # top_annotation = c(abc,efg),
  # left_annotation = hij,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = col_fun,
  name = 'Z-score'
  # row_gap = unit(4,'mm'),
  # row_names_gp = gpar(fontsize = 8)
  )

