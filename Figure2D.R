rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','clusterProfiler','GSVA','ggplot2','ggpubr')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
count = read.csv('count_matrix.csv', row.names = 1)
patient = read.csv('clinical_48.csv', header = T, row.names = 1)
clinic = read.csv('clinical_83.csv', header = T, row.names = 1)
clinic = clinic %>% .[.$sample_id %in% patient$sample_id, ] %>% arrange(sample_id)
clinic = clinic %>% mutate(group = case_when(
  mpr == 'MPR' & sample_time == 'baseline' ~ 'MPR_baseline',
  mpr == 'MPR' & sample_time == 'surgery' ~ 'MPR_surgery',
  mpr == 'No-MPR' & sample_time == 'baseline' ~ 'NMPR_baseline',
  mpr == 'No-MPR' & sample_time == 'surgery' ~ 'NMPR_surgery'))
count = count[, colnames(count) %in% clinic$sample_id]
table(colnames(count) == clinic$sample_id)
cluster = clinic[, colnames(clinic) %in% c('sample_id','group','mpr')]

## akr
akr_set = list(akr_set = c('AKR1C1','AKR1C2','AKR1C3'))
akr_score = gsva(as.matrix(count), akr_set, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
akr_score = akr_score %>% t() %>% as.data.frame() %>% z_score_col()
akr_score = akr_score %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample_id")
akr_score = merge(cluster, akr_score, by = 'sample_id')
akr_score$akr_set = apply(akr_score$akr_set, 2, as.numeric)

df = akr_score
df = akr_score[akr_score$group %in% c('MPR_baseline','NMPR_baseline'), ]
df = akr_score[akr_score$group %in% c('MPR_surgery','NMPR_surgery'), ]

## compare
wilcox.test(akr_set1 ~ group, data = df)
wilcox.test(akr_set2 ~ group, data = df)
wilcox.test(akr_set3 ~ group, data = df)
wilcox.test(akr_set1 ~ mpr, data = df)
wilcox.test(akr_set2 ~ mpr, data = df)
wilcox.test(akr_set3 ~ mpr, data = df)

## plot
ggplot(df, aes(mpr, akr_set1, fill = mpr)) +
  geom_violin(
    width = 0.7,
    alpha = 0.7,
    trim = FALSE) +
  geom_jitter(
    color = 'black',
    width = 0.15,
    size = 1.6,
    alpha = 0.6) +
  scale_y_continuous(
    limits = c(-4, 4),
    breaks = c(-4, -2, 0, 2, 4)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(),
    legend.position = 'none')

ggviolin(
  df, 'group', 'akr_set1',
  draw_quantiles = 0.5,
  orientation = 'horiz',
  add = 'jitter'
  )

