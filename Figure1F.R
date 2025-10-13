rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','GSVA','pheatmap','ggplot2')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
count = read.csv('tpm_matrix.csv', row.names = 1)
sample = read.csv('clinical_48.csv', header = T, row.names = 1)
clinic = read.csv('clinical_83.csv', header = T, row.names = 1)
clinic = clinic %>% .[.$sample_id %in% sample$sample_id, ] %>% arrange(sample_id)
clinic = clinic %>% mutate(group = case_when(
  mpr == 'MPR' & sample_time == 'baseline' ~ 'MPR_baseline',
  mpr == 'MPR' & sample_time == 'surgery' ~ 'MPR_surgery',
  mpr == 'No-MPR' & sample_time == 'baseline' ~ 'NMPR_baseline',
  mpr == 'No-MPR' & sample_time == 'surgery' ~ 'NMPR_surgery'))
cluster = clinic[, colnames(clinic) %in% c('sample_id','group')]
plasma = c('MZB1','JCHAIN','SLAMF7','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHGP','IGLC2','IGLC3','IGLC7')
plasma_score = gsva(as.matrix(count), list(plasma = plasma), method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
plasma_score = plasma_score %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample_id")
df = merge(cluster, plasma_score, by = 'sample_id')

df = df[df$group %in% c('MPR_baseline','NMPR_baseline'), ]
res = wilcox.test(plasma ~ group, data = df)
res = aov(plasma ~ group, data = df)

ggplot(df, aes(group, plasma, fill = group)) +
  geom_boxplot(
    outlier.colour = "grey",
    outlier.size = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(),
    legend.position = 'none')


