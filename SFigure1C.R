rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','openxlsx','dittoSeq','ggplot2','FactoMineR','factoextra')
for (i in ps) {library(i, character.only = TRUE)}

## data
clinic = read.xlsx('clinical_64.xlsx')
count_48 = read.csv('tpm_matrix_48.csv', row.names = 1)
sample = colnames(count_48)
tumor = clinic[clinic$ID %in% sample, ]
lymph = clinic[clinic$t_or.ln == 'ln', ]

final = rbind(tumor, lymph)
colnames(feature) = c('sample_id','sample_time','sample_type','histology','patient_name','CPR','MPR')
sample = feature$sample_id
count = count_64[, colnames(count_64) %in% sample]
sample = colnames(count)
feature = feature[match(sample, feature$sample_id), ]

count_filter = count[rowSums(count) > 0, ]
pca = prcomp(t(count_filter), scale. = TRUE)
summary(pca)

## T vs LN
pca_data = data.frame(pca$x[, 1:30])
pca_data$cluster = as.factor(feature$sample_type)
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.5) +
  stat_ellipse(level = 0.95, show.legend = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black'))
group = factor(feature$sample_type, levels = c('t', 'ln'))
dat = as.data.frame(t(count_filter))
dat = na.omit(dat)
dat$group = group
pca_data = PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(
  pca_data, geom.ind = 'point',
  col.ind = dat$group,
  palette = c('#00AFBB', '#E7B800'),
  addEllipses = T,
  legend.title = "Groups") + theme_bw()

## MPR vs NMPR
feature = feature[feature$sample_type == 't', ]
table(feature$sample_id == colnames(count_48))
count_filter = count_48[rowSums(count_48) > 0, ]
pca = prcomp(t(count_filter), scale. = TRUE)
summary(pca)
group = factor(feature$MPR, levels = c('MPR', 'NMPR'))
dat = as.data.frame(t(count_filter))
dat = na.omit(dat)
dat$group = group
pca_data = PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(
  pca_data, geom.ind = 'point',
  col.ind = dat$group,
  palette = c('#00AFBB', '#E7B800'),
  addEllipses = T,
  legend.title = "Groups") + theme_bw()

## baseline vs surgery
feature = feature[feature$sample_type == 't', ]
table(feature$sample_id == colnames(count_48))
count_filter = count_48[rowSums(count_48) > 0, ]
pca = prcomp(t(count_filter), scale. = TRUE)
summary(pca)
group = factor(feature$sample_time, levels = c('baseline', 'surgery'))
dat = as.data.frame(t(count_filter))
dat = na.omit(dat)
dat$group = group
pca_data = PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(
  pca_data, geom.ind = 'point',
  col.ind = dat$group,
  palette = c('#00AFBB', '#E7B800'),
  addEllipses = T,
  legend.title = "Groups") + theme_bw()

