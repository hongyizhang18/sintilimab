rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','openxlsx','ggplot2','pheatmap','clusterProfiler','GSVA','reshape2','gghalves','ggpubr')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
count = read.csv('tpm_matrix.csv', row.names = 1)
clinic = read.csv('clinical_62.csv', header = T, row.names = 1)
clinic = clinic[clinic$sample_type == 't', ]
count = count[, colnames(count) %in% clinic$sample_id]
table(colnames(count) == clinic$sample_id)

## hallmark
hallmark = read.gmt('Hallmark.v2023.2.Hs.gmt')
hallmark_lists = list()
for (pathway in unique(hallmark$term)) {
  hallmark_lists[[pathway]] = hallmark[hallmark$term == pathway, 2]
}
hallmark_score = gsva(as.matrix(count), hallmark_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## kegg
kegg = read.gmt('C2.KEGG_Medicus.v2023.2.Hs.gmt')
kegg_lists = list()
for (pathway in unique(kegg$term)) {
  kegg_lists[[pathway]] = kegg[kegg$term == pathway, 2]
}
kegg_score = gsva(as.matrix(count), kegg_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## reactome
reactome = read.gmt('C2.Reactome.v2023.2.Hs.gmt')
reactome_lists = list()
for (pathway in unique(reactome$term)) {
  reactome_lists[[pathway]] = reactome[reactome$term == pathway, 2]
}
reactome_score = gsva(as.matrix(count), reactome_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## wiki
wiki = read.gmt('C2.Wikipathways.v2023.2.Hs.gmt')
wiki_lists = list()
for (pathway in unique(wiki$term)) {
  wiki_lists[[pathway]] = wiki[wiki$term == pathway, 2]
}
wiki_score = gsva(as.matrix(count), wiki_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## C4
c4 = read.gmt('C4.All.v2023.2.Hs.gmt')
c4_lists = list()
for (pathway in unique(c4$term)) {
  c4_lists[[pathway]] = c4[c4$term == pathway, 2]
}
c4_score = gsva(as.matrix(count), c4_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## C5
c5 = read.gmt('C5.GOBP.v2023.2.Hs.gmt')
c5_lists = list()
for (pathway in unique(c5$term)) {
  c5_lists[[pathway]] = c5[c5$term == pathway, 2]
}
c5_score = gsva(as.matrix(count), c5_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## C7
c7 = read.gmt('C7.ImmuneSigDB.v2023.2.Hs.gmt')
c7_lists = list()
for (pathway in unique(c7$term)) {
  c7_lists[[pathway]] = c7[c7$term == pathway, 2]
}
c7_score = gsva(as.matrix(count), c7_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

## metabolism
metab = read.gmt('Metabolism.v2023.2.Hs.gmt')
metab_lists = list()
for (pathway in unique(metab$term)) {
  metab_lists[[pathway]] = metab[metab$term == pathway, 2]
}
metab_score = gsva(as.matrix(count), metab_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

sig_hall = hallmark_score %>% t() %>% as.data.frame()
sig_kegg = kegg_score %>% t() %>% as.data.frame()
sig_reac = reactome_score %>% t() %>% as.data.frame()
sig_wiki = wiki_score %>% t() %>% as.data.frame()
sig_meta = metab_score %>% t() %>% as.data.frame()
sig_c4 = c4_score %>% t() %>% as.data.frame()
sig_c4 = sig_c4[, !colnames(sig_c4) %in% grep("MODULE", colnames(sig_c4), value = TRUE)]
sig_c5 = c5_score %>% t() %>% as.data.frame()
sig_c7 = c7_score %>% t() %>% as.data.frame()
sig = cbind(sig_c7, sig_meta)
# sig = cbind(sig_c5[, sample(ncol(sig_c5), 1000)], sig_hall, sig_meta)
sig = sig %>% t() %>% na.omit() %>% t() %>% as.data.frame()

## variance
gini = apply(sig, 2, gini.index)
shannon = apply(sig, 2, shannon.entropy)
metric = apply(sig, 2, metric.entropy)
pearson = apply(sig, 2, pearson.variation)
robust = apply(sig, 2, robust.variation)
entropy = data.frame(
#  gini = gini,
  shannon = shannon,
  metric = metric,
  pearson = pearson,
  robust = robust
  ) %>% t() %>% as.data.frame()

entropy = t(apply(entropy, 1, outliers))
entropy = entropy[, !apply(is.na(entropy), 2, any)]
summary(t(entropy))

entropy = entropy %>% z_score_row() %>% as.data.frame() %>% t()
# entropy[1, ] = -1 * entropy[1, ]

cluster = data.frame(
  id = colnames(entropy),
  metab = factor(c(rep('no', ncol(entropy)-152), rep('yes', 152))))

cluster = tibble::column_to_rownames(cluster, var = 'id')
bk = c(seq(-4,-0.1,by = 0.01), seq(0,4,by = 0.01))

## plot
pheatmap(
  entropy,
  cluster_rows = F, cluster_cols = F, # 是否按行或列聚类
  show_rownames = T, show_colnames = F, # 是否展示行或列名称
  annotation_col = cluster, # 添加分类图例
  scale = 'row', # 对行或列标准化
  # color = colorRampPalette(c('#0a9396','#edf2f4','#bb3e03'))(100),
  color = c(colorRampPalette(colors = c("#0a9396","#f8f9fa"))(length(bk)/2),
            colorRampPalette(colors = c("#f8f9fa","#bb3e03"))(length(bk)/2)),
  breaks = bk,
  annotation_colors = list(metab = c('yes' = '#eb5e28', 'no' = '#e0e1dd')),
  fontsize = 4,
  fontsize_col = 5,
  fontsize_row = 5,
  border = F # 不显示小格子的边框
)


com = cbind(t(as.data.frame(entropy)), cluster)
wilcox.test(com$gini ~ com$metab, exact = T,correct = T, paired = F, conf.level = 0.95)
wilcox.test(com$shannon ~ com$metab, exact = T,correct = T, paired = F, conf.level = 0.95)
wilcox.test(com$metric ~ com$metab, exact = T,correct = T, paired = F, conf.level = 0.95)
wilcox.test(com$pearson ~ com$metab, exact = T,correct = T, paired = F, conf.level = 0.95)
wilcox.test(com$robust ~ com$metab, exact = T,correct = T, paired = F, conf.level = 0.95)

df = com %>% pivot_longer(cols = c(shannon, metric, pearson, robust), names_to = 'name', values_to = 'value')
immun = df %>% filter(metab == 'no')
metab = df %>% filter(metab == 'yes')
ggplot() +
  geom_half_violin(
    data = immun, aes(x = name, y = value),
    position = position_dodge(width = 1),
    side = 'l', scale = 'width',
    colour = NA, fill = '#1ba7b3') +
  geom_half_violin(
    data = metab, aes(x = name, y = value),
    position = position_dodge(width = 1),
    side = 'r', scale = 'width',
    colour = NA, fill = '#dfb424') +
  scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(
    data = df, aes(x = name, y = value, fill = metab),
    stat = 'summary', fun = mean,
    position = position_dodge(width = 0.7)) +
  stat_summary(
    data = df, aes(x = name, y = value, fill = metab),
    fun.min = function(x){quantile(x)[2]},
    fun.max = function(x){quantile(x)[4]},
    geom = 'errorbar', color = 'black',
    width = 0.01, size = 0.5,
    position = position_dodge(width = 0.7)) +
  stat_compare_means(
    data = df, aes(x = name, y = value, fill = metab),
    label = 'p.signif')


## MPR patients
mpr = clinic[clinic$MPR == 'MPR', ]
mpr_tpm = count[, colnames(count) %in% mpr$sample_id]
table(colnames(mpr_tpm) == mpr$sample_id)
# C7
c7_lists = list()
for (pathway in unique(c7$term)) {
  c7_lists[[pathway]] = c7[c7$term == pathway, 2]}
c7_mpr = gsva(as.matrix(mpr_tpm), c7_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
c7_mpr_sig = c7_mpr %>% t() %>% as.data.frame()
# metabolism
metab_lists = list()
for (pathway in unique(metab$term)) {
  metab_lists[[pathway]] = metab[metab$term == pathway, 2]}
metab_mpr = gsva(as.matrix(mpr_tpm), metab_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
metab_mpr_sig = metab_mpr %>% t() %>% as.data.frame()
metab_mpr_sig = metab_mpr_sig %>% t() %>% na.omit() %>% t() %>% as.data.frame()
gini_mpr = apply(metab_mpr_sig, 2, gini.index)
shannon_mpr = apply(metab_mpr_sig, 2, shannon.entropy)
metric_mpr = apply(metab_mpr_sig, 2, metric.entropy)

## non-MPR patients
nmpr = clinic[clinic$MPR == 'NMPR', ]
nmpr_tpm = count[, colnames(count) %in% nmpr$sample_id]
table(colnames(nmpr_tpm) == nmpr$sample_id)
# C7
c7_lists = list()
for (pathway in unique(c7$term)) {
  c7_lists[[pathway]] = c7[c7$term == pathway, 2]}
c7_nmpr = gsva(as.matrix(nmpr_tpm), c7_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
# metabolism
metab_lists = list()
for (pathway in unique(metab$term)) {
  metab_lists[[pathway]] = metab[metab$term == pathway, 2]}
metab_nmpr = gsva(as.matrix(nmpr_tpm), metab_lists, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
metab_nmpr_sig = metab_nmpr %>% t() %>% as.data.frame()
metab_nmpr_sig = metab_nmpr_sig %>% t() %>% na.omit() %>% t() %>% as.data.frame()
gini_nmpr = apply(metab_nmpr_sig, 2, gini.index)
shannon_nmpr = apply(metab_nmpr_sig, 2, shannon.entropy)
metric_nmpr = apply(metab_nmpr_sig, 2, metric.entropy)

## compare
wilcox.test(gini_mpr, gini_nmpr, paired = T, alternative = 'two.sided', exact = F, correct = F)
wilcox.test(shannon_mpr, shannon_nmpr, paired = T, alternative = 'two.sided', exact = F, correct = F)
wilcox.test(metric_mpr, metric_nmpr, paired = T, alternative = 'two.sided', exact = F, correct = F)

## plot
gini_mpr = melt(gini_mpr, value.name = "MPR")
shannon_mpr = melt(shannon_mpr, value.name = "MPR")
metric_mpr = melt(metric_mpr, value.name = "MPR")
gini_nmpr = melt(gini_nmpr, value.name = "NMPR")
shannon_nmpr = melt(shannon_nmpr, value.name = "NMPR")
metric_nmpr = melt(metric_nmpr, value.name = "NMPR")

colnames(gini_mpr) = 'gini'
gini_mpr = mutate(gini_mpr, response = 'MPR')
gini_mpr = tibble::rownames_to_column(gini_mpr, var = 'pathway')
colnames(gini_nmpr) = 'gini'
gini_nmpr = mutate(gini_nmpr, response = 'NMPR')
gini_nmpr = tibble::rownames_to_column(gini_nmpr, var = 'pathway')
df = rbind(gini_mpr, gini_nmpr)
df = subset(df, gini <= 40)
ggplot(df, aes(x = response, y = gini)) +
  geom_point() +
  geom_line(aes(group = pathway)) +
  theme_minimal() +
  labs(x = "Gini Coefficient (Dataset 1)", y = "Gini Coefficient (Dataset 2)", title = "Comparison of Gini Coefficients for Pathways")

colnames(shannon_mpr) = 'shannon'
shannon_mpr = mutate(shannon_mpr, response = 'MPR')
shannon_mpr = tibble::rownames_to_column(shannon_mpr, var = 'pathway')
colnames(shannon_nmpr) = 'shannon'
shannon_nmpr = mutate(shannon_nmpr, response = 'NMPR')
shannon_nmpr = tibble::rownames_to_column(shannon_nmpr, var = 'pathway')
df = rbind(shannon_mpr, shannon_nmpr)
# df = subset(df, gini <= 40)
ggplot(df, aes(x = response, y = shannon)) +
  geom_point() +
  geom_line(aes(group = pathway)) +
  theme_minimal() +
  labs(x = "Gini Coefficient (Dataset 1)",
       y = "Gini Coefficient (Dataset 2)", title =
         "Comparison of Gini Coefficients for Pathways")

colnames(metric_mpr) = 'metric'
metric_mpr = mutate(metric_mpr, response = 'MPR')
metric_mpr = tibble::rownames_to_column(metric_mpr, var = 'pathway')
colnames(metric_nmpr) = 'metric'
metric_nmpr = mutate(metric_nmpr, response = 'NMPR')
metric_nmpr = tibble::rownames_to_column(metric_nmpr, var = 'pathway')
df = rbind(metric_mpr, metric_nmpr)
# df = subset(df, gini <= 40)
ggplot(df, aes(x = response, y = metric)) +
  geom_point() +
  geom_line(aes(group = pathway)) +
  theme_minimal() +
  labs(x = "Gini Coefficient (Dataset 1)",
       y = "Gini Coefficient (Dataset 2)",
       title = "Comparison of Gini Coefficients for Pathways")

