rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','plot1cell','Seurat','dior','RColorBrewer','ggplot2','ggunchull','viridis')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
sce_scissors = dior::read_h5(file = 'subset_data.h5', target.object = 'seurat')
head(sce_scissors@meta.data)

Idents(sce_scissors) = sce_scissors$all_cluster_annotation2
circ = prepare_circlize_data(sce_scissors, scale = 0.8)
print(colnames(sce_scissors@meta.data))
print(colnames(circ))
cluster_colors = rand_color(length(names(table(sce_scissors$all_cluster_annotation2))))
scissor_colors = c('#e5e5e5','#14213d','#fca311')

## plot
df = as.data.frame(sce_scissors[["umap"]]@cell.embeddings)
df$cluster = sce_scissors@meta.data$all_cluster_annotation2
ggplot(df, aes(x = UMAP_1, y = UMAP_2, fill = cluster, color = cluster)) +
  scale_color_manual(values = c(
    'Plasma cell' = '#D1352B',
    'B cell' = '#E887BD',
    'Epithelium' = '#3C76AF',
    'Fibroblast/Endothelium' = '#8FA4AE',
    'Mast cell' = '#F5CEE4',
    'Myeloid cell' = '#68AD57',
    'Neutrophil' = '#B382B9',
    'NK cell' = '#F6E77E',
    'pDC' = '#D0BAD8',
    'T cell' = '#EE924E')) +
  geom_point(size = 0.1) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
  )

## circlize
plot_circlize(
  circ, do.label = F, repel = F,
  col.use = cluster_colors,
  bg.color = 'white',
  label.cex = 1,
  pt.size = 0,
  kde2d.n = 300
)
add_track(
  circ, group = 'scissor',
  colors = scissor_colors, track_num = 2
)

