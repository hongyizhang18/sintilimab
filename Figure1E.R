rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','clusterProfiler','enrichplot','ggplot2','gg.gap','glmnet','GSVA')
for (i in ps) {library(i, character.only = TRUE)}

## signature
adv = c('MZB1','DERL3','JSRP1','TNFRSF17','SLAMF7','IGHG2','IGHGP','IGLV3-1','IGLV6-57','IGHA2','IGKV4-1','IGKV1-12','IGLC7','IGLL5')
neo = c('MZB1','JCHAIN','SLAMF7','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHGP','IGLC2','IGLC3','IGLC7')
cd8 = c('IFNG','CXCL9','CD8A','GZMA','GZMB','CXCL10','PRF1','TBX21')
ddr = c('CXCL10','MX1','IDO1','IFI44L','CD2','GBP5','PRAME','ITGAL','LRP4','APOL3','CDR1','FYB','TSPAN7','RAC2','KLHDC7B','GRB14','AC138128.1','KIF26A','CD274','CD109','ETV7','MFAP5','OLFM4','Pl15','FOSB','FAM19A5','NLRC5','PRICKLE1','EGR1','CLDN10','ADAMTS4','SP140L','ANXA1','RSAD2','ESR1','IKZF3','OR2I1P','EGFR','NAT1','LATS2','CYP2B6','PTPRC','PPP1R1A','AL137218.1')
ifn = c('CD8A','CCL5','CD27','CD274','PDCD1LG2','CD276','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PSMB10','STAT1','TIGIT')
genes = Reduce(union, list(adv, neo, cd8, ddr, ifn))

## data
count_48 = read.csv('tpm_matrix.csv', row.names = 1)
clinic_48 = read.csv('clinical_48.csv', row.names = 1)
count = count_48[, colnames(count_48) %in% clinic$sample_id]
clinic = clinic_48[clinic_48$sample_time == 'baseline', ]

versus = count[genes, ] %>% na.omit()
result = apply(versus, 1, function(x) {t.test(x ~ clinic$MPR)})
pval = sapply(result, function(x) x$p.value)

plasma = list(adv = adv, neo = neo, cd8 = cd8, ddr = ddr, ifn = ifn)
plasma_score = gsva(as.matrix(count), plasma, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
plasma_score = as.data.frame(t(plasma_score))

# immune = read.csv('TIMER2_immune_cell.csv')
# immune = immune %>% filter(grepl("plasma", immune$cell_type)) %>% tibble::column_to_rownames('cell_type') %>% t() %>% as.data.frame()
# immune = immune[rownames(immune) %in% clinic$sample_id, ]

df = plasma_score %>%
  z_score_col() %>%
  cbind(MPR = clinic$MPR) %>%
  as.data.frame()
df[1:8] = apply(df[1:8], 2, as.numeric)

plot = df %>% arrange(desc(neo))
plot$num = c(1:nrow(plot))
p1 = ggplot(plot, aes(num, neo, fill = MPR)) +
  geom_bar(stat = 'identity',
           colour = '#1b263b',
           position = 'dodge',
           width = 0.8) +
  scale_fill_manual(values = c('#D75640','#338EC1')) +
  scale_y_continuous(limits = c(-2.5, 1),
                     breaks = c(-2.5, -2, -1, 0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_line(color = "black",size = 0.5),
        axis.title.y = element_text(color = "black",size = 12,vjust = 1),
        axis.line.y = element_line(color = "black",size = 0.5,lineend = "butt"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
gg.gap(plot = p1,
       segments = c(-2, -1),
       tick_width = 1,
       rel_heights = c(0.2, 0, 1),
       ylim = c(-2.5, 1.0))

result = apply(plot[, 1:8], 2, function(x) {t.test(x ~ plot$MPR)})
pval = sapply(result, function(x) x$p.value)
pval
ggplot(plot, aes(x = MPR, y = plasma8)) + 
  geom_boxplot(size = 0.8, fill = c('#D75640','#338EC1'), outlier.color = 'white') +
  scale_y_continuous(limits = c(-3, 1),
                     breaks = c(-3, -2, -1, 0, 1)) +
  theme_bw() +
  theme(axis.line = element_line(size = 0.5, colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(vjust = 1),
        axis.title.x = element_text(vjust = -1))


gsea_files = list.files(pattern = "NO_MPR.*\\.txt$")
gsea_list = lapply(gsea_files, read.table, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
nmpr = do.call(rbind, gsea_list)
nmpr_filtered = filter(nmpr, abs(NES) > 1.8)
nmpr_filtered = nmpr_filtered[!grepl("CHR", nmpr_filtered$NAME), ]
pathway = nmpr_filtered$NAME

ddr_genesets = read.csv('/Users/jonas/Desktop/senior/2024-03-neoadjuvant/DDR_genesets.csv')
ddr_sig = lapply(ddr_genesets, function(x) x[x != ''])

ddr_score = gsva(as.matrix(count), ddr_sig, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
pval = apply(ddr_score, 1, test_loop)

