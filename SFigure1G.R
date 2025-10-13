rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','openxlsx','survival','survminer','forestplot','GSVA')
for (i in ps) {library(i, character.only = TRUE)}
source('currency.R')

## data
oak_tpm = read.table('OAK_TPM.txt', header = T, row.names = 1)
oak_surv = read.csv('OAK_survival_ORR.csv', header = T)
oak_surv$Sample = gsub('-', '.', oak_surv$Sample)
pop_tpm = read.table('POPLAR_TPM.txt', header = T, row.names = 1)
pop_surv = read.csv('POPLAR_survival_ORR.csv', header = T)
pop_surv$Sample = gsub('-', '.', pop_surv$Sample)

# oak_tide = apply(oak_tpm, 2, as.numeric)
# rownames(oak_tide) = rownames(oak_tpm)
# oak_tide = t(apply(oak_tide, 1, function(x){x-(mean(x))}))
# write.table(oak_tide, file = 'OAK_TIDE.txt', sep = "\t", quote = F, row.names = T)
# pop_tide = apply(pop_tpm, 2, as.numeric)
# rownames(pop_tide) = rownames(pop_tpm)
# pop_tide = t(apply(pop_tide, 1, function(x){x-(mean(x))}))
# write.table(pop_tide, file = 'POPLAR_TIDE.txt', sep = "\t", quote = F, row.names = T)

## signature
# Homologous recombination deficiency (HRD)can predict the therapeutic outcomes of immuno-neoadjuvant therapy in NSCLC patients
neo = c('MZB1','JCHAIN','SLAMF7','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHGP','IGLC2','IGLC3','IGLC7')
adv = c('MZB1','DERL3','JSRP1','TNFRSF17','SLAMF7','IGHG2','IGHGP','IGLV3-1','IGLV6-57','IGHA2','IGKV4-1','IGKV1-12','IGLC7','IGLL5')
cd8 = c('IFNG','CXCL9','CD8A','GZMA','GZMB','CXCL10','PRF1','TBX21')
ddr = c('CXCL10','MX1','IDO1','IFI44L','CD2','GBP5','PRAME','ITGAL','LRP4','APOL3','CDR1','FYB','TSPAN7','RAC2','KLHDC7B','GRB14','AC138128.1','KIF26A','CD274','CD109','ETV7','MFAP5','OLFM4','Pl15','FOSB','FAM19A5','NLRC5','PRICKLE1','EGR1','CLDN10','ADAMTS4','SP140L','ANXA1','RSAD2','ESR1','IKZF3','OR2I1P','EGFR','NAT1','LATS2','CYP2B6','PTPRC','PPP1R1A','AL137218.1')
ifn = c('CD8A','CCL5','CD27','CD274','PDCD1LG2','CD276','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PSMB10','STAT1','TIGIT')
signatures = list(neo = neo, adv = adv, cd8 = cd8, ddr = ddr, ifn = ifn)

oak_immun = oak_surv[oak_surv$ACTARM == 'MPDL3280A', ]
oak_tpm = oak_tpm[, colnames(oak_tpm) %in% oak_immun$Sample]
oak_info = oak_surv[match(colnames(oak_tpm), oak_surv$Sample), ]
table(colnames(oak_tpm) == oak_info$Sample)
pop_immun = pop_surv[pop_surv$ACTARM == 'MPDL3280A', ]
pop_tpm = pop_tpm[, colnames(pop_tpm) %in% pop_immun$Sample]
pop_info = pop_surv[match(colnames(pop_tpm), pop_surv$Sample), ]
table(colnames(pop_tpm) == pop_info$Sample)

## OAK
oak_plasma = gsva(as.matrix(oak_tpm), plasma, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
oak_plasma = t(oak_plasma)
colnames(oak_plasma) = 'plasma'
table(rownames(oak_plasma) == oak_info$Sample)
oak_info = cbind(oak_info, oak_plasma)
oak_info$group = ifelse(oak_info$plasma > median(oak_info$plasma), 'high', 'low')

## OAK-OS
oak_osfit = survfit(Surv(OS_MONTHS, OS_CENSOR)~group, data = oak_info)
ggsurvplot(oak_osfit,
           surv.scale = 'percent', surv.median.line = 'hv',
           size = 1.5, linetype = c(1,1), palette = c("#B3252D","#1C6198"),
           pval.size = 5, pval = TRUE, conf.int = TRUE,
           risk.table.fontsize = 5, risk.table.col = 'strata', risk.table = TRUE, risk.table.y.text.col = TRUE,
           xlim = c(-2, 37), ylim = c(-0.05, 1.05), break.time.by = 7,
           xlab = 'Time in months', ylab = 'Overall survival',
           legend = c(0.15, 0.15),
           legend.title = '',
           legend.labs = c('high', 'low'),
           legend.text = element_text(size = 9),
           axes.offset = FALSE,
           ggtheme = theme_bw() +
             theme(panel.grid = element_blank(),
                   legend.title = element_blank()),
           tables.theme = theme_bw() +
             theme(panel.grid = element_blank(),
                   panel.border = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank()))
## OAK-PFS
oak_pfsfit = survfit(Surv(PFS_MONTHS, PFS_CENSOR)~group, data = oak_info)
ggsurvplot(oak_pfsfit,
           surv.scale = 'percent', surv.median.line = 'hv',
           size = 1.5, linetype = c(1,1), palette = c("#B3252D","#1C6198"),
           pval.size = 5, pval = TRUE, conf.int = TRUE,
           risk.table.fontsize = 5, risk.table.col = 'strata', risk.table = TRUE, risk.table.y.text.col = TRUE,
           xlim = c(-2, 33), ylim = c(-0.05, 1.05), break.time.by = 6,
           xlab = 'Time in months', ylab = 'Overall survival',
           legend = c(0.15, 0.15),
           legend.title = '',
           legend.labs = c('high', 'low'),
           legend.text = element_text(size = 9),
           axes.offset = FALSE,
           ggtheme = theme_bw() +
             theme(panel.grid = element_blank(),
                   legend.title = element_blank()),
           tables.theme = theme_bw() +
             theme(panel.grid = element_blank(),
                   panel.border = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank()))
write.csv(oak_info, file = 'OAK-survival.csv')

## POPLAR
pop_plasma = gsva(as.matrix(pop_tpm), plasma, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
pop_plasma = t(pop_plasma)
colnames(pop_plasma) = 'plasma'
table(rownames(pop_plasma) == pop_info$Sample)
pop_info = cbind(pop_info, pop_plasma)
pop_info$group = ifelse(pop_info$plasma > median(pop_info$plasma), 'high', 'low')

## POPLAR-OS
pop_osfit = survfit(Surv(OS_MONTHS, OS_CENSOR)~group, data = pop_info)
ggsurvplot(pop_osfit,
           surv.scale = 'percent', surv.median.line = 'hv',
           size = 1.5, linetype = c(1,1), palette = c("#B3252D","#1C6198"),
           pval.size = 5, pval = TRUE, conf.int = TRUE,
           risk.table.fontsize = 5, risk.table.col = 'strata', risk.table = TRUE, risk.table.y.text.col = TRUE,
           xlim = c(-2, 50), ylim = c(-0.05, 1.05), break.time.by = 9,
           xlab = 'Time in months', ylab = 'Overall survival',
           legend = c(0.15, 0.15),
           legend.title = '',
           legend.labs = c('high', 'low'),
           legend.text = element_text(size = 9),
           axes.offset = FALSE,
           ggtheme = theme_bw() +
             theme(panel.grid = element_blank(),
                   legend.title = element_blank()),
           tables.theme = theme_bw() +
             theme(panel.grid = element_blank(),
                   panel.border = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank()))
## POPLAR-PFS
pop_pfsfit = survfit(Surv(PFS_MONTHS, PFS_CENSOR)~group, data = pop_info)
ggsurvplot(pop_pfsfit,
           surv.scale = 'percent', surv.median.line = 'hv',
           size = 1.5, linetype = c(1,1), palette = c("#B3252D","#1C6198"),
           pval.size = 5, pval = TRUE, conf.int = TRUE,
           risk.table.fontsize = 5, risk.table.col = 'strata', risk.table = TRUE, risk.table.y.text.col = TRUE,
           xlim = c(-2, 44), ylim = c(-0.05, 1.05), break.time.by = 7,
           xlab = 'Time in months', ylab = 'Overall survival',
           legend = c(0.15, 0.15),
           legend.title = '',
           legend.labs = c('high', 'low'),
           legend.text = element_text(size = 9),
           axes.offset = FALSE,
           ggtheme = theme_bw() +
             theme(panel.grid = element_blank(),
                   legend.title = element_blank()),
           tables.theme = theme_bw() +
             theme(panel.grid = element_blank(),
                   panel.border = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank()))
write.csv(pop_info, file = 'POPLAR-survival.csv')

## OAK multi-factor cox
oak_tide = read.csv('OAK_TIDE.csv', header = T)
oak_tide = oak_tide[match(oak_info$Sample, oak_tide$Patient), ]
table(oak_info$Sample == oak_tide$Patient)
colnames(oak_tide)[2:3] = c('benefit', 'tide')

oak_sig = gsva(as.matrix(oak_tpm), signatures, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
oak_sig = t(oak_sig)
table(rownames(oak_sig) == oak_info$Sample)
oak_pdl1 = t(oak_tpm[rownames(oak_tpm) == 'CD274', ])
table(rownames(oak_pdl1) == oak_info$Sample)
oak_info = cbind(oak_info, oak_sig, oak_pdl1, oak_tide[3])

oak_info$plasma_group = ifelse(oak_info$plasma > median(oak_info$plasma), '1', '0')
oak_info$prev_group = ifelse(oak_info$prev > median(oak_info$prev), '1', '0')
oak_info$cd8t_group = ifelse(oak_info$cd8t > median(oak_info$cd8t), '1', '0')
oak_info$ddir_group = ifelse(oak_info$ddir > median(oak_info$ddir), '1', '0')
oak_info$ifng_group = ifelse(oak_info$ifng > median(oak_info$ifng), '1', '0')
oak_info = oak_info %>% mutate(tps_group = case_when(
  PDL1_TC_22C3 == 'NE' ~ NA,
  PDL1_TC_22C3 == '[0,1)' ~ 0,
  PDL1_TC_22C3 %in% c('[1,50)','[50,100]') ~ 1))
oak_info = oak_info %>% mutate(tmb_group = case_when(
  tTMB == 'NE' ~ NA,
  tTMB == '<16' ~ 0,
  tTMB == '>=16' ~ 1))
oak_info$tide_group = ifelse(oak_info$tide == 'True', '1', '0')

unicox = unicoxph(
  data = oak_info,
  variates = grep('_group', colnames(oak_info), value = TRUE),
  os = 'OS_MONTHS', status = 'OS_CENSOR')

multicox = coxph(
  Surv(OS_MONTHS, OS_CENSOR) ~ 
    plasma_group + cd8t_group + ifng_group + ddir_group + tide_group,
  data = oak_info)
summary(multicox)


## POPLAR multi-factor cox
pop_tide = read.csv('POPLAR_TIDE.csv', header = T)
pop_tide = pop_tide[match(pop_info$Sample, pop_tide$Patient), ]
table(pop_info$Sample == pop_tide$Patient)
colnames(pop_tide)[2:3] = c('benefit', 'tide')

pop_sig = gsva(as.matrix(pop_tpm), signatures, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
pop_sig = t(pop_sig)
table(rownames(pop_sig) == pop_info$Sample)
pop_pdl1 = t(pop_tpm[rownames(pop_tpm) == 'CD274', ])
table(rownames(pop_pdl1) == pop_info$Sample)
pop_info = cbind(pop_info, pop_sig, pop_pdl1, pop_tide[3])

pop_info$plasma_group = ifelse(pop_info$plasma > median(pop_info$plasma), '1', '0')
pop_info$prev_group = ifelse(pop_info$prev > median(pop_info$prev), '1', '0')
pop_info$cd8t_group = ifelse(pop_info$cd8t > median(pop_info$cd8t), '1', '0')
pop_info$ddir_group = ifelse(pop_info$ddir > median(pop_info$ddir), '1', '0')
pop_info$ifng_group = ifelse(pop_info$ifng > median(pop_info$ifng), '1', '0')
pop_info$tide_group = ifelse(pop_info$tide == 'True', '1', '0')

unicox = unicoxph(
  data = pop_info,
  os = 'OS_MONTHS', status = 'OS_CENSOR',
  variates = grep('_group', colnames(pop_info), value = TRUE))

multicox = coxph(
  Surv(OS_MONTHS, OS_CENSOR) ~ 
    plasma_group + cd8t_group + ifng_group + ddir_group + tide_group,
  data = pop_info)
summary(multicox)

## Interaction P
intermodel = glm(
  OS_CENSOR ~ plasma_group + cd8t_group + plasma_group * cd8t_group,
  family = binomial,
  data = oak_info)
summary(intermodel)

## forestplot
forest = read.csv('forest.csv', header = F)
for (i in 4:9) {forest[,i] = as.numeric(forest[,i])}
forestplot(
  forest[,c(1:3)],
  mean = cbind(forest[,7], forest[,4]),
  lower = cbind(forest[,8], forest[,5]),
  upper = cbind(forest[,9], forest[,6]),
  graph.pos = 4, 
  graphwidth = unit(0.3, 'npc'),
  colgap = unit(2, 'mm'),
  lineheight = 'auto',
  zero = 1, clip = c(0, 3), xticks = c(0, 1, 2, 3),
  hrzl_lines = list(
    '1' = gpar(lty = 1, lwd = 2),
    '2' = gpar(lty = 2)),
  col = fpColors(
    zero = 'grey50',
    box = c('#EE5826', '#016BB4'),
    lines = c('#EE5826', '#016BB4')),
  grid = structure(
    c(0, 1, 2, 3),
    gp = gpar(col = 'grey', lty = 2, lwd = 1)),
  lwd.ci = 2,
  boxsize = 0.2,
  line.margin = 0.8,
  ci.vertices = T,
  ci.vertices.height = 0.1
)
