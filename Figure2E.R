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

# Homologous recombination deficiency (HRD)can predict the therapeutic outcomes of immuno-neoadjuvant therapy in NSCLC patients
akr_set = list(c('AKR1C1','AKR1C2','AKR1C3'))

# OAK
oak_immun = oak_surv[oak_surv$ACTARM == 'MPDL3280A', ]
oak_tpm = oak_tpm[, colnames(oak_tpm) %in% oak_immun$Sample]
oak_info = oak_surv[match(colnames(oak_tpm), oak_surv$Sample), ]
table(colnames(oak_tpm) == oak_info$Sample)
# POPLAR
pop_immun = pop_surv[pop_surv$ACTARM == 'MPDL3280A', ]
pop_tpm = pop_tpm[, colnames(pop_tpm) %in% pop_immun$Sample]
pop_info = pop_surv[match(colnames(pop_tpm), pop_surv$Sample), ]
table(colnames(pop_tpm) == pop_info$Sample)

# OAK-signature
oak_akr = gsva(as.matrix(oak_tpm), akr_set, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
oak_akr = t(oak_akr)
colnames(oak_akr) = 'akr'
table(rownames(oak_akr) == oak_info$Sample)
oak_info = cbind(oak_info, oak_akr)
oak_info$group = ifelse(oak_info$akr > median(oak_info$akr), 'high', 'low')

# OAK-OS
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
# OAK-PFS
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

# POPLAR-signature
pop_akr = gsva(as.matrix(pop_tpm), akr_set, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
pop_akr = t(pop_akr)
colnames(pop_akr) = 'akr'
table(rownames(pop_akr) == pop_info$Sample)
pop_info = cbind(pop_info, pop_akr)
pop_info$group = ifelse(pop_info$akr > median(pop_info$akr), 'high', 'low')

# POPLAR-OS
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
# POPLAR-PFS
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
