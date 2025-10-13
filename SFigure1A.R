rm(list = ls())
set.seed(123)
setwd('/Users/jonas/Desktop/')

ps = c('tidyverse','pheatmap','openxlsx')
for (i in ps) {library(i, character.only = TRUE)}

## data
clinic = read.csv('clinical_83.csv', header = T, row.names = 1)
clinic$baseline = ifelse(clinic$sample_time == 'baseline', clinic$sample_type, NA)
clinic$surgery = ifelse(clinic$sample_time == 'surgery', clinic$sample_type, NA)
sample1 = read.xlsx('clinical_64.xlsx')
sample2 = read.csv('clinical_62.csv', header = T, row.names = 1)
sample_dup = setdiff(sample1$ID, sample2$sample_id)
patient_dup = clinic %>% group_by(patient_name) %>% filter(n() > 1)
patient = clinic[!clinic$sample_id %in% patient_dup$sample_id, ]
patient_df = data.frame()
for (id in unique(patient_dup$patient_name)) {
  df = patient_dup[patient_dup$patient_name == id, ]
  row = sapply(colnames(df), function(x) {
    var = unique(df[[x]])
    if (length(var) == 1) {
      return(var)
    } else {
      return(paste(var, collapse = ", "))
    }
  })
  patient_df = rbind(patient_df, as.data.frame(t(row), stringsAsFactors = FALSE))
}
patient = rbind(patient, patient_df)
write.csv(patient, file = 'patient_66.csv', fileEncoding = "GB18030")

df = read.csv('patient_66_rna.csv', header = T, fileEncoding = "GB18030")
# df = mutate(df, age_group = ifelse(age >= 65, '≥65', '<65'))
# anno = select(df, c('patient_name','surgery','baseline','pfs_status','os_status','response','pdl1','smoking','pathology.x','cTNM','age_group'))
anno = dplyr::select(df, c('patient_name','WES','surgery','baseline','pfs_status','os_status','pdl1','response'))
anno = anno %>% tibble::column_to_rownames(var = 'patient_name')
anno_order = anno %>% arrange(response, baseline, surgery)
matrix = matrix(runif(2 * 66), nrow = 2, ncol = 66)
colnames(matrix) = rownames(anno_order)

## plot
pheatmap(
  scale(matrix), annotation_col = anno_order,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  color = colorRampPalette(colors = c("blue","white","red"))(100),
  annotation_colors = list(
    response = c('CPR' = '#ed6a5a', 'MPR' = '#ffd166', 'NO-MPR' = '#5ca4a9'),
    pdl1 = c('0%' = '#dec9e9', '1-49%' = '#b43e8f', '>50%' = '#742294'),
    os_status = c('0' = '#a8dadc', '1' = '#457b9d'),
    pfs_status = c('0' = '#a8dadc', '1' = '#457b9d'),
    baseline = c('LN' = '#00a5cf', 'T' = '#5dd39e', 'T, LN' = '#f3de2c'),
    surgery = c('LN' = '#00a5cf', 'T' = '#5dd39e', 'T, LN' = '#f3de2c'),
    WES = c('1' = '#cf5cc2')),
  fontsize = 4,
  fontsize_col = 5,
  fontsize_row = 5,
  border = T)

