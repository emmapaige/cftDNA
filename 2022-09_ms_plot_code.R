
#### 2022-09-12 Simplified code for creating plots

set.seed(42)

library(here)
library(R.utils)
library(tidyverse)
library(tidytext)
library(readxl)

library(dgof)
library(glmnet)
library(survival)
library(survminer)
# library(UpSetR)

library(grid)
library(gridExtra)
library(reshape2)

###################################################################################################################################################################################################
###################################################################################################################################################################################################

#### my functions

quart <- function(x) {
  x <- sort(x)
  n <- length(x)
  m <- (n+1)/2
  if (floor(m) != m) {
    l <- m-1/2; u <- m+1/2
  } else {
    l <- m-1; u <- m+1
  }
  c(Q1 = median(x[1:l]), Q2 = median(x), Q3 = median(x[u:n]))
}

normalize.vector <- function(vector) {
  
  if ( near(max(vector), 0) & near(min(vector), 0) ) {
    return(vector)
  } else { 
    return((vector - min(vector))/(max(vector) - min(vector)))
  }
}

my_homogeneity_test <- function(mat) {
  
  groups <- rownames(mat)
  
  for (i in 1:(nrow(mat) - 1)) {
    for (j in (i+1):nrow(mat)) {
      
      dummy.x <- mat[i, ]
      dummy.y <- mat[j, ]
      
      #### sometimes we might not have all bins, but the test is still meaningful
      #### for that particular pair
      nonzero.x <- which(dummy.x > 0)
      nonzero.y <- which(dummy.y > 0)
      
      gd.idcs <- base::union(nonzero.x, nonzero.y)
      dof <- length(gd.idcs) - 1
      
      if (dof > 0) {
        dummy.x <- dummy.x[gd.idcs]
        dummy.y <- dummy.y[gd.idcs]
        
        n1 <- sum(dummy.x)
        n2 <- sum(dummy.y)
        N <- n1 + n2
        
        empirical_p_vec <- (dummy.x + dummy.y) / N
        
        P_1 <- n1 * empirical_p_vec
        Q_1 <- sum( (dummy.x - P_1)^2 / P_1 )
        
        P_2 <- n2 * empirical_p_vec
        Q_2 <- sum( (dummy.y - P_2)^2 / P_2 )
        
        Q <- Q_1 + Q_2
        
        chi.out <- list(p.value = pchisq(Q, df = dof, lower.tail = FALSE))
        
        #### I don't think below is working correctly,
        #### possibly because it is testing for independence instead of homogeneity?
        # chi.out <- chisq.test(x = dummy.x, y = dummy.y, simulate.p.value = TRUE, B = 2000)
        
      } else {
        chi.out <- list(p.value = NA)
      }
      
      
      if (i == 1 & j == 2) {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           p_value = chi.out$p.value)
      } else {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           p_value = chi.out$p.value) %>% 
          bind_rows(dummy_tb)
      }
    }
  }
  
  dummy_tb <- dummy_tb %>% 
    arrange(p_value)
  
  return(dummy_tb)
}

find.pairwise.group.mutation.delta <- function(mat) {
  
  groups <- rownames(mat)
  
  for (i in 1:(nrow(mat) - 1)) {
    for (j in (i+1):nrow(mat)) {
      
      dummy.x <- mat[i, ]
      dummy.y <- mat[j, ]
      
      #### focus on mutations where at least one vector has the mutation
      nonzero.x <- which(dummy.x > 0)
      nonzero.y <- which(dummy.y > 0)
      
      gd.idcs <- base::union(nonzero.x, nonzero.y)
      
      delta.vec <- dummy.x[gd.idcs] - dummy.y[gd.idcs]
      
      if (i == 1 & j == 2) {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           mutation = colnames(mat)[gd.idcs],
                           delta = delta.vec)
      } else {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           mutation = colnames(mat)[gd.idcs],
                           delta = delta.vec) %>% 
          bind_rows(dummy_tb)
      }
    }
  }
  
  dummy_tb <- dummy_tb %>% 
    unite(col = 'comparison', group_1:group_2, sep = '_vs_', remove = FALSE) %>% 
    arrange(abs(delta))
  
  return(dummy_tb)
}

###################################################################################################################################################################################################
###################################################################################################################################################################################################

#### load data

full_data <- read_xlsx(here("20190118_GNE01_Cumulative_Results_Report_annotated.xlsx")) %>% 
  separate(Customer.SampleId, 
           c("customer", "sampleId"), 
           sep = "-", 
           remove = FALSE,
           extra = "merge") %>% 
  filter(!(customer %in% c("SAAWZH", "SAAEFC", "SAAFXZ", "SAAGMX", "SAAGVE", "SAAKTZ", "SAANDA", "SAAPNI", "SAAPUL"))) %>% 
  select(-timepoint) %>% 
  rename(PFS_time = "Progression.Free.Survival.Time",
         timepoint = "Timepoint") %>% 
  filter(!is.na(PFS_time) & (PFS_time < 80)) %>% 
  mutate(timepoint = parse_factor(timepoint,
                                  levels = c("PRE", "OT", "EOS")),
         ordered = TRUE)
#### start with 368 unique customers
#### filter out customer based on Federico's data cleaning
#### email (see 2021-04-19 email)
#### 2 customers removed based on PFS_time (declared outliers by Federico and Corbin during EDA)
#### 66 customers have NA for PFS_time
#### end with 295, so some of the customers from Federico's list are already not in the excel data

###################################################################
###################################################################

#### make some of the columns easier for downstream work on strings 
data_tb <- full_data %>% 
  mutate(Mut_aa = if_else(Variant_type == "CNV",
                          "CNV",
                          Mut_aa),
         Mut_aa = if_else(Variant_type == "NONE",
                          "NONE",
                          Mut_aa),
         Mut_aa = if_else(Variant_type == "Indel",
                          "Indel",
                          Mut_aa),
         Mut_aa = if_else((Variant_type == "SNV" & is.na(Mut_aa)),
                          VAR_TYPES,
                          Mut_aa)) %>% 
  mutate(VAR_TYPES = if_else((is.na(VAR_TYPES) & Variant_type == "CNV"),
                             "CNV",
                             VAR_TYPES),
         VAR_TYPES = if_else((is.na(VAR_TYPES) & Gene == "NONE"),
                             "NONE",
                             VAR_TYPES)) %>% 
  unite(gene_mutation, Gene, Mut_aa, 
        sep = "-", remove = FALSE,
        na.rm = FALSE) %>% 
  mutate(gene_mutation = if_else(Variant_type == "NONE",
                                 "NONE",
                                 gene_mutation))

pre_customers <- data_tb %>% 
  filter(timepoint == "PRE") %>% 
  pull(customer) %>% 
  unique()

ot_customers <- data_tb %>% 
  filter(timepoint == "OT") %>% 
  pull(customer) %>% 
  unique()

eos_customers <- data_tb %>% 
  filter(timepoint == "EOS") %>% 
  pull(customer) %>% 
  unique()

all_timepoint_customers <- intersect(pre_customers, intersect(ot_customers, eos_customers))
# length(all_timepoint_customers)##[1] 190

data_tb <- data_tb %>% 
  filter(customer %in% all_timepoint_customers)

quart_progress <- data_tb %>% 
  select(customer, PFS_time) %>% 
  unique() %>% 
  pull(PFS_time) %>% 
  quart()

alt_quart_progress <- data_tb %>% 
  select(customer, PFS_time) %>% 
  unique() %>% 
  mutate(PFS_time = log2(.5 + PFS_time)) %>% 
  pull(PFS_time) %>% 
  quart()

#### do <= for Q3 to make PFS quartile symmetric, i.e., 
#### more in Q2 and Q3 than Q1 and Q4 instead of more in Q2 and Q4 than Q1 and Q3
cox_data_tb <- data_tb %>% 
  mutate(tumor_class = "Q4",
         tumor_class = if_else(PFS_time <= quart_progress[3],
                               "Q3",
                               tumor_class),
         tumor_class = if_else(PFS_time < quart_progress[2],
                               "Q2",
                               tumor_class),
         tumor_class = if_else(PFS_time < quart_progress[1],
                               "Q1",
                               tumor_class),
         gene_mutation = if_else(is.na(gene_mutation),
                                 "NONE",
                                 gene_mutation)) %>% 
  select(customer, sampleId, timepoint, tumor_class, PFS_time, PFS.Status,
         Cancertype, Arm, Variant_type, VAR_TYPES,
         gene_mutation, Gene, Mut_aa, 
         Copy_number, Indel_type, Variant_consequence,
         cfDNA, Percentage, max_maf,
         Chromosome, position, Exon, Mut_nt) %>% 
  rename(PFS_status = PFS.Status, copy_number = Copy_number,
         Cancer_type = Cancertype)

PFS_time_tb <- cox_data_tb %>% 
  select(customer, PFS_time) %>% 
  unique()

mutation_load_tb <- cox_data_tb %>% 
  mutate(load = if_else(Gene != 'NONE',
                        1,
                        0)) %>% 
  group_by(customer, timepoint, PFS_status, PFS_time, tumor_class, Arm, Variant_type) %>% 
  summarize(variant_mutation_load = sum(load)) %>% 
  group_by(customer, timepoint) %>% 
  mutate(mutation_load = sum(variant_mutation_load),
         fraction_variant_mutation_load = variant_mutation_load / mutation_load) %>% 
  ungroup() %>% 
  mutate(count_bin = 1,
         count_bin = if_else(mutation_load >= 1,
                             2,
                             count_bin),
         count_bin = if_else(mutation_load >= 3,
                             3,
                             count_bin),
         count_bin = if_else(mutation_load >= 5,
                             4,
                             count_bin),
         count_bin = if_else(mutation_load >= 7,
                             5,
                             count_bin),
         count_bin = if_else(mutation_load >= 9,
                             6,
                             count_bin),
         count_bin = if_else(mutation_load >= 11,
                             7,
                             count_bin))

load_change_tb <- mutation_load_tb %>% 
  select(customer, timepoint, Arm, PFS_time, PFS_status, tumor_class, mutation_load) %>% 
  unique() %>% 
  pivot_wider(names_from = timepoint, values_from = mutation_load) %>% 
  mutate(delta_1 = OT - PRE,
         delta_2 = EOS - OT,
         G = if_else(delta_2 > 0,
                     delta_2,
                     0),
         L = if_else(delta_2 < 0,
                     -delta_2,
                     0),
         delta_log_1 = log2(1 + OT) - log2(1 + PRE),
         delta_log_2 = log2(1 + EOS) - log2(1 + OT),
         fc_1 = (1 + OT) / (1 + PRE),
         fc_2 = (1 + EOS) / (1 + OT),
         fc_3 = (1 + EOS) / (1 + PRE),
         relative_fc =  (2*fc_3) / (fc_1 + fc_2),
         mean_load = (PRE + OT + EOS)/3,
         mutation_load_score = log2(1 + mean_load * relative_fc),
         scale_factor = if_else(delta_log_1 < 0,
                                (1 + L) / (1 + G),
                                (1 + G) / (1 + L)),
         my_score = if_else(delta_log_1 != 0,
                            scale_factor * delta_log_1,
                            delta_log_1 + delta_log_2))

load_change_tb %>% 
  ggplot(aes(tumor_class, mutation_load_score)) + 
  geom_boxplot(fill = 'purple') + 
  coord_flip() + 
  xlab('PFS quartiles') + 
  ylab('mutation load score')

load_change_tb %>% 
  ggplot(aes(tumor_class, mutation_load_score, fill = Arm)) + 
  geom_boxplot() + 
  coord_flip() + 
  xlab('PFS quartiles') + 
  ylab('mutation load score') + 
  scale_fill_manual(values = c('blue', 'red'))

stats_load_change_tb <- load_change_tb %>% 
  group_by(tumor_class) %>% 
  summarize(mean_1 = mean(delta_log_1),
            median_1 = median(delta_log_1),
            sd_1 = sd(delta_log_1),
            mad_1 = mad(delta_log_1),
            mean_2 = mean(delta_log_2),
            median_2 = median(delta_log_2),
            sd_2 = sd(delta_log_2),
            mad_2 = mad(delta_log_2))

###################################################################################################################################################################################################
###################################################################################################################################################################################################
#### gene mutations

#### time-arm-quartile
timepoint_arm_class_gene_tb <- cox_data_tb %>% 
  count(Gene, timepoint, Arm, tumor_class, customer) %>% 
  mutate(patient_occurrence = 1, 
         Arm = parse_factor(Arm, levels = c("A", "B")),
         tumor_class = parse_factor(tumor_class, levels = str_c("Q", 1:4))) %>% 
  group_by(Gene, timepoint, Arm, tumor_class) %>% 
  summarize(N_total = sum(n),
            N_patient = sum(patient_occurrence)) %>% 
  ungroup() %>% 
  unite(col = "Time-Arm-Quartile", timepoint, Arm, tumor_class, sep = "-", remove = FALSE)

timepoint_arm_class_gene_variant_tb <- cox_data_tb %>% 
  count(Gene, Variant_type, timepoint, Arm, tumor_class, customer) %>% 
  unite(col = 'gene_variant', Gene, Variant_type, sep = '-', remove = TRUE) %>% 
  mutate(patient_occurrence = 1, 
         Arm = parse_factor(Arm, levels = c("A", "B")),
         tumor_class = parse_factor(tumor_class, levels = str_c("Q", 1:4))) %>% 
  group_by(gene_variant, timepoint, Arm, tumor_class) %>% 
  summarize(N_total = sum(n),
            N_patient = sum(patient_occurrence)) %>% 
  ungroup() %>% 
  unite(col = "Time-Arm-Quartile", timepoint, Arm, tumor_class, sep = "-", remove = FALSE)

timepoint_arm_class_gene_allele_tb <- cox_data_tb %>% 
  count(gene_mutation, timepoint, Arm, tumor_class, customer) %>% 
  mutate(patient_occurrence = 1, 
         Arm = parse_factor(Arm, levels = c("A", "B")),
         tumor_class = parse_factor(tumor_class, levels = str_c("Q", 1:4))) %>% 
  group_by(gene_mutation, timepoint, Arm, tumor_class) %>% 
  summarize(N_total = sum(n),
            N_patient = sum(patient_occurrence)) %>% 
  ungroup() %>% 
  unite(col = "Time-Arm-Quartile", timepoint, Arm, tumor_class, sep = "-", remove = FALSE)

#### time-arm
timepoint_arm_gene_tb <- timepoint_arm_class_gene_tb %>% 
  group_by(Gene, timepoint, Arm) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Time-Arm', timepoint, Arm, sep = '-', remove = FALSE)

timepoint_arm_gene_variant_tb <- timepoint_arm_class_gene_variant_tb %>% 
  group_by(gene_variant, timepoint, Arm) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Time-Arm', timepoint, Arm, sep = '-', remove = FALSE)

timepoint_arm_gene_allele_tb <- timepoint_arm_class_gene_allele_tb %>% 
  group_by(gene_mutation, timepoint, Arm) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Time-Arm', timepoint, Arm, sep = '-', remove = FALSE)

#### time-quartile
timepoint_class_gene_tb <- timepoint_arm_class_gene_tb %>% 
  group_by(Gene, timepoint, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Time-Quartile', timepoint, tumor_class, sep = '-', remove = FALSE)

timepoint_class_gene_variant_tb <- timepoint_arm_class_gene_variant_tb %>% 
  group_by(gene_variant, timepoint, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Time-Quartile', timepoint, tumor_class, sep = '-', remove = FALSE)

timepoint_class_gene_allele_tb <- timepoint_arm_class_gene_allele_tb %>% 
  group_by(gene_mutation, timepoint, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Time-Quartile', timepoint, tumor_class, sep = '-', remove = FALSE)

#### arm-quartile
arm_class_gene_tb <- timepoint_arm_class_gene_tb %>% 
  group_by(Gene, Arm, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Arm-Quartile', Arm, tumor_class, sep = '-', remove = FALSE)

arm_class_gene_variant_tb <- timepoint_arm_class_gene_variant_tb %>% 
  group_by(gene_variant, Arm, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Arm-Quartile', Arm, tumor_class, sep = '-', remove = FALSE)

arm_class_gene_allele_tb <- timepoint_arm_class_gene_allele_tb %>% 
  group_by(gene_mutation, Arm, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup() %>% 
  unite(col = 'Arm-Quartile', Arm, tumor_class, sep = '-', remove = FALSE)

#### time
timepoint_gene_tb <- timepoint_arm_class_gene_tb %>% 
  group_by(Gene, timepoint) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

timepoint_gene_variant_tb <- timepoint_arm_class_gene_variant_tb %>% 
  group_by(gene_variant, timepoint) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

timepoint_gene_allele_tb <- timepoint_arm_class_gene_allele_tb %>% 
  group_by(gene_mutation, timepoint) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

#### quartile
class_gene_tb <- timepoint_arm_class_gene_tb %>% 
  group_by(Gene, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

class_gene_variant_tb <- timepoint_arm_class_gene_variant_tb %>% 
  group_by(gene_variant, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

class_gene_allele_tb <- timepoint_arm_class_gene_allele_tb %>% 
  group_by(gene_mutation, tumor_class) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

#### arm
arm_gene_tb <- timepoint_arm_class_gene_tb %>% 
  group_by(Gene, Arm) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

arm_gene_variant_tb <- timepoint_arm_class_gene_variant_tb %>% 
  group_by(gene_variant, Arm) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

arm_gene_allele_tb <- timepoint_arm_class_gene_allele_tb %>% 
  group_by(gene_mutation, Arm) %>% 
  summarize(N_total = sum(N_total),
            N_patient = sum(N_patient)) %>% 
  ungroup()

####################################################################################
####################################################################################

#### scaled log-fold change per mutation

macro_mutation_tb <- cox_data_tb %>% 
  unite(col = 'gene_variant', Gene, Variant_type, sep = '-', remove = FALSE) %>% 
  select(customer, timepoint, tumor_class, PFS_time, PFS_status, Arm, Gene, gene_variant, gene_mutation) %>% 
  mutate(c = 1)

#### gene level mutations only
gene_mutation_feature_tb <- macro_mutation_tb %>% 
  select(-c(gene_variant, gene_mutation)) %>% 
  group_by(customer, PFS_status, PFS_time, tumor_class, Arm, Gene, timepoint) %>% 
  summarize(individual_mutation_load = sum(c)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = timepoint, values_from = individual_mutation_load, values_fill = 0) %>% 
  group_by(customer, PFS_status, PFS_time, tumor_class, Arm, Gene) %>% 
  mutate(fc_1 = (1 + OT) / (1 + PRE),
         fc_2 = (1 + EOS) / (1 + OT),
         fc_3 = (1 + EOS) / (1 + PRE),
         relative_fc =  (2*fc_3) / (fc_1 + fc_2),
         mean_load = (PRE + OT + EOS)/3,
         gene_score = log2(1 + mean_load * relative_fc)) %>% 
  ungroup() %>% 
  select(customer, PFS_time, PFS_status, tumor_class, Arm, Gene, gene_score)

pre_gene_mutation_feature_mat <- gene_mutation_feature_tb %>% 
  group_by(customer) %>% 
  pivot_wider(names_from = Gene, values_from = gene_score, values_fill = 0) %>% 
  ungroup()

gene_mutation_feature_mat <- pre_gene_mutation_feature_mat %>% 
  select(-c(customer, PFS_time, PFS_status, tumor_class, Arm)) %>% 
  as.matrix()
rownames(gene_mutation_feature_mat) <- pre_gene_mutation_feature_mat %>% pull(customer)

#### STARTING WITH THIS BOX PLOT IDEA 
#### LET US MAKE A FEW PLOTS THAT ILLUSTRATE THE UTILITY
#### OF THE GENE SCORE AND THEN USE GENE SCORE IN MODELS
p <- gene_mutation_feature_tb %>% 
  ggplot(aes(Gene, gene_score)) + 
  geom_boxplot(fill = 'purple') + 
  coord_flip()
ggsave(here('ms_plots', 'coarse_gene_mutation_score_box_plots.png'),
       p)

p <- gene_mutation_feature_tb %>% 
  ggplot(aes(Gene, gene_score, fill = Arm)) + 
  geom_boxplot() + 
  coord_flip() + 
  scale_fill_manual(values = c('blue', 'red')) + 
  facet_wrap( ~ Arm)
ggsave(here('ms_plots', 'by_arm_gene_mutation_score_box_plots.png'),
       p)

#### gene-variant level mutations only
gene_variant_mutation_feature_tb <- macro_mutation_tb %>% 
  select(-c(Gene, gene_mutation)) %>% 
  group_by(customer, PFS_status, PFS_time, tumor_class, Arm, gene_variant, timepoint) %>% 
  summarize(individual_mutation_load = sum(c)) %>% 
  ungroup() %>% 
  mutate(gene_variant = if_else(gene_variant == 'NONE-NONE',
                                'NONE',
                                gene_variant)) %>% 
  pivot_wider(names_from = timepoint, values_from = individual_mutation_load, values_fill = 0) %>% 
  group_by(customer, PFS_status, PFS_time, tumor_class, Arm, gene_variant) %>% 
  mutate(fc_1 = (1 + OT) / (1 + PRE),
         fc_2 = (1 + EOS) / (1 + OT),
         fc_3 = (1 + EOS) / (1 + PRE),
         relative_fc =  (2*fc_3) / (fc_1 + fc_2),
         mean_load = (PRE + OT + EOS)/3,
         gene_score = log2(1 + mean_load * relative_fc)) %>% 
  ungroup() %>% 
  select(customer, PFS_time, PFS_status, tumor_class, Arm, gene_variant, gene_score)

pre_gene_variant_mutation_feature_mat <- gene_variant_mutation_feature_tb %>% 
  group_by(customer) %>% 
  pivot_wider(names_from = gene_variant, values_from = gene_score, values_fill = 0) %>% 
  ungroup()

gene_variant_mutation_feature_mat <- pre_gene_variant_mutation_feature_mat %>% 
  select(-c(customer, PFS_time, PFS_status, tumor_class, Arm)) %>% 
  as.matrix()
rownames(gene_variant_mutation_feature_mat) <- pre_gene_variant_mutation_feature_mat %>% pull(customer)

gene_variant_mutation_feature_tb %>% 
  ggplot(aes(gene_variant, gene_score)) + 
  geom_boxplot() + 
  coord_flip()



#### gene-allele level mutations only
gene_allele_mutation_feature_tb <- macro_mutation_tb %>% 
  select(-c(Gene, gene_variant)) %>% 
  group_by(customer, PFS_status, PFS_time, tumor_class, Arm, gene_mutation, timepoint) %>% 
  summarize(individual_mutation_load = sum(c)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = timepoint, values_from = individual_mutation_load, values_fill = 0) %>% 
  group_by(customer, PFS_status, PFS_time, tumor_class, Arm, gene_mutation) %>% 
  mutate(fc_1 = (1 + OT) / (1 + PRE),
         fc_2 = (1 + EOS) / (1 + OT),
         fc_3 = (1 + EOS) / (1 + PRE),
         relative_fc =  (2*fc_3) / (fc_1 + fc_2),
         mean_load = (PRE + OT + EOS)/3,
         gene_score = log2(1 + mean_load * relative_fc)) %>% 
  ungroup() %>% 
  select(customer, PFS_time, PFS_status, tumor_class, Arm, gene_mutation, gene_score)

pre_gene_allele_mutation_feature_mat <- gene_allele_mutation_feature_tb %>% 
  group_by(customer) %>% 
  pivot_wider(names_from = gene_mutation, values_from = gene_score, values_fill = 0) %>% 
  ungroup()

gene_allele_mutation_feature_mat <- pre_gene_allele_mutation_feature_mat %>% 
  select(-c(customer, PFS_time, PFS_status, tumor_class, Arm)) %>% 
  as.matrix()
rownames(gene_allele_mutation_feature_mat) <- pre_gene_allele_mutation_feature_mat %>% pull(customer)

#### NEED TO USE THESE FEATURE MATRICES IN VARIOUS MODELS
#### I THINK WE SHOULD TRY:
#### (1) NAIVE LOGISTIC REGRESSION ON Q1 VS Q4 WITH ELASTIC-NET
#### (2) MULTINOMIAL REGRESSION ON Q1-4 WITH ELASTIC-NET
#### (3) DECISION TREE ON Q1-4
#### (4) RANDOM FOREST ON Q1-4
#### (5) COX PROPORTIONAL HAZZARDS ON PFS TIME WITH ELASTIC-NET

macro_response_tb <- cox_data_tb %>% 
  select(customer, PFS_time, PFS_status, tumor_class, Arm) %>% 
  unique() %>% 
  arrange(customer) %>% 
  mutate(log_PFS_time = log2(.5 + PFS_time))

#### glm.net tutorial recommends standardizing response
#### using 1/N variance, which is not default var() function
raw_log_PFS <- macro_response_tb %>% pull(log_PFS_time)
mean_lPFS <- mean(raw_log_PFS)
sd_lPFS <- sqrt( sum((raw_log_PFS - mean_lPFS)^2) / length(raw_log_PFS) )

macro_response_tb <- macro_response_tb %>% 
  group_by(customer) %>% 
  mutate(log_PFS_time = (log_PFS_time - mean_lPFS)/sd_lPFS) %>% 
  ungroup()

arm_A_response_tb <- macro_response_tb %>% 
  filter(Arm == 'A')

arm_B_response_tb <- macro_response_tb %>% 
  filter(Arm == 'B')

round(nrow(arm_A_response_tb)*.25)#[1] 23
round(nrow(arm_B_response_tb)*.25)#[1] 25

set.seed(23)
train_A_customers <- arm_A_response_tb %>% 
  pull(customer)
test_A_customers <- sort(sample(train_A_customers, 23))
train_A_customers <- sort(setdiff(train_A_customers, test_A_customers))

binary_train_A_customers <- arm_A_response_tb %>% 
  filter(customer %in% train_A_customers, tumor_class %in% c('Q1', 'Q4')) %>% 
  pull(customer)

binary_test_A_customers <- arm_A_response_tb %>% 
  filter(customer %in% test_A_customers, tumor_class %in% c('Q1', 'Q4')) %>% 
  pull(customer)

train_B_customers <- arm_B_response_tb %>% 
  pull(customer)
test_B_customers <- sort(sample(train_B_customers, 25))
train_B_customers <- sort(setdiff(train_B_customers, test_B_customers))

binary_train_B_customers <- arm_B_response_tb %>% 
  filter(customer %in% train_B_customers, tumor_class %in% c('Q1', 'Q4')) %>% 
  pull(customer)

binary_test_B_customers <- arm_B_response_tb %>% 
  filter(customer %in% test_B_customers, tumor_class %in% c('Q1', 'Q4')) %>% 
  pull(customer)

#### train data
train_macro_response_tb <- macro_response_tb %>% 
  filter(customer %in% sort(c(train_A_customers, train_B_customers)))

train_arm_A_response_tb <- arm_A_response_tb %>% 
  filter(customer %in% train_A_customers)

train_arm_B_response_tb <- arm_B_response_tb %>% 
  filter(customer %in% train_B_customers)

binary_train_macro_response_tb <- macro_response_tb %>% 
  filter(customer %in% sort(c(binary_train_A_customers, binary_train_B_customers)))

binary_train_arm_A_response_tb <- arm_A_response_tb %>% 
  filter(customer %in% binary_train_A_customers)

binary_train_arm_B_response_tb <- arm_B_response_tb %>% 
  filter(customer %in% binary_train_B_customers)

feature_mat_train_idcs <- which(rownames(gene_mutation_feature_mat) %in% sort(c(train_A_customers, train_B_customers)))
arm_A_feature_mat_train_idcs <- which(rownames(gene_mutation_feature_mat) %in% train_A_customers)
arm_B_feature_mat_train_idcs <- which(rownames(gene_mutation_feature_mat) %in% train_B_customers)

binary_feature_mat_train_idcs <- which(rownames(gene_mutation_feature_mat) %in% sort(c(binary_train_A_customers, binary_train_B_customers)))
binary_arm_A_feature_mat_train_idcs <- which(rownames(gene_mutation_feature_mat) %in% binary_train_A_customers)
binary_arm_B_feature_mat_train_idcs <- which(rownames(gene_mutation_feature_mat) %in% binary_train_B_customers)


#### test data
test_macro_response_tb <- macro_response_tb %>% 
  filter(customer %in% sort(c(test_A_customers, test_B_customers)))

test_arm_A_response_tb <- arm_A_response_tb %>% 
  filter(customer %in% test_A_customers)

test_arm_B_response_tb <- arm_B_response_tb %>% 
  filter(customer %in% test_B_customers)

binary_test_macro_response_tb <- macro_response_tb %>% 
  filter(customer %in% sort(c(binary_test_A_customers, binary_test_B_customers)))

binary_test_arm_A_response_tb <- arm_A_response_tb %>% 
  filter(customer %in% binary_test_A_customers)

binary_test_arm_B_response_tb <- arm_B_response_tb %>% 
  filter(customer %in% binary_test_B_customers)

feature_mat_test_idcs <- which(rownames(gene_mutation_feature_mat) %in% sort(c(test_A_customers, test_B_customers)))
arm_A_feature_mat_test_idcs <- which(rownames(gene_mutation_feature_mat) %in% test_A_customers)
arm_B_feature_mat_test_idcs <- which(rownames(gene_mutation_feature_mat) %in% test_B_customers)

binary_feature_mat_test_idcs <- which(rownames(gene_mutation_feature_mat) %in% sort(c(binary_test_A_customers, binary_test_B_customers)))
binary_arm_A_feature_mat_test_idcs <- which(rownames(gene_mutation_feature_mat) %in% binary_test_A_customers)
binary_arm_B_feature_mat_test_idcs <- which(rownames(gene_mutation_feature_mat) %in% binary_test_B_customers)


#### linear model
#### here we use log-PFS as response

#### first try gene_mutation_feature_mat
#### then try gene_variant_mutation_feature_mat
#### finally try gene_allele_mutation_feature_mat

dummy_response <- train_macro_response_tb %>% 
  pull(log_PFS_time) %>% 
  as.matrix()
rownames(dummy_response) <- train_macro_response_tb %>% 
  pull(customer)

feature_mat <- gene_allele_mutation_feature_mat
dummy_X <- feature_mat[feature_mat_train_idcs, ]
mean(rownames(dummy_response) == rownames(dummy_X))

test_X <- feature_mat[feature_mat_test_idcs, ]
test_y <- test_macro_response_tb %>% 
  pull(log_PFS_time)

alpha.seq <- seq(0, 1, .1)
mse_vec <- rep(NA, length(alpha.seq))
for (i in seq_along(alpha.seq)) {
  train_out <- cv.glmnet(x = dummy_X, 
                         y = dummy_response, 
                         alpha = alpha.seq[i],
                         family = 'gaussian', 
                         nfolds = 3, 
                         intercept = T)
  print(train_out)
  test_yhat <- predict(train_out, newx = test_X, s = 'lambda.min')
  
  mse_vec[i] <- sqrt(sum((test_y - test_yhat)^2))
}

tibble(alpha = alpha.seq,
       mse = mse_vec) %>% 
  ggplot(aes(alpha, mse)) + 
  geom_line() + 
  ylim(0, 8)

# train_gaussian_lasso_out <- cv.glmnet(x = dummy_X, 
#                                       y = dummy_response, 
#                                       alpha = 1,
#                                       family = 'gaussian', 
#                                       nfolds = 3, 
#                                       intercept = T)
# print(train_gaussian_lasso_out)
# coef(train_gaussian_lasso_out, s = 'lambda.min')
# 
# gaussian_lasso_yhat <- predict(train_gaussian_lasso_out, newx = gene_allele_mutation_feature_mat[feature_mat_test_idcs, ], s = 'lambda.min')
# dummy_y_test <- test_macro_response_tb %>% 
#   pull(log_PFS_time)
# sqrt(sum((dummy_y_test - gaussian_lasso_yhat)^2))##[1] 7.055718 gene ##[1] 7.195984 gene-variant ##[1] 7.195984 gene-allele
# 
# train_gaussian_ridge_out <- cv.glmnet(x = dummy_X, 
#                                       y = dummy_response, 
#                                       alpha = 0,
#                                       family = 'gaussian', 
#                                       nfolds = 3, 
#                                       intercept = T)
# print(train_gaussian_ridge_out)
# coef(train_gaussian_ridge_out, s = 'lambda.min')
# 
# gaussian_ridge_yhat <- predict(train_gaussian_ridge_out, newx = gene_allele_mutation_feature_mat[feature_mat_test_idcs, ], s = 'lambda.min')
# dummy_y_test <- test_macro_response_tb %>% 
#   pull(log_PFS_time)
# sqrt(sum((dummy_y_test - gaussian_ridge_yhat)^2))##[1] 7.072706 gene ##[1] 7.203398 gene-variant ##[1] 7.177026 gene-allele
# 
# train_gaussian_alpha_.5_out <- cv.glmnet(x = dummy_X, 
#                                          y = dummy_response, 
#                                          alpha = .5,
#                                          family = 'gaussian', 
#                                          nfolds = 3, 
#                                          intercept = T)
# print(train_gaussian_alpha_.5_out)
# coef(train_gaussian_alpha_.5_out, s = 'lambda.min')
# 
# gaussian_alpha_.5_yhat <- predict(train_gaussian_alpha_.5_out, newx = gene_allele_mutation_feature_mat[feature_mat_test_idcs, ], s = 'lambda.min')
# dummy_y_test <- test_macro_response_tb %>% 
#   pull(log_PFS_time)
# sqrt(sum((dummy_y_test - gaussian_alpha_.5_yhat)^2))##[1] 7.112753 gene ##[1] 7.195984 gene-variant ##[1] 7.225609 gene-allele


#### multinomial regression

dummy_response <- train_macro_response_tb %>% 
  mutate(tumor_class = parse_factor(tumor_class,
                                    levels = str_c('Q', 1:4))) %>% 
  pull(tumor_class)


feature_mat <- gene_mutation_feature_mat
dummy_X <- feature_mat[feature_mat_train_idcs, ]

test_X <- feature_mat[feature_mat_test_idcs, ]
test_y <- test_macro_response_tb %>% 
  pull(tumor_class)

alpha.seq <- seq(0, 1, .1)
mean_class_error <- rep(NA, length(alpha.seq))
for (i in seq_along(alpha.seq)) {
  train_out <- cv.glmnet(x = dummy_X, 
                         y = dummy_response, 
                         alpha = alpha.seq[i],
                         family = 'multinomial', 
                         nfolds = 3, 
                         intercept = T)
  print(train_out)
  test_yhat <- predict(train_out, newx = test_X, s = 'lambda.min', type='class')
  
  mean_class_error[i] <- mean(test_y == as.vector(test_yhat))
  
}

tibble(alpha = alpha.seq,
       mce = mean_class_error) %>% 
  ggplot(aes(alpha, mce)) + 
  geom_line() + 
  ylim(0, .5)



#### logistic regression

dummy_response <- binary_train_macro_response_tb %>% 
  mutate(tumor_class = parse_factor(tumor_class,
                                    levels = str_c('Q', c(1,4)))) %>% 
  pull(tumor_class)

feature_mat <- gene_mutation_feature_mat
dummy_X <- feature_mat[binary_feature_mat_train_idcs, ]

test_X <- feature_mat[binary_feature_mat_test_idcs, ]
test_y <- binary_test_macro_response_tb %>% 
  pull(tumor_class)

alpha.seq <- seq(0, 1, .1)
mean_class_error <- rep(NA, length(alpha.seq))
for (i in seq_along(alpha.seq)) {
  train_out <- cv.glmnet(x = dummy_X, 
                         y = dummy_response, 
                         alpha = alpha.seq[i],
                         family = 'binomial', 
                         nfolds = 3, 
                         intercept = T)
  test_yhat <- predict(train_out, newx = test_X, s = 'lambda.min', type='class')
  
  mean_class_error[i] <- mean(test_y == as.vector(test_yhat))
  
}

tibble(alpha = alpha.seq,
       mce = mean_class_error) %>% 
  ggplot(aes(alpha, mce)) + 
  geom_line() + 
  ylim(0, .75)





#### cox proportional hazzards

train_cox_model_y <- Surv(time = train_macro_response_tb$PFS_time, 
                          event = if_else(train_macro_response_tb$PFS_status == 'Censor',
                                          0, 
                                          1)) 

feature_mat <- gene_variant_mutation_feature_mat
dummy_X <- feature_mat[feature_mat_train_idcs, ]

test_X <- feature_mat[feature_mat_test_idcs, ]
test_cox_model_y <- Surv(time = test_macro_response_tb$PFS_time, 
                         event = if_else(test_macro_response_tb$PFS_status == 'Censor',
                                         0, 
                                         1))
alpha.seq <- seq(0, 1, .1)
mean_class_error <- rep(NA, length(alpha.seq))
for (i in seq_along(alpha.seq)) {
  # i <- 5
  train_out <- cv.glmnet(x = dummy_X, 
                         y = train_cox_model_y, 
                         alpha = alpha.seq[i],
                         family = 'cox', 
                         nfolds = 3,
                         type.measure = 'C')
  test_yhat <- predict(train_out, newx = test_X, s = 'lambda.min')
  
  survfit_obj <- survival::survfit(train_out, 
                                   s = "lambda.min", 
                                   x = dummy_X, 
                                   y = train_cox_model_y, newx = test_X)
  plot(survfit_obj)
  
  # mean_class_error[i] <- mean(test_y == as.vector(test_yhat))
  
}

dummy_coef <- coef(train_out, s = 'lambda.min')
rownames(dummy_coef)[which(!near(dummy_coef, 0))]
dummy_coef[which(!near(dummy_coef, 0))]

# tibble(alpha = alpha.seq,
#        mce = mean_class_error) %>% 
#   ggplot(aes(alpha, mce)) + 
#   geom_line() + 
#   ylim(0, .75)

#### decision tree

#### random forest
library(randomForest)
dummy_response <- train_macro_response_tb %>% 
  mutate(tumor_class = parse_factor(tumor_class,
                                    levels = str_c('Q', 1:4))) %>% 
  pull(tumor_class)

feature_mat <- gene_allele_mutation_feature_mat
dummy_X <- feature_mat[feature_mat_train_idcs, ]

test_X <- feature_mat[feature_mat_test_idcs, ]
test_y <- test_macro_response_tb %>% 
  mutate(tumor_class = parse_factor(tumor_class,
                                    levels = str_c('Q', 1:4))) %>% 
  pull(tumor_class)

rF_out <- randomForest(x = dummy_X, y = dummy_response, xtest = test_X, ytest = test_y, ntree = 2000)

#### gene
#   randomForest(x = dummy_X, y = dummy_response, xtest = test_X,      ytest = test_y, ntree = 2000) 
# Type of random forest: classification
# Number of trees: 2000
# No. of variables tried at each split: 8
# 
# OOB estimate of  error rate: 70.42%
# Confusion matrix:
#   Q1 Q2 Q3 Q4 class.error
# Q1 12  9  3 10   0.6470588
# Q2  7 13 10  8   0.6578947
# Q3  6  9  9 14   0.7631579
# Q4  5  7 12  8   0.7500000
# Test set error rate: 54.17%
# Confusion matrix:
#   Q1 Q2 Q3 Q4 class.error
# Q1  8  0  1  4   0.3846154
# Q2  2  4  1  3   0.6000000
# Q3  0  2  4  4   0.6000000
# Q4  1  2  6  6   0.6000000

# gene-variant
#   randomForest(x = dummy_X, y = dummy_response, xtest = test_X,      ytest = test_y, ntree = 2000) 
# Type of random forest: classification
# Number of trees: 2000
# No. of variables tried at each split: 10
# 
# OOB estimate of  error rate: 69.72%
# Confusion matrix:
#   Q1 Q2 Q3 Q4 class.error
# Q1 13  8  7  6   0.6176471
# Q2  8 11  9 10   0.7105263
# Q3  9 10  6 13   0.8421053
# Q4  5  5  9 13   0.5937500
# Test set error rate: 45.83%
# Confusion matrix:
#   Q1 Q2 Q3 Q4 class.error
# Q1  9  1  2  1   0.3076923
# Q2  2  6  0  2   0.4000000
# Q3  2  1  4  3   0.6000000
# Q4  2  2  4  7   0.5333333

# gene-allele
#   randomForest(x = dummy_X, y = dummy_response, xtest = test_X,      ytest = test_y, ntree = 2000) 
# Type of random forest: classification
# Number of trees: 2000
# No. of variables tried at each split: 25
# 
# OOB estimate of  error rate: 77.46%
# Confusion matrix:
#   Q1 Q2 Q3 Q4 class.error
# Q1  6  9  9 10   0.8235294
# Q2  3  8 15 12   0.7894737
# Q3  4 13 10 11   0.7368421
# Q4  3 11 10  8   0.7500000
# Test set error rate: 64.58%
# Confusion matrix:
#   Q1 Q2 Q3 Q4 class.error
# Q1  4  3  5  1   0.6923077
# Q2  1  5  2  2   0.5000000
# Q3  1  1  4  4   0.6000000
# Q4  3  3  5  4   0.7333333

###################################################################################################################################################################################################
###################################################################################################################################################################################################
#### ms p-values

#### QUESTION: Are there statistically significant differences in mutation load per patient?
chisq_mutation_load_tb <- mutation_load_tb %>% 
  select(customer, timepoint, tumor_class, Arm, count_bin) %>% 
  unique() %>% 
  mutate(p = 1,
         count_bin = str_c('bin_', count_bin)) %>% 
  group_by(timepoint, tumor_class, Arm, count_bin) %>% 
  summarize(n = sum(p)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = count_bin, values_from = n, values_fill = 0) %>% 
  unite(col = group, timepoint:Arm, sep = '-')

coarse_chisq_mutation_load_tb <- mutation_load_tb %>% 
  select(customer, timepoint, tumor_class, count_bin) %>% 
  unique() %>% 
  mutate(p = 1,
         count_bin = str_c('bin_', count_bin)) %>% 
  group_by(timepoint, tumor_class, count_bin) %>% 
  summarize(n = sum(p)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = count_bin, values_from = n, values_fill = 0) %>% 
  unite(col = group, timepoint:tumor_class, sep = '-')

mat_chisq_mutation_load <- chisq_mutation_load_tb %>% 
  select(-group) %>% 
  as.matrix()
mat_chisq_mutation_load <- mat_chisq_mutation_load[ , c(7, 1:6)]
rownames(mat_chisq_mutation_load) <- chisq_mutation_load_tb %>% 
  pull(group)

mat_coarse_chisq_mutation_load <- coarse_chisq_mutation_load_tb %>% 
  select(-group) %>% 
  as.matrix()
mat_coarse_chisq_mutation_load <- mat_coarse_chisq_mutation_load[ , c(7, 1:6)]
rownames(mat_coarse_chisq_mutation_load) <- coarse_chisq_mutation_load_tb %>% 
  pull(group)

coarse_chisq_pvalue_tb <- my_homogeneity_test(mat_coarse_chisq_mutation_load) %>% 
  separate(group_1, into = c('time_1', 'class_1')) %>% 
  separate(group_2, into = c('time_2', 'class_2'))

coarse_chisq_pvalue_tb %>% 
  filter(time_1 == time_2)

# QUESTION: What is the difference between patients in Q2 and Q4 at the OT timepoint?

coarse_chisq_pvalue_tb %>% 
  filter(class_1 == class_2)

arm_chisq_pvalue_tb <- my_homogeneity_test(mat_chisq_mutation_load) %>% 
  separate(group_1, into = c('time_1', 'class_1', 'Arm_1')) %>% 
  separate(group_2, into = c('time_2', 'class_2', 'Arm_2'))

arm_chisq_pvalue_tb %>% 
  filter(time_1 == time_2, class_1 == class_2, Arm_1 != Arm_2)

arm_chisq_pvalue_tb %>% 
  filter(time_1 != time_2, class_1 == class_2, Arm_1 == Arm_2)

###############################################################################
#### QUESTION: Are there statistically significant different mutations per patient?

#### N.b.: NONE is included below, which I believe makes sense bc we would like to include a
#### "no mutation" event in the comparisons below.

#######################################
#### Q: between the two treatment arms?

#### genes
pre_arm_gene_mat <- arm_gene_tb %>% 
  # select(-N_total) %>%
  # pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)
  select(-N_patient) %>%
  pivot_wider(names_from = Gene, values_from = N_total, values_fill = 0)

arm_gene_mat <- pre_arm_gene_mat %>%  
  select(-Arm) %>% 
  as.matrix()
rownames(arm_gene_mat) <- pre_arm_gene_mat %>% pull(Arm)

#### using N_patient
my_homogeneity_test(arm_gene_mat)$p_value#[1] 0.003483149

#### using N_total
my_homogeneity_test(arm_gene_mat)$p_value#[1] 0.0002885396

#### gene-variant
pre_arm_gene_variant_mat <- arm_gene_variant_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = gene_variant, values_from = N_total, values_fill = 0)

arm_gene_variant_mat <- pre_arm_gene_variant_mat %>%  
  select(-Arm) %>% 
  as.matrix()
rownames(arm_gene_variant_mat) <- pre_arm_gene_variant_mat %>% pull(Arm)

#### using N_patient
my_homogeneity_test(arm_gene_variant_mat)$p_value#[1] 0.0008191125

#### using N_total
my_homogeneity_test(arm_gene_variant_mat)$p_value#[1] 1.088052e-05

#### gene-allele
pre_arm_gene_allele_mat <- arm_gene_allele_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = gene_mutation, values_from = N_total, values_fill = 0)

arm_gene_allele_mat <- pre_arm_gene_allele_mat %>%  
  select(-Arm) %>% 
  as.matrix()
rownames(arm_gene_allele_mat) <- pre_arm_gene_allele_mat %>% pull(Arm)

#### using N_patient
my_homogeneity_test(arm_gene_allele_mat)$p_value#[1] 8.080541e-25

#### using N_total
my_homogeneity_test(arm_gene_allele_mat)$p_value#[1] 3.829418e-25

#######################################
#### Q: between the timepoints?

#### genes
pre_timepoint_gene_mat <- timepoint_gene_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = Gene, values_from = N_total, values_fill = 0)

timepoint_gene_mat <- pre_timepoint_gene_mat %>%  
  select(-timepoint) %>% 
  as.matrix()
rownames(timepoint_gene_mat) <- pre_timepoint_gene_mat %>% pull(timepoint)

#### using N_patient
my_homogeneity_test(timepoint_gene_mat)
# group_1 group_2  p_value
# <chr>   <chr>      <dbl>
# 1 PRE     OT      5.56e-20
# 2 PRE     EOS     1.46e- 3
# 3 OT      EOS     3.15e- 1

#### using N_total
my_homogeneity_test(timepoint_gene_mat)
# group_1 group_2  p_value
# <chr>   <chr>      <dbl>
# 1 PRE     OT      5.24e-21
# 2 PRE     EOS     1.77e- 3
# 3 OT      EOS     1.63e- 1

#### gene-variant
pre_timepoint_gene_variant_mat <- timepoint_gene_variant_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = gene_variant, values_from = N_total, values_fill = 0)

timepoint_gene_variant_mat <- pre_timepoint_gene_variant_mat %>%  
  select(-timepoint) %>% 
  as.matrix()
rownames(timepoint_gene_variant_mat) <- pre_timepoint_gene_variant_mat %>% pull(timepoint)

#### using N_patient
my_homogeneity_test(timepoint_gene_variant_mat)
# group_1 group_2  p_value
# <chr>   <chr>      <dbl>
# 1 PRE     OT      2.46e-23
# 2 PRE     EOS     2.38e- 4
# 3 OT      EOS     2.83e- 1

#### using N_total
my_homogeneity_test(timepoint_gene_variant_mat)
# group_1 group_2  p_value
# <chr>   <chr>      <dbl>
# 1 PRE     OT      9.63e-25
# 2 PRE     EOS     3.38e- 5
# 3 OT      EOS     1.86e- 1

#### gene-allele
pre_timepoint_gene_allele_mat <- timepoint_gene_allele_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = gene_mutation, values_from = N_total, values_fill = 0)

timepoint_gene_allele_mat <- pre_timepoint_gene_allele_mat %>%  
  select(-timepoint) %>% 
  as.matrix()
rownames(timepoint_gene_allele_mat) <- pre_timepoint_gene_allele_mat %>% pull(timepoint)

#### using N_patient
my_homogeneity_test(timepoint_gene_allele_mat)
# group_1 group_2 p_value
# <chr>   <chr>     <dbl>
# 1 PRE     OT       0.0159
# 2 PRE     EOS      1.00  
# 3 OT      EOS      1.00  

#### using N_total
my_homogeneity_test(timepoint_gene_allele_mat)
# group_1 group_2 p_value
# <chr>   <chr>     <dbl>
# 1 PRE     OT       0.0120
# 2 PRE     EOS      1.00  
# 3 OT      EOS      1.00  

#######################################
#### Q: between the PFS quartiles?

#### genes
pre_class_gene_mat <- class_gene_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = Gene, values_from = N_total, values_fill = 0)

class_gene_mat <- pre_class_gene_mat %>%  
  select(-tumor_class) %>% 
  as.matrix()
rownames(class_gene_mat) <- pre_class_gene_mat %>% pull(tumor_class)

#### using N_patient
my_homogeneity_test(class_gene_mat)
# group_1 group_2   p_value
# <chr>   <chr>       <dbl>
# 1 Q1      Q2      0.0000105
# 2 Q2      Q4      0.0000290
# 3 Q1      Q4      0.000521 
# 4 Q1      Q3      0.000842 
# 5 Q2      Q3      0.000951 
# 6 Q3      Q4      0.0179  

#### using N_total
my_homogeneity_test(class_gene_mat)
# group_1 group_2      p_value
# <chr>   <chr>          <dbl>
# 1 Q1      Q2      0.0000000753
# 2 Q1      Q4      0.000000587 
# 3 Q2      Q4      0.00000339  
# 4 Q1      Q3      0.00000837  
# 5 Q2      Q3      0.000392    
# 6 Q3      Q4      0.00992 

#### gene-variant
pre_class_gene_variant_mat <- class_gene_variant_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = gene_variant, values_from = N_total, values_fill = 0)

class_gene_variant_mat <- pre_class_gene_variant_mat %>%  
  select(-tumor_class) %>% 
  as.matrix()
rownames(class_gene_variant_mat) <- pre_class_gene_variant_mat %>% pull(tumor_class)

#### using N_patient
my_homogeneity_test(class_gene_variant_mat)
# group_1 group_2    p_value
# <chr>   <chr>        <dbl>
# 1 Q1      Q2      0.00000357
# 2 Q1      Q4      0.00000778
# 3 Q2      Q4      0.0000219 
# 4 Q1      Q3      0.000558  
# 5 Q2      Q3      0.00122   
# 6 Q3      Q4      0.0330

#### using N_total
my_homogeneity_test(class_gene_variant_mat)
# group_1 group_2       p_value
# <chr>   <chr>           <dbl>
# 1 Q1      Q4      0.00000000123
# 2 Q1      Q2      0.00000000530
# 3 Q1      Q3      0.000000873  
# 4 Q2      Q4      0.00000961   
# 5 Q2      Q3      0.000306     
# 6 Q3      Q4      0.0100  

#### gene-allele
pre_class_gene_allele_mat <- class_gene_allele_tb %>% 
  select(-N_total) %>% 
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)
# select(-N_patient) %>% 
# pivot_wider(names_from = gene_mutation, values_from = N_total, values_fill = 0)

class_gene_allele_mat <- pre_class_gene_allele_mat %>%  
  select(-tumor_class) %>% 
  as.matrix()
rownames(class_gene_allele_mat) <- pre_class_gene_allele_mat %>% pull(tumor_class)

#### using N_patient
my_homogeneity_test(class_gene_allele_mat)
# group_1 group_2  p_value
# <chr>   <chr>      <dbl>
# 1 Q1      Q2      4.10e-20
# 2 Q1      Q3      8.89e-18
# 3 Q1      Q4      1.42e-14
# 4 Q2      Q4      2.38e-11
# 5 Q2      Q3      5.82e-11
# 6 Q4      Q3      1.68e- 6

#### using N_total
my_homogeneity_test(class_gene_allele_mat)
# group_1 group_2  p_value
# <chr>   <chr>      <dbl>
# 1 Q1      Q2      1.88e-20
# 2 Q1      Q3      1.68e-18
# 3 Q1      Q4      2.34e-15
# 4 Q2      Q4      2.07e-11
# 5 Q2      Q3      3.68e-11
# 6 Q4      Q3      1.01e- 6

#######################################
#### Q: between all possible patient groups?

#### genes
pre_gene_mat <- timepoint_arm_class_gene_tb %>% 
  rename(group = `Time-Arm-Quartile`) %>% 
  select(-c(timepoint, Arm, tumor_class, N_total)) %>% 
  pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)
# select(-c(timepoint, Arm, tumor_class, N_patient)) %>% 
# pivot_wider(names_from = Gene, values_from = N_total, values_fill = 0)

gene_mat <- pre_gene_mat %>%  
  select(-group) %>%  
  as.matrix()
rownames(gene_mat) <- pre_gene_mat %>% pull(group)

gene_pvalue_tb <- my_homogeneity_test(gene_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'Arm', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'Arm', 'tumor_class'), 2, sep = '_'), sep = '-')

gene_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 == tumor_class_2, Arm_1 != Arm_2)

gene_pvalue_tb %>% 
  filter(timepoint_1 != timepoint_2, tumor_class_1 == tumor_class_2, Arm_1 == Arm_2)

gene_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 != tumor_class_2, Arm_1 == Arm_2)

#### gene-variant
pre_gene_variant_mat <- timepoint_arm_class_gene_variant_tb %>% 
  rename(group = `Time-Arm-Quartile`) %>%
  select(-c(timepoint, Arm, tumor_class, N_total)) %>% 
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)
# select(-c(timepoint, Arm, tumor_class, N_patient)) %>% 
# pivot_wider(names_from = gene_variant, values_from = N_total, values_fill = 0)

gene_variant_mat <- pre_gene_variant_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(gene_variant_mat) <- pre_gene_variant_mat %>% pull(group)

gene_variant_pvalue_tb <- my_homogeneity_test(gene_variant_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'Arm', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'Arm', 'tumor_class'), 2, sep = '_'), sep = '-')

gene_variant_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 == tumor_class_2, Arm_1 != Arm_2)

gene_variant_pvalue_tb %>% 
  filter(timepoint_1 != timepoint_2, tumor_class_1 == tumor_class_2, Arm_1 == Arm_2)

gene_variant_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 != tumor_class_2, Arm_1 == Arm_2)

#### gene-allele
pre_gene_allele_mat <- timepoint_arm_class_gene_allele_tb %>% 
  rename(group = `Time-Arm-Quartile`) %>% 
  select(-c(timepoint, Arm, tumor_class, N_total)) %>% 
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)
# select(-c(timepoint, Arm, tumor_class, N_patient)) %>% 
# pivot_wider(names_from = gene_mutation, values_from = N_total, values_fill = 0)

gene_allele_mat <- pre_gene_allele_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(gene_allele_mat) <- pre_gene_allele_mat %>% pull(group)

gene_allele_pvalue_tb <- my_homogeneity_test(gene_allele_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'Arm', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'Arm', 'tumor_class'), 2, sep = '_'), sep = '-')

gene_allele_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 == tumor_class_2, Arm_1 != Arm_2)

gene_allele_pvalue_tb %>% 
  filter(timepoint_1 != timepoint_2, tumor_class_1 == tumor_class_2, Arm_1 == Arm_2)

gene_allele_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 != tumor_class_2, Arm_1 == Arm_2)

#######################################
#### Q: between PFS quartile and timepoint patient groups?

#### genes
pre_time_class_gene_mat <- timepoint_class_gene_tb %>% 
  rename(group = `Time-Quartile`) %>% 
  select(-c(timepoint, tumor_class, N_total)) %>% 
  pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)

time_class_gene_mat <- pre_time_class_gene_mat %>%  
  select(-group) %>%  
  as.matrix()
rownames(time_class_gene_mat) <- pre_time_class_gene_mat %>% pull(group)

time_class_gene_pvalue_tb <- my_homogeneity_test(time_class_gene_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'tumor_class'), 2, sep = '_'), sep = '-')

time_class_gene_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 != tumor_class_2)

time_class_gene_pvalue_tb %>% 
  filter(timepoint_1 != timepoint_2, tumor_class_1 == tumor_class_2)

time_class_gene_delta_tb <- find.pairwise.group.mutation.delta(time_class_gene_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'tumor_class'), 1, sep = '_'), sep = '-', remove = FALSE) %>% 
  separate(group_2, into = str_c(c('timepoint', 'tumor_class'), 2, sep = '_'), sep = '-', remove = FALSE)

p <- time_class_gene_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(timepoint_1 == timepoint_2, timepoint_1 == 'OT', tumor_class_1 != tumor_class_2) %>%
  ggplot(aes(mutation, delta)) + 
  geom_col(fill = 'purple') + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')
ggsave(here('ms_plots', 'PFS_gene_mutation_delta_at_OT_timepoint.png'),
       p)

#### gene-variant
pre_time_class_gene_variant_mat <- timepoint_class_gene_variant_tb %>% 
  rename(group = `Time-Quartile`) %>% 
  select(-c(timepoint, tumor_class, N_total)) %>%  
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)

time_class_gene_variant_mat <- pre_time_class_gene_variant_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(time_class_gene_variant_mat) <- pre_time_class_gene_variant_mat %>% pull(group)

time_class_gene_variant_pvalue_tb <- my_homogeneity_test(time_class_gene_variant_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'tumor_class'), 2, sep = '_'), sep = '-')

time_class_gene_variant_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 != tumor_class_2)

time_class_gene_variant_pvalue_tb %>% 
  filter(timepoint_1 != timepoint_2, tumor_class_1 == tumor_class_2)

time_class_gene_variant_delta_tb <- find.pairwise.group.mutation.delta(time_class_gene_variant_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'tumor_class'), 1, sep = '_'), sep = '-', remove = FALSE) %>% 
  separate(group_2, into = str_c(c('timepoint', 'tumor_class'), 2, sep = '_'), sep = '-', remove = FALSE)

p <- time_class_gene_variant_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(timepoint_1 == timepoint_2, timepoint_1 == 'OT', tumor_class_1 != tumor_class_2) %>%
  ggplot(aes(mutation, delta)) + 
  geom_col(fill = 'purple') + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')
ggsave(here('ms_plots', 'PFS_gene_variant_mutation_delta_at_OT_timepoint.png'),
       p)

#### gene-allele
pre_time_class_gene_allele_mat <- timepoint_class_gene_allele_tb %>% 
  rename(group = `Time-Quartile`) %>% 
  select(-c(timepoint, tumor_class, N_total)) %>%  
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)

time_class_gene_allele_mat <- pre_time_class_gene_allele_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(time_class_gene_allele_mat) <- pre_time_class_gene_allele_mat %>% pull(group)

time_class_gene_allele_pvalue_tb <- my_homogeneity_test(time_class_gene_allele_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'tumor_class'), 2, sep = '_'), sep = '-')

time_class_gene_allele_pvalue_tb %>% 
  filter(timepoint_1 == timepoint_2, tumor_class_1 != tumor_class_2)

time_class_gene_allele_pvalue_tb %>% 
  filter(timepoint_1 != timepoint_2, tumor_class_1 == tumor_class_2)


time_class_gene_allele_delta_tb <- find.pairwise.group.mutation.delta(time_class_gene_allele_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'tumor_class'), 1, sep = '_'), sep = '-', remove = FALSE) %>% 
  separate(group_2, into = str_c(c('timepoint', 'tumor_class'), 2, sep = '_'), sep = '-', remove = FALSE)

time_class_gene_allele_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(timepoint_1 == timepoint_2, timepoint_1 == 'PRE', tumor_class_1 != tumor_class_2) %>%
  # mutate(comparison = parse_factor(comparison,
  #                                  levels = c('B-Q1_vs_A-Q1', 'A-Q2_vs_B-Q2', 'A-Q3_vs_B-Q3', 'B-Q4_vs_A-Q4'))) %>% 
  ggplot(aes(mutation, delta)) + 
  geom_col(fill = 'purple') + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')


p <- time_class_gene_allele_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(timepoint_1 == timepoint_2, timepoint_1 == 'OT', 
         tumor_class_1 != tumor_class_2, tumor_class_1 == 'Q1') %>%
  # mutate(comparison = parse_factor(comparison,
  #                                  levels = c('B-Q1_vs_A-Q1', 'A-Q2_vs_B-Q2', 'A-Q3_vs_B-Q3', 'B-Q4_vs_A-Q4'))) %>% 
  ggplot(aes(mutation, delta)) + 
  geom_col(fill = 'purple') + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 8)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')
ggsave(here('ms_plots', 'Q1_PFS_gene_allele_mutation_delta_at_OT_timepoint.png'),
       p)

p <- time_class_gene_allele_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(timepoint_1 == timepoint_2, timepoint_1 == 'OT', 
         tumor_class_1 != tumor_class_2, tumor_class_1 != 'Q1') %>%
  ggplot(aes(mutation, delta)) + 
  geom_col(fill = 'purple') + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 6)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')
ggsave(here('ms_plots', 'non-Q1_PFS_gene_allele_mutation_delta_at_OT_timepoint.png'),
       p)

#######################################
#### Q: between arm and PFS quartile patient groups?

#### genes
pre_arm_class_gene_mat <- arm_class_gene_tb %>% 
  rename(group = `Arm-Quartile`) %>% 
  select(-c(Arm, tumor_class, N_total)) %>% 
  pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)

arm_class_gene_mat <- pre_arm_class_gene_mat %>%  
  select(-group) %>%  
  as.matrix()
rownames(arm_class_gene_mat) <- pre_arm_class_gene_mat %>% pull(group)

arm_class_gene_pvalue_tb <- my_homogeneity_test(arm_class_gene_mat) %>% 
  separate(group_1, into = str_c(c('Arm', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('Arm', 'tumor_class'), 2, sep = '_'), sep = '-')

arm_class_gene_pvalue_tb %>% 
  filter(Arm_1 == Arm_2, tumor_class_1 != tumor_class_2)

arm_class_gene_pvalue_tb %>% 
  filter(Arm_1 != Arm_2, tumor_class_1 == tumor_class_2)

arm_class_gene_delta_tb <- find.pairwise.group.mutation.delta(arm_class_gene_mat) %>% 
  separate(group_1, into = str_c(c('Arm', 'tumor_class'), 1, sep = '_'), sep = '-', remove = FALSE) %>% 
  separate(group_2, into = str_c(c('Arm', 'tumor_class'), 2, sep = '_'), sep = '-', remove = FALSE)

arm_class_gene_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(Arm_1 != Arm_2, tumor_class_1 == tumor_class_2) %>%
  mutate(comparison = parse_factor(comparison,
                                   levels = c('B-Q1_vs_A-Q1', 'A-Q2_vs_B-Q2', 'A-Q3_vs_B-Q3', 'B-Q4_vs_A-Q4'))) %>% 
  ggplot(aes(mutation, delta, fill = tumor_class_1)) + 
  geom_col() + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')


#### gene-variant
pre_arm_class_gene_variant_mat <- arm_class_gene_variant_tb %>% 
  rename(group = `Arm-Quartile`) %>% 
  select(-c(Arm, tumor_class, N_total)) %>% 
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)

arm_class_gene_variant_mat <- pre_arm_class_gene_variant_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(arm_class_gene_variant_mat) <- pre_arm_class_gene_variant_mat %>% pull(group)

arm_class_gene_variant_pvalue_tb <- my_homogeneity_test(arm_class_gene_variant_mat) %>% 
  separate(group_1, into = str_c(c('Arm', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('Arm', 'tumor_class'), 2, sep = '_'), sep = '-')

arm_class_gene_variant_pvalue_tb %>% 
  filter(Arm_1 == Arm_2, tumor_class_1 != tumor_class_2)

arm_class_gene_variant_pvalue_tb %>% 
  filter(Arm_1 != Arm_2, tumor_class_1 == tumor_class_2)

arm_class_gene_variant_delta_tb <- find.pairwise.group.mutation.delta(arm_class_gene_variant_mat) %>% 
  separate(group_1, into = str_c(c('Arm', 'tumor_class'), 1, sep = '_'), sep = '-', remove = FALSE) %>% 
  separate(group_2, into = str_c(c('Arm', 'tumor_class'), 2, sep = '_'), sep = '-', remove = FALSE)

arm_class_gene_variant_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(Arm_1 != Arm_2, tumor_class_1 == tumor_class_2) %>%
  mutate(comparison = parse_factor(comparison,
                                   levels = c('B-Q1_vs_A-Q1', 'B-Q2_vs_A-Q2', 'A-Q3_vs_B-Q3', 'B-Q4_vs_A-Q4')),
         mutation = if_else(mutation == 'NONE-NONE',
                            'NONE',
                            mutation)) %>% 
  ggplot(aes(mutation, delta, fill = tumor_class_1)) + 
  geom_col() + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')


#### gene-allele
pre_arm_class_gene_allele_mat <- arm_class_gene_allele_tb %>% 
  rename(group = `Arm-Quartile`) %>% 
  select(-c(Arm, tumor_class, N_total)) %>% 
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)

arm_class_gene_allele_mat <- pre_arm_class_gene_allele_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(arm_class_gene_allele_mat) <- pre_arm_class_gene_allele_mat %>% pull(group)

arm_class_gene_allele_pvalue_tb <- my_homogeneity_test(arm_class_gene_allele_mat) %>% 
  separate(group_1, into = str_c(c('Arm', 'tumor_class'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('Arm', 'tumor_class'), 2, sep = '_'), sep = '-')

arm_class_gene_allele_pvalue_tb %>% 
  filter(Arm_1 == Arm_2, tumor_class_1 != tumor_class_2)

arm_class_gene_allele_pvalue_tb %>% 
  filter(Arm_1 != Arm_2, tumor_class_1 == tumor_class_2)

arm_class_gene_allele_delta_tb <- find.pairwise.group.mutation.delta(arm_class_gene_allele_mat) %>% 
  separate(group_1, into = str_c(c('Arm', 'tumor_class'), 1, sep = '_'), sep = '-', remove = FALSE) %>% 
  separate(group_2, into = str_c(c('Arm', 'tumor_class'), 2, sep = '_'), sep = '-', remove = FALSE)

arm_class_gene_allele_delta_tb %>% 
  group_by(comparison) %>% 
  slice_max(order_by = abs(delta), n = 15) %>% 
  filter(Arm_1 != Arm_2, tumor_class_1 == tumor_class_2) %>% 
  mutate(comparison = parse_factor(comparison,
                                   levels = c('B-Q1_vs_A-Q1', 'A-Q2_vs_B-Q2', 'A-Q3_vs_B-Q3', 'B-Q4_vs_A-Q4'))) %>% 
  ggplot(aes(mutation, delta, fill = tumor_class_1)) + 
  geom_col() + 
  facet_wrap( ~ comparison, scales = 'free') + 
  coord_flip() + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7)) + 
  ylab('delta [(N_patient in group 1) - (N_patient in group 2)] per gene')

#######################################
#### Q: between arm and timepoint patient groups?

#### genes
pre_time_arm_gene_mat <- timepoint_arm_gene_tb %>% 
  rename(group = `Time-Arm`) %>% 
  select(-c(Arm, timepoint, N_total)) %>% 
  pivot_wider(names_from = Gene, values_from = N_patient, values_fill = 0)

time_arm_gene_mat <- pre_time_arm_gene_mat %>%  
  select(-group) %>%  
  as.matrix()
rownames(time_arm_gene_mat) <- pre_time_arm_gene_mat %>% pull(group)

time_arm_gene_pvalue_tb <- my_homogeneity_test(time_arm_gene_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'Arm'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'Arm'), 2, sep = '_'), sep = '-')

time_arm_gene_pvalue_tb %>% 
  filter(Arm_1 == Arm_2, timepoint_1 != timepoint_2)

time_arm_gene_pvalue_tb %>% 
  filter(Arm_1 != Arm_2, timepoint_1 == timepoint_2)

#### gene-variant
pre_time_arm_gene_variant_mat <- timepoint_arm_gene_variant_tb %>% 
  rename(group = `Time-Arm`) %>% 
  select(-c(Arm, timepoint, N_total)) %>% 
  pivot_wider(names_from = gene_variant, values_from = N_patient, values_fill = 0)

time_arm_gene_variant_mat <- pre_time_arm_gene_variant_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(time_arm_gene_variant_mat) <- pre_time_arm_gene_variant_mat %>% pull(group)

time_arm_gene_variant_pvalue_tb <- my_homogeneity_test(time_arm_gene_variant_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'Arm'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'Arm'), 2, sep = '_'), sep = '-')

time_arm_gene_variant_pvalue_tb %>% 
  filter(Arm_1 == Arm_2, timepoint_1 != timepoint_2)

time_arm_gene_variant_pvalue_tb %>% 
  filter(Arm_1 != Arm_2, timepoint_1 == timepoint_2)

#### gene-allele
pre_time_arm_gene_allele_mat <- timepoint_arm_gene_allele_tb %>% 
  rename(group = `Time-Arm`) %>% 
  select(-c(Arm, timepoint, N_total)) %>% 
  pivot_wider(names_from = gene_mutation, values_from = N_patient, values_fill = 0)

time_arm_gene_allele_mat <- pre_time_arm_gene_allele_mat %>%  
  select(-group) %>% 
  as.matrix()
rownames(time_arm_gene_allele_mat) <- pre_time_arm_gene_allele_mat %>% pull(group)

time_arm_gene_allele_pvalue_tb <- my_homogeneity_test(time_arm_gene_allele_mat) %>% 
  separate(group_1, into = str_c(c('timepoint', 'Arm'), 1, sep = '_'), sep = '-') %>% 
  separate(group_2, into = str_c(c('timepoint', 'Arm'), 2, sep = '_'), sep = '-')

time_arm_gene_allele_pvalue_tb %>% 
  filter(Arm_1 == Arm_2, timepoint_1 != timepoint_2)

time_arm_gene_allele_pvalue_tb %>% 
  filter(Arm_1 != Arm_2, timepoint_1 == timepoint_2)

###################################################################################################################################################################################################
###################################################################################################################################################################################################
#### ms plots

p <- cox_data_tb %>% 
  select(customer, PFS_time, Arm) %>% 
  unique() %>% 
  ggplot(aes(PFS_time, fill = Arm)) + 
  geom_density(alpha = .5) +
  geom_vline(xintercept = quart_progress[1],
             color = "black") +
  geom_vline(xintercept = quart_progress[2],
             color = "black") +
  geom_vline(xintercept = quart_progress[3],
             color = "black") +
  scale_fill_manual(values = c("blue", "red")) + 
  # ggtitle("Distribution of progression-free survival time in each treatment [n = 190]",
  #         subtitle = "black verticle lines mark where PFS quartiles split") + 
  xlab("Progression Free Survival Time")

ggsave(here('ms_plots', 'pfs_time_density.png'),
       p,
       dpi = 800)

p <- cox_data_tb %>% 
  select(customer, PFS_time, Arm) %>% 
  unique() %>% 
  ggplot(aes(log2(.5 + PFS_time), fill = Arm)) + 
  geom_density(alpha = .5) +
  geom_vline(xintercept = log2(.5 + quart_progress[1]),
             color = "black") +
  geom_vline(xintercept = log2(.5 + quart_progress[2]),
             color = "black") +
  geom_vline(xintercept = log2(.5 + quart_progress[3]),
             color = "black") +
  scale_fill_manual(values = c("blue", "red")) + 
  # ggtitle("Distribution of progression-free survival time in each treatment [n = 190]",
  #         subtitle = "black verticle lines mark where PFS quartiles split") + 
  xlab("Logarithmic Progression Free Survival Time [log2(.5 + PFS time)]")

ggsave(here('ms_plots', 'log_pfs_time_density.png'),
       p,
       dpi = 800)

p <- mutation_load_tb %>%
  select(customer, timepoint, Arm, mutation_load) %>% 
  unique() %>% 
  group_by(timepoint, Arm) %>% 
  summarize(group_mutation_load = sum(mutation_load)) %>% 
  ungroup() %>%  
  ggplot(aes(timepoint, group_mutation_load, fill = Arm)) + 
  geom_col() +
  scale_fill_manual(values = c('blue', 'red')) +
  xlab('time point') + 
  ylab('total number of gene mutations [CNV + SNV + Indel]')

ggsave(here('ms_plots', 'total_mutations_per_timepoint.png'),
       p,
       dpi = 800)


p <- mutation_load_tb %>% 
  filter(mutation_load != 0) %>% 
  mutate(Variant_type = parse_factor(Variant_type, levels = c('SNV', 'CNV', 'Indel'))) %>% 
  ggplot(aes(Variant_type, fraction_variant_mutation_load, fill = Variant_type)) + 
  geom_violin() + 
  facet_wrap( ~ timepoint) +
  theme(legend.position = 'none') + 
  xlab('variant') + 
  ylab('fraction of mutation load per patient')

ggsave(here('ms_plots', 'fraction_variant_mutation_load_per_timepoint.png'),
       p,
       dpi = 800)


p <- mutation_load_tb %>% 
  filter(mutation_load != 0) %>% 
  mutate(Variant_type = parse_factor(Variant_type, levels = c('SNV', 'CNV', 'Indel'))) %>% 
  ggplot(aes(Variant_type, fraction_variant_mutation_load, fill = Variant_type)) + 
  geom_violin() + 
  facet_wrap(Arm ~ timepoint) +
  theme(legend.position = 'none') + 
  xlab('variant') + 
  ylab('fraction of mutation load per patient')

ggsave(here('ms_plots', 'fraction_variant_mutation_load_per_timepoint_and_arm.png'),
       p,
       dpi = 800)


p <- mutation_load_tb %>% 
  select(-c(Variant_type, variant_mutation_load, fraction_variant_mutation_load)) %>% 
  unique() %>% 
  ggplot(aes(log2(1 + mutation_load))) + 
  geom_density(fill = 'purple') + 
  facet_wrap(timepoint ~ tumor_class) +
  xlab('log2(1 + mutation load)') + 
  ylab('density')

ggsave(here('ms_plots', 'logarithmic_mutation_load_per_timepoint_and_quartile.png'),
       p,
       dpi = 800)

p <- mutation_load_tb %>% 
  mutate(p = 1) %>% 
  group_by(timepoint, tumor_class, count_bin) %>% 
  summarize(N = sum(p)) %>% 
  ungroup() %>% 
  ggplot(aes(count_bin, N)) + 
  geom_col(fill = 'purple', alpha = .85) + 
  facet_wrap(timepoint ~ tumor_class) +
  xlab('binned mutation load') + 
  ylab('number of patients')

ggsave(here('ms_plots', 'binned_mutation_load_per_timepoint_and_quartile.png'),
       p,
       dpi = 800)


p <- mutation_load_tb %>% 
  select(-c(Variant_type, variant_mutation_load, fraction_variant_mutation_load)) %>% 
  unique() %>% 
  filter(Arm == 'A') %>% 
  ggplot(aes(log2(1 + mutation_load))) + 
  geom_density(fill = 'blue') + 
  facet_wrap(timepoint ~ tumor_class) +
  xlab('log2(1 + mutation load)') + 
  ylab('density')

ggsave(here('ms_plots', 'arm_A_logarithmic_mutation_load_per_timepoint_and_quartile.png'),
       p,
       dpi = 800)

p <- mutation_load_tb %>% 
  filter(Arm == 'A') %>% 
  mutate(p = 1) %>% 
  group_by(timepoint, tumor_class, count_bin) %>% 
  summarize(N = sum(p)) %>% 
  ungroup() %>% 
  ggplot(aes(count_bin, N)) + 
  geom_col(fill = 'blue', alpha = .95) + 
  facet_wrap(timepoint ~ tumor_class) +
  xlab('binned mutation load') + 
  ylab('number of patients')

ggsave(here('ms_plots', 'arm_A_binned_mutation_load_per_timepoint_and_quartile.png'),
       p,
       dpi = 800)

p <- mutation_load_tb %>% 
  select(-c(Variant_type, variant_mutation_load, fraction_variant_mutation_load)) %>% 
  unique() %>% 
  filter(Arm == 'B') %>% 
  ggplot(aes(log2(1 + mutation_load))) + 
  geom_density(fill = 'red') + 
  facet_wrap(timepoint ~ tumor_class) +
  xlab('log2(1 + mutation load)') + 
  ylab('density')

ggsave(here('ms_plots', 'arm_B_logarithmic_mutation_load_per_timepoint_and_quartile.png'),
       p,
       dpi = 800)

p <- mutation_load_tb %>% 
  filter(Arm == 'B') %>% 
  mutate(p = 1) %>% 
  group_by(timepoint, tumor_class, count_bin) %>% 
  summarize(N = sum(p)) %>% 
  ungroup() %>% 
  ggplot(aes(count_bin, N)) + 
  geom_col(fill = 'red', alpha = .95) + 
  facet_wrap(timepoint ~ tumor_class) +
  xlab('binned mutation load') + 
  ylab('number of patients')

ggsave(here('ms_plots', 'arm_B_binned_mutation_load_per_timepoint_and_quartile.png'),
       p,
       dpi = 800)



load_change_tb %>% 
  ggplot(aes(tumor_class, log2(fc_1 + fc_2) - log2(fc_1*fc_2))) + 
  geom_boxplot(fill = 'purple') +
  coord_flip() +
  xlab('PFS quartile') + 
  ylab('log2[(PRE to OT fold change) + (OT to EOS fold change)] - log2[ PRE to EOS fold change ]')

load_change_tb %>% 
  ggplot(aes(tumor_class, log2(fc_1 + fc_2) - log2(fc_1*fc_2), fill = Arm)) + 
  geom_boxplot() +
  scale_fill_manual(values = c('blue', 'red')) +
  coord_flip() +
  xlab('PFS quartile') + 
  ylab('log2[(PRE to OT fold change) + (OT to EOS fold change)] - log2[ PRE to EOS fold change ]')

load_change_tb %>% 
  ggplot(aes(tumor_class, my_score)) + 
  geom_boxplot(fill = 'purple') +
  coord_flip() +
  xlab('PFS quartile') + 
  ylab('scaled fold-change of mutation load per patient')

load_change_tb %>% 
  ggplot(aes(tumor_class, my_score, fill = Arm)) + 
  geom_boxplot() +
  scale_fill_manual(values = c('blue', 'red')) +
  coord_flip() +
  xlab('PFS quartile') + 
  ylab('scaled fold-change of mutation load per patient')



















