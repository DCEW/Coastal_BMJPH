# Sensitivity analyses for coastal paper 1

#S2.1 - Count of England LSOAs defined as 'coastal' based on each derived 
#candidate variable. Total counts are those available for modelling 
#(i.e. after excluding outliers).

# found at line 1166 in 2.Test_candidates.R

#s3.1 LSOA11s by G_25_5 definition

#s3.2 LSOA21s by G_25_5 definition

# Start : add mean_dist---------------------------------------
rm(list=ls())
pacman::p_load(tidyverse, janitor, readxl, Hmisc, modelr, broom, purrr, here)
dat <- read.csv(paste0(Sys.getenv("coastal_clean"),"/modeldata.csv"))
source("functions.R")

#association - correlation matrices (age, sex, ethnicity, IMDquintile)
candidates <- c("GR_10_5", "G_25_5")

# for each candidate, test its correlation with the covariates
covars <- c("meanage", "pcfem", "pcwhite", "quint", "rural", "mean_dist")


# Modelling ---------------------------------------------------------------

#geography only 
g_models <- c("g_25_5")

#rural
gr_models <- c("gr_10_5")

modelsets <- c(g_models, gr_models)

#'outcome' = combination of cancer site and metric e.g. lung prevalence
dat <- dat |>
  select(-c(X.1, X))
idcols <- names(dat)[names(dat) %in% c("LSOA11","RGN","rural","quint","coastline","allpcd","LAD20CD", "LAD20NM","n","YEAR", "total")]
counts <- names(dat)[grepl("Cases", names(dat))]
ncands <- names(dat)[grepl("G", names(dat))]
outcome_cols <- names(dat)[!names(dat) %in% c(counts, covars, candidates, idcols, ncands)]

dat_long <- dat |>
  pivot_longer(cols=all_of(outcome_cols), names_to = "outcome", values_to = "value") |>
  mutate(cancer = case_when(grepl("Lung",outcome) ~ "Lung",
                            grepl("All", outcome) ~ "All",
                            grepl("Breast", outcome) ~ "Breast",
                            grepl("Prostate", outcome) ~ "Prostate",
                            grepl("Melanoma", outcome) ~ "Melanoma",
                            grepl("chd", outcome) ~ "CHD",
                            .default = "Bowel"),
         out_type = case_when(grepl("inc_Cases", outcome) ~ "IncidenceCount",
                              grepl("inc_CR", outcome) ~ "CrudeIncidence",
                              grepl("inc_ASR", outcome) ~ "Incidence_ASR",
                              grepl("mort_CR", outcome) ~ "Mortaltiy_Crude",
                              grepl("mort_Cases", outcome) ~ "Mortality_Count",
                              grepl("mort_ASR", outcome) ~ "Mortality_ASR",
                              grepl("prev_Cases", outcome) ~ "PrevalenceCases",
                              grepl("prev_CR|_prev", outcome) ~ "Prevalence_Crude"),
         sex = case_when(grepl("^P_|^p", outcome) ~ "Persons",
                         grepl("^M_", outcome) ~ "Males",
                         grepl("^F_", outcome) ~ "Females")) |>
  filter(!out_type %in% c("CrudeIncidence", "DeathsCrude"))


# Exclude LSOA outliers -----
dat_long <- dat_long |>
  filter(!LSOA11 %in% c("E01033739","E01032703"))

dat <- dat |>
  filter(!LSOA11 %in% c("E01033739","E01032703"))

by_out <- dat_long |>
  group_by(outcome) |>
  nest()

g_25_5_model <- function(df){
  lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite + mean_dist, data = df)
}

gr_10_5_model <- function(df){
  lm(value ~ GR_10_5 + quint + meanage + pcfem + pcwhite + mean_dist, data = df)
}

models <- c(gr_10_5_model, g_25_5_model)

names(models) <- modelsets

b <- map(modelsets, addmodels, by_out)
bo <- b |>
  reduce(inner_join, by=c("outcome","data"))

by_out <- bo |>
  mutate(resids_gr_10_5 = map2(data, gr_10_5, add_residuals),
         resids_g_25_5 = map2(data, g_25_5, add_residuals),
         pred_gr_10_5 = map2(data, gr_10_5, add_predictions),
         pred_g_25_5 = map2(data, g_25_5, add_predictions))

glances <- bind_rows(map(modelsets, glance_fun, by_out))

by_out <- by_out |>
  mutate(summary_gr_10_5 = map(gr_10_5, summary),
         summary_g_25_5 = map(g_25_5, summary),
         tidy_g_25_5 = map(summary_g_25_5, tidy),
         tidy_gr_10_5 = map(summary_gr_10_5, tidy))

ests <- data.frame()

for ( i in modelsets){
  print(i)
  est <- by_out |>
    select(outcome, paste0("tidy_",i)) |>
    unnest(paste0("tidy_",i)) |>
    filter(term %in% candidates) |>
    pivot_longer(-c(outcome, term), names_to = "varpart", values_to = "value")
  ests <- bind_rows(ests,est)
}

ests <- ests |>
  mutate(cancer = case_when(grepl("Lung",outcome) ~ "Lung",
                            grepl("All", outcome) ~ "All",
                            grepl("Breast", outcome) ~ "Breast",
                            grepl("Prostate", outcome) ~ "Prostate",
                            grepl("Melanoma", outcome) ~ "Melanoma",
                            grepl("chd", outcome) ~ "CHD",
                            .default = "Bowel"),
         out_type = case_when(grepl("inc_Cases", outcome) ~ "IncidenceCount",
                              grepl("inc_CR", outcome) ~ "CrudeIncidence",
                              grepl("inc_ASR", outcome) ~ "Incidence_ASR",
                              grepl("mort_CR", outcome) ~ "Mortaltiy_Crude",
                              grepl("mort_Cases", outcome) ~ "Mortality_Count",
                              grepl("mort_ASR", outcome) ~ "Mortality_ASR",
                              grepl("prev_Cases", outcome) ~ "PrevalenceCases",
                              grepl("prev_CR|_prev", outcome) ~ "Prevalence_Crude"),
         sex = case_when(grepl("^P_|^p", outcome) ~ "Persons",
                         grepl("^M_", outcome) ~ "Males",
                         grepl("^F_", outcome) ~ "Females"),
         modelform = case_when(grepl("^G1|^G_", term) ~ "Geographic",
                               grepl("^GD1|^GD_", term) ~ "Geog/Deprivation",
                               grepl("^GR", term) ~ "Geog/Rural",
                               .default = "Geog/Rural/Deprivation"))  |>
  mutate(model = tolower(term))

all <- inner_join(glances, ests[ests$varpart == "p.value"|ests$varpart == "estimate",], by=c("outcome","model"))

all$sig <- ifelse(all$varpart == "p.value" & all$value < 0.05, 1,0)

library(lmtest)
w <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite + mean_dist, data = dat_long[dat_long$outcome == "P_inc_ASR_All",])
wo <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite , data = dat_long[dat_long$outcome == "P_inc_ASR_All",])
lrtest(w, wo)
BIC(w)
BIC(wo)
# with mean_dist is better
summary(w)
# all variables sig

w <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite + mean_dist, data = dat_long[dat_long$outcome == "prev_CR_All",])
wo <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite , data = dat_long[dat_long$outcome == "prev_CR_All",])
lrtest(w, wo)
BIC(w)
BIC(wo)

summary(w)
# all variables sig
## The coastal variable explains variation in all cancer incidence and 
#prevalence even after mean GP distance accounted for 


#-------------------------------------------------------------------------------------
#SA 2 - overfit model
rm(list=ls())
pacman::p_load(tidyverse, janitor, readxl, Hmisc, modelr, broom, purrr, here)
dat <- read.csv(paste0(Sys.getenv("coastal_clean"),"/modeldata.csv"))
source("functions.R")

#add extra lsoa descriptive data
#https://data.cdrc.ac.uk/dataset/access-healthy-assets-hazards-ahah
ahah <- read.csv(paste0(Sys.getenv("coastal_clean"), "AHAH_V3_0.csv")) |>
  select(lsoa11, ah3gp, ah3hosp, ah3e) |>
  filter(grepl("E", lsoa11))
#https://pldr.org/dataset/2noyv/small-area-mental-health-index-samhi
samhi <- read.csv(paste0(Sys.getenv("coastal_clean"), "samhi_21_01_v4.00_2011_2019_LSOA_tall.csv")) |>
  filter(year == 2011) |>
  select(lsoa11, samhi_index)

#add these to the existing data
#association - correlation matrices (age, sex, ethnicity, IMDquintile)
candidates <- c("G1", "G_10_1", "G_10_3", "G_10_5", 
                "GD1", "GD_10_1", "GD_10_3", "GD_10_5", 
                "GR1", "GR_10_1", "GR_10_3", "GR_10_5", 
                "GDR1", "GDR_10_1", "GDR_10_3", "GDR_10_5",
                "G_25_1","G_25_3","G_25_5",
                "GD_25_1","GD_25_3","GD_25_5",
                "GR_25_1","GR_25_3","GR_25_5",
                "GDR_25_1","GDR_25_3","GDR_25_5",
                "G_50_1","G_50_3","G_50_5",
                "GD_50_1","GD_50_3","GD_50_5",
                "GR_50_1","GR_50_3","GR_50_5",
                "GDR_50_1","GDR_50_3","GDR_50_5")

# for each candidate, test its correlation with the covariates
covars <- c("meanage", "pcfem", "pcwhite", "quint", "rural", "mean_dist", "ah3hosp", "ah3e", "samhi_index")


# Modelling ---------------------------------------------------------------

#geography only (has coastline or has 10 postcodes within 1,3,5 km of coast)
g_models <- c("g1", "g_10_1", "g_10_3", "g_10_5", 
              "g_25_1", "g_25_3", "g_25_5",
              "g_50_1", "g_50_3", "g_50_5")
#prevalence ~ coastal + imd + rural + mean age + %female + %pcwhite
#incidence ~ coastal + imd + rural + mean age + %female + %pcwhite
#mortality ~ coastal + imd + rural + mean age + %female + %pcwhite

#imd
gd_models <- c("gd1", "gd_10_1", "gd_10_3", "gd_10_5", 
               "gd_25_1", "gd_25_3", "gd_25_5",
               "gd_50_1", "gd_50_3", "gd_50_5")
#prevalence ~ coastal + rural + mean age + %female + %pcwhite
#incidence ~ coastal + rural + mean age + %female + %pcwhite
#mortality ~ coastal + rural + mean age + %female + %pcwhite

#rural
gr_models <- c("gr1", "gr_10_1", "gr_10_3", "gr_10_5", 
               "gr_25_1", "gr_25_3", "gr_25_5",
               "gr_50_1", "gr_50_3", "gr_50_5")
#prevalence ~ coastal + imd + mean age + %female + %pcwhite
#incidence ~ coastal + imd + mean age + %female + %pcwhite
#mortality ~ coastal + imd + mean age + %female + %pcwhite

#combination models
gdr_models <- c("gdr1", "gdr_10_1", "gdr_10_3", "gdr_10_5", 
                "gdr_25_1", "gdr_25_3", "gdr_25_5",
                "gdr_50_1", "gdr_50_3", "gdr_50_5")

modelsets <- c(g_models, gd_models, gr_models, gdr_models)


#'outcome' = combination of cancer site and metric e.g. lung prevalence
dat <- dat |>
  select(-c(X.1, X)) |>
  left_join(ahah, by = c("LSOA11" = "lsoa11"))
dat <- dat |>
  left_join(samhi, by= c("LSOA11" = "lsoa11"))
idcols <- names(dat)[names(dat) %in% c("LSOA11","RGN","rural","quint","coastline","allpcd","LAD20CD", "LAD20NM","n","YEAR", "total", "ah3gp", "ah3e","ah3hposp","samhi_index")]
counts <- names(dat)[grepl("Cases", names(dat))]
ncands <- names(dat)[grepl("^G", names(dat))]
outcome_cols <- names(dat)[!names(dat) %in% c(counts, covars, candidates, idcols, ncands)]

dat_long <- dat |>
  pivot_longer(cols=all_of(outcome_cols), names_to = "outcome", values_to = "value") |>
  mutate(cancer = case_when(grepl("Lung",outcome) ~ "Lung",
                            grepl("All", outcome) ~ "All",
                            grepl("Breast", outcome) ~ "Breast",
                            grepl("Prostate", outcome) ~ "Prostate",
                            grepl("Melanoma", outcome) ~ "Melanoma",
                            grepl("chd", outcome) ~ "CHD",
                            .default = "Bowel"),
         out_type = case_when(grepl("inc_Cases", outcome) ~ "IncidenceCount",
                              grepl("inc_CR", outcome) ~ "CrudeIncidence",
                              grepl("inc_ASR", outcome) ~ "Incidence_ASR",
                              grepl("mort_CR", outcome) ~ "Mortaltiy_Crude",
                              grepl("mort_Cases", outcome) ~ "Mortality_Count",
                              grepl("mort_ASR", outcome) ~ "Mortality_ASR",
                              grepl("prev_Cases", outcome) ~ "PrevalenceCases",
                              grepl("prev_CR|_prev", outcome) ~ "Prevalence_Crude"),
         sex = case_when(grepl("^P_|^p", outcome) ~ "Persons",
                         grepl("^M_", outcome) ~ "Males",
                         grepl("^F_", outcome) ~ "Females")) |>
  filter(!out_type %in% c("CrudeIncidence", "DeathsCrude"))


# Exclude LSOA outliers -----
dat_long <- dat_long |>
  filter(!LSOA11 %in% c("E01033739","E01032703"))

#correlations
cors <- data.frame(var = c(NA, covars))
for(cand in candidates){
  data <- dat[,c(cand, covars)] 
  data[[cand]][data[[cand]] == 2] <- 0 #set the 'comparator' to 0 so all are binary
  out <- bind_cols(rcorr(as.matrix(data), type = "spearman")$r[1,], rcorr(as.matrix(data), type = "spearman")$P[1,]) |>
    mutate(sig = ifelse(`...2`<0.05, `...1`,NA)) |> #select only the significant ones
    select(sig)
  names(out)[1] <- names(data)[1]
  out <- out |>
    mutate(var = names(data)) |>
    select(var, 1)
  cors <- bind_cols(cors, out) 
}
odd <- seq_len(ncol(cors)) %% 2
corsodd <- cors[,odd==1]

longcors <- corsodd |>
  pivot_longer(cols = -c("var...1"), names_to = "cand",values_to = "Pearson") |>
  filter(!is.na(Pearson))

cor_heatmap <- longcors |>
  mutate(cov = case_when(`var...1` == "meanage" ~ "Mean Age (years)",
                         `var...1` == "pcfem" ~ "Proportion Female",
                         `var...1` == "pcwhite" ~ "Proportion White Ethnicity",
                         `var...1` == "rural" ~ "Rural RUC11 code",
                         `var...1` == "total" ~ "Population in 2016",
                         `var...1` == "mean_dist" ~ "Mean distance to GP",
                         `var...1` == "ah3hosp" ~ "Distance to hospital",
                         `var...1` == "ah3e" ~ "Air quality domain score",
                         `var...1` == "samhi_index" ~ "SAMHI Index",
                         .default = "IMD Quintile"),
         bin = ntile(abs(Pearson), n=5)) |>
  ggplot(aes(x=cov, y=cand, fill = as.factor(bin))) +
  geom_tile() +
  geom_text(aes(label = round(Pearson,2)), colour = "black", size = 2.5) +
  labs(x="Covariate", y="Candidate variable", 
       title = "Candidate:Covariate Correlation Coefficients",
       caption = "Statistically significant (P<0.05) Spearman correlation coefficients shown") +
  guides(fill = "none") +
  theme_minimal() +
  scale_fill_brewer() 
cor_heatmap

#remodel with these added
# Define models -----

by_out <- dat_long |>
  group_by(outcome) |>
  nest()

g1_model <- function(df){
  lm(value ~ G1 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_10_1_model <- function(df){
  lm(value ~ G_10_1 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_10_3_model <- function(df){
  lm(value ~ G_10_3 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_10_5_model <- function(df){
  lm(value ~ G_10_5 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_25_1_model <- function(df){
  lm(value ~ G_25_1 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_25_3_model <- function(df){
  lm(value ~ G_25_3 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_25_5_model <- function(df){
  lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_50_1_model <- function(df){
  lm(value ~ G_50_1 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_50_3_model <- function(df){
  lm(value ~ G_50_3 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
g_50_5_model <- function(df){
  lm(value ~ G_50_5 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}


gd1_model <- function(df){
  lm(value ~ GD1 +rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_10_1_model <- function(df){
  lm(value ~ GD_10_1 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_10_3_model <- function(df){
  lm(value ~ GD_10_3 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_10_5_model <- function(df){
  lm(value ~ GD_10_5 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_25_1_model <- function(df){
  lm(value ~ GD_25_1 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_25_3_model <- function(df){
  lm(value ~ GD_25_3 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_25_5_model <- function(df){
  lm(value ~ GD_25_5 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_50_1_model <- function(df){
  lm(value ~ GD_50_1 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_50_3_model <- function(df){
  lm(value ~ GD_50_3 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gd_50_5_model <- function(df){
  lm(value ~ GD_50_5 + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}



gr1_model <- function(df){
  lm(value ~ GR1 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_10_1_model <- function(df){
  lm(value ~ GR_10_1 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_10_3_model <- function(df){
  lm(value ~ GR_10_3 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_10_5_model <- function(df){
  lm(value ~ GR_10_5 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_25_1_model <- function(df){
  lm(value ~ GR_25_1 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_25_3_model <- function(df){
  lm(value ~ GR_25_3 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_25_5_model <- function(df){
  lm(value ~ GR_25_5 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_50_1_model <- function(df){
  lm(value ~ GR_50_1 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_50_3_model <- function(df){
  lm(value ~ GR_50_3 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gr_50_5_model <- function(df){
  lm(value ~ GR_50_5 + quint + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}

gdr1_model <- function(df){
  lm(value ~ GDR1 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_10_1_model <- function(df){
  lm(value ~ GDR_10_1 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_10_3_model <- function(df){
  lm(value ~ GDR_10_3 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_10_5_model <- function(df){
  lm(value ~ GDR_10_5 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_25_1_model <- function(df){
  lm(value ~ GDR_25_1 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_25_3_model <- function(df){
  lm(value ~ GDR_25_3 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_25_5_model <- function(df){
  lm(value ~ GDR_25_5 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_50_1_model <- function(df){
  lm(value ~ GDR_50_1 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_50_3_model <- function(df){
  lm(value ~ GDR_50_3 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}
gdr_50_5_model <- function(df){
  lm(value ~ GDR_50_5 + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = df)
}

models <- c(g1_model, g_10_1_model, g_10_3_model, g_10_5_model,
            g_25_1_model, g_25_3_model, g_25_5_model,
            g_50_1_model, g_50_3_model, g_50_5_model,
            gd1_model, gd_10_1_model, gd_10_3_model, gd_10_5_model,
            gd_25_1_model, gd_25_3_model, gd_25_5_model,
            gd_50_1_model, gd_50_3_model, gd_50_5_model,
            gr1_model, gr_10_1_model, gr_10_3_model, gr_10_5_model,
            gr_25_1_model, gr_25_3_model, gr_25_5_model,
            gr_50_1_model, gr_50_3_model, gr_50_5_model,
            gdr1_model, gdr_10_1_model, gdr_10_3_model, gdr_10_5_model,
            gdr_25_1_model, gdr_25_3_model, gdr_25_5_model,
            gdr_50_1_model, gdr_50_3_model, gdr_50_5_model)

names(models) <- modelsets

b <- map(modelsets, addmodels, by_out)
bo <- b |>
  reduce(inner_join, by=c("outcome","data"))

by_out <- bo |>
  mutate(resids_g1 = map2(data, g1, add_residuals),
         resids_g_10_1 = map2(data, g_10_1, add_residuals),
         resids_g_10_3 = map2(data, g_10_3, add_residuals),
         resids_g_10_5 = map2(data, g_10_5, add_residuals),
         resids_g_25_1 = map2(data, g_25_1, add_residuals),
         resids_g_25_3 = map2(data, g_25_3, add_residuals),
         resids_g_25_5 = map2(data, g_25_5, add_residuals),
         resids_g_50_1 = map2(data, g_50_1, add_residuals),
         resids_g_50_3 = map2(data, g_50_3, add_residuals),
         resids_g_50_5 = map2(data, g_50_5, add_residuals),
         resids_gd1 = map2(data, gd1, add_residuals),
         resids_gd_10_1 = map2(data, gd_10_1, add_residuals),
         resids_gd_10_3 = map2(data, gd_10_3, add_residuals),
         resids_gd_10_5 = map2(data, gd_10_5, add_residuals),
         resids_gd_25_1 = map2(data, gd_25_1, add_residuals),
         resids_gd_25_3 = map2(data, gd_25_3, add_residuals),
         resids_gd_25_5 = map2(data, gd_25_5, add_residuals),
         resids_gd_50_1 = map2(data, gd_50_1, add_residuals),
         resids_gd_50_3 = map2(data, gd_50_3, add_residuals),
         resids_gd_50_5 = map2(data, gd_50_5, add_residuals),
         resids_gr1 = map2(data, gr1, add_residuals),
         resids_gr_10_1 = map2(data, gr_10_1, add_residuals),
         resids_gr_10_3 = map2(data, gr_10_3, add_residuals),
         resids_gr_10_5 = map2(data, gr_10_5, add_residuals),
         resids_gr_25_1 = map2(data, gr_25_1, add_residuals),
         resids_gr_25_3 = map2(data, gr_25_3, add_residuals),
         resids_gr_25_5 = map2(data, gr_25_5, add_residuals),
         resids_gr_50_1 = map2(data, gr_50_1, add_residuals),
         resids_gr_50_3 = map2(data, gr_50_3, add_residuals),
         resids_gr_50_5 = map2(data, gr_50_5, add_residuals),
         resids_gdr1 = map2(data, gdr1, add_residuals),
         resids_gdr_10_1 = map2(data, gdr_10_1, add_residuals),
         resids_gdr_10_3 = map2(data, gdr_10_3, add_residuals),
         resids_gdr_10_5 = map2(data, gdr_10_5, add_residuals),
         resids_gdr_25_1 = map2(data, gdr_25_1, add_residuals),
         resids_gdr_25_3 = map2(data, gdr_25_3, add_residuals),
         resids_gdr_25_5 = map2(data, gdr_25_5, add_residuals),
         resids_gdr_50_1 = map2(data, gdr_50_1, add_residuals),
         resids_gdr_50_3 = map2(data, gdr_50_3, add_residuals),
         resids_gdr_50_5 = map2(data, gdr_50_5, add_residuals),
         pred_g1 = map2(data, g1, add_predictions),
         pred_g_10_1 = map2(data, g_10_1, add_predictions),
         pred_g_10_3 = map2(data, g_10_3, add_predictions),
         pred_g_10_5 = map2(data, g_10_5, add_predictions),
         pred_g_25_1 = map2(data, g_25_1, add_predictions),
         pred_g_25_3 = map2(data, g_25_3, add_predictions),
         pred_g_25_5 = map2(data, g_25_5, add_predictions),
         pred_g_50_1 = map2(data, g_50_1, add_predictions),
         pred_g_50_3 = map2(data, g_50_3, add_predictions),
         pred_g_50_5 = map2(data, g_50_5, add_predictions),
         pred_gd1 = map2(data, gd1, add_predictions),
         pred_gd_10_1 = map2(data, gd_10_1, add_predictions),
         pred_gd_10_3 = map2(data, gd_10_3, add_predictions),
         pred_gd_10_5 = map2(data, gd_10_5, add_predictions),
         pred_gd_25_1 = map2(data, gd_25_1, add_predictions),
         pred_gd_25_3 = map2(data, gd_25_3, add_predictions),
         pred_gd_25_5 = map2(data, gd_25_5, add_predictions),
         pred_gd_50_1 = map2(data, gd_50_1, add_predictions),
         pred_gd_50_3 = map2(data, gd_50_3, add_predictions),
         pred_gd_50_5 = map2(data, gd_50_5, add_predictions),
         pred_gr1 = map2(data, gr1, add_predictions),
         pred_gr_10_1 = map2(data, gr_10_1, add_predictions),
         pred_gr_10_3 = map2(data, gr_10_3, add_predictions),
         pred_gr_10_5 = map2(data, gr_10_5, add_predictions),
         pred_gr_25_1 = map2(data, gr_25_1, add_predictions),
         pred_gr_25_3 = map2(data, gr_25_3, add_predictions),
         pred_gr_25_5 = map2(data, gr_25_5, add_predictions),
         pred_gr_50_1 = map2(data, gr_50_1, add_predictions),
         pred_gr_50_3 = map2(data, gr_50_3, add_predictions),
         pred_gr_50_5 = map2(data, gr_50_5, add_predictions),
         pred_gdr1 = map2(data, gdr1, add_predictions),
         pred_gdr_10_1 = map2(data, gdr_10_1, add_predictions),
         pred_gdr_10_3 = map2(data, gdr_10_3, add_predictions),
         pred_gdr_10_5 = map2(data, gdr_10_5, add_predictions),
         pred_gdr_25_1 = map2(data, gdr_25_1, add_predictions),
         pred_gdr_25_3 = map2(data, gdr_25_3, add_predictions),
         pred_gdr_25_5 = map2(data, gdr_25_5, add_predictions),
         pred_gdr_50_1 = map2(data, gdr_50_1, add_predictions),
         pred_gdr_50_3 = map2(data, gdr_50_3, add_predictions),
         pred_gdr_50_5 = map2(data, gdr_50_5, add_predictions))

glances <- bind_rows(map(modelsets, glance_fun, by_out))


by_out <- by_out |>
  mutate(summary_g1 = map(g1, summary),
         summary_g101 = map(g_10_1, summary),
         summary_g103 = map(g_10_3, summary),
         summary_g105 = map(g_10_5, summary),
         summary_g251 = map(g_25_1, summary),
         summary_g253 = map(g_25_3, summary),
         summary_g255 = map(g_25_5, summary),
         summary_g501 = map(g_50_1, summary),
         summary_g503 = map(g_50_3, summary),
         summary_g505 = map(g_50_5, summary),
         summary_gd1 = map(gd1, summary),
         summary_gd101 = map(gd_10_1, summary),
         summary_gd103 = map(gd_10_3, summary),
         summary_gd105 = map(gd_10_5, summary),
         summary_gd251 = map(gd_25_1, summary),
         summary_gd253 = map(gd_25_3, summary),
         summary_gd255 = map(gd_25_5, summary),
         summary_gd501 = map(gd_50_1, summary),
         summary_gd503 = map(gd_50_3, summary),
         summary_gd505 = map(gd_50_5, summary),
         summary_gr1 = map(gr1, summary),
         summary_gr101 = map(gr_10_1, summary),
         summary_gr103 = map(gr_10_3, summary),
         summary_gr105 = map(gr_10_5, summary),
         summary_gr251 = map(gr_25_1, summary),
         summary_gr253 = map(gr_25_3, summary),
         summary_gr255 = map(gr_25_5, summary),
         summary_gr501 = map(gr_50_1, summary),
         summary_gr503 = map(gr_50_3, summary),
         summary_gr505 = map(gr_50_5, summary),
         summary_gdr1 = map(gdr1, summary),
         summary_gdr101 = map(gdr_10_1, summary),
         summary_gdr103 = map(gdr_10_3, summary),
         summary_gdr105 = map(gdr_10_5, summary),
         summary_gdr251 = map(gdr_25_1, summary),
         summary_gdr253 = map(gdr_25_3, summary),
         summary_gdr255 = map(gdr_25_5, summary),
         summary_gdr501 = map(gdr_50_1, summary),
         summary_gdr503 = map(gdr_50_3, summary),
         summary_gdr505 = map(gdr_50_5, summary),
         tidy_g1 = map(summary_g1, tidy),
         tidy_g_10_1 = map(summary_g101, tidy),
         tidy_g_10_3 = map(summary_g103, tidy),
         tidy_g_10_5 = map(summary_g105, tidy),
         tidy_g_25_1 = map(summary_g251, tidy),
         tidy_g_25_3 = map(summary_g253, tidy),
         tidy_g_25_5 = map(summary_g255, tidy),
         tidy_g_50_1 = map(summary_g501, tidy),
         tidy_g_50_3 = map(summary_g503, tidy),
         tidy_g_50_5 = map(summary_g505, tidy),
         tidy_gd1 = map(summary_gd1, tidy),
         tidy_gd_10_1 = map(summary_gd101, tidy),
         tidy_gd_10_3 = map(summary_gd103, tidy),
         tidy_gd_10_5 = map(summary_gd105, tidy),
         tidy_gd_25_1 = map(summary_gd251, tidy),
         tidy_gd_25_3 = map(summary_gd253, tidy),
         tidy_gd_25_5 = map(summary_gd255, tidy),
         tidy_gd_50_1 = map(summary_gd501, tidy),
         tidy_gd_50_3 = map(summary_gd503, tidy),
         tidy_gd_50_5 = map(summary_gd505, tidy),
         tidy_gr1 = map(summary_gr1, tidy),
         tidy_gr_10_1 = map(summary_gr101, tidy),
         tidy_gr_10_3 = map(summary_gr103, tidy),
         tidy_gr_10_5 = map(summary_gr105, tidy),
         tidy_gr_25_1 = map(summary_gr251, tidy),
         tidy_gr_25_3 = map(summary_gr253, tidy),
         tidy_gr_25_5 = map(summary_gr255, tidy),
         tidy_gr_50_1 = map(summary_gr501, tidy),
         tidy_gr_50_3 = map(summary_gr503, tidy),
         tidy_gr_50_5 = map(summary_gr505, tidy),
         tidy_gdr1 = map(summary_gdr1, tidy),
         tidy_gdr_10_1 = map(summary_gdr101, tidy),
         tidy_gdr_10_3 = map(summary_gdr103, tidy),
         tidy_gdr_10_5 = map(summary_gdr105, tidy),
         tidy_gdr_25_1 = map(summary_gdr251, tidy),
         tidy_gdr_25_3 = map(summary_gdr253, tidy),
         tidy_gdr_25_5 = map(summary_gdr255, tidy),
         tidy_gdr_50_1 = map(summary_gdr501, tidy),
         tidy_gdr_50_3 = map(summary_gdr503, tidy),
         tidy_gdr_50_5 = map(summary_gdr505, tidy))


ests <- data.frame()

for ( i in modelsets){
  print(i)
  est <- by_out |>
    select(outcome, paste0("tidy_",i)) |>
    unnest(paste0("tidy_",i)) |>
    filter(term %in% candidates) |>
    pivot_longer(-c(outcome, term), names_to = "varpart", values_to = "value")
  ests <- bind_rows(ests,est)
}

ests <- ests |>
  mutate(cancer = case_when(grepl("Lung",outcome) ~ "Lung",
                            grepl("All", outcome) ~ "All",
                            grepl("Breast", outcome) ~ "Breast",
                            grepl("Prostate", outcome) ~ "Prostate",
                            grepl("Melanoma", outcome) ~ "Melanoma",
                            grepl("chd", outcome) ~ "CHD",
                            .default = "Bowel"),
         out_type = case_when(grepl("inc_Cases", outcome) ~ "IncidenceCount",
                              grepl("inc_CR", outcome) ~ "CrudeIncidence",
                              grepl("inc_ASR", outcome) ~ "Incidence_ASR",
                              grepl("mort_CR", outcome) ~ "Mortaltiy_Crude",
                              grepl("mort_Cases", outcome) ~ "Mortality_Count",
                              grepl("mort_ASR", outcome) ~ "Mortality_ASR",
                              grepl("prev_Cases", outcome) ~ "PrevalenceCases",
                              grepl("prev_CR|_prev", outcome) ~ "Prevalence_Crude"),
         sex = case_when(grepl("^P_|^p", outcome) ~ "Persons",
                         grepl("^M_", outcome) ~ "Males",
                         grepl("^F_", outcome) ~ "Females"),
         modelform = case_when(grepl("^G1|^G_", term) ~ "Geographic",
                               grepl("^GD1|^GD_", term) ~ "Geog/Deprivation",
                               grepl("^GR", term) ~ "Geog/Rural",
                               .default = "Geog/Rural/Deprivation"))  |>
  mutate(model = tolower(term))

all <- inner_join(glances, ests[ests$varpart == "p.value"|ests$varpart == "estimate",], by=c("outcome","model"))

all$sig <- ifelse(all$varpart == "p.value" & all$value < 0.05, 1,0)

allcan <- all |> 
  filter(cancer == "All", varpart == "p.value") |> 
  mutate(sig = ifelse(value < 0.05, 1,0)) |>
  group_by(out_type) |>
  mutate(scaledbic = ifelse(sig ==1, scale(BIC) |>as.vector(), NA))

g255<- all |>
  filter(model == "g_25_5",
         varpart == "p.value")

best <- all |>
  select(-sig) |>
  filter(outcome %in% c("P_inc_ASR_All","prev_CR_All")) |>
  pivot_wider(names_from = varpart, values_from = value) |>
  mutate(sig = ifelse(p.value < 0.05, 1,0)) |>
  group_by(out_type) |>
  mutate(scaledbic = ifelse(sig ==1, scale(BIC) |> as.vector(), NA)) |>
  group_by(outcome) |>
  mutate(star = ifelse(scaledbic == min(scaledbic, na.rm=T), 1, 0))

b <- best |> select(outcome, model, star)

allcan <- left_join(allcan, b, by=c("outcome","model")) |>
  arrange(outcome, BIC)

library(NDRSAfunctions)
allcan |> 
  ggplot(aes(x=out_type, y=term)) +
  geom_tile(aes(fill = scaledbic), colour= "white") +
  scale_fill_gradientn(colours = get_ndrs_palette("NCRASBlues")(20), 
                       na.value = "white") +
  labs(x = "Metric", y = "Candidate", fill = "Scaled BIC",
       caption = "'Best' Indicates the significant (Wald p-value <0.05) candidate with
the lowest scaled Bayesian Information Criterion (BIC) per outcome.
BIC values are scaled within each cancer metric.
ASR age-standardised rate.") +
  geom_text(data = allcan[allcan$star == 1,], aes(label = "Best")) +
  theme_bw() 
# G_25_5 still best for incidence and prevalence

bestout <- best |>
  filter(p.value < 0.05) |>
  group_by(outcome) |>
  filter(scaledbic == min(scaledbic, na.rm=T))

inc <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = dat_long[dat_long$outcome == "P_inc_ASR_All",])

prev2 <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite + ah3hosp + ah3e + samhi_index, data = dat_long[dat_long$outcome == "prev_CR_All",])

library(jtools)
mod_outputs <- jtools::export_summs(inc, prev, scale = TRUE, to.file = "docx", file.name = paste0(Sys.getenv("coastal_clean"),"/ModelsTable2.docx"), error_format = "[{conf.low}, {conf.high}]")


summary(inc)

summary(prev)


#we want ALL LSOAs not only the ones used in modelling
## s3.1 LSOA11.......................................................................
cands <-read.csv(paste0(Sys.getenv("coastal_clean"), "/lsoa_candidates.csv"))

out <- cands |>
  select(LSOA11, G_25_5) |>
  rename(`Coastal (as per G_25_5)` = G_25_5)

write.csv(out, paste0(Sys.getenv("coastal_clean"), "/LSOA11_G255.csv"), row.names = F)


## s3.2 LSOA21 -------------------------------------------------------------------------
rm(list=ls())
pacman::p_load(tidyverse, sf, data.table, httr, jsonlite, broom, tidyr, purrr, 
               remotes, readxl, here, scales, assertr)
library(NDRSAfunctions)
source("functions.R")

gb <- read_sf(paste0(Sys.getenv("coastal_clean"),"/Countries_December_2021_GB_BUC_2022_3618480958259313072.geojson"))
#make the shapefile a simple line
gb_coast <- st_cast(st_union(gb),"MULTILINESTRING")
rm(gb)


#CASREF01 connection
x <- NDRSAfunctions::createConnection(username = Sys.getenv("analysis_username"))

#collect the postcodes and keep only the useful bits
pcd <- NDRSAfunctions::dbGetQueryOracle(x, "select * from ANALYSISNCR.NSPLC21_202405") |>
  select(PCD, PCD2, PCDS, DOINTR, DOTERM, LATITUDE, LONGITUDE, LSOA21, LAUA, RGN, 
         NHSER, CTRY, RU11IND) |>
  filter(CTRY == "E92000001", LATITUDE > 40 & LATITUDE < 70) |> 
  filter(DOTERM == ""|is.na(DOTERM))
#32844 LSOAs

# Set CRS to match shapefile (WGS)  
pcd <- pcd |>
  st_as_sf(coords = c("LONGITUDE", "LATITUDE")) |>
  st_set_crs(4326)

coast_pcds_1km <- st_is_within_distance(gb_coast, pcd, dist = 1000) # Is within 1km (1000 meters) of the UK boundary
coast_pcds_3km <- st_is_within_distance(gb_coast, pcd, dist = 3000) # Repeat but for 3km (3000 meters)
coast_pcds_5km <- st_is_within_distance(gb_coast, pcd, dist = 5000) # Repeat but for 5km (5000 meters)

pcd$coastal_1km <- 0 # Create blank variable (i.e., not coastal postcode)
pcd$coastal_1km[coast_pcds_1km[[1]]] <- 1 # Define postcode as being located by the coast
pcd$coastal_3km <- 0 # Repeat same process but for 3 km
pcd$coastal_3km[coast_pcds_3km[[1]]] <- 1
pcd$coastal_5km <- 0 # Repeat same process but for 5 km
pcd$coastal_5km[coast_pcds_5km[[1]]] <- 1


#save and tidy up
st_write(pcd, paste0(Sys.getenv("coastal_clean"),"/postcodes_coastal_definition21.shp"))

pcd <- st_read(paste0(Sys.getenv("coastal_clean"),"/postcodes_coastal_definition21.shp")) |>
  st_as_sf(coords = c("long", "lat")) 
pcd <- pcd |>
  st_drop_geometry()


lsoa_coasts <- read_sf(paste0(Sys.getenv("coastal_clean"),"/Lower_layer_Super_Output_Areas_2021_EW_BGC_V3_-6823567593069184824.geojson")) |>
  filter(stringr::str_detect(LSOA21CD, "^E")) 
sf_use_s2(FALSE)
int <- st_intersects(gb_coast, lsoa_coasts)

lsoa_coasts$intersects <-0
lsoa_coasts$intersects[int[[1]]] <- 1

lsoa_coasts <- lsoa_coasts |>
  select(LSOA21CD, intersects)

lsoas <- pcd |>
  st_drop_geometry() |>
  select(PCD, LSOA21, RGN, NHSER, cstl_1k, cstl_3k, cstl_5k) |>
  inner_join(lsoa_coasts, by=c("LSOA21" = "LSOA21CD")) |>
  group_by(LSOA21) |>
  mutate(coastline = ifelse(sum(intersects > 0,na.rm=T), 1,0), # does the LSOA have any coastline?
         allpcd = n(),
         G_25_5 = ifelse((sum(cstl_5k,na.rm=T)/max(allpcd,na.rm=T))>=0.25, 1, 0)) |>
  select(LSOA21, G_25_5) |>
  distinct() |>
  droplevels() |>
  ungroup()

write.csv(lsoas, file = paste0(Sys.getenv("coastal_clean"), "/LSOA21_G255.csv"))
