## COASTAL ANALYSIS

# This script takes the candidate coastal variables generated (and mapped) in 
#"1. Derive_variables.R" and runs various analyses to understand their 
#relationships with each other, other important covariates and with cancer 
#outcomes.


## 1. Combine all relevant data sources ----------------------------------------
options("install.lock" = FALSE)
pacman::p_load(tidyverse, janitor, readxl, Hmisc, modelr, broom, purrr, here,
               assert, jtools, officer, huxtable, table1)
library(NDRSAfunctions)
source("functions.R")

# Candidates per LSOA (created in '1.Derive_variables.R') ---
cands <-read.csv(paste0(Sys.getenv("coastal_clean"), "/lsoa_candidates.csv"))

#obtain persons-level cancer data from CAS
prev_p <- read.csv(paste0(Sys.getenv("coastal_clean"), 
                          "/inc_mort_prev_rates_by_lsoa_persons_FINAL.csv"), 
                   stringsAsFactors = F)
prev_s <- read.csv(paste0(Sys.getenv("coastal_clean"), 
                          "/inc_mort_prev_rates_by_lsoa_sex_specific_FINAL.csv"), 
                   stringsAsFactors = F)

# MEAN AGE 

# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/
#populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates
#LSOA-level counts of people per single year of age in 2016
age <- read_excel(paste0(Sys.getenv("coastal_clean"),
                         "/SAPE19DT1-mid-2016-lsoa-syoa-estimates.xls"), 
                  sheet = "Mid-2016 Persons", range="A5:CQ35106")

agelong <- age |>
  pivot_longer(`0`:`90+`, names_to = "age", values_to = "count") |>
  filter(!is.na(`...3`)) |>
  select(-c(`Area Names`)) |>
  mutate(age2 = gsub("\\+$","",age))  |>
  mutate(mult = (as.numeric(age2)*count)) |>
  select(-age) |>
  rename(age = age2)

#calculate the mean age (years) per LSOA
meanage<- agelong |>
  group_by(`Area Codes`) |>
  summarise(mean = sum(mult)/`All Ages`[1]) |>
  filter(grepl("^E",`Area Codes`))

# GENDER (% female) 
sex <- read_excel(paste0(Sys.getenv("coastal_clean"),
                         "/SAPE19DT1-mid-2016-lsoa-syoa-estimates.xls"), 
                  sheet = "Mid-2016 Females", range="A5:CQ35106") |>
  select(`Area Codes`, `...3`, `All Ages`) |>
  inner_join(age[,c("Area Codes","...3", "All Ages")], by="Area Codes") |>
  mutate(fem = as.numeric(`All Ages.x`),
         all = as.numeric(`All Ages.y`)) |>
  mutate(pcfem = fem/all) |>
  select(`Area Codes`, pcfem)

# ETHNICITY (% white) 


eth <- read_excel(paste0(Sys.getenv("coastal_clean"),"/Nomis 2016 eth by LA.xlsx")) |>
  slice(-c(1L:5L)) |>
  row_to_names(row_number=1) |>
  clean_names() |>
  select(!starts_with("na")) |>
  slice(-1L) |>
  mutate_at(-1, as.numeric) |> 
  mutate(pcwhite = `t18_4_all_ages_white_all_people` / `t18_1_all_ages_all_all_people`) 
#n=368, includes scotland and wales

eth <- eth |>
  select(1,ncol(eth))
names(eth)[1] <- "LA"

# add CHD prevalence and distance to GP
gp <- read.csv(paste0(Sys.getenv("coastal_clean"),"/GPP_3_04_GPP_network_distances_LSOA11.csv")) |>
  filter(year == "2016") |>
  select(lsoa11, mean_dist)

chd <- read.csv(paste0(Sys.getenv("coastal_clean"),"/QOF_CHD_LSOA.csv")) |>
  filter(year == "2016") |>
  select(lsoa11, p_rate) |>
  rename(chd_prev = p_rate)

#Use the CAS database to access LA to LSOA lookup
con <- createConnection(username = Sys.getenv("analysis_username"))

lu <- dbGetQueryOracle(con, "select * from ANALYSISNCR.LSOA_CCG_LA_STP_CA_202004") |>
  select(LSOA11CD, LAD20CD, LAD20NM)


#modify lu to match names in eth
# 2021 Corby merged with East Northamptonshire, Kettering, and Wellingborough to 
#make North Northamptonshire
# 2021 Daventry, Northampton and South Northamptonshire merged to become West 
#Northamptonshire

lu$LAD20NM[lu$LAD20NM %in% 
             c("Corby","East Northamptonshire","Kettering","Wellingborough")] <- 
  "North Northamptonshire"
lu$LAD20NM[lu$LAD20NM %in% 
             c("Daventry","South Northamptonshire","Northampton")] <- 
  "West Northamptonshire"

eth <- inner_join(eth, lu, by=c("LA" = "LAD20NM"))

# Add all the covriates into a single file
covs <- inner_join(meanage, sex, by="Area Codes")
covs <- inner_join(covs, eth, by=c("Area Codes" = "LSOA11CD")) |>
  select(-c(LA, LAD20CD))
names(covs)[c(1,2)] = c("lsoa","meanage")
covs <- inner_join(covs, gp, by=c("lsoa" = "lsoa11"))
covs <- inner_join(covs, chd, by=c("lsoa" ="lsoa11"))

#Cancer data 

#prev_p is the LSOA level nonsex-specific rates (all, bowel, lung)
#prev_s is sex-specific cases

#cases and ASRs of incidence, mortality, plus CR prevalence
#all persons for all bowel lung
#f for breast
#m for prostate

#keep the relevant rows
dat_p <- prev_p |>
  select(Geography, 
         P_inc_ASR_All, P_inc_ASR_Bowel, P_inc_ASR_Lung, P_inc_ASR_Melanoma,
         P_mort_ASR_All, P_mort_ASR_Bowel, P_mort_ASR_Lung, P_mort_ASR_Melanoma,
         P_inc_Cases_All, P_inc_Cases_Bowel, P_inc_Cases_Lung,P_inc_Cases_Melanoma,
         P_mort_Cases_All,P_mort_Cases_Bowel, P_mort_Cases_Lung, P_mort_Cases_Melanoma,
         prev_CR_All, prev_CR_Bowel, prev_CR_Lung,P_prev_CR_Melanoma,
         prev_Cases_All, prev_Cases_Bowel, prev_Cases_Lung, P_prev_Cases_Melanoma)
#sex-specific rates (prostate and breast)
dat_s <- prev_s |>
  select(Geography, 
         F_inc_Cases_Breast, M_inc_Cases_Prostate,
         F_inc_ASR_Breast, M_inc_ASR_Prostate,
         F_mort_Cases_Breast, M_mort_Cases_Prostate,
         F_mort_ASR_Breast, M_mort_ASR_Prostate,
         F_prev_Cases_Breast, M_prev_Cases_Prostate,
         F_prev_CR_Breast, M_prev_CR_Prostate)

prev <- inner_join(dat_p, dat_s, by="Geography")

cand_covs <- inner_join(cands, covs, by=c("LSOA11" = "lsoa")) 

dat <- inner_join(cand_covs, prev, by=c("LSOA11" = "Geography"))


# Case Number checks 


#check that there are enough cases of inc/prev/mort in each coastal definition to 
#use for analysis - within each region/LA

#cut off is 30 cases per area

dat <- inner_join(dat, lu, by=c("LSOA11" = "LSOA11CD"))

regioncounts <- dat |>
  group_by(RGN) |>
  summarise(counti_all = sum(`P_inc_Cases_All`, na.rm = T),
            counti_breast = sum(`F_inc_Cases_Breast`, na.rm = T),
            counti_prost = sum(`M_inc_Cases_Prostate`, na.rm = T),
            counti_colo = sum(`P_inc_Cases_Bowel`, na.rm = T),
            counti_lung = sum(`P_inc_Cases_Lung`, na.rm = T),
            counti_mel = sum(`P_inc_Cases_Melanoma`, na.rm = T),
            #counts of prevalent cases,
            countp_all = sum(`prev_Cases_All`, na.rm = T),
            countp_breast = sum(`F_prev_Cases_Breast`, na.rm = T),
            countp_prost = sum(`M_prev_Cases_Prostate`, na.rm = T),
            countp_colo = sum(`prev_Cases_Bowel`, na.rm = T),
            countp_lung = sum(`prev_Cases_Lung`, na.rm = T),
            countp_mel = sum(`P_prev_Cases_Melanoma`, na.rm = T),
            #count deaths,
            countd_all = sum(`P_mort_Cases_All`, na.rm = T),
            countd_breast = sum(`F_mort_Cases_Breast`, na.rm = T),
            countd_prost = sum(`M_mort_Cases_Prostate`, na.rm = T),
            countd_colo = sum(`P_mort_Cases_Bowel`, na.rm = T),
            countd_lung = sum(`P_mort_Cases_Lung`, na.rm = T),
            countd_mel = sum(`P_mort_Cases_Melanoma`, na.rm = T))

la<-dat |>
  group_by(LAD20NM) |>
  summarise(counti_all = sum(`P_inc_Cases_All`, na.rm = T),
            counti_breast = sum(`F_inc_Cases_Breast`, na.rm = T),
            counti_prost = sum(`M_inc_Cases_Prostate`, na.rm = T),
            counti_colo = sum(`P_inc_Cases_Bowel`, na.rm = T),
            counti_lung = sum(`P_inc_Cases_Lung`, na.rm = T),
            counti_mel = sum(`P_inc_Cases_Melanoma`, na.rm = T),
            #counts of prevalent cases,
            countp_all = sum(`prev_Cases_All`, na.rm = T),
            countp_breast = sum(`F_prev_Cases_Breast`, na.rm = T),
            countp_prost = sum(`M_prev_Cases_Prostate`, na.rm = T),
            countp_colo = sum(`prev_Cases_Bowel`, na.rm = T),
            countp_lung = sum(`prev_Cases_Lung`, na.rm = T),
            countp_mel = sum(`P_prev_Cases_Melanoma`, na.rm = T),
            #count deaths,
            countd_all = sum(`P_mort_Cases_All`, na.rm = T),
            countd_breast = sum(`F_mort_Cases_Breast`, na.rm = T),
            countd_prost = sum(`M_mort_Cases_Prostate`, na.rm = T),
            countd_colo = sum(`P_mort_Cases_Bowel`, na.rm = T),
            countd_lung = sum(`P_mort_Cases_Lung`, na.rm = T),
            countd_mel = sum(`P_mort_Cases_Melanoma`, na.rm = T))

#minimum values
reg_long <-regioncounts |>
  pivot_longer(-c(RGN), names_to = "Region", values_to = "count")
min(reg_long$count)
#499

la_Long <- la |>
  pivot_longer(-c(LAD20NM), names_to = "LA", values_to = "count")
min(la_Long$count)
#1
#City of London and Isles of Scilly obviously small


# Remove Local Authorities with too few cases per LA 
dat <- dat |>
  filter(!LAD20NM %in% c("Isles of Scilly","City of London"))
n_distinct(dat$LSOA11)
#32837

#check it's one row per LSOA
dat <- dat |>
  group_by(LSOA11) |>
  mutate(n=n())

write.csv(dat, paste0(Sys.getenv("coastal_clean"),"/modeldata.csv"))

# 2. MODELLING-------------------------------------------------------------------
rm(list=ls())
pacman::p_load(tidyverse, janitor, readxl, Hmisc, modelr, broom, purrr, here)
dat <- read.csv(paste0(Sys.getenv("coastal_clean"),"/modeldata.csv"))
source("functions.R")

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
covars <- c("meanage", "pcfem", "pcwhite", "quint", "rural", "mean_dist")

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
  select(-c(X.1, X))
idcols <- names(dat)[names(dat) %in% c("LSOA11","RGN","rural","quint","coastline","allpcd","LAD20CD", "LAD20NM","n","YEAR", "total")]
counts <- names(dat)[grepl("Cases", names(dat))]
outcome_cols <- names(dat)[!names(dat) %in% c(counts, covars, candidates, idcols)]

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

#investigate outliers 

# summary(dat_long$value[dat_long$outcome == "P_inc_ASR_Bowel"])
# boxplot(dat_long$value[dat_long$outcome == "P_inc_ASR_Bowel"])
# hist(dat_long$value[dat_long$outcome == "P_inc_ASR_Bowel"])
# 
# bit <- dat_long |> filter(outcome == "P_inc_ASR_Bowel", value >2000)
# 
# boxplot(dat_long$value[dat_long$outcome == "P_mort_ASR_Bowel"])
# bit <- dat_long |> filter(outcome == "P_mort_ASR_Bowel", value >2000)
# 
# #remove this grenwich LSOA from all analyses
# boxplot(dat_long$value[dat_long$outcome == "M_mort_ASR_Prostate"])
# bit <- dat_long |> filter(outcome == "M_mort_ASR_Prostate", value >1500)
# 
# boxplot(dat_long$value[dat_long$outcome == "P_mort_ASR_All"])
# bit <- dat_long |> filter(outcome == "P_mort_ASR_All", value > 2500) 
# 

# Exclude LSOA outliers 
dat_long <- dat_long |>
  filter(!LSOA11 %in% c("E01033739","E01032703"))

dat <- dat |>
  filter(!LSOA11 %in% c("E01033739","E01032703"))

saveRDS(dat_long, file = paste0(Sys.getenv("coastal_clean"), "/dat_long.RDS"))
## 3. CORRELATIONS --------------------------------------
# Figure 1 ----

cors <- data.frame(var = c(NA, covars))
for(cand in candidates){
  data <- dat[,c(cand, covars)] 
  data[[cand]][data[[cand]] == 2] <- 0 #set the 'comparator' to 0 so all are binary
  out <- bind_cols(rcorr(as.matrix(data), type = "spearman")$r[1,], 
                   rcorr(as.matrix(data), type = "spearman")$P[1,]) |>
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
write.csv(corsodd, paste0(Sys.getenv("coastal_clean"), "/Correlations_no_outliers.csv"))

#find large correlations
longcors <- corsodd |>
  pivot_longer(cols = -c("var...1"), names_to = "cand",values_to = "Pearson") |>
  filter(!is.na(Pearson))

cor_heatmap <- longcors |>
  filter(`var...1` != "mean_dist") |>
  droplevels() |>
  mutate(cov = case_when(`var...1` == "meanage" ~ "Mean Age (years)",
                         `var...1` == "pcfem" ~ "Proportion Female",
                         `var...1` == "pcwhite" ~ "Proportion White Ethnicity",
                         `var...1` == "rural" ~ "Rural RUC11 code",
                         `var...1` == "total" ~ "Population in 2016",
                         `var...1` == "mean_dist" ~ "Mean distance to GP",
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
ggsave(paste0(Sys.getenv("coastal_clean"),"/Corr_heatmap_no_outliers.png"), 
       dpi = 600, width = 210, height = 297, units = "mm")

# Define models 

by_out <- dat_long |>
  group_by(outcome) |>
  nest()

g1_model <- function(df){
  lm(value ~ G1 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_10_1_model <- function(df){
  lm(value ~ G_10_1 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_10_3_model <- function(df){
  lm(value ~ G_10_3 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_10_5_model <- function(df){
  lm(value ~ G_10_5 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_25_1_model <- function(df){
  lm(value ~ G_25_1 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_25_3_model <- function(df){
  lm(value ~ G_25_3 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_25_5_model <- function(df){
  lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_50_1_model <- function(df){
  lm(value ~ G_50_1 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_50_3_model <- function(df){
  lm(value ~ G_50_3 + quint + rural + meanage + pcfem + pcwhite , data = df)
}
g_50_5_model <- function(df){
  lm(value ~ G_50_5 + quint + rural + meanage + pcfem + pcwhite , data = df)
}


gd1_model <- function(df){
  lm(value ~ GD1 +rural + meanage + pcfem + pcwhite , data = df)
}
gd_10_1_model <- function(df){
  lm(value ~ GD_10_1 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_10_3_model <- function(df){
  lm(value ~ GD_10_3 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_10_5_model <- function(df){
  lm(value ~ GD_10_5 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_25_1_model <- function(df){
  lm(value ~ GD_25_1 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_25_3_model <- function(df){
  lm(value ~ GD_25_3 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_25_5_model <- function(df){
  lm(value ~ GD_25_5 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_50_1_model <- function(df){
  lm(value ~ GD_50_1 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_50_3_model <- function(df){
  lm(value ~ GD_50_3 + rural + meanage + pcfem + pcwhite , data = df)
}
gd_50_5_model <- function(df){
  lm(value ~ GD_50_5 + rural + meanage + pcfem + pcwhite , data = df)
}



gr1_model <- function(df){
  lm(value ~ GR1 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_10_1_model <- function(df){
  lm(value ~ GR_10_1 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_10_3_model <- function(df){
  lm(value ~ GR_10_3 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_10_5_model <- function(df){
  lm(value ~ GR_10_5 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_25_1_model <- function(df){
  lm(value ~ GR_25_1 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_25_3_model <- function(df){
  lm(value ~ GR_25_3 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_25_5_model <- function(df){
  lm(value ~ GR_25_5 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_50_1_model <- function(df){
  lm(value ~ GR_50_1 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_50_3_model <- function(df){
  lm(value ~ GR_50_3 + quint + meanage + pcfem + pcwhite , data = df)
}
gr_50_5_model <- function(df){
  lm(value ~ GR_50_5 + quint + meanage + pcfem + pcwhite , data = df)
}

gdr1_model <- function(df){
  lm(value ~ GDR1 + meanage + pcfem + pcwhite , data = df)
}
gdr_10_1_model <- function(df){
  lm(value ~ GDR_10_1 + meanage + pcfem + pcwhite , data = df)
}
gdr_10_3_model <- function(df){
  lm(value ~ GDR_10_3 + meanage + pcfem + pcwhite , data = df)
}
gdr_10_5_model <- function(df){
  lm(value ~ GDR_10_5 + meanage + pcfem + pcwhite , data = df)
}
gdr_25_1_model <- function(df){
  lm(value ~ GDR_25_1 + meanage + pcfem + pcwhite , data = df)
}
gdr_25_3_model <- function(df){
  lm(value ~ GDR_25_3 + meanage + pcfem + pcwhite , data = df)
}
gdr_25_5_model <- function(df){
  lm(value ~ GDR_25_5 + meanage + pcfem + pcwhite , data = df)
}
gdr_50_1_model <- function(df){
  lm(value ~ GDR_50_1 + meanage + pcfem + pcwhite , data = df)
}
gdr_50_3_model <- function(df){
  lm(value ~ GDR_50_3 + meanage + pcfem + pcwhite , data = df)
}
gdr_50_5_model <- function(df){
  lm(value ~ GDR_50_5 + meanage + pcfem + pcwhite , data = df)
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

# obtain residuals and predictions for all models
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

#plot the residuals plots for each model/outcome
plots <- map(modelsets, plot_fun, by_out)

#get the overall model info for every model/outcome
glances <- bind_rows(map(modelsets, glance_fun, by_out))

#obtain model outputs
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

#select the results for the candidate of each model and combine
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

all <- inner_join(glances, 
                  ests[ests$varpart == "p.value"|ests$varpart == "estimate",], 
                  by=c("outcome","model"))

saveRDS(all, paste0(Sys.getenv("coastal_clean"),"/all.RDS"))

# 4. Model outputs -----------------------------------------------------------

# is every candidate a significant addition to every model?
all$sig <- ifelse(all$varpart == "p.value" & all$value < 0.05, 1,0)

allcan <- all |> 
  filter(cancer == "All", varpart == "p.value") |> 
  mutate(sig = ifelse(value < 0.05, 1,0)) |>
  group_by(out_type) |>
  mutate(scaledbic = ifelse(sig ==1, scale(BIC) |>as.vector(), NA))
#is G_25_5 significant in all outcomes?
g255<- all |>
  filter(model == "g_25_5",
         varpart == "p.value")

#does every candidate have at least one significant model?
mods <- allcan |> group_by(model, out_type, sig) |> summarise( n=n())
table(mods$n[mods$sig == 1], useNA = "always")

#best model per cancer based on BIC
#combine BIC with candidate p-value info
best <- all |>
  select(-sig) |>
  pivot_wider(names_from = varpart, values_from = value) |>
  mutate(sig = ifelse(p.value < 0.05, 1,0)) |>
  group_by(out_type) |>
  mutate(scaledbic = ifelse(sig ==1, scale(BIC) |> as.vector(), NA)) |>
  group_by(outcome) |>
  mutate(star = ifelse(scaledbic == min(scaledbic, na.rm=T), 1, 0))

b <- best |> select(outcome, model, star)

##redo this bit so that each is programatically chosen as the 'best'
## Figure 2.----------------------
allcan <- left_join(allcan, b, by=c("outcome","model")) |>
  arrange(outcome, BIC)


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
ggsave(paste0(Sys.getenv("coastal_clean"), "/Allcancer_models_tile_corrected.png"), 
       dpi = 600, width = 210, height = 297, units = "mm")

# same but for CHD
chd_best <- all |> 
  filter(cancer == "CHD", varpart == "p.value") |> 
  mutate(sig = ifelse(value < 0.05, 1,0)) |>
  group_by(out_type) |>
  mutate(scaledbic = ifelse(sig ==1, scale(BIC) |>as.vector(), NA)) 
chd_best <- left_join(chd_best, b, by=c("outcome","model")) |>
  arrange(outcome, BIC)

chd_best |> 
  ggplot(aes(x=out_type, y=term)) +
  geom_tile(aes(fill = scaledbic), colour= "white") +
  scale_fill_gradientn(colours = get_ndrs_palette("NCRASBlues")(20), 
                       na.value = "white") +
  labs(x = "Coronary Heart Disease Prevalence", y = "Candidate", fill = "Scaled BIC",
       caption = "'Best' Indicates the significant (Wald p-value <0.05) candidate with
the lowest scaled Bayesian Information Criterion (BIC).") +
  geom_text(data = chd_best[chd_best$star == 1,], aes(label = "Best")) +
  theme_bw() 
ggsave(paste0(Sys.getenv("coastal_clean"), "/CHD_models_tile_corrected.png"), 
       dpi = 600, width = 210, height = 297, units = "mm")


#plot all the 'bests'
best |>
  filter(p.value < 0.05, outcome != "chd_prev") |>
  group_by(outcome) |>
  filter(scaledbic == min(scaledbic, na.rm=T)) |>
  ggplot(aes(x=out_type, y=term)) +
  geom_tile(aes(fill = adj.r.squared), colour= "white") +
  scale_fill_distiller(palette = "RdPu") +
  facet_wrap(~ cancer) + 
  theme_light() +
  labs(title = "Selected candidate per model",
       y="", x = "Outcome Type", fill = "Adjusted R-squared",
       caption = "Candidates with Wald P-values >0.05 are omitted") +
  theme(axis.text.x = element_text(angle=90, vjust =0.5))
ggsave(paste0(Sys.getenv("coastal_clean"), "/Best_models.png"))


bestout <- best |>
  filter(p.value < 0.05) |>
  group_by(outcome) |>
  filter(scaledbic == min(scaledbic, na.rm=T))
write.csv(bestout, paste0(Sys.getenv("coastal_clean"), "/Best_models.csv"))

#check whether G_25_5 is significant for all cancers & CHD
g255 <- best |>
  filter(model == "g_25_5")
# G_25_5 is also significant for CHD


# All Cancer: inc/prev = G_25_5, mort = GR_10_5

all |>
  filter(term == "G_25_5") |>
  mutate(mapcol = case_when(value < 0.05 ~ 1,
                            .default = 0.8)) |>
  ggplot(aes(x=cancer, y=adj.r.squared, group = out_type, fill = out_type,
             alpha = mapcol)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Adjusted R-squared", x = "", fill = "Metric type",
       title = "Adjusted R-squared of the G_25_5 model",
       caption = "Adjusted for candidate, mean age, proportion female, proportion white ethnicity, rurality and IMD quintile.\n
       Wald significance of the candidate variable is denoted by opacity.") +
  theme_bw() +
  guides(alpha = FALSE)
ggsave(paste0(Sys.getenv("coastal_clean"), "/Candidate_pvalues.png"))

#check a lrtest with All mortality to use G_25_5 anyway
with <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite , 
           data = dat_long[dat_long$outcome == "P_mort_ASR_All",])
without <- lm(value ~  quint + rural + meanage + pcfem + pcwhite , 
              data = dat_long[dat_long$outcome == "P_mort_ASR_All",])
library(lmtest)
lrtest(without, with)
#not sig
withgr <- lm(value ~  G_25_5 + quint + meanage + pcfem + pcwhite , 
             data = dat_long[dat_long$outcome == "P_mort_ASR_All",])
nogr <- lm(value ~  rural + quint + meanage + pcfem + pcwhite , 
           data = dat_long[dat_long$outcome == "P_mort_ASR_All",])
both <- lm(value ~  G_25_5 + rural + quint + meanage + pcfem + pcwhite , 
           data = dat_long[dat_long$outcome == "P_mort_ASR_All",])
summary(both)
BIC(withgr)
BIC(nogr)
# so G_25_5 is of no benefit in mortality of all cancer models


# Repeat LA counts ----
# Now count cases per region and per LA based on G_10_5
regioncounts <- dat |>
  filter(G_25_5 == 1)|>
  group_by(RGN) |>
  summarise(counti_all = sum(`P_inc_Cases_All`, na.rm = T),
            counti_breast = sum(`F_inc_Cases_Breast`, na.rm = T),
            counti_prost = sum(`M_inc_Cases_Prostate`, na.rm = T),
            counti_colo = sum(`P_inc_Cases_Bowel`, na.rm = T),
            counti_lung = sum(`P_inc_Cases_Lung`, na.rm = T),
            counti_mel = sum(`P_inc_Cases_Melanoma`, na.rm = T),
            #counts of prevalent cases,
            countp_all = sum(`prev_Cases_All`, na.rm = T),
            countp_breast = sum(`F_prev_Cases_Breast`, na.rm = T),
            countp_prost = sum(`M_prev_Cases_Prostate`, na.rm = T),
            countp_colo = sum(`prev_Cases_Bowel`, na.rm = T),
            countp_lung = sum(`prev_Cases_Lung`, na.rm = T),
            countp_mel = sum(`P_prev_Cases_Melanoma`, na.rm = T),
            #count deaths,
            countd_all = sum(`P_mort_Cases_All`, na.rm = T),
            countd_breast = sum(`F_mort_Cases_Breast`, na.rm = T),
            countd_prost = sum(`M_mort_Cases_Prostate`, na.rm = T),
            countd_colo = sum(`P_mort_Cases_Bowel`, na.rm = T),
            countd_lung = sum(`P_mort_Cases_Lung`, na.rm = T),
            countd_mel = sum(`P_mort_Cases_Melanoma`, na.rm = T))

la<-dat |>
  filter(G_25_5 == 1)|>
  group_by(LAD20NM) |>
  summarise(counti_all = sum(`P_inc_Cases_All`, na.rm = T),
            counti_breast = sum(`F_inc_Cases_Breast`, na.rm = T),
            counti_prost = sum(`M_inc_Cases_Prostate`, na.rm = T),
            counti_colo = sum(`P_inc_Cases_Bowel`, na.rm = T),
            counti_lung = sum(`P_inc_Cases_Lung`, na.rm = T),
            counti_mel = sum(`P_inc_Cases_Melanoma`, na.rm = T),
            #counts of prevalent cases,
            countp_all = sum(`prev_Cases_All`, na.rm = T),
            countp_breast = sum(`F_prev_Cases_Breast`, na.rm = T),
            countp_prost = sum(`M_prev_Cases_Prostate`, na.rm = T),
            countp_colo = sum(`prev_Cases_Bowel`, na.rm = T),
            countp_lung = sum(`prev_Cases_Lung`, na.rm = T),
            countp_mel = sum(`P_prev_Cases_Melanoma`, na.rm = T),
            #count deaths,
            countd_all = sum(`P_mort_Cases_All`, na.rm = T),
            countd_breast = sum(`F_mort_Cases_Breast`, na.rm = T),
            countd_prost = sum(`M_mort_Cases_Prostate`, na.rm = T),
            countd_colo = sum(`P_mort_Cases_Bowel`, na.rm = T),
            countd_lung = sum(`P_mort_Cases_Lung`, na.rm = T),
            countd_mel = sum(`P_mort_Cases_Melanoma`, na.rm = T))

#minimum values
reg_long <-regioncounts |>
  pivot_longer(-c(RGN), names_to = "Region", values_to = "count")
min(reg_long$count[!grepl("mel", reg_long$Region)])


la_long <- la |>
  pivot_longer(-c(LAD20NM), names_to = "LA", values_to = "count")

min(la_long$count[!grepl("mel", la_long$LA)])

la_long <- la_long |>
  mutate(cancer = case_when(grepl("breast", LA) ~ "Breast",
                            grepl("lung", LA) ~ "Lung",
                            grepl("colo", LA) ~ "Bowel",
                            grepl("prost", LA) ~ "Prostate",
                            grepl("mel", LA) ~ "Melanoma",
                            grepl("chd", LA) ~ "CHD",
                            .default = "All"),
         metric = case_when(grepl("i_",LA) ~ "Incidence",
                            grepl("d_",LA) ~ "Mortality",
                            grepl("p_|_prev", LA) ~ "Prevalence"))



metrics <- c("Incidence", "Prevalence","Mortality")


G255 <- data.frame()

#For g_25_5, if summed to LA level, how many LA would have < 10 cases?
for (m in metrics) {
  xs <- la_long |>
    filter(metric == m, count < 10) |>
    group_by(cancer) |>
    summarise(nlad = n_distinct(LAD20NM)) |>
    mutate(proplad = round((nlad/114)*100, 2),
           metric = m)
  G255 <- bind_rows(G255, xs)
}        


# G_25_5 is the best model (inc/prev)
# model output
inc <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite, 
          data = dat_long[dat_long$outcome == "P_inc_ASR_All",])

mort <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite, 
           data = dat_long[dat_long$outcome == "P_mort_ASR_All",])

prev <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite, 
           data = dat_long[dat_long$outcome == "prev_CR_All",])

chd <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite, 
          data = dat_long[dat_long$outcome == "chd_prev",])


mod_outputs <- jtools::export_summs(inc, prev, mort, chd, scale = TRUE, 
                                    to.file = "docx", 
                                    file.name = paste0(Sys.getenv("coastal_clean"),"/ModelsTable.docx"), error_format = "[{conf.low}, {conf.high}]")

# Descriptive stats Table 1-----
# describe cancer stats and covs by the most 'inclusive' geog only candidate


dat$rural <- factor(dat$rural, levels = 0:1, 
                    labels = c("Not rural (RUC 2011 codes A-C)", 
                               "Rural (RUC 2011 codes D, E or F)"))
dat$G_25_5 <- factor(dat$G_25_5, levels = 0:1, 
                     labels = c("Not coastal", "Coastal"))
dat$GR_10_5 <- factor(dat$GR_10_5, levels = 0:1, 
                      labels = c("Not coastal", "Coastal"))

label(dat$rural) <- "Rural"
label(dat$quint) <- "Index of Multiple Deprivation Quintile (1 = most deprived)"
label(dat$meanage) <- "Mean population age in 2016 (years)"
label(dat$pcfem) <- "Proportion of female residents in 2016"
label(dat$pcwhite) <- "Proportion of white ethnicity residents in 2016"
label(dat$P_inc_ASR_All) <- "Age standardised incidence rate of all cancers"
label(dat$P_inc_ASR_Bowel) <- "Age standardised incidence rate of colorectal cancer"
label(dat$P_inc_ASR_Lung) <- "Age standardised incidence rate of lung cancer"
label(dat$F_inc_ASR_Breast) <- "Age standardised incidence rate of breast cancer in women"
label(dat$M_inc_ASR_Prostate) <- "Age standardised incidence rate of prostate cancer in men"
label(dat$P_mort_ASR_All) <- "Age standardised mortality rate of all cancers"
label(dat$P_mort_ASR_Bowel) <- "Age standardised mortality rate of colorectal cancer"
label(dat$P_mort_ASR_Lung) <- "Age standardised mortality rate of lung cancer"
label(dat$F_mort_ASR_Breast) <- "Age standardised mortality rate of breast cancer in women"
label(dat$M_mort_ASR_Prostate) <- "Age standardised mortality rate of prostate cancer in men"
label(dat$prev_CR_All) <- "Crude prevalence of all cancers"
label(dat$prev_CR_Bowel) <- "Crude prevalence of colorectal cancer"
label(dat$prev_CR_Lung) <- "Crude prevalence of lung cancer"
label(dat$F_prev_CR_Breast) <- "Crude prevalence of breast cancer in women"
label(dat$M_prev_CR_Prostate) <- "Crude prevalence of prostate cancer in men"

table1( ~ quint + meanage + rural + pcfem + pcwhite + 
          P_inc_ASR_All + P_inc_ASR_Bowel + P_inc_ASR_Lung + F_inc_ASR_Breast + M_inc_ASR_Prostate +
          P_mort_ASR_All + P_mort_ASR_Bowel + P_mort_ASR_Lung + F_mort_ASR_Breast +M_mort_ASR_Prostate +
          prev_CR_All + prev_CR_Bowel + prev_CR_Lung + F_prev_CR_Breast + M_prev_CR_Prostate | G_25_5, data = dat)
# Table 1 with only ALL cancers
mt1 <- table1( ~ quint + meanage + rural + pcfem + pcwhite + 
                 P_inc_ASR_All + 
                 P_mort_ASR_All + 
                 prev_CR_All | G_25_5, data = dat)
write.table (mt1 , paste0(Sys.getenv("coastal_clean"),"/Table1_All.csv"), 
             col.names = T, row.names=F, append= T, sep=',')

#gr_10_5 table1 of mortality
mt2 <- table1( ~ quint + meanage + rural + pcfem + pcwhite + 
                 P_inc_ASR_All + 
                 P_mort_ASR_All + 
                 prev_CR_All | GR_10_5, data = dat)

## table of count of coastal LSOAs per model
dat <- read.csv(paste0(Sys.getenv("coastal_clean"),"/modeldata.csv"))
definitions <- data.frame(var = c("G1", "G_10_1", "G_10_3", "G_10_5", 
                                  "GD1", "GD_10_1", "GD_10_3", "GD_10_5", 
                                  "GR1", "GR_10_1", "GR_10_3", "GR_10_5", 
                                  "GDR1", "GDR_10_1", "GDR_10_3", "GDR_10_5",
                                  "G_25_1", "G_25_3", "G_25_5",
                                  "GD_25_1", "GD_25_3", "GD_25_5",
                                  "GR_25_1", "GR_25_3", "GR_25_5",
                                  "GDR_25_1", "GDR_25_3", "GDR_25_5",
                                  "G_50_1", "G_50_3", "G_50_5",
                                  "GD_50_1", "GD_50_3", "GD_50_5",
                                  "GR_50_1", "GR_50_3", "GR_50_5",
                                  "GDR_50_1", "GDR_50_3", "GDR_50_5"),
                          definition = c("Intersects coastline only",
                                         "10 or more postcodes within 1km of coast",
                                         "10 or more postcodes within 3km of coast",
                                         "10 or more postcodes within 5km of coast",
                                         "Coastline and IMD quintile 1 or 2",
                                         "10 or more postcodes within 1km of coast and IMD quintile 1 or 2",
                                         "10 or more postcodes within 3km of coast and IMD quintile 1 or 2",
                                         "10 or more postcodes within 5km of coast and IMD quintile 1 or 2",
                                         "Coastline and rural (RUC11 codes D/E/F)",
                                         "10 or more postcodes within 1km of coast and rural (RUC11 codes D/E/F)",
                                         "10 or more postcodes within 3km of coast and rural (RUC11 codes D/E/F)",
                                         "10 or more postcodes within 5km of coast and rural (RUC11 codes D/E/F)",
                                         "Coastline and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "10 or more postcodes within 1km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "10 or more postcodes within 3km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "10 or more postcodes within 5km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "25% of postcodes within 1km of coast",
                                         "25% of postcodes within 3km of coast",
                                         "25% of postcodes within 5km of coast",
                                         "25% of postcodes within 1km of coast and IMD quintile 1 or 2",
                                         "25% of postcodes within 3km of coast and IMD quintile 1 or 2",
                                         "25% of postcodes within 5km of coast and IMD quintile 1 or 2",
                                         "25% of postcodes within 1km of coast and rural (RUC11 codes D/E/F)",
                                         "25% of postcodes within 3km of coast and rural (RUC11 codes D/E/F)",
                                         "25% of postcodes within 5km of coast and rural (RUC11 codes D/E/F)",
                                         "25% of postcodes within 1km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "25% of postcodes within 3km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "25% of postcodes within 5km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "50% of postcodes within 1km of coast",
                                         "50% of postcodes within 3km of coast",
                                         "50% of postcodes within 5km of coast",
                                         "50% of postcodes within 1km of coast and IMD quintile 1 or 2",
                                         "50% of postcodes within 3km of coast and IMD quintile 1 or 2",
                                         "50% of postcodes within 5km of coast and IMD quintile 1 or 2",
                                         "50% of postcodes within 1km of coast and rural (RUC11 codes D/E/F)",
                                         "50% of postcodes within 3km of coast and rural (RUC11 codes D/E/F)",
                                         "50% of postcodes within 5km of coast and rural (RUC11 codes D/E/F)",
                                         "50% of postcodes within 1km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "50% of postcodes within 3km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)",
                                         "50% of postcodes within 5km of coast\n and IMD quintile 1 or 2 and rural (RUC11 codes D/E/F)")
)

lsoas_per_cand <- dat |>
  select(`G1`:`GDR_25_5`) |>
  summarise(across(`G1`:`GDR_25_5`, sum)) |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = "var") 

totalpop <- sum(dat$total)

pop_per_cand <- dat |>
  select(G1:GDR_25_5, total) |>
  pivot_longer(cols = c(G1:GDR_25_5), names_to = "cand", values_to = "yn") |>
  group_by(cand) |>
  summarise(tpop = sum(total[yn == 1])) |>
  mutate(prop = round(tpop/totalpop*100, 1))


tab <- inner_join(lsoas_per_cand, definitions, by="var") |>
  rename(Candidate = var, `Count of LSOAs` = V1, Definition = definition) 

tab <- inner_join(tab, pop_per_cand, by=c("Candidate" = "cand")) |>
  select(-tpop) |>
  rename("Proportion of England Population" = prop) |>
  select(1,2,4,3)
write.csv(tab, paste0(Sys.getenv("coastal_clean"),"/TableS21.csv"))

#plot scatter of r2 versus BIC for all models, facet_wrapped on metric
pd <- all |> filter(cancer == "All", sex == "Persons",
                    varpart == "p.value", value < 0.05)
library(ggrepel)
pd |>
  ggplot(aes(x = adj.r.squared, y = BIC, colour = term )) +
  geom_point() +
  facet_wrap(~ out_type) +
  theme_bw() +
  geom_label_repel(data = pd[pd$term== "G_25_5",],
                   aes(label = term))

#what % of excess variation in inc/prev does G_25_5 explain? R2
summary(lm(value ~ quint + rural + meanage + pcfem + pcwhite , data = dat_long[dat_long$outcome == "P_inc_ASR_All",]))

summary(lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite , data = dat_long[dat_long$outcome == "P_inc_ASR_All",]))


#is chd model lrt improved bu G_25_5
x <- lm(value ~ quint + rural + meanage + pcfem + pcwhite , data = dat_long[dat_long$outcome == "chd_prev",])
y <- lm(value ~ G_25_5 + quint + rural + meanage + pcfem + pcwhite , data = dat_long[dat_long$outcome == "chd_prev",])
library(lmtest)
lrtest(x, y)
