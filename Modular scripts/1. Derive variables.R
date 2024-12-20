## Coastal variable derivation

# AIMS: what version of 'coastal' is (most) associated with poor cancer outcomes at LSOA level?
# compared to inland (a.k.a. 'not coastal')
#1. geographic? proximity to coast 
#2. deprived? geog + low IMD, vs geog+not low IMD 
#3. rural? geog+ rural vs geog+not rural 
#4. geog + rural+ imd 

#all 'geog' variables are to have a version where its not just 'has some coastline'
#but is also defined by certain % population living within x km of coast,
#and a certain count of population living within x km of coast.


# 1. LOAD GEOGRAPHIC DATA ------------------------------------------------------------

#install and/or load all relevant packages 

pacman::p_load(tidyverse, sf, data.table, httr, jsonlite, broom, tidyr, purrr, 
               remotes, readxl, here, scales, assertr)

source("functions.R")
# Load country boundaries for Great Britain from ONS API. 
# Ultra-generalised to 500m - this should be ok for our distance metrics 

gb <- read_sf(paste0(Sys.getenv("coastal_clean"),"/Countries_December_2021_GB_BUC_2022_3618480958259313072.geojson"))

#make the shapefile a simple line
gb_coast <- st_cast(st_union(gb),"MULTILINESTRING")
rm(gb)

# Create a connection to CAS to get all England postcodes and their locations
library(NDRSAfunctions)

#CASREF01 connection
x <- createConnection(username = Sys.getenv("analysis_username"))

#collect the postcodes and keep only the useful bits
pcd <- dbGetQueryOracle(x, "select * from ANALYSISNCR.NSPL_202205") |>
  select(PCD, PCD2, PCDS, DOINTR, DOTERM, LATITUDE, LONGITUDE, LSOA11, LAUA, 
         CALNCV, RGN, NHSER, CTRY, RU11IND, IMD) |>
  filter(CTRY == "E92000001", LATITUDE > 40 & LATITUDE < 70) |> 
  filter(DOTERM == ""|is.na(DOTERM))

# Set CRS to match shapefile (WGS)  
pcd <- pcd |>
  st_as_sf(coords = c("LONGITUDE", "LATITUDE")) |>
  st_set_crs(4326)         

#identify whether each individual postcode is within x km of coastline
#TAKES AGES so saved below

coast_pcds_1km <- st_is_within_distance(gb_coast, pcd, dist = 1000) # Is within 
#1km (1000 meters) of the UK boundary
coast_pcds_3km <- st_is_within_distance(gb_coast, pcd, dist = 3000) 
coast_pcds_5km <- st_is_within_distance(gb_coast, pcd, dist = 5000) 

pcd$coastal_1km <- 0 # Create blank variable (i.e., not coastal postcode)
pcd$coastal_1km[coast_pcds_1km[[1]]] <- 1 
pcd$coastal_3km <- 0 # Repeat same process but for 3 km
pcd$coastal_3km[coast_pcds_3km[[1]]] <- 1
pcd$coastal_5km <- 0 # Repeat same process but for 5 km
pcd$coastal_5km[coast_pcds_5km[[1]]] <- 1


#save and tidy up

st_write(pcd, paste0(Sys.getenv("coastal_clean"),"/postcodes_coastal_definition.shp"))

# coast_pcds_1km <- read.csv(paste0(Sys.getenv("spatial"),"/coast_pcds_1km.csv"))
# coast_pcds_3km <- read.csv(paste0(Sys.getenv("spatial"),"/coast_pcds_3km.csv"))
# coast_pcds_5km <- read.csv(paste0(Sys.getenv("spatial"),"/coast_pcds_5km.csv"))

pcd <- st_read(paste0(Sys.getenv("coastal_clean"),"/postcodes_coastal_definition.shp")) |>
  st_as_sf(coords = c("long", "lat"))

pcd <- pcd |> st_drop_geometry() |>
  select(-IMD)
saveRDS(pcd, file = paste0(Sys.getenv("coastal_clean"), "/pcd.RDS"))

#define LSOAs that have any coastline 
lsoa_coasts <- 
  read_sf(paste0(Sys.getenv("coastal_clean"),"/LSOA_(Dec_2011)_Boundaries_Generalised_Clipped_(BGC)_EW_V3.geojson")) |>
  filter(stringr::str_detect(LSOA11CD, "^E")) 
int <- st_intersects(gb_coast, lsoa_coasts)

# for each LSOA, does it intersect the England coastline
lsoa_coasts$intersects <-0
lsoa_coasts$intersects[int[[1]]] <- 1

lsoa_coasts <- lsoa_coasts |>
  select(LSOA11CD, intersects)

# 2. DEFINE CANDIDATES --------------------------------------------------------------

#a PCD is coastal if it is within x km of coastline/estuary

#an LSOA (based only on geography) is coastal if
#1. has any coastline/estuary
#2. has >=10 coastal postcodes
#3. has >=50% coastal postcodes
#4. has >=25% of coastal postcodes

#We derive 40 candidate variables in total, split into: 
#just geography (n=7)
#any coastline "G"
#>10 postcodes within x km: 1,3,5 ("G_10_x")
#>=50% of postcodes: 1,3,5 ("G_50_x")
#>=25% of poscodes: 1,3,5 ("G_25_x")

#geography + deprivation  
# has any coastline & lsoa IMD==1/2 ("GD")
# >10: 1,3,5, & IMD=1/2 ("GD_10_x")
#>=50%: 1,3,5, & IMD=1/2 ("GD_50_x")
#>=25%: 1,3,5, & IMD=1/2 ("GD_25_x") 

#geography + rural (based on rural urban index from ONS) 
# any coast and RU D1/2 E1/2 F1/2 ("GR")
# >10: 1,3,5, & rural ("GR_10_x")
# >=50%: 1,3,5 & rural ("GR_50_x")
# >=25%: 1,3,5 & rural ("GR_25_x")

#geography + deprived + rural 
# any coast & imd 1/2 & rural ("GDR")
#>10: 1,3,5, & imd 1/2 & rural ("GDR_10_x")
#>=50%: 1,3,5, & imd 1/2 & rural ("GDR_50_x")
# >=25%: 1,3,5, & imd 1/2 & rural ("GDR_25_x")

#lsoa level rural and IMD info

rural <- read.csv(paste0(Sys.getenv("coastal_clean"),"/Rural_Urban_Classification_(2011)_of_Lower_Layer_Super_Output_Areas_in_England_and_Wales.csv")) |>
  select(2,4) |>
  mutate(rural = ifelse(grepl("^D|^E|^F", RUC11CD), 1,0)) |>
  select(LSOA11CD, rural)
imd <- read.csv(paste0(Sys.getenv("coastal_clean"),"/Index_of_Multiple_Deprivation_(Dec_2019)_Lookup_in_England.csv")) |>
  select(c("LSOA11CD","IMD19")) |>
  rename(lsoa = "LSOA11CD") |>
  mutate(quint = ntile(IMD19, 5)) |>
  select(lsoa, quint)

covs <- inner_join(imd, rural, by=c("lsoa" = "LSOA11CD"))

#Create the set of 40 candidates of LSOA binary.
lsoas <- pcd |>
  select(PCD, LSOA11, RGN, NHSER, cstl_1k, cstl_3k, cstl_5k) |>
  inner_join(lsoa_coasts, by=c("LSOA11" = "LSOA11CD")) |>
  inner_join(covs, by = c("LSOA11" = "lsoa")) |>
  group_by(LSOA11) |>
  mutate(coastline = ifelse(sum(intersects > 0,na.rm=T), 1,0), 
         # does the LSOA have any coastline?
         allpcd = n()) |> #how many postcodes does each lsoa have?
  mutate(G1 = ifelse(coastline == 1, 1, 0),
         G_10_1 = ifelse(sum(cstl_1k,na.rm=T) >= 10, 1, 0),
         G_10_3 = ifelse(sum(cstl_3k,na.rm=T) >= 10, 1, 0),
         G_10_5 = ifelse(sum(cstl_5k,na.rm=T) >= 10, 1, 0),
         G_25_1 = ifelse((sum(cstl_1k,na.rm=T)/
                            max(allpcd,na.rm=T))>=0.25, 1, 0),
         G_25_3 = ifelse((sum(cstl_3k,na.rm=T)/
                            max(allpcd,na.rm=T))>=0.25, 1, 0),
         G_25_5 = ifelse((sum(cstl_5k,na.rm=T)/
                            max(allpcd,na.rm=T))>=0.25, 1, 0),
         G_50_1 = ifelse((sum(cstl_1k,na.rm=T)/
                            max(allpcd,na.rm=T))>=0.5, 1, 0),
         G_50_3 = ifelse((sum(cstl_3k,na.rm=T)/
                            max(allpcd,na.rm=T))>=0.5, 1, 0),
         G_50_5 = ifelse((sum(cstl_5k,na.rm=T)/
                            max(allpcd,na.rm=T))>=0.5, 1, 0),
         GD1 = case_when(coastline ==1 & quint <= 2 ~ 1,
                         .default = 0),
         GD_10_1 = case_when(sum(cstl_1k,na.rm=T) >= 10 & 
                               max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_10_3 = case_when(sum(cstl_3k,na.rm=T) >= 10 & 
                               max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_10_5 = case_when(sum(cstl_5k,na.rm=T) >= 10 & 
                               max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_50_1 = case_when(sum(cstl_1k,na.rm=T)/max(allpcd,na.rm=T) >= 
                               0.5 & max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_50_3 = case_when(sum(cstl_3k,na.rm=T)/max(allpcd,na.rm=T) >= 
                               0.5 & max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_50_5 = case_when(sum(cstl_5k,na.rm=T)/max(allpcd,na.rm=T) >= 
                               0.5 & max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_25_1 = case_when(sum(cstl_1k,na.rm=T)/max(allpcd,na.rm=T) >= 
                               0.25 & max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_25_3 = case_when(sum(cstl_3k,na.rm=T)/max(allpcd,na.rm=T) >= 
                               0.25 & max(quint) <= 2  ~ 1,
                             .default = 0),
         GD_25_5 = case_when(sum(cstl_5k,na.rm=T)/max(allpcd,na.rm=T) >= 
                               0.25 & max(quint) <= 2  ~ 1,
                             .default = 0),
         GR1 = case_when(coastline == 1 & max(rural) == 1 ~ 1,
                         .default = 0),
         GR_10_1 = case_when(sum(cstl_1k,na.rm=T) >= 10 & 
                               max(rural) == 1 ~ 1,
                             .default = 0),
         GR_10_3 = case_when(sum(cstl_3k,na.rm=T) >= 10 & 
                               max(rural) == 1 ~ 1,
                             .default = 0),
         GR_10_5 = case_when(sum(cstl_5k,na.rm=T) >= 10 & 
                               max(rural) == 1 ~ 1,
                             .default = 0),
         GR_50_1 = case_when(sum(cstl_1k,na.rm=T)/
                               max(allpcd,na.rm=T) >= 0.5 & max(rural) == 1 ~ 1,
                             .default = 0),
         GR_50_3 = case_when(sum(cstl_3k,na.rm=T)/
                               max(allpcd,na.rm=T) >= 0.5 & max(rural) == 1 ~ 1,
                             .default = 0),
         GR_50_5 = case_when(sum(cstl_5k,na.rm=T)/
                               max(allpcd,na.rm=T) >= 0.5 & max(rural) == 1 ~ 1,
                             .default = 0),
         GR_25_1 = case_when(sum(cstl_1k,na.rm=T)/
                               max(allpcd,na.rm=T) >= 0.25 & max(rural) == 1 ~ 1,
                             .default = 0),
         GR_25_3 = case_when(sum(cstl_3k,na.rm=T)/
                               max(allpcd,na.rm=T) >= 0.25 & max(rural) == 1 ~ 1,
                             .default = 0),
         GR_25_5 = case_when(sum(cstl_5k,na.rm=T)/
                               max(allpcd,na.rm=T) >= 0.25 & max(rural) == 1 ~ 1,
                             .default = 0),
         GDR1 = ifelse(coastline == 1 & quint <= 2 & max(rural) == 1, 1, 0),
         GDR_10_1 = ifelse(sum(cstl_1k,na.rm=T) >= 10 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_10_3 = ifelse(sum(cstl_3k,na.rm=T) >= 10 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_10_5 = ifelse(sum(cstl_5k,na.rm=T) >= 10 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_50_1 = ifelse((sum(cstl_1k,na.rm=T)/
                              max(allpcd,na.rm=T)) >=0.5 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_50_3 = ifelse((sum(cstl_3k,na.rm=T)/
                              max(allpcd,na.rm=T)) >=0.5 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_50_5 = ifelse((sum(cstl_5k,na.rm=T)/
                              max(allpcd,na.rm=T)) >=0.5 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_25_1 = ifelse((sum(cstl_1k,na.rm=T)/
                              max(allpcd,na.rm=T)) >=0.25 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_25_3 = ifelse((sum(cstl_3k,na.rm=T)/
                              max(allpcd,na.rm=T)) >=0.25 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0),
         GDR_25_5 = ifelse((sum(cstl_5k,na.rm=T)/
                              max(allpcd,na.rm=T)) >=0.25 & 
                             max(quint) <= 2 & max(rural) == 1, 1, 0)) |>
  select(LSOA11, RGN, rural, quint, coastline, allpcd, 
         G1, G_10_1, G_10_3, G_10_5, G_25_1, G_25_3, G_25_5, G_50_1, G_50_3, 
         G_50_5,GD1, GD_10_1, GD_10_3, GD_10_5, GD_50_1, GD_50_3, GD_50_5, 
         GD_25_1, GD_25_3, GD_25_5, GR1, GR_10_1, GR_10_3, GR_10_5, GR_50_1, 
         GR_50_3, GR_50_5, GR_25_1, GR_25_3, GR_25_5, GDR1, GDR_10_1, GDR_10_3, 
         GDR_10_5, GDR_50_1, GDR_50_3, GDR_50_5, GDR_25_1, GDR_25_3, GDR_25_5) |>
  distinct() |>
  droplevels() |>
  ungroup()


#write out the LSOA-level definitions data
#lsoa populations to use as labels for maps - 
#use 2016 as this is start of cancer data period used
pop <- dbGetQueryOracle(x, "select * from ONS2020.POPULATIONS") |>
  pivot_longer(cols=-c(YEAR, SEX, LSOA11), names_to = "Age", values_to = "count") |>
  group_by(YEAR, LSOA11) |>
  summarise(total = sum(count, na.rm=T)) |>
  filter(YEAR == 2016)

lsoas <- inner_join(lsoas, pop, by="LSOA11")
write.csv(lsoas, paste0(Sys.getenv("coastal_clean"), "/lsoa_candidates.csv"))

## 3. MAPPING ----
mapdat <- left_join(lsoas, lsoa_coasts, by=c("LSOA11" = "LSOA11CD"))

#get sensible definitions for each candidate to add to the maps
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


vars <- c("G1", "G_10_1", "G_10_3", "G_10_5", 
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



map.pal <- RColorBrewer::brewer.pal(12, "Paired")[c(10,11)]

map(vars, mapplot)
