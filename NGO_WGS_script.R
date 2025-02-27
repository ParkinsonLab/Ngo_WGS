#Sept 26th 2023 Duncan Carruthers-Lay

# Load libraries and set working directory 
library(tidyverse)
library(cluster)
library(factoextra)
library(countrycode)
library(varhandle)
library(ggpubr)
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(ggnewscale)
library(janitor)
library(ggmosaic)
library(gplots) 
library(RColorBrewer)
library(maps)
library(heatmaply)
library(jsonlite)
library(phytools)
library(mclust)
#library(clusterProfiler)
library(ggrepel)
library(vcfR)
library(qqman)
library(manhattanly)
library(gggenes)

setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC")

# Filtering of initial set of runs from the NCBI SRA  ---------------------
NCBI_SRA_fastq <- read.csv("SRA/SRA_RunTable_nofastq.txt", na.strings = c("", " ", NA)) %>%
  select(!isolate_name_alias)
NCBI_SRA_nofastq <- read.csv("SRA/SRA_RunTable_Sept26_2023.txt", na.strings = c("", " ", NA))

NCBI_SRA_raw <- bind_rows(NCBI_SRA_nofastq, NCBI_SRA_fastq)

# Date of isolation -------------------------------------------------------

NCBI_SRA_date <- NCBI_SRA_raw %>%
  select(Run, Collection_Date, Isolation_date, Collection.Date, collection_date, year.month) %>%
  mutate(Collection_Date = case_when(grepl(x=Collection_Date, 'missing', ignore.case = TRUE) ~ NA_character_,
                                     grepl(x=Collection_Date, 'not', ignore.case = TRUE) ~ NA_character_,
                                     TRUE ~ Collection_Date)) %>%
  mutate(Collection.Date = as.character(Collection.Date)) %>%
  rowwise() %>%
  mutate(Collection_Date = case_when(is.na(Collection_Date) ~ Isolation_date,
                                     TRUE ~ Collection_Date)) %>%
  mutate(Collection_Date = case_when(is.na(Collection_Date) ~ Collection.Date,
                                     TRUE ~ Collection_Date)) %>%
  mutate(Collection_Date = case_when(is.na(Collection_Date) ~ collection_date,
                                     TRUE ~ Collection_Date)) %>%
  mutate(Collection_Date = case_when(is.na(Collection_Date) ~ year.month,
                                     TRUE ~ Collection_Date)) %>%
  mutate(Collection_Date = case_when(nchar(Collection_Date) == 4 ~ Collection_Date,
                                     nchar(Collection_Date) == 10 & 
                                       str_detect(Collection_Date, "^[0-9]{4}") ~ substr(Collection_Date, 1,4),
                                     nchar(Collection_Date) == 10 &
                                       str_detect(Collection_Date, "^[0-9]{2}\\-") ~ substr(Collection_Date, 7,10),
                                     nchar(Collection_Date) == 21 ~ substr(Collection_Date, 1,4),
                                     nchar(Collection_Date) == 7 ~ substr(Collection_Date, 1,4),
                                     TRUE ~ Collection_Date)) %>%
  mutate(Collection_Date = case_when(grepl(x=Collection_Date, 'miss', ignore.case = TRUE) ~ NA_character_,
                                     grepl(x=Collection_Date, 'not', ignore.case = TRUE) ~ NA_character_,
                                     TRUE ~ Collection_Date)) %>%
  select(Run, Collection_Date) %>%
  filter(!is.na(Collection_Date))

no_date <- NCBI_SRA_raw %>%
  filter(!Run %in% NCBI_SRA_date$Run)

colnames(NCBI_SRA_date)[2] <- "Year"


# Country -----------------------------------------------------------------

NCBI_SRA_Country <- NCBI_SRA_raw %>%
  select(Run, Isolation_source, BioProject, geo_loc_name_country, geo_loc_name_country_continent,
         collected_by, geo_loc_name, Lat_Lon, geographic_location_.country_and.or_sea.,
         COUNTRY, Region, lat_lon, Location, isolation.source) %>%
  mutate(across(where(is.character), 
                ~if_else(. %in% c("missing", "Missing", "not applicable", 
                                  "Not Given", "not collected", "Not applicable",
                                  "Not Recorded", "not provided"), NA_character_, .))) %>%
  mutate(geo_loc_name = str_replace(geo_loc_name, "\\:.*", "")) %>%
  mutate(COUNTRY = str_replace(COUNTRY, "\\:.*", "")) %>%
  mutate(Country = case_when(is.na(geo_loc_name_country) ~ geo_loc_name,
                             geo_loc_name_country == "Viet Nam" ~ "Vietnam",
                             Isolation_source == "Western Australia" ~ "Australia",
                             geo_loc_name_country == "uncalculated" ~ NA,
                             TRUE ~ geo_loc_name_country)) %>%
  mutate(Country = case_when(is.na(Country) ~ COUNTRY,
                             TRUE ~ Country)) %>%
  mutate(Country = case_when(collected_by == "PHO" ~ "Canada",
                             TRUE ~ Country)) %>%
  mutate(Country = case_when(BioProject == "PRJNA992923" ~ "Japan",
                             BioProject == "PRJNA317462" ~ "USA",
                             BioProject == "PRJEB8940" ~ "Singapore",
                             BioProject == "PRJNA329501" ~ "USA",
                             BioProject == "PRJEB17615" ~ "Phillipines",
                             BioProject == "PRJEB25118" ~ "United Kingdom",
                             BioProject == "PRJNA926517" ~ "Belgium",
                             BioProject == "PRJEB34287" ~ "Ireland",
                             BioProject == "PRJDB12873" ~ "Japan",
                             BioProject == "PRJEB10016" ~ "USA",
                             BioProject == "PRJEB29480" ~ "Australia",
                             BioProject == "PRJDB10182" ~ "Japan",
                             BioProject == "PRJDB10572" ~ "Japan",
                             BioProject == "PRJDB6504" ~ "Japan",
                             BioProject == "PRJNA266539" ~ "Canada",
                             BioProject == "PRJNA322254" ~ "Korea",
                             BioProject == "PRJEB10104" ~ "Kenya",
                             BioProject == "PRJEB14933" ~ "United Kingdom",
                             BioProject == "PRJEB19989" ~ "United Kingdom",
                             BioProject == "PRJEB2124" ~ "United Kingdom",
                             BioProject == "PRJEB23008" ~ "United Kingdom",
                             BioProject == "PRJEB29738" ~ "Phillipines",
                             BioProject == "PRJEB7904" ~ "USA",
                             BioProject == "PRJEB58139" ~ NA,
                             BioProject == "PRJNA315363" ~ "United Kingdom",
                             BioProject == "PRJEB2090" ~ "USA",
                             BioProject == "PRJEB2999" ~ "USA",
                             BioProject == "PRJEB32746" ~ "United Kingdom",
                             TRUE ~ Country)) %>%
  mutate(Country = case_when(str_detect(Country, "Azerbaidjan") ~ "Azerbaidjan",
                             Country == "Viet Nam" ~ "Vietnam",
                             Country == "not provided" ~ NA,
                             TRUE ~ Country)) %>%
  select(Run, Country) %>%
  filter(!is.na(Country))

unique(NCBI_SRA_Country$Country)

# Sex ---------------------------------------------------------------------

NCBI_SRA_sex <- NCBI_SRA_raw %>%
  select(Run, sex, host_sex) %>%
  mutate(Sex = case_when(sex == "male" ~ "Male",
                         sex == "female" ~ "Female",
                         TRUE ~ NA)) %>%
  mutate(Sex = case_when(host_sex == "male" ~ "Male",
                         host_sex == "female" ~ "Female",
                         TRUE ~ Sex)) %>%
  select(Run, Sex) %>%
  filter(!is.na(Sex))
# Site of Isolation -------------------------------------------------------

NCBI_SRA_site <- NCBI_SRA_raw %>%
  select(Run, Isolation_source, host_tissue_sampled, isolation_source_host_associated, Source,
         isolation.source, isolation_source) %>%
  left_join(NCBI_SRA_sex) %>%
  mutate(Site = case_when(str_detect(Isolation_source, regex('vag', ignore_case = T)) |
                            str_detect(Isolation_source, regex('Cerv', ignore_case = T)) |
                            str_detect(Isolation_source, regex('endometrium', ignore_case = T)) ~ "Cervicovaginal",
                          str_detect(Isolation_source, regex('rect', ignore_case = T)) |
                            str_detect(Isolation_source, regex('ano', ignore_case = T)) |
                            str_detect(Isolation_source, regex('anal', ignore_case = T)) |
                            str_detect(Isolation_source, regex('anus', ignore_case = T)) |
                            str_detect(Isolation_source, regex('proctum', ignore_case = T)) ~ "Rectal",
                          str_detect(Isolation_source, regex('penis', ignore_case = T)) |
                            str_detect(Isolation_source, regex('urethra', ignore_case = T)) |
                            str_detect(Isolation_source, regex('foreskin', ignore_case = T)) |
                            str_detect(Isolation_source, regex('scrotum', ignore_case = T)) |
                            str_detect(Isolation_source, regex('penile', ignore_case = T)) |
                            str_detect(Isolation_source, regex("semen", ignore_case = T)) |
                            str_detect(Isolation_source, regex("ureter", ignore_case = T)) ~ "Penile",
                          str_detect(Isolation_source, regex('throat', ignore_case = T)) |
                            str_detect(Isolation_source, regex('phary', ignore_case = T)) |
                            str_detect(Isolation_source, regex('tonsil', ignore_case = T)) |
                            str_detect(Isolation_source, regex('oral', ignore_case = T)) ~ "Pharynx",
                          str_detect(Isolation_source, regex('blood', ignore_case = T)) |
                            str_detect(Isolation_source, regex('synovial', ignore_case = T)) |
                            str_detect(Isolation_source, regex('knee', ignore_case = T)) |
                            str_detect(Isolation_source, regex('hand', ignore_case = T)) |
                            str_detect(Isolation_source, regex("joint", ignore_case = T)) |
                            str_detect(Isolation_source, regex('articular', ignore_case = T)) |
                            str_detect(Isolation_source, regex("periton", ignore_case = T))~ "Disseminated",
                          str_detect(Isolation_source, regex('eye', ignore_case = T)) |
                            str_detect(Isolation_source, regex('conjunctiva', ignore_case = T)) ~ "Ocular",
                          str_detect(Isolation_source, regex('urin', ignore_case = T)) |
                            str_detect(Isolation_source, regex('genital', ignore_case = T)) |
                            str_detect(Isolation_source, regex('reproductive', ignore_case = T))~ "Urogenital", 
                          TRUE ~ NA)) %>%
  mutate(Site = case_when(str_detect(host_tissue_sampled, regex('vag', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('Cerv', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('endometrium', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("uterus", ignore_case = T))~ "Cervicovaginal",
                          str_detect(host_tissue_sampled, regex('rect', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('ano', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('anal', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('anus', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('proctum', ignore_case = T)) ~ "Rectal",
                          str_detect(host_tissue_sampled, regex('penis', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('urethra', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('foreskin', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('scrotum', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('penile', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("semen", ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("ureter", ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("uret", ignore_case = T)) ~ "Penile",
                          str_detect(host_tissue_sampled, regex('throat', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('phary', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('tonsil', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('oral', ignore_case = T)) ~ "Pharynx",
                          str_detect(host_tissue_sampled, regex('blood', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('synovial', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('knee', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('hand', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("joint", ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('articular', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("ankle", ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex("periton", ignore_case = T))~ "Disseminated",
                          str_detect(host_tissue_sampled, regex('eye', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('conjunctiva', ignore_case = T)) ~ "Ocular",
                          str_detect(host_tissue_sampled, regex('urin', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('genital', ignore_case = T)) |
                            str_detect(host_tissue_sampled, regex('reproductive', ignore_case = T))~ "Urogenital", 
                          TRUE ~ Site)) %>%
  mutate(Site = case_when(str_detect(isolation_source_host_associated, regex('urine', ignore_case = T)) ~ "Urogenital",
                          str_detect(isolation_source_host_associated, regex('uret', ignore_case = T)) ~ "Penile",
                          TRUE ~ Site)) %>%
  mutate(Site = case_when(Site == "Urogenital" & Sex == "Male" ~ "Urogenital_Male",
                          Site == "Urogenital" & Sex == "Female" ~ "Urogenital_Female",
                          TRUE ~ Site)) %>%
  select(Run, Site) %>%
  filter(!is.na(Site))


site_count <- NCBI_SRA_site %>%
  count(Site) 

# Final list of accessions ------------------------------------------------
#Add in WHO reads from Daniel
WHO_strains_list <- c("ERR349898", "ERR363606", "ERR363595", "ERR388421",
                      "ERR388420", "ERR363586", "ERR349899", "ERR349905",
                      "ERR449479", "ERR1448254", "ERR439282", "ERR449539",
                      "ERR449540", "ERR1448255")
NCBI_SRA <- NCBI_SRA_raw %>%
  filter(Run %in% NCBI_SRA_Country$Run &
           Run %in% NCBI_SRA_date$Run |
           Run %in% WHO_strains_list)
SRA_metadata <- NCBI_SRA %>%
  select(Run, BioProject, BioSample) %>%
  left_join(NCBI_SRA_Country) %>%
  left_join(NCBI_SRA_date) %>%
  left_join(NCBI_SRA_sex) %>%
  left_join(NCBI_SRA_site) %>%
  distinct(BioSample, .keep_all = T) %>%
  mutate(WHO = case_when(Run %in% WHO_strains_list ~ "WHO"))

#Make metadata file and list of accessions
#write.csv(SRA_metadata, "SRA_metadata_Sep2023.csv", quote = F, row.names = F)
#write.table(SRA_metadata$Run, "SRA_accessions_Sep2023.txt", quote = F, row.names = F, col.names = F)

# Defining reference set --------------------------------------------------
#Going to start with the manual curation of sequences for the reference collection.
#Focus is going to be on metadata diversity among country and date.
all_combos <- SRA_metadata %>%
  count(Country, Year)
metadata_CY <- SRA_metadata %>%
  unite(c("Country", "Year"), remove = T, col = "CY") %>%
  mutate(CY = case_when(Site == "Ocular" ~ "Ocular",
                        Site == "Disseminated" ~ "DGI",
                        TRUE ~ CY))
combos_grouped <- all_combos %>%
  unite(c("Country", "Year"), remove = T, col = "CY")

#separate out the sets for jaccard calculations
rare_combos <- combos_grouped %>%
  filter(n < 3)
rare_runs <- metadata_CY %>%
  filter(CY %in% rare_combos$CY)
Multi3_combos <- combos_grouped %>%
  filter(n <= 25 & n >=3)
#Add in ocular and disseminated strains
OC_DGI <- metadata_CY %>%
  filter(CY == "Ocular" | CY == "DGI")
Multi3 <- metadata_CY %>%
  filter(CY %in% Multi3_combos$CY) %>%
  filter(!Run %in% rare_runs$Run) %>%
  rbind(OC_DGI) %>%
  select(Run) %>%
  distinct()
USA1 <- metadata_CY %>%
  filter(CY == "USA_2018" | CY == "USA_2019") %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run)
USA2 <- metadata_CY %>%
  filter(CY == "USA_2021" | CY == "USA_2022") %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run)
AUS_UK <- SRA_metadata %>%
  filter(Country == "Australia" | Country == "United Kingdom") %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run)
NOR_CAN <- SRA_metadata %>%
  filter(Country == "Norway" | Country == "Canada") %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run)
USA3 <- SRA_metadata %>%
  filter(Country == "USA") %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(!Run %in% USA1$Run) %>%
  filter(!Run %in% USA2$Run) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run)
Multi_1_countries <- c("Austria", "Netherlands", "New Zealand", "China", "Portugal", "Brazil")
Multi1 <- SRA_metadata %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(Country %in% Multi_1_countries) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run)
Multi2_countries <- c("Vietnam", "Japan", "Argentina", "Hong Kong", "Spain", "Uganda", "Switzerland", 
                      "South Africa", "Belgium", "Madagascar", "Dominican Republic", "Malawi")
Multi2 <- SRA_metadata %>%
  filter(!Run %in% rare_runs$Run) %>%
  filter(Country %in% Multi2_countries) %>%
  filter(!Run %in% Multi3$Run) %>%
  select(Run) 

# Processing jaccard results ----------------------------------------------
#Make a function that does the k-medoids PAM clustering. First select a k by silhouette method. For sets with >1k sequences; make max k = 15.
#For between 16-999 sequence sets, find optimal k with k max of 15. For 5-14 find optimal k with kmax of n-1.  
set.seed(127405675)
ref_set <- data.frame(Run = NULL, Cluster_size = NULL, CY = NULL)
ref_set_clust <- data.frame(Run = NULL, Cluster = NULL, CY = NULL)

reference_selection <- function(CoYe) {who
  CY_meta <- metadata_CY %>%
    filter(CY == CoYe)
  CY_jaccard <- jaccard %>%
    filter(V1 %in% CY_meta$Run & V2 %in% CY_meta$Run)
  CY_matrix <- CY_jaccard %>%
    pivot_wider(names_from = V2, values_from = V3) %>%
    column_to_rownames(var = "V1")
  if (nrow(CY_meta) >= 1000) {
    CY_cluster <- pam(CY_matrix, k = 15, diss = T, nstart = 100, pamonce = 6)
  } else if (nrow(CY_meta) < 1000 & nrow(CY_meta) > 15) {
    optimal_k <- fviz_nbclust(CY_matrix, cluster::pam, method = "silhouette", k.max = 15, print.summary = T)  
    CY_cluster <- pam(CY_matrix, k = optimal_k$layers[[3]]$data$xintercept, diss = T, nstart = 100, pamonce = 6) 
  } else if (nrow(CY_meta) <= 15) {
    optimal_k <- fviz_nbclust(CY_matrix, cluster::pam, method = "silhouette", k.max = (nrow(CY_meta) - 1), print.summary = T)  
    CY_cluster <- pam(CY_matrix, k = optimal_k$layers[[3]]$data$xintercept, diss = T, nstart = 100, pamonce = 6)
  }
  ref_seqs <- data.frame(Run = CY_cluster$medoids, Cluster_size = summary(as.factor(CY_cluster$clustering)))
  ref_seqs$CY <- CoYe
  ref_set <<- ref_set %>%
    rbind(ref_seqs)
  ref_clust <- data.frame(CY_cluster$clustering) %>%
    rownames_to_column("Run")
  colnames(ref_clust)[2] <- "Cluster"
  ref_clust$CY <- CoYe
  ref_set_clust <<- ref_set_clust %>%
    rbind(ref_clust)
}

#USA3
#List of CY to test
USA3_CY <- metadata_CY %>%
  filter(Run %in% USA3$Run) %>%
  select(CY) %>%
  unique()
jaccard <- read.delim("jaccard/USA3_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(USA3_CY$CY,reference_selection)

#Multi1
Multi1_CY <-metadata_CY %>%
  filter(Run %in% Multi1$Run) %>%
  select(CY) %>%
  unique() 
jaccard <- read.delim("jaccard/Multi1_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(Multi1_CY$CY,reference_selection)

#NORCAN
NORCAN_CY <-metadata_CY %>%
  filter(Run %in% NOR_CAN$Run) %>%
  select(CY) %>%
  unique() 
jaccard <- read.delim("jaccard/NORCAN_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(NORCAN_CY$CY,reference_selection)

#USA1
USA1_CY <-metadata_CY %>%
  filter(Run %in% USA1$Run) %>%
  select(CY) %>%
  unique() 
jaccard <- read.delim("jaccard/USA1_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(USA1_CY$CY,reference_selection)

#USA2
USA2_CY <-metadata_CY %>%
  filter(Run %in% USA2$Run) %>%
  select(CY) %>%
  unique() 
jaccard <- read.delim("jaccard/USA2_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(USA2_CY$CY,reference_selection)

#AUS_UK
AUS_UK_CY <-metadata_CY %>%
  filter(Run %in% AUS_UK$Run) %>%
  select(CY) %>%
  unique() 
jaccard <- read.delim("jaccard/AUS_UK_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(AUS_UK_CY$CY,reference_selection)

#Multi2
Multi2_CY <-metadata_CY %>%
  filter(Run %in% Multi2$Run) %>%
  select(CY) %>%
  unique() 
jaccard <- read.delim("jaccard/Multi2_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(Multi2_CY$CY,reference_selection)

#Multi3
Multi3_CY <- metadata_CY %>%
  filter(Run %in% Multi3$Run) %>%
  select(CY) %>%
  unique()
jaccard <- read.delim("jaccard/Multi3_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(Multi3_CY$CY,reference_selection)

write.csv(ref_set_clust, "reference_set_clusters.csv", quote = F, row.names = F)
# Final reference set -----------------------------------------------------
# Add rare set
final_ref_set <- SRA_metadata %>%
  filter(Run %in% ref_set$Run |
           Run %in% rare_runs$Run) %>%
  filter(is.na(WHO)) %>%
  distinct()
#Add WHO
WHO_meta <- read.csv("WHO_metadata.csv")
final_ref_set <- final_ref_set %>%
  rbind(WHO_meta) %>%
  distinct()

#write.csv(final_ref_set, "final_ref_set.csv", row.names = F, quote = F)

#final_ref_acc <- final_ref_set %>%
  #select(Run)
#write.table(final_ref_acc, "final_ref_acc.txt", row.names = F, col.names = F, quote = F)


# Initial reference set data ----------------------------------------------
final_ref_set <- read.csv("final_ref_set.csv")

final_site_count <- final_ref_set %>%
  count(Site)
final_country_count <- final_ref_set %>%
  count(Country)
final_sex_count <- final_ref_set %>%
  count(Sex)
final_year_count <- final_ref_set %>%
  count(Year)
final_bioproject_count <- final_ref_set %>%
  count(BioProject)
# Continent clustering ----------------------------------------------------

final_ref_set <- read.csv("final_ref_set.csv")

#Cluster the 1307 by continent and year again.
set.seed(127405675)
ref_cont_set <- data.frame(Run = NULL, Cluster_size = NULL, CY = NULL)
ref_cont_set_clust <- data.frame(Run = NULL, Cluster = NULL, CY = NULL)

reference_selection_continent <- function(ConYe) {
  CY_meta <- final_metadata_CY %>%
    filter(ConY == ConYe)
  CY_jaccard <- jaccard %>%
    filter(V1 %in% CY_meta$Run & V2 %in% CY_meta$Run)
  CY_matrix <- CY_jaccard %>%
    pivot_wider(names_from = V2, values_from = V3) %>%
    column_to_rownames(var = "V1")
  optimal_k <- fviz_nbclust(CY_matrix, cluster::pam, method = "silhouette", k.max = (nrow(CY_meta) - 1), print.summary = T)  
  CY_cluster <- pam(CY_matrix, k = optimal_k$layers[[3]]$data$xintercept, diss = T, nstart = 100, pamonce = 6)
  ref_seqs <- data.frame(Run = CY_cluster$medoids, Cluster_size = summary(as.factor(CY_cluster$clustering)))
  ref_seqs$CY <- ConYe
  ref_cont_set <<- ref_cont_set %>%
    rbind(ref_seqs)
  ref_cont_clust <- data.frame(CY_cluster$clustering) %>%
    rownames_to_column("Run")
  colnames(ref_cont_clust)[2] <- "Cluster"
  ref_cont_clust$CY <- ConYe
  ref_cont_set_clust <<- ref_cont_set_clust %>%
    rbind(ref_cont_clust)
}

#Reference set inception
final_ref_set$Continent <- countrycode(sourcevar = final_ref_set$Country, origin = "country.name",
                                       destination = "region")
final_ref_set <- final_ref_set %>%
  mutate(Continent = case_when(Country == "Australia" ~ "Oceania",
                               Country == "New Zealand" ~ "Oceania",
                               Continent == "East Asia & Pacific" ~ "Asia",
                               Country == "Azerbaijan" ~ "Asia",
                               Continent == "Europe & Central Asia" ~ "Europe",
                               Continent == "Latin America & Caribbean" ~ "South & Central America",
                               Continent == "Sub-Saharan Africa" ~ "Africa",
                               Country == "Israel" ~ "Asia",
                               WHO == "WHO_L" ~ "Asia",
                               TRUE ~ Continent))
continent_mapping <- final_ref_set %>%
  count(Continent, Country)
final_metadata_CY <- final_ref_set %>%
  unite(c("Continent", "Year"), remove = T, col = "ConY")
ConY_count <- final_metadata_CY %>%
  count(ConY)
cluster_ConY <- ConY_count %>%
  filter(n >= 3)
no_clust_ConY <- ConY_count %>%
  filter(n < 3)
final_CY <- final_metadata_CY %>%
  filter(ConY %in% cluster_ConY$ConY) %>%
  select(ConY) %>%
  distinct()
jaccard <- read.delim("final_ref_jaccard.txt", header = F, sep = " ") %>%
  distinct()
lapply(final_CY$ConY, reference_selection_continent)

no_clust_set <- final_metadata_CY %>%
  filter(ConY %in% no_clust_ConY$ConY)
continent_clust_set <- final_metadata_CY %>%
  filter(Run %in% ref_cont_set$Run)

continent_clust_set <- continent_clust_set %>%
  rbind(no_clust_set)
#add in all DGI and ocular strains from ref set as well as the WHO strains
DGI_WHO_ref <- final_metadata_CY %>%
  filter(Site == "Ocular" | Site == "Disseminated" | !is.na(WHO))
final_continent_set <- continent_clust_set %>%
  rbind(DGI_WHO_ref) %>%
  separate(ConY, into = c("Continent", "Year"), sep = "_") %>%
  distinct()
#write.csv(final_continent_set, "final_continent_set.csv", quote = F, row.names = F)
final_continent_acc <- final_continent_set %>%
  select(Run)
write.table(final_continent_acc, "final_continent_acc.txt", quote = F, row.names = F, col.names = F)
# MIC/MLST metadata -------------------------------------------------------
#Get metadata from MLST to supplement
isolate_metadata <- read.csv("isolate_metadata.csv")
pubMLST_1 <- readxl::read_excel("MIC_data/BIGSdb_3631124_6817505180_90348.xlsx")
pubMLST_2 <- readxl::read_excel("MIC_data/BIGSdb_3631164_2308297437_23309.xlsx")
pubMLST_3 <- readxl::read_excel("MIC_data/BIGSdb_3631505_4988794147_06338.xlsx")

pubMLST <- rbind(pubMLST_1, pubMLST_2) %>%
  rbind(pubMLST_3)

pub_meta <- pubMLST %>%
  filter(run_accession %in% isolate_metadata$Run) %>%
  select(!c("epidemiology", "disease", "serogroup", "genogroup", "capsule_group", "serotype", "sero_subtype", 
            "NCBI_assembly_accession", "clonal_complex (MLST)"))
no_pub <- isolate_metadata %>%
  filter(!Run %in% pub_meta$run_accession) %>%
  count(BioProject)
missing_info <- isolate_metadata %>%
  filter(Run %in% pub_meta$run_accession) %>%
  filter(is.na(Sex) | is.na(Site))
sex_site_mlst <- pub_meta %>%
  select(run_accession, sex, source) %>%
  filter(!is.na(sex) | !is.na(source))

MIC_pub_data <- pub_meta %>%
  mutate(AZM_R = case_when(azithromycin_mic <= 1 ~ "S",
                           azithromycin_mic >= 2 ~ "R",
                           between(azithromycin_mic, 1, 2) ~ "DS",
                           TRUE ~ NA)) %>%
  mutate(PEN_R = case_when(penicillin_mic <= 0.06 ~ "S",
                           penicillin_mic >= 2 ~ "R",
                           between(penicillin_mic, 0.06, 0.12) ~ "DS",
                           between(penicillin_mic, 0.12, 2) ~ "I",
                           TRUE ~ NA)) %>%
  mutate(TET_R = case_when(tetracycline_mic <= 0.25 ~ "S",
                           tetracycline_mic >= 2 ~ "R",
                           between(tetracycline_mic, 0.25, 0.5) ~ "DS",
                           between(tetracycline_mic, 0.5, 2) ~ "I",
                           TRUE ~ NA)) %>%
  mutate(SPT_R = case_when(spectinomycin_mic <= 32 ~ "S",
                           spectinomycin_mic >= 128 ~ "R",
                           between(spectinomycin_mic, 32, 63) ~ "DS",
                           between(spectinomycin_mic, 64, 128) ~ "I",
                           TRUE ~ NA)) %>%
  mutate(CIP_R = case_when(ciprofloxacin_mic <= 0.06 ~ "S",
                           ciprofloxacin_mic >= 1 ~ "R",
                           between(ciprofloxacin_mic, 0.06, 0.11) ~ "DS",
                           between(ciprofloxacin_mic, 0.12, 1) ~ "I",
                           TRUE ~ NA)) %>%
  mutate(CRO_R = case_when(ceftriaxone_mic <= 0.124 ~ "S",
                           ceftriaxone_mic >= 0.125 ~ "DS",
                           #between(ciprofloxacin_mic, 0.06, 0.11) ~ "DS",
                           #between(ciprofloxacin_mic, 0.12, 1) ~ "I",
                           TRUE ~ NA)) %>%
  mutate(CFM_R = case_when(cefixime_mic <= 0.24 ~ "S",
                           cefixime_mic >= 0.25 ~ "DS",
                           #between(ciprofloxacin_mic, 0.06, 0.11) ~ "DS",
                           #between(ciprofloxacin_mic, 0.12, 1) ~ "I",
                           TRUE ~ NA)) %>%
  select(run_accession, CFM_R, CRO_R, CIP_R, SPT_R, TET_R, PEN_R, AZM_R)
colnames(MIC_pub_data)[1] <- "name"

write.csv(MIC_pub_data, "MIC_data/pubmlst_MIC.csv", quote = F, row.names = F)

# PopNet results ----------------------------------------------------------
PopNet_nodes <- read.csv("popnet_meta.csv")
isolate_metadata <- read.csv("isolate_metadata.csv")
colnames(PopNet_nodes)[3] <- "Run"

isolate_metadata <- left_join(isolate_metadata, PopNet_nodes)
#write.csv(isolate_metadata, "isolate_metadata.csv", quote = F, row.names = F)

#chromosome painting - not sure if will include
chromosome_painting <- read.delim("PopNet/chromosome_paintings.tsv") %>%
  unite(Segment, c(CHR, POS))
isolate_cluster <- isolate_metadata %>%
  select(Run, group) %>%
  arrange(group)
cp_trans <- as.data.frame(t(chromosome_painting)) %>%
  rownames_to_column() %>%
  row_to_names(row_number = 1) %>%
  mutate(across(-Segment, as.numeric))
colnames(cp_trans)[1] <- "Strain"

cp_long <- pivot_longer(data = cp_trans,
                        cols = -Strain,
                        names_to = "Segment",
                        values_to = "PopNet_group") %>%
  filter(Strain != "None")
cp_long$PopNet_group <- as.factor(cp_long$PopNet_group)

#order by popnet group and chromosome segment
cp_long$Strain <- factor(cp_long$Strain,
                         levels = isolate_cluster$Run,
                         ordered = T)
cp_long$Segment <- factor(x=cp_long$Segment,
                          levels = chromosome_painting$Segment,
                          ordered = T)
#new color palette
PopNet_colours <- isolate_metadata %>%
  select(group, color) %>%
  distinct() %>%
  filter(!is.na(group)) %>%
  arrange(group)
group_segments <- cp_long %>%
  count(PopNet_group)
#plot
ggplot(data = cp_long, aes(x=Segment, y=Strain)) +
  geom_tile(aes(fill = PopNet_group)) +
  theme(axis.text.x=element_text(size = 1, angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 0.5),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = PopNet_colours$color)
cp_heatmap_plot

#zoomed in chromosome painting
IKE_isolates <-isolate_metadata %>%
  filter(group == "K" | group == "E" | group == "I")
IKE_long <- cp_long %>%
  filter(Strain %in% IKE_isolates$Run)

ggplot(data = IKE_long, aes(x=Segment, y=Strain)) +
  geom_tile(aes(fill = PopNet_group)) +
  theme(axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = PopNet_colours$color)

# Pangenome ---------------------------------------------------------------
panaroo_gc <- read.delim("panaroo_results/gene_presence_absence.rtab")
panaroo_gpa <- read.csv("panaroo_results/gene_presence_absence.csv")
roary_gpa <- read.csv("panaroo_results/gene_presence_absence_roary.csv") %>% 
  select(Gene, Isolates)
gene_meta <- panaroo_gpa %>%
  select(Gene, Annotation)
eggnog_annotations <- read.delim("IPNC_annotations.tsv")

gene_meta <- left_join(gene_meta, eggnog_annotations)
gene_meta <- left_join(gene_meta, roary_gpa)
gene_meta <- gene_meta %>%
  mutate(Pangenome = case_when(Isolates > 463 ~ "Core",
                               Isolates < 73 ~ "Cloud",
                               TRUE ~ "Shell")) %>%
  mutate(across(everything(), str_replace_all, ",", ";"))

#gene_meta <- write.csv(gene_meta, "gene_metadata.csv", quote = F, row.names = F)

isolates_per_genomes <- panaroo_gc %>%
  column_to_rownames("Gene")
summary(colSums(isolates_per_genomes))

#heatmap
isolate_metadata <- read.csv("isolate_metadata.csv")

gene_meta <- read.csv("gene_metadata.csv")

Blac <- gene_meta %>%
  filter(Blac_p == "Yes")

gene_hmap_df <- panaroo_gc %>%
  column_to_rownames("Gene") %>%
  t() 

my_group <- as.numeric(as.factor(substr(rownames(gene_hmap_df), 1 , 1)))

isolate_popnet_colours <- data.frame(Run = row.names(gene_hmap_df)) %>%
  left_join(isolate_metadata)
popnet_groups <- data.frame(group = as.factor(isolate_popnet_colours$group)) %>%
  left_join(PopNet_colours)

pal <- colorRampPalette(c("grey95", "#023958"))(2)

heatmap.2(gene_hmap_df, Colv = FALSE, labRow = NA, labCol = NA, dendrogram = 'row',
        xlab = "Genes", ylab = "Isolates", RowSideColors = isolate_popnet_colours$color,
        col = pal, key = FALSE)

#interactive heatmap with heatmaply of only accessory genes to speed up
gene_meta <- read.csv("gene_metadata.csv")
accessory_genes <- gene_meta %>%
  filter(Pangenome != "Core")

acces_hmap_df <- panaroo_gc %>% 
  filter(Gene %in% accessory_genes$Gene) %>%
  column_to_rownames("Gene") %>%
  t()

hmap_gradient <- ggplot2::scale_fill_gradient2(low = "lightgrey", high = "navy",
                                               midpoint = 0.5, limits =c(0,1))
heatmaply(acces_hmap_df, Colv = NA,
          scale_fill_gradient_fun = hmap_gradient, hide_colorbar = TRUE,
          showticklabels = FALSE)


# Mobile elements ---------------------------------------------------------
isolate_metadata <- read.csv("isolate_metadata.csv")
gene_meta <- read.csv("gene_metadata.csv")

#GGI 
GGI_blast <- read.delim("mobile_elements/results/GGI_blast.tsv", header = F) %>%
  select(V1, V2)
colnames(GGI_blast) <- c("Gene", "GGI_gene")

#cryptic 
cryptic_blast <- read.delim("mobile_elements/results/cryptic_blast.tsv", header = F) %>%
  select(V1, V2)
colnames(cryptic_blast) <- c("Gene", "cryptic_gene")

#Blac
Blac_blast <- read.delim("mobile_elements/results/Blac_blast.tsv", header = F) %>%
  select(V1, V2)
colnames(Blac_blast) <- c("Gene", "Blac_gene")

#Conj
pConj_blast <- read.delim("mobile_elements/results/pconj_blast.tsv", header = F) %>%
  select(V1, V2)
colnames(pConj_blast) <- c("Gene", "Conjugative_gene")

pConj_D_tet_blast <- read.delim("mobile_elements/results/pconj_d_tet_blast.tsv", header = F) #%>%
  select(V1, V2)
colnames(pConj_D_tet_blast) <- c("Gene", "Conjugative_gene")

pConj_US_tet_blast <- read.delim("mobile_elements/results/pconj_US_tet_blast.tsv", header = F) %>%
  select(V1, V2)
colnames(pConj_US_tet_blast) <- c("Gene", "Conjugative_gene")

unified_pConj <- rbind(pConj_blast, pConj_D_tet_blast, pConj_US_tet_blast) %>%
  group_by(Gene) %>%
  summarise(Conjugative_gene = toString(Conjugative_gene))
unified_pConj$Conjugative_gene <- str_replace_all(unified_pConj$Conjugative_gene, ",", ";")

#merge with gene metadata
gene_meta <- gene_meta %>%
  left_join(GGI_blast) %>%
  left_join(cryptic_blast) %>%
  left_join(Blac_blast) %>%
  left_join(unified_pConj)

#identify mobile elements for genes
gene_meta <- gene_meta %>%
  mutate(GGI = case_when(!is.na(GGI_gene) ~ "Yes",
                         is.na(GGI_gene) ~ "No")) %>%
  mutate(cryptic_p = case_when(!is.na(cryptic_gene) ~ "Yes",
                               is.na(cryptic_gene) ~ "No")) %>%
  mutate(Blac_p = case_when(!is.na(Blac_gene) ~ "Yes",
                            is.na(Blac_gene) ~ "No")) %>%
  mutate(pConj = case_when(!is.na(Conjugative_gene) ~ "Yes",
                           is.na(Conjugative_gene) ~ "No")) 
gene_meta[gene_meta == "-"] <- NA

#add which mobile element each gene belongs to
gene_meta <- gene_meta %>%
  mutate(mobile_element = case_when(pConj == "Yes" ~ "pConj",
                                    GGI == "Yes" ~ "GGI",
                                    cryptic_p == "Yes" ~ "cryptic",
                                    Blac_p == "Yes" ~ "Blac",
                                    TRUE ~ NA))

write.csv(gene_meta, "gene_metadata.csv", quote = F, row.names = F)

#identify which isolates have which mobile elements
isolate_metadata <- read.csv("isolate_metadata.csv")

#isolates with GGI genes
panaroo_gc <- read.delim("panaroo_results/gene_presence_absence.rtab")
gene_meta <- read.csv("gene_metadata.csv")
GGI_genes <- gene_meta %>%
  filter(!is.na(GGI_gene))
panaroo_GGI <- panaroo_gc %>%
  filter(Gene %in% GGI_genes$Gene)
GGI_count <- colSums(panaroo_GGI[,2:ncol(panaroo_GGI)])
GGI_count <- GGI_count[GGI_count > 55]
GGI_isolates <- names(GGI_count)

isolate_metadata <- isolate_metadata %>%
  mutate(GGI = case_when(Run %in% GGI_isolates ~ "Yes",
                         TRUE ~ "No"))
#isolates with Blac plasmid
panaroo_gc <- read.delim("panaroo_results/gene_presence_absence.rtab")
gene_meta <- read.csv("gene_metadata.csv")
Blac_genes <- gene_meta %>%
  filter(!is.na(Blac_gene))
panaroo_Blac <- panaroo_gc %>%
  filter(Gene %in% Blac_genes$Gene)
Blac_isolates <- as.vector((panaroo_Blac[3,] == 1))
names(Blac_isolates) <- colnames(panaroo_Blac)
Blac_isolates <- Blac_isolates[Blac_isolates == TRUE]

isolate_metadata <- isolate_metadata %>%
  mutate(Blac = case_when(Run %in% names(Blac_isolates) ~ "Yes",
                         TRUE ~ "No"))
#isolates with cryptic plasmid
panaroo_gc <- read.delim("panaroo_results/gene_presence_absence.rtab")
gene_meta <- read.csv("gene_metadata.csv")
cryptic_genes <- gene_meta %>%
  filter(!is.na(cryptic_gene))
panaroo_cryptic <- panaroo_gc %>%
  filter(Gene %in% cryptic_genes$Gene)
cryptic_count <- colSums(panaroo_cryptic[,2:ncol(panaroo_cryptic)])
cryptic_count <- cryptic_count[cryptic_count >= 8]
cryptic_isolates <- names(cryptic_count)

isolate_metadata <- isolate_metadata %>%
  mutate(cryptic = case_when(Run %in% cryptic_isolates ~ "Yes",
                         TRUE ~ "No"))

#isolates with conjugative plasmids; need to distinguish between Tet and no Tet
panaroo_gc <- read.delim("panaroo_results/gene_presence_absence.rtab")
gene_meta <- read.csv("gene_metadata.csv")
conj_genes <- gene_meta %>%
  filter(!is.na(Conjugative_gene))
panaroo_conj <- panaroo_gc %>%
  filter(Gene %in% conj_genes$Gene)
conj_count <- colSums(panaroo_conj[,2:ncol(panaroo_conj)])
conj_count <- conj_count[conj_count >= 40]
conj_isolates <- names(conj_count)

tet_isolates <- as.vector((panaroo_conj[44,] == 1))
names(tet_isolates) <- colnames(panaroo_conj)
tet_isolates <- tet_isolates[tet_isolates == TRUE]

conj_tet_isolates <- intersect(names(tet_isolates), conj_isolates)
conj_isolates <- setdiff(conj_isolates, names(tet_isolates))

isolate_metadata <- isolate_metadata %>%
  mutate(pConj = case_when(Run %in% conj_isolates ~ "Yes",
                             TRUE ~ "No")) %>%
  mutate(pConj_Tet = case_when(Run %in% conj_tet_isolates ~ "Yes",
                               TRUE ~ "No"))
#overlap of mobile elements
isolate_metadata <- isolate_metadata %>%
  mutate(GGI_Blac = case_when(Blac == "Yes" & GGI == "Yes" ~ "Both",
                              Blac == "Yes" & GGI == "No" ~ "Blac",
                              Blac == "No" & GGI == "Yes" ~ "GGI",
                              TRUE ~ "None")) %>%
  mutate(GGI_pConj = case_when(pConj == "Yes" & GGI == "Yes" ~ "Both",
                               pConj_Tet == "Yes" & GGI == "Yes" ~ "Both",
                               pConj == "Yes" & GGI == "No" ~ "PConj",
                               pConj_Tet == "Yes" & GGI == "No" ~ "PConj",
                               pConj == "No" & GGI == "Yes" ~ "GGI",
                               pConj_Tet == "No" & GGI == "Yes" ~ "GGI",
                               TRUE ~ "None")) %>%
  mutate(Blac_pConj = case_when(pConj == "Yes" & Blac == "Yes" ~ "Both",
                                pConj_Tet == "Yes" & Blac == "Yes" ~ "Both",
                                pConj == "Yes" & Blac == "No" ~ "PConj",
                                pConj_Tet == "Yes" & Blac == "No" ~ "PConj",
                                pConj == "No" & Blac == "Yes" ~ "Blac",
                                pConj_Tet == "No" & Blac == "Yes" ~ "Blac",
                                TRUE ~ "None")) %>%
  mutate(Mobile_element = case_when(pConj == "No" & Blac == "No" & pConj == "No" 
                                    & pConj_Tet == "No" & GGI == "No" ~ "No",
                                    TRUE ~ "Yes"))
write.csv(isolate_metadata, "isolate_metadata.csv", quote = F, row.names = F)

# AMR ---------------------------------------------------------------------

#Pathogenwatch predictions
pwatch_AMR <- read.csv("pwatch/amr_predict.csv")
AMR_metadata <- pwatch_AMR %>%
  select(Run = Genome.Name, Azithromycin, Cefixime, Ceftriaxone, Penicillin, Ciprofloxacin, Tetracycline)

isolate_metadata <- read.csv("isolate_metadata.csv")
isolate_metadata <- isolate_metadata %>%
  left_join(AMR_metadata) %>%
  mutate(MDR_drug = case_when(Azithromycin == "RESISTANT" & Cefixime == "NOT_FOUND" ~ "AZM",
                              Azithromycin == "NOT_FOUND" & Cefixime == "RESISTANT" ~ "CFM",
                              Azithromycin == "NOT_FOUND" & Cefixime == "NOT_FOUND" ~ "None",
                              Azithromycin == "RESISTANT" & Cefixime == "RESISTANT" ~ "Both",
                              TRUE ~ "None"))

write.csv(isolate_metadata, "isolate_metadata.csv", quote = F, row.names = F)

#Figure displaying predicted resistance by antimicrobial
pwatch_AMR_df <- pwatch_AMR %>%
  select(!c("Genome.ID", "Version", "Library.Version", "Spectinomycin", "Sulfonamides")) %>%
  pivot_longer(!Genome.Name, names_to = "Antibiotic", values_to = "Resistance") %>%
  count(Antibiotic, Resistance)
pwatch_AMR_df$Resistance <- str_replace_all(pwatch_AMR_df$Resistance, "NOT_FOUND", "SUSCEPTIBLE")
pwatch_AMR_df$Resistance <- factor(pwatch_AMR_df$Resistance, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))

pwatch_AMR_plot <- ggplot(data = pwatch_AMR_df, aes(x=Antibiotic, y=n, fill = Resistance)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#F7BF5E", "#B75EF7")) +
  theme_bw(base_size = 18) +
  ylab("# of resistant isolates") +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.75))
pwatch_AMR_plot

ggsave("paper assembly/Figures/pwatch_amr.png", pwatch_AMR_plot, height = 8, width = 11, units = "in", dpi = 300)

#Comparison of predictions and MIC
MIC_pub_data <- read.csv("MIC_data/pubmlst_MIC.csv")

#Ciprofloxacin
MIC_CIP <- MIC_pub_data %>%
  select(name, CIP_R) %>%
  filter(!is.na(CIP_R)) %>%
  filter(CIP_R != "I")

pwatch_CIP <- pwatch_AMR %>%
  filter(Genome.Name %in% MIC_CIP$name) %>%
  select(Genome.Name, Ciprofloxacin)
colnames(pwatch_CIP) <- c("name", "Predicted")

CIP_df <- MIC_CIP %>%
  left_join(pwatch_CIP) %>%
  filter(!is.na(Predicted)) %>%
  count(CIP_R, Predicted)
colnames(CIP_df)[1] <- "CIP_MIC"
CIP_df$Predicted <- str_replace_all(CIP_df$Predicted, "NOT_FOUND", "SUSCEPTIBLE")
CIP_df$Predicted <- factor(CIP_df$Predicted, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))
CIP_df[nrow(CIP_df) + 1,] = c("R","SUSCEPTIBLE", 0)
CIP_df$n <- as.numeric(CIP_df$n)

CIP_plot <- ggplot(data = CIP_df, aes(x=CIP_MIC, y=n, fill = Predicted)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#B75EF7")) +
  theme_bw() +
  ylab("# of isolates") +
  xlab("Ciprofloxacin")
CIP_plot
#Tetracycline
MIC_TET <- MIC_pub_data %>%
  select(name, TET_R) %>%
  filter(!is.na(TET_R))
MIC_TET$TET_R <- str_replace_all(MIC_TET$TET_R, "DS", "I")

pwatch_TET <- pwatch_AMR %>%
  filter(Genome.Name %in% MIC_TET$name) %>%
  select(Genome.Name, Tetracycline)
colnames(pwatch_TET) <- c("name", "Predicted")

TET_df <- MIC_TET %>%
  left_join(pwatch_TET) %>%
  filter(!is.na(Predicted)) %>%
  count(TET_R, Predicted)
colnames(TET_df)[1] <- "TET_MIC"
TET_df$Predicted <- str_replace_all(TET_df$Predicted, "NOT_FOUND", "SUSCEPTIBLE")
TET_df$Predicted <- factor(TET_df$Predicted, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))
TET_df[nrow(TET_df) + 1,] = c("S","RESISTANT", 0)
TET_df$n <- as.numeric(TET_df$n)
TET_df$TET_MIC <- factor(TET_df$TET_MIC, levels = c("S", "I", "R"))

TET_plot <- ggplot(data = TET_df, aes(x=TET_MIC, y=n, fill = Predicted)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#F7BF5E", "#B75EF7")) +
  theme_bw() +
  ylab("# of isolates") +
  xlab("Tetracycline")
TET_plot

#Penicillin
MIC_PEN <- MIC_pub_data %>%
  select(name, PEN_R) %>%
  filter(!is.na(PEN_R))
MIC_PEN$PEN_R <- str_replace_all(MIC_PEN$PEN_R, "DS", "I")

pwatch_PEN <- pwatch_AMR %>%
  filter(Genome.Name %in% MIC_PEN$name) %>%
  select(Genome.Name, Penicillin)
colnames(pwatch_PEN) <- c("name", "Predicted")

PEN_df <- MIC_PEN %>%
  left_join(pwatch_PEN) %>%
  filter(!is.na(Predicted)) %>%
  count(PEN_R, Predicted)
colnames(PEN_df)[1] <- "PEN_MIC"
PEN_df$Predicted <- str_replace_all(PEN_df$Predicted, "NOT_FOUND", "SUSCEPTIBLE")
PEN_df$Predicted <- factor(PEN_df$Predicted, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))
PEN_df[nrow(PEN_df) + 1,] = c("I","SUSCEPTIBLE", 0)
PEN_df[nrow(PEN_df) + 1,] = c("R","SUSCEPTIBLE", 0)
PEN_df[nrow(PEN_df) + 1,] = c("S","RESISTANT", 0)
PEN_df$n <- as.numeric(PEN_df$n)
PEN_df$PEN_MIC <- factor(PEN_df$PEN_MIC, levels = c("S", "I", "R"))

PEN_plot <- ggplot(data = PEN_df, aes(x=PEN_MIC, y=n, fill = Predicted)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#F7BF5E", "#B75EF7")) +
  theme_bw() +
  ylab("# of isolates") +
  xlab("Penicillin")
PEN_plot

#Azithromycin
MIC_AZM <- MIC_pub_data %>%
  select(name, AZM_R) %>%
  filter(!is.na(AZM_R))
MIC_AZM$AZM_R <- str_replace_all(MIC_AZM$AZM_R, "DS", "R")

pwatch_AZM <- pwatch_AMR %>%
  filter(Genome.Name %in% MIC_AZM$name) %>%
  select(Genome.Name, Azithromycin)
colnames(pwatch_AZM) <- c("name", "Predicted")

AZM_df <- MIC_AZM %>%
  left_join(pwatch_AZM) %>%
  filter(!is.na(Predicted)) %>%
  count(AZM_R, Predicted)
colnames(AZM_df)[1] <- "AZM_MIC"
AZM_df$Predicted <- str_replace_all(AZM_df$Predicted, "NOT_FOUND", "SUSCEPTIBLE")
AZM_df$Predicted <- factor(AZM_df$Predicted, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))
AZM_df[nrow(AZM_df) + 1,] = c("S","RESISTANT", 0)
AZM_df$n <- as.numeric(AZM_df$n)

AZM_plot <- ggplot(data = AZM_df, aes(x=AZM_MIC, y=n, fill = Predicted)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#B75EF7")) +
  theme_bw() +
  ylab("# of isolates") +
  xlab("Azithromycin")
AZM_plot

#Cefixime
MIC_CFM <- MIC_pub_data %>%
  select(name, CFM_R) %>%
  filter(!is.na(CFM_R))
MIC_CFM$CFM_R <- str_replace_all(MIC_CFM$CFM_R, "DS", "R")

pwatch_CFM <- pwatch_AMR %>%
  filter(Genome.Name %in% MIC_CFM$name) %>%
  select(Genome.Name, Cefixime)
colnames(pwatch_CFM) <- c("name", "Predicted")

CFM_df <- MIC_CFM %>%
  left_join(pwatch_CFM) %>%
  filter(!is.na(Predicted)) %>%
  count(CFM_R, Predicted)
colnames(CFM_df)[1] <- "CFM_MIC"
CFM_df$Predicted <- str_replace_all(CFM_df$Predicted, "NOT_FOUND", "SUSCEPTIBLE")
CFM_df$Predicted <- factor(CFM_df$Predicted, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))
CFM_df[nrow(CFM_df) + 1,] = c("R","SUSCEPTIBLE", 0)
CFM_df$n <- as.numeric(CFM_df$n)

CFM_plot <- ggplot(data = CFM_df, aes(x=CFM_MIC, y=n, fill = Predicted)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#B75EF7")) +
  theme_bw() +
  ylab("# of isolates") +
  xlab("Cefixime")
CFM_plot

#Ceftriaxone
MIC_CRO <- MIC_pub_data %>%
  select(name, CRO_R) %>%
  filter(!is.na(CRO_R))
MIC_CRO$CRO_R <- str_replace_all(MIC_CRO$CRO_R, "DS", "R")

pwatch_CRO <- pwatch_AMR %>%
  filter(Genome.Name %in% MIC_CRO$name) %>%
  select(Genome.Name, Ceftriaxone)
colnames(pwatch_CRO) <- c("name", "Predicted")

CRO_df <- MIC_CRO %>%
  left_join(pwatch_CRO) %>%
  filter(!is.na(Predicted)) %>%
  count(CRO_R, Predicted)
colnames(CRO_df)[1] <- "CRO_MIC"
CRO_df$Predicted <- str_replace_all(CRO_df$Predicted, "NOT_FOUND", "SUSCEPTIBLE")
CRO_df$Predicted <- factor(CRO_df$Predicted, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"))
CRO_df[nrow(CRO_df) + 1,] = c("S","RESISTANT", 0)
CRO_df$n <- as.numeric(CRO_df$n)

CRO_plot <- ggplot(data = CRO_df, aes(x=CRO_MIC, y=n, fill = Predicted)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#5EF79E", "#B75EF7")) +
  theme_bw() +
  ylab("# of isolates") +
  xlab("Ceftriaxone")
CRO_plot

#bring the plots together
ggarrange(PEN_plot, TET_plot, CIP_plot, AZM_plot, CFM_plot, CRO_plot, 
          ncol = 2, nrow = 3, common.legend = T, legend = "right")

#Tetracycline resistance in group B
TET_cladeB <- isolate_metadata %>%
  filter(group == "B") %>%
  count(Tetracycline)

#PathogenWatch resistance mechanisms for clades + AMR
isolate_metadata <- read.csv("isolate_metadata.csv")
pwatch_snps <- read.csv("pwatch/amr_snps_genes.csv") %>% 
  rename(Run = Genome.Name)
popnet_clades <- isolate_metadata %>%
  select(Run, group, Azithromycin, Cefixime, Ceftriaxone, Penicillin, Ciprofloxacin, Tetracycline)
popnet_pwatch_snps <- left_join(popnet_clades, pwatch_snps)

#Azithromycin related SNPs
AZM_snps <- popnet_pwatch_snps %>%
  filter(Azithromycin == "RESISTANT")


#missing pathogenwatch annotations
pwatch_na <- popnet_pwatch_snps %>%
  filter(is.na(Tetracycline))
write.table(pwatch_na$Run, "missing_pwatch.txt", quote = F, row.names = F)

# Isolate stats and map -------------------------------------------------------------
isolate_metadata <- read.csv("isolate_metadata.csv")

#Isolates by year
isolate_metadata %>%
  count(Year_bin)

#isolates per country
country_isols <- isolate_metadata %>%
  count(Country) %>%
  filter(!is.na(Country)) %>%
  select(region = Country, Isolates = n) %>%
  mutate(region = recode(str_trim(region), "Hong Kong" = "China",
                         "United Kingdom" = "UK"))

world <- map_data("world")

map_df <- left_join(world, country_isols, relationship = "many-to-many")

map_ngo <- map_df %>%
  filter(!is.na(Isolates))

ngo_map <- ggplot() + 
  coord_fixed(1.3) +
  geom_polygon(data = map_df,aes(x = long, y = lat, group = group),
               fill = "lightgrey") +
  geom_polygon(data = map_ngo, aes(x = long, y = lat, group = group,
                                            fill = Isolates),
               color = "black") +
  scale_fill_gradient(low = "thistle1", high = "purple4",
                      na.value = "lightgrey") +
  theme(axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "white"))

ggsave("paper assembly/Figures/ngo_map.pdf", ngo_map, width = 8, height = 5, units = "in")

#Year and continent
cont_year <- isolate_metadata %>%
  count(Continent, Year) %>%
  filter(Year != 1905)
cont_year_bubble <-ggplot(data = cont_year, aes(x=Year, y=Continent, size=n, color = Continent)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none",
         size=guide_legend(title = "Isolates")) +
  theme_bw(base_size = 18)
  
cont_year_bubble

ggsave("paper assembly/Figures/cont_year_bubble.pdf", cont_year_bubble, height = 6, width = 12, units = "in")
ggsave("paper assembly/Figures/cont_year_bubble.png", cont_year_bubble, dpi = 300, height = 6, width = 12, units = "in")
#sex per continent
continent_sex <- isolate_metadata %>%
  count(Continent, Sex)
#site per continent
continent_site <- isolate_metadata %>%
  count(Continent, Site)
#isolates per continent
continent_count <- isolate_metadata %>%
  count(Continent)

#Let's do the same for the initial set of 30k
initial_metadata <- read.csv("SRA_metadata_final.csv")

#get continent from country
initial_metadata$Continent <- countrycode(sourcevar = initial_metadata$Country, origin = "country.name",
                                          destination = "region")
initial_metadata <- initial_metadata %>%
  mutate(Continent = case_when(Country == "Australia" ~ "Oceania",
                               Country == "New Zealand" ~ "Oceania",
                               Continent == "East Asia & Pacific" ~ "Asia",
                               Country == "Azerbaijan" ~ "Asia",
                               Continent == "Europe & Central Asia" ~ "Europe",
                               Continent == "Latin America & Caribbean" ~ "South & Central America",
                               Continent == "Sub-Saharan Africa" ~ "Africa",
                               Country == "Israel" ~ "Asia",
                               WHO == "WHO_L" ~ "Asia",
                               TRUE ~ Continent))
#isolates per country from intial set
init_country_isols <- initial_metadata %>%
  count(Country) %>%
  filter(!is.na(Country)) %>%
  select(region = Country, Isolates = n) %>%
  mutate(region = recode(str_trim(region), "Hong Kong" = "China",
                         "United Kingdom" = "UK"))

world <- map_data("world")

init_map_df <- left_join(world, init_country_isols, relationship = "many-to-many")

init_map_ngo <- init_map_df %>%
  filter(!is.na(Isolates))

ggplot() + 
  coord_fixed(1.3) +
  geom_polygon(data = init_map_df,aes(x = long, y = lat, group = group),
               fill = "lightgrey") +
  geom_polygon(data = init_map_ngo, aes(x = long, y = lat, group = group,
                                   fill = Isolates),
               color = "black") +
  scale_fill_gradient(low = "lightcyan", high = "royalblue",
                      na.value = "lightgrey") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"))

#Year 
init_year_hist <- ggplot(data = initial_metadata) +
  geom_histogram(aes(x=Year), binwidth = 1,
                 color = "black", fill = "grey") +
  theme_bw() +
  theme(axis.title.y = element_blank())
init_year_hist

#Year and continent
init_cont_year <- initial_metadata %>%
  count(Continent, Year) %>%
  filter(Year != 1905) %>%
  filter(!is.na(Continent))
init_cont_year_bubble <-ggplot(data = init_cont_year, aes(x=Year, y=Continent, size=n, color = Continent)) +
  geom_point() +
  scale_fill_brewer(palette = "Dark2") +
  guides(color = "none",
         size=guide_legend(title = "Isolates")) +
  theme_bw() 
init_cont_year_bubble

#initial continent count
init_continent_count <- initial_metadata %>%
  count(Continent)

#compare continent distribution between original and reference set
init_continent_share <- init_continent_count %>%
  mutate(init_share = n/sum(n)*100) %>%
  filter(!is.na(Continent)) %>%
  select(Continent, Initial = init_share)

continent_share <- continent_count %>%
  mutate(Reference = n/sum(n)*100) %>%
  select(Continent, Reference)

continent_compare_df <- continent_share %>%
  left_join(init_continent_share) %>%
  pivot_longer(!Continent, names_to = "Set", values_to = "Share")

cont_compare <- ggplot(data = continent_compare_df, aes(x=Continent, y=Share, fill=Set)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() +
  scale_fill_manual(values = c("darkorange", "lightblue3")) +
  ylab("% of isolates") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave("paper assembly/Figures/cont_compare.pdf", height = 6, width = 8, units = "in")

#Compare year distribution for initial vs ref
initial_year <- initial_metadata %>%
  select(Year)
initial_year$Set <- "Initial"
ref_year <- isolate_metadata %>%
  select(Year)
ref_year$Set <- "Reference"
year_comp <- rbind(initial_year, ref_year) %>%
  filter(Year != 1905)
year_comp$Set <- factor(year_comp$Set)

year_comp_violin <- ggplot(year_comp, aes(x=Set, y=Year, fill = Set)) +
  geom_violin() +
  scale_fill_manual(values = c("darkorange", "lightblue3")) +
  theme_bw()
year_comp_box <- ggplot(year_comp, aes(x=Set, y=Year, fill = Set)) +
  geom_boxplot() +
  scale_fill_manual(values = c("darkorange", "lightblue3")) +
  theme_bw()
ggsave("paper assembly/Figures/year_comp_violin.pdf", year_comp_violin, height = 6, width = 6, units = "in")
ggsave("paper assembly/Figures/year_comp_box.pdf", year_comp_box, height = 6, width = 6, units = "in")

# Phase variable genes ----------------------------------------------------

#Read in perf predictions and match to gene IDs
perf_predictions <- read.delim("phase/perf_predictions.tsv", header = F) %>%
  dplyr::select(Gene = V9, Repeat_class = V4, Repeat_length = V5, Num_repeats = V7, Actual_repeat = V8, Region = V13, Type = V14) %>%
  distinct()
ssr_genes <- read.csv("panaroo_results/gene_presence_absence.csv") %>%
  dplyr::select(!Non.unique.Gene.name) %>%
  pivot_longer(cols = -c("Gene", "Annotation"), names_to = "Isolate") %>%
  dplyr::select(Gene_name = Gene, Annotation, Isolate, Gene = value) %>%
  separate_longer_delim(Gene, ";") %>%
  inner_join(perf_predictions)

#let's try and figure out which ones might actually be PV
gene_meta <- read.csv("gene_metadata.csv")
isolate_metadata <- read.csv("isolate_metadata.csv")
PopNet_clades <- isolate_metadata %>%
  select(Isolate = Run, Clade = group)

no_rare_repeats <- ssr_genes %>%
  select(!c("Region", "Type")) %>%
  group_by(Gene_name) %>%
  count(Gene_name, Repeat_class) %>%
  filter(n >= 25) %>%
  distinct(Gene_name, Repeat_class) %>%
  mutate(Combo = paste0(Gene_name, "_", Repeat_class))

varied_repeats <- ssr_genes %>%
  select(!c("Region", "Type")) %>%
  mutate(Combo = paste0(Gene_name, "_", Repeat_class)) %>%
  filter(Combo %in% no_rare_repeats$Combo) %>%
  select(!Combo) %>%
  distinct(Gene_name, Repeat_class, Num_repeats) %>%
  count(Gene_name, Repeat_class)

potential_pv <- varied_repeats %>%
  filter(n >= 2) %>%
  rename(Gene = Gene_name) %>%
  left_join(gene_meta)

pv_orthologs <- potential_pv %>%
  ungroup() %>%
  count(seed_ortholog)

pv_combos <- potential_pv %>%
  unite("pv_combo", c("Gene", "Repeat_class")) %>%
  select(pv_combo)

potential_pv_genes <- ssr_genes %>%
  unite("pv_combo", c("Gene_name", "Repeat_class"), remove = F) %>%
  filter(pv_combo %in% pv_combos$pv_combo) %>%
  select(!pv_combo) %>%
  left_join(PopNet_clades)

pv_genes_repeats <- potential_pv_genes %>%
  select(Gene_name, Annotation, Repeat_class, Repeat_length, Num_repeats, Actual_repeat, Region) %>%
  distinct()

pv_repeat_count <- potential_pv_genes %>%
  count(Gene_name, Annotation, Repeat_class, Num_repeats)

write.csv(potential_pv, "potential_pv_genes.csv", quote = F, row.names = F)

write.csv(pv_genes_repeats, "pv_repeat_distributions.csv", quote = F, row.names = F)

#Take a look at the MS11 passage experiment
MS11_predictions <- read.delim("phase/MS11_passaging/MS11_perf.tsv", header = F) %>%
  dplyr::select(Gene = V9, Repeat_class = V4, Repeat_length = V5, Num_repeats = V7, Actual_repeat = V8, Region = V13, Type = V14) %>%
  distinct()

MS11_ssr_genes <- read.csv("phase/MS11_passaging/panaroo_results/gene_presence_absence.csv") %>%
  dplyr::select(!Non.unique.Gene.name) %>%
  pivot_longer(cols = -c("Gene", "Annotation"), names_to = "Isolate") %>%
  dplyr::select(Gene_name = Gene, Annotation, Isolate, Gene = value) %>%
  separate_longer_delim(Gene, ";") %>%
  inner_join(MS11_predictions)

MS11_varied_repeats <- MS11_ssr_genes %>%
  select(!c("Region", "Type", "Actual_repeat")) %>%
  group_by(Gene_name) %>%
  distinct(Repeat_class, Repeat_length, Num_repeats) %>%
  count(Gene_name, Repeat_class) %>%
  filter(n >= 2)

Ms11_potential_PV <- MS11_ssr_genes %>%
  filter(Gene_name %in% MS11_varied_repeats$Gene_name) %>%
  filter(Annotation != "hypothetical protein")


#Results of capillary gel electrophoresis
setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC/phase/MS11_passaging/cpge_10/Results")

sample_names <- list.files(pattern="\\.csv$") %>%
  str_remove("\\-\\w*.csv") %>%
  unique()

process_cpge <- function(sample_id) {
  A <- read.csv(paste0(sample_id,"-A.csv")) %>%
    filter(Dye == "Blue") %>%
    select(Size, Height) %>%
    mutate(Sample = paste0(sample_id, "-A")) %>%
    filter(Height >= 0.1*sum(Height)) %>%
    filter(Size > 200) %>%
    mutate(Ratio = Height/sum(Height))
  B <- read.csv(paste0(sample_id,"-B.csv")) %>%
    filter(Dye == "Blue") %>%
    select(Size, Height) %>%
    mutate(Sample = paste0(sample_id, "-B")) %>%
    filter(Height >= 0.1*sum(Height)) %>%
    filter(Size >= 200) %>%
    mutate(Ratio = Height/sum(Height))
  ratio <- rbind(A,B)
  return(ratio)
}

cpge_ratio <- do.call(rbind.data.frame, lapply(sample_names, process_cpge)) %>%
  mutate(Size = trunc(Size)) %>%
  mutate(Sample = str_remove(Sample, "\\-\\w*")) %>%
  group_by(Sample, Size) %>%
  summarise(mean = mean(Ratio))

write.csv(cpge_ratio, "../cpge_ratios.csv", quote = F, row.names = F)

setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC")

#plot CFUs from Jess
MS11_CFU <- read.csv("phase/MS11_passaging_CFU.csv") %>%
  pivot_longer(cols = c("P1", "P2", "P3"), names_to = "Passage", values_to = "CFU") %>%
  mutate(Passage = factor(Passage)) %>%
  filter(!is.na(CFU))

MS11_CFU_plot <- ggplot(MS11_CFU, aes(x=Passage, y= CFU, fill = Passage, group = Passage)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  scale_y_continuous(trans = "log10") +
  theme_bw(base_size = 18) +
  stat_summary(fun = "mean", geom = "crossbar") +
  scale_fill_brewer(palette = "Dark2")

ggsave("paper assembly/Figures/MS11_CFU_plot.pdf", MS11_CFU_plot, height = 8, width = 8, units = "in")
ggsave("paper assembly/Figures/MS11_CFU_plot.png", MS11_CFU_plot, dpi = 300, height = 8, width = 10, units = "in")

# Tree --------------------------------------------------------------------

#Let's try and make use of itol
isolate_metadata <- read.csv("isolate_metadata.csv")

#PopNet
color_PopNet <- isolate_metadata %>%
  select(Run, color) %>%
  mutate(Run = paste0(Run, "_0"))
write.csv(color_PopNet, "color_strip_PopNet.csv", row.names = F, quote = F)

#MDR
color_MDR <- isolate_metadata %>% 
  select(Run, MDR) %>%
  filter(MDR == "Yes") %>%
  mutate(Run = paste0(Run, "_0")) %>%
  select(Run)
color_MDR$color <- "#005500"
write.csv(color_MDR, "color_strip_MDR.csv", row.names = F, quote = F)  

#ggtree
tree_df <- read.csv("isolate_metadata.csv")
#make matrix for heatmap of tree
add_heatmap <- function(ggtree, tree, mt.df, metadata.col, pal, width, offset, guide='legend',
                        guide_title='NA', ...){
  
  labels <- tree$tip.label
  mt.filt <- mt.df %>% filter(Run %in% labels)
  
  # Make heatmap matrix
  mat <- data.frame(val = mt.filt[[metadata.col]])
  rownames(mat) <- mt.filt$Run
  colnames(mat) <- guide_title
  
  ggtree.scale <- ggtree + ggnewscale::new_scale_fill()
  
  ggtree_with_hm <-
    ggtree::gheatmap(ggtree.scale, mat, width=width,
                     colnames=TRUE, offset=offset, colnames_angle = 90,
                     font.size = 2, colnames_position = 'top', ...) +
    scale_fill_manual(values=pal, guide=guide, name=guide_title)
  
  return(ggtree_with_hm)
  
}

PopNet_colours <- isolate_metadata %>%
  select(group, color) %>%
  distinct() %>%
  filter(!is.na(group)) %>%
  arrange(group)
PopNet_pal <- as.character(PopNet_colours$color)
names(PopNet_pal) <- PopNet_colours$group

continent_pal <- scales::brewer_pal(palette = 'Dark2')(length(unique(tree_df$Continent)))
names(continent_pal) <- unique(tree_df$Continent)

#Let's make the tree with resistance to AMR
MDR_tree_df <- read.csv("isolate_metadata.csv") %>%
  mutate(MDR_drug = case_when(MDR_drug == "Both" ~ "XDR",
                              MDR_drug == "AZM" ~ "MDR-AZM",
                              MDR_drug == "CFM" ~ "MDR-CFM",
                              TRUE ~ MDR_drug)) %>%
  mutate(MDR_drug = factor(MDR_drug, levels = c("XDR", "MDR-AZM", "MDR-CFM", "None")))

MDR_pal <- c("red", "blue", "green", "grey70")
names(MDR_pal) <- levels(MDR_tree_df$MDR_drug)
MDR_pal["None"] <- "grey70"

MDR_tree <- read.newick("cfml_continent.output.labelled_tree.newick") %>%
  drop.tip(c("GB3_FA19", "GB2_FA1090", "GB1_MS11"))
MDR_circ_tree <- ggtree(MDR_tree, size = 0.5, layout = "fan", open.angle = 10) +
  geom_tiplab(align=TRUE, linetype = '99191919', color='grey60',
              linesize = 0.1, size=0)
MDR.tree.plot <- add_heatmap(tree = MDR_tree, mt.df = MDR_tree_df, metadata.col = "group", ggtree = MDR_circ_tree,
                         width = 0.075, offset = 0, guide_title = "PopNet", pal = PopNet_pal)
MDR.tree.plot <- add_heatmap(tree = MDR_tree, mt.df = MDR_tree_df, metadata.col = "MDR_drug", ggtree = MDR.tree.plot,
                         width = 0.075, offset = 0.00002, guide_title = "MDR", pal = MDR_pal)
MDR.tree.plot <- add_heatmap(tree = MDR_tree, mt.df = MDR_tree_df, metadata.col = "Continent", ggtree = MDR.tree.plot,
                             width = 0.075, offset = 0.00001, guide_title = "Continent", pal = continent_pal)
MDR.tree.plot
ggsave(plot = MDR.tree.plot, filename = "paper assembly/Figures/ggtree_MDR_cont.pdf",
       width=8.5, height = 10, units = 'in')


#Let's make a tree with the top PopNet clades, cgMLST, and poppunk
top_cgMLST <- isolate_metadata %>%
  count(Ng_cgc_300) %>%
  filter(!is.na(Ng_cgc_300)) %>%
  filter(n >= 10)
top_ppunk <- isolate_metadata %>%
  count(ppunk_cluster) %>%
  filter(!is.na(ppunk_cluster)) %>%
  filter(n >= 10)
  
clade_tree_df <- isolate_metadata %>%
  select(Run, group, Ng_cgc_300, ppunk_cluster) %>%
  mutate(Ng_cgc_300 = factor(Ng_cgc_300)) %>%
  mutate(ppunk_cluster = factor(ppunk_cluster)) %>%
  mutate(Ng_cgc_300 = case_when(!Ng_cgc_300 %in% top_cgMLST$Ng_cgc_300 ~ "Other",
                                TRUE ~ Ng_cgc_300)) %>%
  mutate(ppunk_cluster = case_when(!ppunk_cluster %in% top_ppunk$ppunk_cluster ~ "Other",
                                TRUE ~ ppunk_cluster)) %>%
  mutate(Ng_cgc_300 = factor(Ng_cgc_300)) %>%
  mutate(ppunk_cluster = factor(ppunk_cluster))
  

cgMLST_pal <- scales::brewer_pal(palette = 'Paired')(length(unique(clade_tree_df$Ng_cgc_300)))
names(cgMLST_pal) <- levels(clade_tree_df$Ng_cgc_300)
cgMLST_pal["Other"] <- "grey95"

ppunk_pal <- scales::brewer_pal(palette = 'Set3')(length(unique(clade_tree_df$ppunk_cluster)))
names(ppunk_pal) <- levels(clade_tree_df$ppunk_cluster)
ppunk_pal["Other"] <- "grey95"

tree <- read.newick("cfml_continent.output.labelled_tree.newick") %>%
  drop.tip(c("GB3_FA19", "GB2_FA1090", "GB1_MS11"))
circ_tree <- ggtree(tree, size = 0.5, layout = "fan", open.angle = 10) +
  geom_tiplab(align=TRUE, linetype = '99191919', color='grey60',
              linesize = 0.1, size=0)
tree.plot <- add_heatmap(tree = tree, mt.df = tree_df, metadata.col = "group", ggtree = circ_tree,
                         width = 0.075, offset = 0, guide_title = "PopNet", pal = PopNet_pal)
tree.plot <- add_heatmap(tree = tree, mt.df = clade_tree_df, metadata.col = "Ng_cgc_300", ggtree = tree.plot,
                         width = 0.075, offset = 0.00001, guide_title = "cgMLST", pal = cgMLST_pal)
tree.plot <- add_heatmap(tree = tree, mt.df = clade_tree_df, metadata.col = "ppunk_cluster", ggtree = tree.plot,
                         width = 0.075, offset = 0.00002, guide_title = "PopPunk", pal = ppunk_pal)
tree.plot

ggsave(plot = tree.plot, filename = "paper assembly/Figures/ggtree_clusters.pdf",
       width=8.5, height = 10, units = 'in')
# pyseer GWAS -------------------------------------------------------------
#make trait file (either binary or continuous) for GWAS testing, lets start with AMR and mobile elements
isolate_metadata <- read.csv("isolate_metadata.csv")

traits <- isolate_metadata %>%
  select(Run, Azithromycin, Cefixime, Ceftriaxone, Penicillin, Ciprofloxacin, 
         Tetracycline, GGI, Blac, cryptic, pConj, pConj_Tet) %>%
  mutate(across(Azithromycin:Tetracycline, ~ case_when(. == "NOT_FOUND" ~ 0,
                                                        . == "RESISTANT" ~ 1,
                                                        . == "INTERMEDIATE" ~ 2))) %>%
  mutate(across(GGI:pConj_Tet, ~ case_when(. == "No" ~ 0,
                                           . == "Yes" ~ 1)))
write.table(traits, "pyseer_traits.tsv", sep = "\t", quote = F, row.names = F)

#process output for combining with snpEff annotations
#Start by reading in and checking out the vcf file
pyseer_vcf <- read.vcfR("gwas/pyseer_fixed_vcf.orf.anno.vcf")

#get annotations 
vcf_anno <- INFO2df(pyseer_vcf) %>%
  select(ANN)
vcfr_tidy <- vcfR2tidy(pyseer_vcf)

vcf_gt <- extract.gt(pyseer_vcf)

vcf_df <- data.frame(vcfr_tidy$fix$CHROM)
vcf_df$POS <- as.character(vcfr_tidy$fix$POS)
vcf_df$REF <- vcfr_tidy$fix$REF
vcf_df$ALT <- vcfr_tidy$fix$ALT
vcf_df$ANN <- vcf_anno$ANN

vcf_df <- vcf_df %>%
  separate(ANN, into = c("Allele", "Variant_ann", "Impact", "Gene_name", "Locus_tag", "Feature", 
                         "Feature_ID", "Transcript", "Rank_total", "Ntd_variant", "Prot_variant",
                         "cDNA_pos", "cDS_pos", "Prot_pos"), sep = "\\|")
colnames(vcf_df)[1] <- "CHROM"

write.csv(vcf_df, "gwas/snpEff_annotations.csv", quote = F, row.names = F)

#match gene metadata with loci
ref_headers <- read.delim("gwas/ref_headers.txt", sep = " ", header = F) %>%
  filter(V2 != "") %>%
  select(Gene = V1, Locus_tag = V2)
blast_results <- read.delim(("gwas/ref_mapping.tsv"), header = F) %>%
  select(Gene_name = V1, Gene = V2)
ref_mapping <- blast_results %>%
  left_join(ref_headers) %>%
  select(!Gene)
colnames(ref_mapping)[1] <- "Gene"
gene_meta <- read.csv("gene_metadata.csv") %>%
  left_join(ref_mapping)

write.csv(gene_meta, "gene_metadata.csv", quote = F, row.names = F)

#integrate gwas results with snpEff annotations and gene metadata
vcf_df <- read.csv("gwas/snpEff_annotations.csv")
gene_meta <- read.csv("gene_metadata.csv")

process_gwas <- function(in_file) {
  gwas_df <- read.delim(file = paste0("gwas/results/", in_file)) %>%
    mutate(log10pval = -log10(lrt.pvalue)) %>%
    filter(variant != "") %>%
    separate(variant, into = c("CHROM1", "CHROM2","POS", "REF", "ALT"), sep = "_") %>%
    unite(c("CHROM1", "CHROM2"), sep = "_", col = "CHROM") %>%
    mutate(POS = as.numeric(POS))
  gwas_ann <- gwas_df %>%
    left_join(vcf_df) %>%
    left_join(gene_meta, relationship = "many-to-many") %>%
    distinct()
  write.csv(gwas_ann, paste0("gwas/results/", str_replace(in_file, ".txt", ""),"_ann.csv"), quote = F, row.names = F)
}

gwas_files <- list.files("gwas/results/")

sapply(gwas_files, process_gwas)

#Plot gwas results for CIP
CIP_gwas <- read.csv("gwas/results/sig_Ciprofloxacin_gwas_ann.csv")

ggplot(CIP_gwas, aes(x=beta, y=log10pval, colour = variant_h2, size = af, label = Gene_name)) +
  geom_point(alpha = 0.5) + 
  geom_text_repel(aes(size=60), show.legend = FALSE, colour='black')

#Look for genes in common across mobile elements or AMR
#Mobile elements
pConj_Tet_gwas <- read.csv("gwas/results/sig_pConj_Tet_gwas_ann.csv")
pConj_gwas <- read.csv("gwas/results/sig_pConj_gwas_ann.csv")
cryptic_gwas <- read.csv("gwas/results/sig_cryptic_gwas_ann.csv")
Blac_gwas <- read.csv("gwas/results/sig_Blac_gwas_ann.csv")

intersect(pConj_Tet_gwas$Gene_name, pConj_gwas$Gene_name)
intersect(intersect(pConj_gwas$Gene_name, Blac_gwas$Gene_name), pConj_Tet_gwas$Gene_name)

intersect(Blac_gwas$Gene_name, cryptic_gwas$Gene_name)

#Manhattan plots for visualization

#Let's make some plots for publication highlighting regions
#AZM
AZM_gwas_raw <- read.delim("gwas/raw/Azithromycin_gwas.txt") %>%
  mutate(log10pval = -log10(lrt.pvalue)) %>%
  filter(variant != "") %>%
  separate(variant, into = c("CHROM1", "CHROM2","POS", "REF", "ALT"), sep = "_", remove = F) %>%
  mutate(CHROM1 = 1) %>%
  select(!CHROM2) %>%
  mutate(POS = as.numeric(POS)) %>%
  left_join(vcf_df) %>%
  left_join(gene_meta, relationship = "many-to-many") %>%
  distinct() %>%
  filter(lrt.pvalue <= 0.05) %>%
  mutate(is_highlight = case_when(POS > 1260000 & POS < 1280000 ~ "mtr",
                                  POS > 1810000 & POS < 1820000 ~ "mqo",
                                  TRUE ~ "NO")) %>%
  mutate(is_annotate = ifelse(log10pval>8 & Impact != "LOW" & is_highlight != "NO", "YES", "NO")) %>%
  distinct(POS, Prot_variant, .keep_all = T)

AZM_manhattan <- ggplot(AZM_gwas_raw, aes(x = POS, y = log10pval)) +
  geom_point(aes(color=is_highlight), alpha = 0.8, size = 1.3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(1.3, 17)) +
  geom_label_repel(data = subset(AZM_gwas_raw, is_annotate=="YES"), aes(label=Prot_variant), size = 2) +
  scale_color_manual(values = c("red", "yellow", "darkgrey"), name = "Gene") +
  geom_hline(yintercept = -log10(7.19E-7), color = "darkred", linewidth = 1.1) +
  theme_classic(base_size = 18) +
  xlab("Base") +
  ylab("-log10pval")
AZM_manhattan
ggsave(AZM_manhattan, file = "paper assembly/Figures/AZM_manhattan.pdf", height = 4, width = 10,
       units = "in", dpi = 300)
ggsave(AZM_manhattan, file = "paper assembly/Figures/AZM_manhattan.png", height = 4, width = 10,
       units = "in", dpi = 300)

#TET
TET_gwas_raw <- read.delim("gwas/raw/Tetracycline_gwas.txt") %>%
  mutate(log10pval = -log10(lrt.pvalue)) %>%
  filter(variant != "") %>%
  separate(variant, into = c("CHROM1", "CHROM2","POS", "REF", "ALT"), sep = "_", remove = F) %>%
  mutate(CHROM1 = 1) %>%
  select(!CHROM2) %>%
  mutate(POS = as.numeric(POS)) %>%
  left_join(vcf_df) %>%
  left_join(gene_meta, relationship = "many-to-many") %>%
  distinct() %>%
  filter(lrt.pvalue <= 0.05) %>%
  mutate(is_highlight = case_when(POS > 1967000 & POS < 1968000 ~ "rpsJ",
                                  POS > 52700 & POS < 58500 ~ "pilC-apbC",
                                  POS > 1920000 & POS < 1940000 ~ "pyrB-rsmB",
                                  TRUE ~ "NO")) %>%
  mutate(is_annotate = ifelse(log10pval>8 & Impact != "LOW" & is_highlight != "NO", "YES", "NO")) %>%
  distinct(POS, Prot_variant, .keep_all = T)

TET_manhattan <- ggplot(TET_gwas_raw, aes(x = POS, y = log10pval)) +
  geom_point(aes(color=is_highlight), alpha = 0.8, size = 1.3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(1.3, 17)) +
  geom_label_repel(data = subset(TET_gwas_raw, is_annotate=="YES"), aes(label=Prot_variant), size = 2) +
  scale_color_manual(values = c("darkgrey", "navyblue", "purple4", "darkgreen"), name = "Gene") +
  geom_hline(yintercept = -log10(7.19E-7), color = "darkred", linewidth = 1.1) +
  theme_classic() +
  xlab("Base") +
  ylab("-log10pval")

TET_manhattan

ggsave(TET_manhattan, file = "paper assembly/Figures/TET_manhattan.pdf", height = 4, width = 10,
       units = "in", dpi = 300)
#PEN
PEN_gwas_raw <- read.delim("gwas/raw/Penicillin_gwas.txt") %>%
  mutate(log10pval = -log10(lrt.pvalue)) %>%
  filter(variant != "") %>%
  separate(variant, into = c("CHROM1", "CHROM2","POS", "REF", "ALT"), sep = "_", remove = F) %>%
  mutate(CHROM1 = 1) %>%
  select(!CHROM2) %>%
  mutate(POS = as.numeric(POS)) %>%
  left_join(vcf_df) %>%
  left_join(gene_meta, relationship = "many-to-many") %>%
  distinct() %>%
  filter(lrt.pvalue <= 0.05) %>%
  mutate(is_highlight = case_when(POS > 1463000 & POS < 1465000 ~ "penA",
                                  POS > 1462700 & POS < 1462900 ~ "murE",
                                  POS > 64000 & POS < 65000 ~ "ychF",
                                  TRUE ~ "NO")) %>%
  mutate(is_annotate = ifelse(log10pval>8 & Impact != "LOW" & is_highlight != "NO", "YES", "NO")) %>%
  distinct(POS, Prot_variant, .keep_all = T)

PEN_manhattan <- ggplot(PEN_gwas_raw, aes(x = POS, y = log10pval)) +
  geom_point(aes(color=is_highlight), alpha = 0.8, size = 1.3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(1.3, 20)) +
  geom_label_repel(data = subset(PEN_gwas_raw, is_annotate=="YES"), aes(label=Prot_variant), size = 2) +
  scale_color_manual(values = c("navyblue", "darkgrey", "purple4", "darkgreen"), name = "Gene") +
  geom_hline(yintercept = -log10(7.19E-7), color = "darkred", linewidth = 1.1) +
  theme_classic() +
  xlab("Base") +
  ylab("-log10pval")

ggsave(PEN_manhattan, file = "paper assembly/Figures/PEN_manhattan.pdf", height = 4, width = 10,
       units = "in", dpi = 300)

#CIP
CIP_gwas_raw <- read.delim("gwas/raw/Ciprofloxacin_gwas.txt") %>%
  mutate(log10pval = -log10(lrt.pvalue)) %>%
  filter(variant != "") %>%
  separate(variant, into = c("CHROM1", "CHROM2","POS", "REF", "ALT"), sep = "_", remove = F) %>%
  mutate(CHROM1 = 1) %>%
  select(!CHROM2) %>%
  mutate(POS = as.numeric(POS)) %>%
  left_join(vcf_df) %>%
  left_join(gene_meta, relationship = "many-to-many") %>%
  distinct() %>%
  filter(lrt.pvalue <= 0.05) %>%
  mutate(is_highlight = case_when(POS > 690000 & POS < 700000 ~ "gyrA",
                                  TRUE ~ "NO")) %>%
  mutate(is_annotate = ifelse(log10pval>8 & Impact != "LOW" & is_highlight != "NO", "YES", "NO")) %>%
  distinct(POS, Prot_variant, .keep_all = T)

CIP_manhattan <- ggplot(CIP_gwas_raw, aes(x = POS, y = log10pval)) +
  geom_point(aes(color=is_highlight), alpha = 0.8, size = 1.3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(1.3, 170)) +
  scale_color_manual(values = c("navyblue", "darkgrey"), name = "Gene") +
  geom_hline(yintercept = -log10(7.19E-7), color = "darkred", linewidth = 1.1) +
  theme_classic(base_size = 18) +
  xlab("Base") +
  ylab("-log10pval")

CIP_manhattan

ggsave(CIP_manhattan, file = "paper assembly/Figures/CIP_manhattan.pdf", height = 4, width = 10,
       units = "in", dpi = 300)

# dNdS --------------------------------------------------------------------
#process json files from BUSTED hyphy
setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC/dNdS/BUSTED_json")
files <- list.files(pattern = "BUSTED.json")
selected_genes <- data.frame(NULL)

test_json <- read_json("dNdS/BUSTED_json/group_204_clean.fasta.BUSTED.json")

BUSTED_process <- function(in_file) {
  BUSTED_json <- read_json(in_file)
  BUSTED_df <- data.frame(t(matrix(unlist(BUSTED_json$`test results`))))
  colnames(BUSTED_df) <- c("LRT", "pval")
  BUSTED_df$dNdS <- unlist(BUSTED_json$fits["MG94xREV with separate rates for branch sets"]$`MG94xREV with separate rates for branch sets`$`Rate Distributions`$`non-synonymous/synonymous rate ratio for *test*`)[1]
  BUSTED_df$Gene <- str_extract(basename(as.character(BUSTED_json$input[1])), ".+(?=_clean.fasta)")
  BUSTED_df <- BUSTED_df %>%
    mutate(Significant = case_when(pval < 0.05 ~ TRUE,
                                   pval >= 0.05 ~ FALSE)) %>%
    relocate(Gene)
  selected_genes <<- selected_genes %>%
    rbind(BUSTED_df)
}

lapply(files, BUSTED_process)

setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC")
write.csv(selected_genes, "dNdS/BUSTED_selected_genes.csv", quote = F, row.names = F)

#merge BUSTED selection with metadata
selected_genes <- read.csv("dNdS/BUSTED_selected_genes.csv")
gene_metadata <- read.csv("gene_metadata.csv")
COG_categories <- read.csv("COG_categories.csv")


#subcellular localizations
pSort <- read.delim("dNdS/psort.txt") %>%
  mutate(Final_Localization = case_when(Final_Localization == "Unknown" & Cytoplasmic_Score >= 5 ~ "Likely_Cytoplasmic",
                                        Final_Localization == "Unknown" & CytoplasmicMembrane_Score >= 5 ~ "Likely_CytoMembrane",
                                        Final_Localization == "Unknown" & Periplasmic_Score >= 5 ~ "Likely_Periplasmic",
                                        Final_Localization == "Unknown" & OuterMembrane_Score >= 5 ~ "Likely_OM",
                                        Final_Localization == "Unknown" & Extracellular_Score >= 5 ~ "Likely_EC",
                                        TRUE ~ Final_Localization)) %>%
  select(Gene=SeqID, Location=Final_Localization) %>%
  mutate(Gene = str_trim(Gene))
pSort_locations <- pSort %>%
  count(Location)

gene_selection <- left_join(gene_metadata, selected_genes) %>%
  left_join(COG_categories) %>%
  left_join(pSort) %>%
  filter(Pangenome == "Core") %>%
  mutate(Significant = case_when(Pangenome == "Core" & is.na(Significant) ~ FALSE,
                                 TRUE ~ Significant)) %>%
  distinct()

core_genes <- gene_metadata %>%
  count(Pangenome)

positive_selection <- gene_selection %>%
  filter(Significant == TRUE) %>%
  distinct()

COG_selection <- gene_selection %>%
  filter(!is.na(COG_category)) %>%
  count(COG_category, Significant)

pSort_selection <- gene_selection %>%
  filter(Location != "Unknown") %>%
  count(Location, Significant)
  
pathway_selection <- gene_selection %>%
  count(KEGG_Pathway, Significant)

module_selection <- gene_selection %>%
  count(KEGG_Module, Significant)

#pSort selection vis
#percent of genes per category
pSort_selection_vis <- pSort_selection %>%
  filter(str_detect(Location, "Likely", negate = T)) %>%
  group_by(Location) %>%
  mutate(percent = 100*n/sum(n)) %>%
  filter(Significant == TRUE)
ggplot(data = pSort_selection_vis) +
  geom_bar(aes(x=Location, y=percent), stat = "identity")

#Clusterprofiler enrichment of locations
test_genes <- positive_selection$Gene
loc2gene <- gene_selection[,c("Location", "Gene")]
enricher(test_genes, TERM2GENE = loc2gene)
psort_enriched <- enricher(test_genes, TERM2GENE = loc2gene)

psort_enriched

#ClusterProfiler enrichment of COG categories

test_genes <- positive_selection$Gene
COG2gene <- gene_selection[,c("COG_category", "Gene")]
enricher(test_genes, TERM2GENE = COG2gene)
COG_enriched <- enricher(test_genes, TERM2GENE = COG2gene)

COG_enriched

#Let's plot some values from busted on a dotplot or scatter plot
dnds_plot_df <- gene_selection %>%
  select(Gene, Significant, LRT, pval, dNdS, COG_category, Location) %>%
  filter(!is.na(pval)) %>%
  filter(!is.na(Location)) %>%
  filter(Location %in% c("Cytoplasmic", "CytoplasmicMembrane", "Extracellular", "Periplasmic", "OuterMembrane"))

#Let's start with boxplot for Location

dnds_box_location <- ggplot(dnds_plot_df, aes(x=Location, y=dNdS, fill = Location)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette = "Set3")

dnds_box_location
ggsave("paper assembly/Figures/dnds_box_location.pdf", dnds_box_location, height = 8, width = 8, units = "in")

#Let's take a look at AMR genes
AMR_selection <- gene_selection %>%
  filter(AMR_name != "") %>%
  select(Gene, AMR_name, AMR_class, LRT, pval, Significant)

#Compare to phase variable gene list
PV_genes <- read.csv("potential_pv_genes.csv") %>%
  distinct()
PV_pos <- positive_selection %>%
  filter(Gene %in% PV_genes$Gene) %>%
  distinct()

no_PV_select <- gene_metadata %>%
  filter(Pangenome == "Core") %>%
  filter(!Gene %in% PV_genes$Gene) %>%
  filter(!Gene %in% positive_selection$Gene)

#process json files from fubar hyphy
setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC/dNdS/FUBAR_json")

files <- list.files(pattern = "FUBAR.json")
selected_sites <- data.frame(NULL)

fubar_process <- function(in_file) {
  fubar_json <- read_json(in_file)
  fubar_df <- data.frame(matrix(unlist(fubar_json$MLE$content$`0`), nrow = length(fubar_json$MLE$content$`0`), byrow = T)) %>%
    select(!(c("X7", "X8")))
  colnames(fubar_df) <- unlist(fubar_json$MLE$headers)[c(1,3,5,7,9,11)]
  fubar_df <- fubar_df %>%
    rownames_to_column("Site")
  fubar_df$Gene <- str_extract(basename(as.character(fubar_json$input[1])), ".+(?=_clean.fasta)")
  gene_sites <- fubar_df %>%
    filter(`Prob[alpha<beta]` > 0.9 | `Prob[alpha>beta]` > 0.9) %>%
    mutate(Selection = case_when(`Prob[alpha>beta]` > 0.9 ~ "Negative",
                                 `Prob[alpha<beta]` > 0.9 ~ "Positive")) %>%
    relocate(Gene)
  selected_sites <<- selected_sites %>%
    rbind(gene_sites)
}

lapply(files, fubar_process)
setwd("C:/Users/dunca/OneDrive - University of Toronto/Documents/Projects/Ngoclustering/New_test_IPNC")
write.csv(selected_sites, "dNdS/FUBAR_selected_sites.csv", quote = F, row.names = F)

#selected sites
selected_sites <- read.csv("dNdS/FUBAR_selected_sites.csv")
#sites per gene
selected_genes <- selected_sites %>%
  count(Gene, Selection)
positive_sites <- selected_sites %>%
  filter(Selection == "Positive")

positive_genes <- selected_genes %>%
  filter(Selection == "Positive")
negative_genes <- selected_genes %>%
  filter(Selection == "Negative")

neg_only_genes <- negative_genes %>%
  filter(!Gene %in% positive_genes$Gene)
quantile(positive_genes$n, probs = seq(0.9,1,0.01))

#check distribution
descdist(positve_genes$n, discrete = F)

gene_meta <- read.csv("gene_metadata.csv")

plot(density(selected_sites$beta.alpha))

#let's take a look at sites under positive selection with posterior probabilities > 0.9 and beta-alpha of >95th percentile (30.775958)
highly_selected <- selected_sites %>%
  filter(beta.alpha >= 30.775958) %>%
  filter(Prob.alpha.beta..1 >= 0.9)

length(unique(highly_selected$Gene))
gene_metadata <- read.csv("gene_metadata.csv")

highly_selected <- highly_selected %>%
  left_join(gene_metadata)

#FUBAR sites for AMR genes
selected_AMR <- positive_selection %>%
  filter(AMR_name != "")
selected_AMR_sites <- selected_sites %>%
  filter(Gene %in% selected_AMR$Gene)

porB_pos_select <- selected_AMR_sites %>%
  filter(Gene == "group_340") %>%
  filter(Selection == "Positive")

# PopNet + Graphia visualization ------------------------------------------
#load in cytoscape nodes and metadata
cyto_node <- read.csv("PopNet/node_file.csv")
isolate_metadata <- read.csv("isolate_metadata.csv")
colnames(isolate_metadata)[1] <- "name"

#read in and parse graphia layout
graphia_layout <- read_json("panaroo_results/graphia/PopNet-node-positions.json")
graphia_layout <- read_json("panaroo_results/graphia/PopNet-node-positions.json") %>%
  data.frame(matrix(unlist(graphia_layout), 
                          nrow = length(graphia_layout), byrow = T)) %>%
  select(name = X2, x_pos = X3, y_pos = X4)

popnet_graphia <- left_join(cyto_node, graphia_layout) %>%
  mutate(x_pos = as.numeric(x_pos)) %>%
  mutate(x_pos = x_pos*6) %>%
  mutate(y_pos = as.numeric(y_pos)) %>%
  mutate(y_pos = y_pos*6) %>%
  left_join(isolate_metadata, by = "name") %>%
  select(!c("color.y", "group.y")) %>%
  rename(color = color.x, group = group.x)
  

write.csv(popnet_graphia, "popnet_graphia.csv", quote = F, row.names = F)
write.table(popnet_graphia, "popnet_graphia.tsv", quote = F, row.names = F, sep = "\t")

#Pangenome
#read in and parse graphia layout for pangenome
graphia_pg_layout <- read_json("panaroo_results/graphia/isolates_pangenome-node-positions.json")
graphia_pg_layout <- read_json("panaroo_results/graphia/isolates_pangenome-node-positions.json") %>%
  data.frame(matrix(unlist(graphia_pg_layout), 
                    nrow = length(graphia_pg_layout), byrow = T)) %>%
  select(name = X2, x_pos = X3, y_pos = X4)

pangenome_cyto <- left_join(cyto_node, graphia_pg_layout) %>%
  mutate(x_pos = as.numeric(x_pos)) %>%
  mutate(x_pos = x_pos*7) %>%
  mutate(y_pos = as.numeric(y_pos)) %>%
  mutate(y_pos = y_pos*7) %>%
  left_join(isolate_metadata, by = "name") %>%
  select(!c("color.y", "group.y")) %>%
  rename(color = color.x, group = group.x)


write.csv(pangenome_cyto, "pangenome_cytoscape.csv", quote = F, row.names = F)

#Accessory genes
acc_gene_attr <- read.csv("panaroo_results/accessory_genes-attributes.csv") %>%
  select(name = Node.Name, mobile_element, MCL.Cluster)
graphia_acc_layout <- read_json("panaroo_results/graphia/accessory_genes-node-positions.json")
graphia_acc_layout <- read_json("panaroo_results/graphia/accessory_genes-node-positions.json") %>%
  data.frame(matrix(unlist(graphia_acc_layout), 
                    nrow = length(graphia_acc_layout), byrow = T)) %>%
  select(name = X2, x_pos = X3, y_pos = X4)
graphia_acc_layout <- graphia_acc_layout %>%
  filter(name %in% acc_gene_attr$name)
gene_gene_pw <- read.delim("panaroo_results/graphia/pw_genes_genes_pw_sim.txt", header = F)
colnames(gene_gene_pw) <- c("name", "Gene2", "Edge")
acc_genes_pw <- gene_gene_pw %>%
  filter(name %in% graphia_acc_layout$name) %>%
  filter(Gene2 %in% graphia_acc_layout$name)

acc_gene_attr <- read.csv("panaroo_results/accessory_genes-attributes.csv") %>%
  select(name = Node.Name, mobile_element, MCL.Cluster)

acc_genes_cyto <- acc_genes_pw %>%
  left_join(graphia_acc_layout, relationship = "many-to-many") %>%
  mutate(x_pos = as.numeric(x_pos)) %>%
  mutate(y_pos = as.numeric(y_pos)) %>%
  left_join(acc_gene_attr, relationship = "many-to-many")

write.csv(acc_genes_cyto, "panaroo_results/acc_genes_cyto.csv", quote = F,
          row.names = F)

# pubMLST submission ------------------------------------------------------
isolate_metadata <- read.csv("isolate_metadata.csv")
#need to make this pubMLST compatible
pubMLST <- isolate_metadata %>%
  select(isolate = Run, country = Country, year = Year, sex = Sex, source = Site, 
         bioproject_accession = BioProject, biosample_accession = BioSample, run_accession = Run) %>%
  mutate(source = case_when(source == "Cervicovaginal" ~ "female reproductive tract",
                            source == "Disseminated" ~ "blood",
                            source == "Ocular" ~ "eye",
                            source == "Penile" ~ "urethral swab",
                            source == "Pharynx" ~ "throat swab",
                            source == "Rectal" ~ "rectal swab",
                            source == "Urogenital" ~ "",
                            source == "Urogenital_Male" ~ "urethral swab",
                            TRUE ~ NA)) %>%
  mutate(sex = case_when(sex == "Male" ~ "male",
                         sex == "Female" ~ "female",
                         TRUE ~ NA)) %>%
  mutate(across(c(isolate, country, sex, source), ~ case_when(is.na(.) ~ "",
                                        TRUE ~ .))) %>%
  mutate(species = "Neisseria gonorrhoeae")


write.table(pubMLST, "pubMLST_metadata.txt", sep = "\t", quote = F, row.names = F)

# pubMLST clustering ------------------------------------------------------
#isolate_metatada 
isolate_metadata <- read.csv("isolate_metadata.csv")

#read in pubMLST cgMLST clusters
pubMLST_cgMLST <- readxl::read_excel("pubMLST_cgMLST.xlsx")
Ng_300_count <- pubMLST_cgMLST %>%
  count(Ng_cgc_300)

#merge and count
isolate_metadata <- isolate_metadata %>%
  left_join(pubMLST_cgMLST)

write.csv(isolate_metadata, "isolate_metadata.csv", quote = F, row.names = F)
cg300 <- isolate_metadata %>%
  count(Ng_cgc_300, group)

# PopPunk clustering ------------------------------------------------------
#Read in clusters
poppunk_clusters <- read.csv("poppunk/poppunk_out_clusters.csv")
colnames(poppunk_clusters) <- c("Run", "ppunk_cluster")

#Cluster counts
poppunk_counts <- poppunk_clusters %>%
  count(ppunk_cluster)

#add poppunk clusters to metadata
isolate_metadata <- read.csv("isolate_metadata.csv")

isolate_metadata <- left_join(isolate_metadata, poppunk_clusters)

write.csv(isolate_metadata, "isolate_metadata.csv", quote = F, row.names = F)


# Clustering comparisons --------------------------------------------------

isolate_metadata <- read.csv("isolate_metadata.csv")

#raw adjusted Rand index for comparison of clustering
adjustedRandIndex(isolate_metadata$group, isolate_metadata$ppunk_cluster)
adjustedRandIndex(isolate_metadata$group, isolate_metadata$Ng_cgc_300)
adjustedRandIndex(isolate_metadata$Ng_cgc_300, isolate_metadata$ppunk_cluster)

#remove clusters with less than 5 isolates
poppunk_counts <- isolate_metadata %>%
  count(ppunk_cluster)
mlst_count <- isolate_metadata %>%
  count(Ng_cgc_300)

ppunk_singles <- poppunk_counts %>%
  filter(n >= 5)
mlst_singles <- mlst_count %>%
  filter(n >= 5)
remove_singletons <- isolate_metadata %>%
  filter(ppunk_cluster %in% ppunk_singles$ppunk_cluster) %>%
  filter(Ng_cgc_300 %in% mlst_singles$Ng_cgc_300)

adjustedRandIndex(remove_singletons$group, remove_singletons$ppunk_cluster)
adjustedRandIndex(remove_singletons$group, remove_singletons$Ng_cgc_300)
adjustedRandIndex(remove_singletons$Ng_cgc_300, remove_singletons$ppunk_cluster)

#chisq test for association
#mlst
popnet_mlst <- isolate_metadata %>%
  filter(!is.na(group)) %>%
  count(group,Ng_cgc_300)
popnet_mlst_table <- popnet_mlst %>%
  pivot_wider(names_from = Ng_cgc_300, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "group")
popnet_mlst_table <- as.matrix(popnet_mlst_table)

chisq.test(popnet_mlst_table, simulate.p.value = TRUE)

#ppunk
popnet_ppunk <- isolate_metadata %>%
  filter(!is.na(group)) %>%
  count(group,ppunk_cluster)
popnet_ppunk_table <- popnet_ppunk %>%
  pivot_wider(names_from = ppunk_cluster, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "group")
popnet_ppunk_table <- as.matrix(popnet_ppunk_table)

chisq.test(popnet_ppunk_table, simulate.p.value = TRUE)

popnet_count <- isolate_metadata %>%
  count(group)

summary(popnet_count$n)


gene_metadata <- read.csv("gene_metadata.csv")

# Correlation plots from graphia ------------------------------------------
#Let's test how these might look with associations between continent of orign and PopNet clades
PopNet_continent_fisher <- read.csv("graphia_fisher/PopNet-enrichment-group-vs-continent.csv") %>%
  mutate(Representation = case_when(Bonferroni.Adjusted > 0.05 ~ NA,
                                   TRUE ~ Representation)) %>%
  mutate(Bonferroni.Adjusted = case_when(Bonferroni.Adjusted > 0.05 ~ NA,
         TRUE ~ Bonferroni.Adjusted)) %>%
  mutate(Continent = case_when(Continent == "South & Central America" ~ paste0("South &", '\n', "Central America"),
                               TRUE ~ Continent))

geo_fisher_plot <- ggplot(PopNet_continent_fisher, aes(x=group, y=Continent)) +
  geom_point(aes(size = Representation, colour = Bonferroni.Adjusted), ) +
  geom_tile(fill=NA, color="grey50") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic(base_size = 14) +
  scale_color_gradient(low = "#023858", high = "#CCDEFC", name = paste0("BF-Adjusted", '\n', "p-value"), 
                       limits = c(0,0.05), breaks = c(0,0.025,0.05)) +
  scale_size_continuous(range = c(3,15), limits = c(1,9), breaks = c(1,5,9)) +
  xlab("PopNet Clade") +
  ylab("Continent")
geo_fisher_plot

#Let's expand this to AMR and mobile elements
#AMR
PopNet_AMR_fisher <- read.csv("graphia_fisher/AMR_group_enrichment.csv") %>%
  mutate(Representation = case_when(Bonferroni.Adjusted > 0.05 ~ NA,
                                    TRUE ~ Representation)) %>%
  mutate(Bonferroni.Adjusted = case_when(Bonferroni.Adjusted > 0.05 ~ NA,
                                         TRUE ~ Bonferroni.Adjusted))

AMR_fisher_plot <- ggplot(PopNet_AMR_fisher, aes(x=group, y=Antimicrobial)) +
  geom_point(aes(size = Representation, colour = Bonferroni.Adjusted)) +
  geom_tile(fill=NA, color="grey50") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic(base_size = 14) +
  scale_size_continuous(range = c(3,15), limits = c(1,9), breaks = c(1,5,9)) +
  scale_color_gradient(low = "#023858", high = "#CCDEFC", name = paste0("BF-Adjusted", '\n', "p-value"), 
                       limits = c(0,0.05), breaks = c(0,0.025,0.05)) +
  xlab("PopNet Clade") +
  ylab("Antimicrobial")

AMR_fisher_plot

ggsave("paper assembly/Figures/AMR_fisher_plot.pdf", test_fisher_plot, height = 4, width = 10,
       units = "in", dpi = 300)
ggsave("paper assembly/Figures/AMR_fisher_plot.png", test_fisher_plot, height = 4, width = 10,
       units = "in", dpi = 300)

#Mobile elements
PopNet_mobile_fisher <- read.csv("graphia_fisher/mobile_PopNet_fisher.csv") %>%
  mutate(Representation = case_when(Bonferroni.Adjusted > 0.05 ~ NA,
                                    TRUE ~ Representation)) %>%
  mutate(Bonferroni.Adjusted = case_when(Bonferroni.Adjusted > 0.05 ~ NA,
                                         TRUE ~ Bonferroni.Adjusted))

mobile_fisher_plot <- ggplot(PopNet_mobile_fisher, aes(x=group, y=Mobile.Element)) +
  geom_point(aes(size = Representation, colour = Bonferroni.Adjusted)) +
  geom_tile(fill=NA, color="grey50") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic(base_size = 14) +
  scale_color_gradient(low = "#023858", high = "#CCDEFC", name = paste0("BF-Adjusted", '\n', "p-value"), 
                       limits = c(0,0.05), breaks = c(0,0.025,0.05)) +
  scale_size_continuous(range = c(3,15), limits = c(1,9), breaks = c(1,5,9)) +
  xlab("PopNet Clade") +
  ylab("Mobile Element")

mobile_fisher_plot

ggsave("paper assembly/Figures/mobile_fisher_plot.pdf", mobile_fisher_plot, height = 8, width = 10,
       units = "in", dpi = 300)

ggsave("paper assembly/Figures/mobile_fisher_plot.png", mobile_fisher_plot, height = 4, width = 10,
       units = "in", dpi = 300)

#Combine fisher plots
graphia_fisher <- ggarrange(geo_fisher_plot, AMR_fisher_plot, mobile_fisher_plot,
          heights = c(6,4),
          labels = c("i)", "ii)", "iii)"),
          ncol = 2, nrow = 2, align = "v", common.legend = T, legend = "right")
ggsave("paper assembly/Figures/graphia_fisher.pdf", graphia_fisher, height = 6, width = 15,
       units = "in", dpi = 300)

#mosaic plots for association with Africa/Asia w/ clade B and plasmids
isolate_metadata <- read.csv("isolate_metadata.csv")

PopNet_colours <- c("#ffaa7f","#e3ff7f","#7fff8c","#7ffffc",
                    "#7f90ff","#df7fff","#7f2a00","#637f00","#007f0c",
                    "#007f7d", "#00107f", "#5f007f")

ggplot(data = isolate_metadata) +
  geom_mosaic(aes(x=product(group, Continent), fill = group), offset = 0.001) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = PopNet_colours) +
  facet_grid(~Continent)

#overlap of SNPs with penicillin plasmid and group B

blac_snps <- read.csv("gwas/results/sig_Blac_gwas_ann.csv")


isolate_metadata <- read.csv("isolate_metadata.csv")

PopNet_res <- isolate_metadata %>%
  select(Run, group, Penicillin) %>%
  mutate(Penicillin = case_when(Penicillin == "NOT_FOUND" ~ "SUSCEPTIBLE",
                                TRUE ~ Penicillin)) %>%
  arrange(group, factor(Penicillin, levels = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT")))
  
PopNet_resistance <- PopNet_res %>%
  select(!group) %>%
  pivot_longer(!Run, names_to = "x", values_to = "fill")
PopNet_resistance$type <- "Penicillin Resistance"

PopNet_ID <- PopNet_res %>%
  select(Run, group) %>%
  mutate(group = case_when(group == "A" ~ 0,
                           group == "B" ~ 1,
                           group == "C" ~ 2,
                           group == "D" ~ 3,
                           group == "E" ~ 4,
                           group == "F" ~ 5,
                           group == "G" ~ 6,
                           group == "H" ~ 7,
                           group == "I" ~ 8,
                           group == "J" ~ 9,
                           group == "K" ~ 10,
                           group == "L" ~ 11))
colnames(PopNet_ID)[2] <- "fill" 
PopNet_ID$x <- "PopNet"
PopNet_ID$type <- "PopNet group"

PopNet_mobile <- isolate_metadata %>%
  select(Run, Blac) %>%
  pivot_longer(!Run, names_to = "x", values_to = "fill") %>%
  distinct()
PopNet_mobile$type <- "Blac plasmid"


PopNet_colours_AMR <- c("lightblue", "royalblue","darkorchid4",
                        "#ffaa7f","#e3ff7f","#7fff8c","#7ffffc",
                        "#7f90ff","#df7fff","#7f2a00","#637f00","#007f0c",
                        "#007f7d", "#00107f", "#5f007f",
                        "forestgreen", "black")
names(PopNet_colours_AMR) <- unique(PopNet_AMR_df$fill)

PEN_Blac_plot <- ggplot(data = PopNet_AMR_df, aes(x=x, y= factor(Run, levels = unique(PopNet_res$Run)))) +
  geom_tile(aes(fill = fill)) + 
  facet_grid(cols = vars(type), scales = "free") +
  scale_fill_manual(name = "Resistance", breaks = c("SUSCEPTIBLE", "INTERMEDIATE", "RESISTANT"), values = PopNet_colours_AMR) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 10))

ggsave(PEN_Blac_plot, file = "paper assembly/Figures/PEN_Blac_plot.pdf", dpi = 300, height = 8, width = 6,
       units = "in")

sex_count <- isolate_metadata %>%
  count(Sex, Site)

# NGO/NGK mapping ---------------------------------------------------------
NGK_mapping <- read.delim("NGK_locus_mapping.txt", header = F) %>%
  select(Gene=V1, NGK_locus=V2, id=V3, eval=V11) %>%
  group_by(Gene) %>%
  filter(eval < 0.0001) %>%
  summarise(NGK_locus = toString(NGK_locus), .groups = 'drop') %>%
  mutate(NGK_locus = str_replace_all(NGK_locus, ",", ";"))

NGO_mapping <- read.delim("NGO_locus_mapping.txt", header = F) %>%
  select(Gene=V1, NGO_locus=V2, id=V3, eval=V11) %>%
  group_by(Gene) %>%
  filter(eval < 0.0001) %>%
  summarise(NGO_locus = toString(NGO_locus), .groups = 'drop') %>%
  mutate(NGO_locus = str_replace_all(NGO_locus, ",", ";"))

gene_meta <- read.csv("gene_metadata.csv") %>%
  left_join(NGK_mapping) %>%
  left_join(NGO_mapping)

write.csv(gene_meta, "gene_metadata.csv", quote = F, row.names = F)

# Accessory genes --------------------------------------------------------
#Get only genes/genomes from accessory gene clusters
acc_gene_cluster <- read.csv("panaroo_results/acc_gene_clusters.csv")
roary_gpa <- read.csv("panaroo_results/gene_presence_absence_roary.csv") %>%
  filter(Gene %in% acc_gene_cluster$Node.Name) %>%
  select_if(function(x) !(all(is.na(x)) | all(x=="")))

#Isolate clade D cluster (MCL cluster 5)
clusterD <- acc_gene_cluster %>%
  filter(MCL.Cluster == "Cluster 5")
Dgenes <- read.csv("panaroo_results/gene_presence_absence_roary.csv") %>%
  filter(Gene %in% clusterD$Node.Name) %>%
  select_if(function(x) !(all(is.na(x)) | all(x=="")))
Disolate_names <- colnames(Dgenes)[12:32]
Disolates <- read.csv("isolate_metadata.csv") %>%
  filter(Run %in% Disolate_names)

#read gff files
cladeD_contigs <- read.gff("panaroo_results/acc_gene_clusters/cluster5/gff/ERR026590.gff") %>%
  filter(type == "CDS") %>%
  separate_wider_delim(attributes, delim = ";", names = "ID", too_many = "drop") %>%
  mutate(ID = str_remove_all(ID, "ID=")) %>%
  filter(seqid == "ERR026590_25" | seqid == "ERR026590_31") %>%
  mutate(CladeD = case_when(ID %in% Dgenes$ERR026590 ~ "CladeD",
                            TRUE ~ NA))

write.csv(cladeD_contigs, "panaroo_results/acc_gene_clusters/cluster5/clade_D_gggenes.csv", quote = F, row.names = F)

gene_pa <- read.csv("panaroo_results/gene_presence_absence_roary.csv")
gene_meta <- read.csv("gene_metadata.csv") 

contig_genes <- gene_pa %>%
  select(Gene, Isolates, Annotation, ERR026590) %>%
  filter(ERR026590 %in% cladeD_contigs$ID) %>%
  left_join(gene_meta)

cladeD_gene_df <- read.csv("panaroo_results/acc_gene_clusters/cluster5/clade_D_gggenes.csv")

ggplot(cladeD_gene_df, aes(xmin = start, xmax = end, y = seqid, fill = CladeD, label = Gene,
                           forward = strand)) +
  geom_gene_arrow() +
  geom_gene_label(align = "left") +
  facet_wrap(~ seqid, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

#Isolate cluster 11 (Possible microvirus)
cluster11 <- acc_gene_cluster %>%
  filter(MCL.Cluster == "Cluster 11")
cluster11_genes <- read.csv("panaroo_results/gene_presence_absence_roary.csv") %>%
  filter(Gene %in% cluster11$Node.Name) %>%
  select_if(function(x) !(all(is.na(x)) | all(x=="")))
cluster11_names <- colnames(cluster11_genes)[11:17]
cluster11_isolates <- read.csv("isolate_metadata.csv") %>%
  filter(Run %in% cluster11_names)

isolate_metadata <- read.csv("isolate_metadata.csv")

#Foldseek results
foldseek <- read.delim("panaroo_results/acc_gene_clusters/foldseek_pred.tsv") %>%
  mutate(DB = case_when(str_detect(Output, "^AF-") ~ "AF",
                        TRUE ~ "PDB")) %>%
  group_by(Query, DB) %>%
  slice_max(order_by = TM_score, n = 3)
write.csv(foldseek, "panaroo_results/acc_gene_clusters/foldseek_filter.csv", quote = F, row.names = F)

gene_meta <- read.csv("gene_metadata.csv")
eggnog_COG <- read.delim("panaroo_results/eggnog_anno.tsv") %>%
  select(Gene = query, COG_category)

COG_plot_df <- gene_meta %>%
  left_join(eggnog_COG) %>%
  select(Gene, COG_category, Pangenome) %>%
  filter(Pangenome == "Core") %>%
  mutate(COG_category = case_when(COG_category == "-" ~ "S",
                                  is.na(COG_category) ~ "S",
                                  TRUE ~ COG_category)) %>%
  count(COG_category, Pangenome) %>%
  separate_longer_position(COG_category, width = 1) %>%
  filter(n > 3) %>%
  aggregate(. ~ COG_category + Pangenome, sum) %>%
  mutate(Pangenome = factor(Pangenome, levels = c("Core", "Accessory"))) %>%
  arrange(-n) %>%
  mutate(COG_category = factor(COG_category, levels = unique(COG_plot_df$COG_category))) %>%
  mutate(Percent = round(n/sum(n)*100, digits = 2))
  
ggplot(data = COG_plot_df, aes(x = COG_category, y = n, fill = COG_category, label = Percent)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("COG category") +
  ylab("# of genes") +
  guides(fill = "none") +
  geom_text(aes(label=Percent), vjust=-0.3, size=3.5)
