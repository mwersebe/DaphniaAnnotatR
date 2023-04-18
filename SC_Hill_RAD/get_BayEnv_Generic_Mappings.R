################################################################################
# Matthew Wersebe
# University of Oklahoma
# April 09, 2023 
# GO Enrichment for BayeEnv: SC Hill RAD
################################################################################
suppressMessages(library(DBI))
suppressMessages(library(RSQLite))

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "SC_RAD_BayScan.sqlite")

library(tidyverse)

outliers <- readr::read_csv("/home/giovanni/SC_Hill_RAD/GOEnrichment/BayEnv/SoutHCenter/SC_BayEnv_outliers.csv")


BayEnv_OrgC <- outliers %>% filter(Factor == "OrgC") %>%
  left_join(., as_tibble(con %>% tbl("VEP")), by = c("Chrom" = "CHROM", "Position" = "POS")) %>%
  select(Chrom, Position, Gene, Feature) %>% distinct() %>% 
  left_join(., as_tibble(con %>% tbl("GFF")), by = c("Gene" = "ID", "Feature" = "parent")) %>%
  select(Chrom, Position, Gene, Feature, name) %>% distinct() %>% na.exclude() %>%
  left_join(., as_tibble(con %>% tbl("Mappings")), by = c("name" = "Gene")) %>% na.exclude() %>%
  group_by(Chrom, Position, Gene) %>% slice(1) %>% ungroup() %>% select(5:10)


BayEnv_OrgC %>% write_tsv(., "SC_outliers_panther.tsv", col_names = F)


#con %>% tbl("Proteins") %>% filter(XP_name == "XP_046631016.1")



con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "Hill_RAD_BayScan.sqlite")


outliers <- readr::read_csv("/home/giovanni/SC_Hill_RAD/GOEnrichment/BayEnv/Hill/Hill_outliers_BayEnv.csv")


outliers <- outliers %>% mutate(Chrom = case_when(
  Chrom == 1 ~ "CHR01",
  Chrom == 2 ~ "CHR02",
  Chrom == 3 ~ "CHR03",
  Chrom == 4 ~ "CHR04",
  Chrom == 5 ~ "CHR05",
  Chrom == 6 ~ "CHR06",
  Chrom == 7 ~ "CHR07",
  Chrom == 8 ~ "CHR08",
  Chrom == 9 ~ "CHR09",
  Chrom == 10 ~ "CHR10",
  Chrom == 11 ~ "CHR11",
  Chrom == 12 ~ "CHR12"
))

#VEP <- as_tibble(con %>% tbl("VEP")) %>% mutate(POS = as.numeric(POS))
#DBI::dbWriteTable(con, "VEP", VEP, overwrite = T)

outliers %>% filter(Factor == "OrgC") %>%
  left_join(., as_tibble(con %>% tbl("VEP")), by = c("Chrom" = "CHROM", "Position" = "POS")) %>%
  select(Chrom, Position, Gene, Feature) %>% distinct() %>% 
  left_join(., as_tibble(con %>% tbl("GFF")), by = c("Gene" = "ID", "Feature" = "parent")) %>%
  select(Chrom, Position, Gene, Feature, name) %>% distinct() %>% na.exclude() %>% 
  left_join(., as_tibble(con %>% tbl("Mappings")), by = c("name" = "Gene")) %>% na.exclude() %>%
  group_by(Chrom, Position, LongName) %>% slice(1) %>% ungroup() %>% select(5:10) %>% 
  group_by(name, LongName) %>% slice(1) %>%
  write_tsv("Hill_OrgC_outliers_panther.tsv", col_names = F)
  



outliers %>% filter(Factor == "Age") %>%
  left_join(., as_tibble(con %>% tbl("VEP")), by = c("Chrom" = "CHROM", "Position" = "POS")) %>%
  select(Chrom, Position, Gene, Feature) %>% distinct() %>% 
  left_join(., as_tibble(con %>% tbl("GFF")), by = c("Gene" = "ID", "Feature" = "parent")) %>%
  select(Chrom, Position, Gene, Feature, name) %>% distinct() %>% na.exclude() %>% 
  left_join(., as_tibble(con %>% tbl("Mappings")), by = c("name" = "Gene")) %>% na.exclude() %>%
  group_by(Chrom, Position, LongName) %>% slice(1) %>% ungroup() %>% select(5:10) %>%
  group_by(name) %>% slice(1) %>%
  write_tsv("Hill_Age_outliers_panther.tsv", col_names = F)



outliers %>% filter(Factor == "P") %>%
  left_join(., as_tibble(con %>% tbl("VEP")), by = c("Chrom" = "CHROM", "Position" = "POS")) %>%
  select(Chrom, Position, Gene, Feature) %>% distinct() %>% 
  left_join(., as_tibble(con %>% tbl("GFF")), by = c("Gene" = "ID", "Feature" = "parent")) %>%
  select(Chrom, Position, Gene, Feature, name) %>% distinct() %>% na.exclude() %>% 
  left_join(., as_tibble(con %>% tbl("Mappings")), by = c("name" = "Gene")) %>% na.exclude() %>%
  group_by(Chrom, Position, LongName) %>% slice(1) %>% ungroup() %>% select(5:10) %>%
  group_by(name) %>% slice(1) %>%
  write_tsv("Hill_P_outliers_panther.tsv", col_names = F)


outliers %>% filter(Factor == "CaCO3") %>%
  left_join(., as_tibble(con %>% tbl("VEP")), by = c("Chrom" = "CHROM", "Position" = "POS")) %>%
  select(Chrom, Position, Gene, Feature) %>% distinct() %>% 
  left_join(., as_tibble(con %>% tbl("GFF")), by = c("Gene" = "ID", "Feature" = "parent")) %>%
  select(Chrom, Position, Gene, Feature, name) %>% distinct() %>% na.exclude() %>% 
  left_join(., as_tibble(con %>% tbl("Mappings")), by = c("name" = "Gene")) %>% na.exclude() %>%
  group_by(Chrom, Position, Gene) %>% slice(1) %>%
  write_tsv("Hill_CaCO3_outliers_panther.tsv", col_names = F)


