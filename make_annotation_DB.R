#!/usr/bin/Rscript
################################################################################
# Matthew Wersebe
# University of Oklahoma
# April 09, 2023 
# Tidy_pop_genomics
################################################################################
## Purpose: Build an SQLite DB for use in Tidy Population Genomics
## Usage: make_genome_DB.R < options >
## 
## 
##
################################################################################
suppressMessages(library(optparse))

# Create Options List
option_list = list(
  make_option(c("-f", "--fasta"), type = "character", default=NA, action = "store",
              help="Protein Sequences in fasta format"),
  make_option(c("-o", "--output"), type = "character", default=NA, action = "store",
              help="Output SQLite Formatted DB"),
  make_option(c("-g", "--gff"), type = "character", default=NA, action = "store",
              help="Gff3 formatted genome annotation"),
  make_option(c("-v", "--vcf"), type = "character", default=NA, action = "store",
              help="vcf file with SNPs"),
  make_option(c("-m", "--mappings"), type = "character", default=8, action = "store",
              help="Pathern DB generic mappings file for Protein Sequences"),
  make_option(c("-e", "--vep"), type = "character", default=8, action = "store",
              help="SNP effect annotation created with VEP"))

opts = parse_args(OptionParser(option_list=option_list))

################################################################################
# Check for requires inputs:

if(is.na(opts$fasta)){
  stop("Protein Fasta is Required")
}else if(is.na(opts$gff)){
  stop("GFF3 File is Required")
}else if(is.na(opts$vcf)){
  stop("VCF file is Required")
}else if(is.na(opts$mappings)){
  stop("Panther DB mappings file is required")
}else if(is.na(opts$vep)){
  stop("VEP SNP effect annotation is required")
}else if(is.na(opts$output)){
  stop("Give Output DB a name")
}else{
  print(paste0("Creating a DB named ", opts$output))
}
################################################################################
# Create a DB to write to:
suppressMessages(library(DBI))
suppressMessages(library(RSQLite))

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = opts$output)
#con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "SC_RAD_BayScan.sqlite")
################################################################################
# Read in the Fasta file and add it as a table:
suppressMessages(library(tidyverse))
suppressMessages(library(tidysq))

print(paste0("Working on: ", opts$fasta))

prots <- tidysq::read_fasta(opts$fasta)
prots <- prots %>% tidyr::separate(name, c("XP_name", "Full_name"), " ", extra = "merge")

DBI::dbWriteTable(con, "Proteins", prots)

rm(prots)
################################################################################
# Read in GFF3 file for parsing:
print(paste0("Working on: ", opts$gff))

suppressMessages(library(annotatr))

## Modified Function from Annotatr
extract_attr <- function(x, key)  {
  has_key <- regexpr(key, x)
  value <- gsub(sprintf('.*%s=([^;]+).*', key), '\\1', x)
  value[has_key == -1] <- NA
  out <- strsplit(value, ',')
  multi_value <- any(sapply(out, length) > 1)
  if (!multi_value) return(unlist(out))
  out
}
## Modified Function from Annotatr 
tidy_gff <- function(x, chroms=NULL,
                     features=c('CDS', 'lnc_RNA')) {
  out <- dplyr::select(dplyr::mutate(dplyr::filter(x, feature %in% features), 
                                     name=extract_attr(attr,'Name'),
                                     parent=extract_attr(attr,'Parent'),
                                     gene=extract_attr(attr,'gene')),
                       chrom, start, end, feature, name, parent, gene)
  if (!is.null(chroms)) {
    out[out$chrom %in% chroms, ]
  }
  out
}


gff <- annotatr::read_gff(opts$gff)

gff <- gff %>% tidy_gff() %>% mutate(ID = stringr::str_remove(gene, "LOC")) %>%
  tidyr::separate_wider_delim(parent, names = c("temp", "parent"), delim = "-") %>%
  select(-c("temp", "gene"))

DBI::dbWriteTable(con, "GFF", gff)

rm(gff)
################################################################################
# Read in VEP output:

print(paste0("Working on: ", opts$vep))

VEP <- readr::read_tsv(opts$vep, comment = "##", col_names = T)

VEP <- VEP %>% dplyr::mutate(Impact = extract_attr(Extra, 'IMPACT')) %>% 
  dplyr::select(1:13, Impact) %>%
  tidyr::separate(Location, c("CHROM", "POS"), sep = ":", remove = F) %>%
  mutate(POS = as.numeric(POS))

DBI::dbWriteTable(con, "VEP", VEP)

rm(VEP)
################################################################################
# Read in VCF 
suppressMessages(library(radiator))

print(paste0("Working on: ", opts$vcf))

VCF <- radiator::tidy_vcf(opts$vcf)

VCF <- VCF %>% mutate(CHROM = case_when(
  CHROM == 1 ~ "CHR01",
  CHROM == 2 ~ "CHR02",
  CHROM ==3 ~ "CHR03",
  CHROM == 4 ~ "CHR04",
  CHROM == 5 ~ "CHR05",
  CHROM == 6 ~ "CHR06",
  CHROM == 7 ~ "CHR07",
  CHROM == 8 ~ "CHR08",
  CHROM == 9 ~ "CHR09",
  CHROM == 10 ~ "CHR10",
  CHROM == 11 ~ "CHR11",
  CHROM == 12 ~ "CHR12"))

DBI::dbWriteTable(con, "VCF", VCF)

rm(VCF)
################################################################################
# Read in Mappings File:

print(paste0("Working on: ", opts$mappings))

Mappings <- readr::read_tsv(opts$mappings, col_names = FALSE)

#Mappings <- readr::read_tsv("SC_F0-13Bv2_protein.pantherDB.tsv", col_names = FALSE)

names(Mappings) <- c("Gene", "PantherFam", "LongName", "Evalue", "Bitscore", "MapRange")

DBI::dbWriteTable(con, "Mappings", Mappings)
