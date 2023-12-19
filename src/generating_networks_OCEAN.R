################################################################################
## Description: generating organism-specific PKNs for OCEAN
## Author: Diego Mañanes
## Date: 23/11/30
################################################################################

## set of dependencies
suppressMessages(library("biomaRt"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))
suppressMessages(library("metaboliteIDmapping"))
suppressMessages(library("R.matlab"))
suppressMessages(library("dplyr"))

if (!requireNamespace("here")) {
  install.packages("here")
}

projectPath <- here::here()
## source
source(file.path(projectPath, "src/final_functions_PKN_OCEAN.R"))

## KEGG compounds: not neeeded
# kegg.all.compunds <- readRDS(
#   file.path(projectPath, "data", "231020_KEGG_All_Compunds", "ALL_COMPUNDS_KEGG.rds")
# )
# KEGG.compounds <- unlist(kegg.all.compunds, recursive = FALSE)
# names(KEGG.compounds) <- sapply(KEGG.compounds, \(x) x$ENTRY)

## Homo sapiens
reactions.map.hs <- read.delim(file.path(projectPath, "data/Human-PKN/reactions.tsv")) 
metabolites.map.hs <- read.delim(file.path(projectPath, "data/Human-PKN/metabolites.tsv"))

OCEAN.pkn.hs <- create_PKN_OCEAN( 
  GSMM.matlab.path = file.path(projectPath, "data/Human-PKN/Human-GEM.mat"),
  GSMM.reactions.map = reactions.map.hs,
  GSMM.metabolites.map = metabolites.map.hs,
  KEGG.compounds = NULL,
  translate.genes = TRUE,
  organism = 9606
)
saveRDS(OCEAN.pkn.hs, file.path(projectPath, "output/OCEAN.PKN.hs.9606.rds"))

## Mus musculus
reactions.map.mm <- read.delim(file.path(projectPath, "data/Mouse-PKN/reactions.tsv")) 
metabolites.map.mm <- read.delim(file.path(projectPath, "data/Mouse-PKN/metabolites.tsv"))
OCEAN.pkn.mm <- create_PKN_OCEAN( 
  GSMM.matlab.path = file.path(projectPath, "data/Mouse-PKN/Mouse-GEM.mat"),
  GSMM.reactions.map = reactions.map.mm,
  GSMM.metabolites.map = metabolites.map.mm,
  KEGG.compounds = NULL,
  translate.genes = FALSE,
  organism = 10090
)
saveRDS(OCEAN.pkn.mm, file.path(projectPath, "output/OCEAN.PKN.mm.10090.rds"))

## Rat novergicus
reactions.map.rn <- read.delim(file.path(projectPath, "data/Rat-PKN/reactions.tsv")) 
metabolites.map.rn <- read.delim(file.path(projectPath, "data/Rat-PKN/metabolites.tsv"))
OCEAN.pkn.rn <- create_PKN_OCEAN( 
  GSMM.matlab.path = file.path(projectPath, "data/Rat-PKN/Rat-GEM.mat"),
  GSMM.reactions.map = reactions.map.rn,
  GSMM.metabolites.map = metabolites.map.rn,
  KEGG.compounds = NULL,
  translate.genes = FALSE,
  organism = 10116
)
saveRDS(OCEAN.pkn.rn, file.path(projectPath, "output/OCEAN.PKN.rn.10116.rds"))
