################################################################################
## Description: generating organism-specific PKNs for COSMOS
## Author: Diego Mañanes
## Date: 23/11/16
################################################################################

## set of dependencies
suppressMessages(library("biomaRt"))
suppressMessages(library("OmnipathR"))
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
source(file.path(projectPath, "src/final_functions_PKN_COSMOS.R"))

## Homo sapiens
reactions.map.hs <- read.delim(
  file.path(projectPath, "data/Human-PKN/reactions.tsv")
) 
metabolites.map.hs <- read.delim(
  file.path(projectPath, "data/Human-PKN/metabolites.tsv")
)
cosmos.pkn.hs <- create_PKN_COSMOS( 
  organism = 9606,
  GSMM.matlab.path = file.path(projectPath, "data/Human-PKN/Human-GEM.mat"),
  GSMM.reactions.map = reactions.map.hs,
  GSMM.metabolites.map = metabolites.map.hs,
  translate.genes = TRUE,
  stitch.actions.file = file.path(
    projectPath, "data/Human-PKN/9606.actions.v5.0.tsv"
  ), 
  stitch.links.file = file.path(
    projectPath, "data/Human-PKN/9606.protein_chemical.links.detailed.v5.0.tsv"
  ), 
  clear_omnipath_cache = T
)
saveRDS(cosmos.pkn.hs, file.path(projectPath, "output/COSMOS.PKN.hs.9606.rds"))

## Mus musculus
reactions.map.mm <- read.delim(
  file.path(projectPath, "data/Mouse-PKN/reactions.tsv")
) 
metabolites.map.mm <- read.delim(
  file.path(projectPath, "data/Mouse-PKN/metabolites.tsv")
)
cosmos.pkn.mm <- create_PKN_COSMOS( 
  organism = 10090,
  GSMM.matlab.path = file.path(projectPath, "data/Mouse-PKN/Mouse-GEM.mat"),
  GSMM.reactions.map = reactions.map.mm,
  GSMM.metabolites.map = metabolites.map.mm,
  stitch.actions.file = file.path(
    projectPath, "data/Mouse-PKN/10090.actions.v5.0.tsv"
  ), 
  stitch.links.file = file.path(
    projectPath, "data/Mouse-PKN/10090.protein_chemical.links.detailed.v5.0.tsv"
  )
)
saveRDS(cosmos.pkn.mm, file.path(projectPath, "output/COSMOS.PKN.mm.10090.rds"))

## Rat novergicus
reactions.map.rn <- read.delim(
  file.path(projectPath, "data/Rat-PKN/reactions.tsv")
) 
metabolites.map.rn <- read.delim(
  file.path(projectPath, "data/Rat-PKN/metabolites.tsv")
)
cosmos.pkn.rn <- create_PKN_COSMOS( 
  organism = 10116,
  GSMM.matlab.path = file.path(projectPath, "data/Rat-PKN/Rat-GEM.mat"),
  GSMM.reactions.map = reactions.map.rn,
  GSMM.metabolites.map = metabolites.map.rn,
  stitch.actions.file = file.path(
    projectPath, "data/Rat-PKN/10116.actions.v5.0.tsv"
  ), 
  stitch.links.file = file.path(
    projectPath, "data/Rat-PKN/10116.protein_chemical.links.detailed.v5.0.tsv"
  )
)
saveRDS(cosmos.pkn.rn, file.path(projectPath, "output/COSMOS.PKN.rn.10116.rds"))
