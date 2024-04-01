################################################################################
## Description: set of functions to generate COSMOS PKN for different organisms
## OmnipathR, GSMM and STITCH resources. 
## Author: Diego Mañanes
## Date: 24/03/28
################################################################################

## set of dependencies
# suppressMessages(library("biomaRt"))
# suppressMessages(library("OmnipathR"))
# suppressMessages(library("readr"))
# suppressMessages(library("stringr"))
# suppressMessages(library("metaboliteIDmapping"))
# suppressMessages(library("R.matlab"))
# suppressMessages(library("dplyr"))

## exported
create_PKN_COSMOS <- function(
  organism,
  GSMM.matlab.path,
  GSMM.reactions.map,
  GSMM.metabolites.map,
  stitch.actions.file, 
  stitch.links.file, 
  translate.genes = FALSE,
  biomart.use.omnipath = TRUE,
  GSMM.reactions.map.col = "rxns", 
  GSMM.metabolites.map.col = "mets",
  clear_omnipath_cache = FALSE,
  GSMM.list.params = list(
    stoich.name = "S",
    reaction.name = "grRules",
    lb.name = "lb",
    ub.name = "ub",
    rev.name = "rev",
    reaction.ID.name = "rxns",
    metabolites.ID.name = "mets",
    metabolites.names.name = "metNames",
    metabolites.fomulas.name = "metFormulas",
    metabolites.inchi.name = "inchis"
  ),
  GSMM.degree.mets.threshold = 400,
  stitch.threshold = 700,
  verbose = TRUE
) {
  ## check organisms
  dataset.biomart <- switch(
    as.character(organism), 
    "9606" = "hsapiens_gene_ensembl",
    "10090" = "mmusculus_gene_ensembl",
    "10116" = "rnorvegicus_gene_ensembl",
    "7955" = "drerio_gene_ensembl",
    "7227" = "dmelanogaster_gene_ensembl",
    "6239" = "celegants_gene_ensembl"
  )
  if (is.null(dataset.biomart)) 
    stop(
      "Chosen organism is not recognizable Available options are: ", 
      paste(c(9606, 10090, 10116, 7955, 7227, 6239), collapse = ", ")
    )
  ## checking if files exist
  if (!file.exists(GSMM.matlab.path)) {
    stop("GSMM.matlab.path file does not exist")
  } else if (!file.exists(stitch.actions.file)) {
    stop("stitch.actions.file file does not exist")
  } else if (!file.exists(stitch.links.file)) {
    stop("stitch.links.file file does not exist")
  }
  if (biomart.use.omnipath == TRUE) {
    if (verbose) message(">>> Using the OmnipathR to retrieve biomart information\n")
    ## get info from BiomartR using OmnipathR
    mapping.biomart <- OmnipathR::biomart_query(
      attrs = c(
        "ensembl_peptide_id",'ensembl_gene_id', 'external_gene_name'
      ),
      dataset = dataset.biomart
    ) %>% as.data.frame()
  } else {
    if (verbose) message("\n>>> Using the BiomaRt R package when needed")
    
    mapping.biomart <- NULL
  }
  ## Omnipath data
  if (verbose) message("\n>>> Getting Omnipath PKN...\n")
  omnipath.PKN <- .retrievingOmnipath(
    organism, clear_omnipath_cache = clear_omnipath_cache
  )
  ## Getting GSSM PKN
  if (verbose) message("\n>>> Getting GSMM PKN...\n")
  gsmm.PKN.list <- .create_GSMM_basal_PKN(
    matlab.path = GSMM.matlab.path,
    reactions.map = GSMM.reactions.map, 
    reactions.map.col = GSMM.reactions.map.col, 
    metabolites.map = GSMM.metabolites.map,
    metabolites.map.col = GSMM.metabolites.map.col,
    list.params.GSMM = GSMM.list.params,
    degree.mets.cutoff = GSMM.degree.mets.threshold,
    verbose = verbose
  )
  if (verbose) message("\n>>> Formatting GSMM PKN for COSMOS...\n")
  gsmm.PKN.list <- .mets_to_HMDB(gsmm.PKN.list)
  if (translate.genes){
    gsmm.PKN.list <- .genes_to_symbol(
      gsmm.PKN.list, organism = organism, mapping.biomart = mapping.biomart
    )
  }
  gsmm.PKN.list <- .format_GSMM_COSMOS(gsmm.PKN.list, verbose = TRUE)
  
  if (verbose) message("\n>>> Getting STITCH PKN...\n")
  stitch.PKN <- .formatSTITCH(
    actions.file = stitch.actions.file, 
    links.file = stitch.links.file, 
    organism = organism,
    omnipath.PKN = omnipath.PKN,
    mapping.biomart = mapping.biomart,
    threshold = stitch.threshold,
    verbose = verbose
  )
  output.final <- .mixing_resources(
    GSMM.network = gsmm.PKN.list[[1]], 
    omnipath.PKN = omnipath.PKN, 
    stitch.PKN = stitch.PKN
  )
  
  if (verbose) message("\nDONE")
  
  return(
    list(
      COSMOS.PKN = output.final, 
      GSSM.mets.map = gsmm.PKN.list[[2]], 
      GSSM.reac.to.gene = gsmm.PKN.list[[3]], 
      reac.map = gsmm.PKN.list[[4]]
    )
  )
}


## retrieving PKN info from omnipathR
.retrievingOmnipath <- function(
    organism = 9606,
    clear_omnipath_cache = F
) {
  if(clear_omnipath_cache)
  {
    OmnipathR::omnipath_cache_remove()
  }
  full_pkn_mm <- as.data.frame(import_omnipath_interactions(organism = organism))
  full_pkn_mm <- full_pkn_mm[!is.na(full_pkn_mm$references),]
  clean_PKN_mm <- full_pkn_mm[
    full_pkn_mm$consensus_stimulation == 1 | 
      full_pkn_mm$consensus_inhibition == 1,
  ]
  clean_PKN_mm$sign <- clean_PKN_mm$consensus_stimulation - 
    clean_PKN_mm$consensus_inhibition
  clean_PKN_mm <- clean_PKN_mm[,c(3 ,4, 16)]
  clean_PKN_supp_mm <- clean_PKN_mm[clean_PKN_mm$sign == 0,]
  clean_PKN_supp_mm$sign <- -1
  clean_PKN_mm[clean_PKN_mm$sign == 0, "sign"] <- 1
  clean_PKN_mm <- as.data.frame(rbind(clean_PKN_mm, clean_PKN_supp_mm))
  names(clean_PKN_mm) <- c("source", "target", "sign")
  
  return(clean_PKN_mm)
}

## helper function
.vecInfoMetab <- function(mat.object, attribs.mat, name) {
  unlist(
    sapply(
      X = mat.object[[which(attribs.mat == name)]],
      FUN = \(elem) {
        if (length(unlist(elem) != 0)) {
          return(elem)
        } else {
          return(NA)
        }
      }
    )
  )
}

## this will be exported, as it can be used in OCEAN as well
.create_GSMM_basal_PKN <- function(
  matlab.path,
  reactions.map,
  metabolites.map,
  reactions.map.col = "rxns",
  metabolites.map.col = "mets",
  list.params.GSMM = list(
    stoich.name = "S",
    reaction.name = "grRules",
    lb.name = "lb",
    ub.name = "ub",
    rev.name = "rev",
    reaction.ID.name = "rxns",
    metabolites.ID.name = "mets",
    metabolites.names.name = "metNames",
    metabolites.fomulas.name = "metFormulas",
    metabolites.inchi.name = "inchis" 
  ),
  degree.mets.cutoff = 400,
  verbose = TRUE
) {
  ## check parameters
  if (!file.exists(matlab.path)) {
    stop("matlab.path file does not exist")
  } else if (!reactions.map.col %in% colnames(reactions.map)) {
    stop("reactions.map.col cannot be found in reactions.map data.frame")
  } else if (!metabolites.map.col %in% colnames(metabolites.map)) {
    stop("metabolites.map.col cannot be found in metabolites.map data.frame")
  } else if (degree.mets.cutoff < 1) {
    stop("degree.mets.cutoff cannot be less than 1")
  }
  
  if (verbose) message("\t>>> Reading matlab file")
  
  mat.object <- readMat(matlab.path)[[1]]
  attribs.mat <- rownames(mat.object)
  ## check elements are in the object
  invisible(
    sapply(
      names(list.params.GSMM), \(idx) {
        if (!list.params.GSMM[[idx]] %in% attribs.mat) {
          stop(
            paste0(
              x, "element in list.params.GSMM (", 
              list.params.GSMM[[idx]], ") is not in matlab object"
            )
          )
        }
      }
    )
  )
  ## obtaining data
  s.matrix <- mat.object[[which(attribs.mat == list.params.GSMM$stoich.name)]]
  reaction.list <- mat.object[[which(attribs.mat == list.params.GSMM$reaction.name)]]
  ##############################################################################
  ## reactions
  # direction reactions
  lbs <- as.data.frame(
    cbind(
      mat.object[[which(attribs.mat == list.params.GSMM$lb.name)]],
      mat.object[[which(attribs.mat == list.params.GSMM$ub.name)]],
      mat.object[[which(attribs.mat == list.params.GSMM$rev.name)]]
    )
  )
  ## this could be done with mutate
  lbs$direction <- ifelse(
    (mat.object[[which(attribs.mat == list.params.GSMM$ub.name)]] + 
       mat.object[[which(attribs.mat == list.params.GSMM$lb.name)]]) >= 0,
    "forward", "backward"
  )
  reversible <- ifelse(
    mat.object[[which(attribs.mat == list.params.GSMM$rev.name)]] == 1, 
    TRUE, FALSE
  ) %>% as.logical()
  reaction.ids <- unlist(
    mat.object[[which(attribs.mat == list.params.GSMM$reaction.ID.name)]]
  )
  ## reaction to genes df
  reaction.to.genes.df <- lapply(
    seq_along(reaction.list),
    \(idx) {
      genes.reac <- unlist(reaction.list[[idx]], recursive = FALSE)
      if (length(genes.reac) != 0) {
        genes <- unique(
          gsub(
            " and ", "_",
            ## TODO: ask about these two gsub calls (_AT[0-9] and [()]): what 
            # are they intendeed to?? I guess it only works for SYMBOL genes
            gsub(
              "[()]", "", 
              gsub("_AT[0-9]+","", strsplit(genes.reac, split = " or ")[[1]])
            )
          )
        )
        return(
          data.frame(
            Gene = genes, Reaction = rep(idx, length(genes)), 
            Reaction.ID = rep(reaction.ids[idx], length(genes))
          )
        )
      } else {
        return(
          data.frame(
            Gene = idx, Reaction = idx, Reaction.ID = reaction.ids[idx]
          )
        )
      }
    }
  ) %>% do.call(rbind, .)
  orphan.reacts <- grepl(pattern = "^\\d+$", reaction.to.genes.df$Gene)
  reaction.to.genes.df[orphan.reacts, "Reaction.ID"] <- paste0(
    "orphanReac.", reaction.to.genes.df[orphan.reacts, "Reaction.ID"]
  )
  reaction.to.genes.df[orphan.reacts, "Gene"] <- paste0(
    reaction.to.genes.df[orphan.reacts, "Gene"], ".",
    reaction.to.genes.df[orphan.reacts, "Reaction.ID"]
  )
  reaction.to.genes.df <- unique(reaction.to.genes.df)
  
  ##############################################################################
  ## metabolites
  metabolites.IDs <- unlist(
    mat.object[[which(attribs.mat == list.params.GSMM$metabolites.ID.name)]]
  )
  metabolites.names <- .vecInfoMetab(
    mat.object = mat.object, attribs.mat = attribs.mat, 
    name = list.params.GSMM$metabolites.names.name
  )
  ## check if IDs are the same and show number of lost metabolites
  # metabolites.map[[metabolites.map.col]]
  inter.metab <- intersect(
    metabolites.map[[metabolites.map.col]], metabolites.IDs
  )
  rownames(metabolites.map) <- metabolites.map[[metabolites.map.col]]
  metabolites.map <- metabolites.map[metabolites.IDs, ]
  ## adding additional information 
  metabolites.formulas <- .vecInfoMetab(
    mat.object = mat.object, attribs.mat = attribs.mat, 
    name = list.params.GSMM$metabolites.fomulas.name
  )
  metabolites.inchi <- .vecInfoMetab(
    mat.object = mat.object, attribs.mat = attribs.mat, 
    name = list.params.GSMM$metabolites.inchi.name
  )
  metabolites.map <- cbind(
    metabolites.map,
    Metabolite.Name = metabolites.IDs,
    Metabolite.Formula = metabolites.formulas,
    Metabolite.Inchi = metabolites.inchi
  )
  metabolites.map[metabolites.map == ""] <- NA
  
  ##############################################################################
  ## SIF file: PKN
  if (verbose) message("\n>>> Generating PKN")
  
  reaction.to.genes.df.reac <- reaction.to.genes.df
  reaction.list <- list()
  for (reac.idx in seq(ncol(s.matrix))) {
    reaction <- s.matrix[, reac.idx]
    #modify gene name so reactions that are catalised by same enzyme stay separated
    reaction.to.genes.df.reac[reaction.to.genes.df$Reaction == reac.idx, 1] <- paste(
      paste0("Gene", reac.idx), 
      reaction.to.genes.df[reaction.to.genes.df$Reaction == reac.idx, 1], 
      sep = "__"
    )
    # get the enzymes associated with reaction
    genes <- reaction.to.genes.df.reac[reaction.to.genes.df.reac$Reaction == reac.idx, 1]
    if (as.character(lbs[reac.idx, 4]) == "forward") {
      reactants <- metabolites.IDs[reaction == -1]
      products <- metabolites.IDs[reaction == 1]
    } else {
      reactants <- metabolites.IDs[reaction == 1]
      products <- metabolites.IDs[reaction == -1]
    }
    reactants <- paste0("Metab__", reactants)
    products <- paste0("Metab__", products)
    number_of_interations <- length(reactants) + length(products)
    # now for each enzyme, we create a two column dataframe recapitulating the 
    # interactions between the metabolites and this enzyme
    reaction.df <- lapply(
      X = as.list(genes), 
      FUN = \(gene) {
        gene.df <- data.frame(
          # reactants followed by the enzyme (the enzyme is repeated as many time as they are products)
          source = c(reactants, rep(gene, number_of_interations - length(reactants))), 
          # enzyme(repeated as many time as they are reactants) followed by products  
          target = c(rep(gene, number_of_interations - length(products)), products) 
        )
        if (reversible[reac.idx]) {
          gene.df.reverse <- data.frame(
            source = c(
              rep(
                paste(gene, "_reverse", sep = ""), 
                number_of_interations - length(products)
              ), 
              products
            ),
            target = c(
              reactants,
              rep(
                paste(gene, "_reverse", sep = ""), 
                number_of_interations - length(reactants)
              )
            )
          )
          gene.df <- rbind(gene.df, gene.df.reverse)
        }
        return(gene.df)
      }
    ) %>% do.call(rbind, .)
    reaction.list[[reac.idx]] <- reaction.df
  }
  reaction.df.all <- do.call(rbind, reaction.list)
  ## removing those reactions with no metab <--> gene
  reaction.df.all <- reaction.df.all[reaction.df.all$source != "Metab__" & 
                                         reaction.df.all$target != "Metab__",]
  ## only complete cases
  reaction.df.all <- reaction.df.all[complete.cases(reaction.df.all),]
  
  ##############################################################################
  ## removing cofactors (metabolites with a high degree)
  metabs.degree <- sort(
    table(
      grep(
        "^Metab__", c(reaction.df.all$source, reaction.df.all$target), 
        value = TRUE
      )
    ), 
    decreasing = TRUE
  )
  if (verbose) 
    message(
      "\t>>> Number of metabolites removed after degree >", 
      degree.mets.cutoff,  ": ", sum(metabs.degree >= degree.mets.cutoff)
    )
  metabs.degree.f <- metabs.degree[metabs.degree < degree.mets.cutoff]
  reactions.df.no.cofac <- reaction.df.all[
    reaction.df.all$source %in% names(metabs.degree.f) | 
      reaction.df.all$target %in% names(metabs.degree.f),
  ]
  mets <- grep(
    pattern = "Metab__",
    x = unique(c(reactions.df.no.cofac[[1]], reactions.df.no.cofac[[1]])),
    value = TRUE
  ) %>% gsub("Metab__", "", .)
  metabolites.map <- metabolites.map[mets, ]
  if (verbose) {
    message(
      "\t>>> Final number of connections: ", nrow(reactions.df.no.cofac)
    )
    message("\n\tDONE")
  }
  
  return(
    list(
      GSMM.PKN = reactions.df.no.cofac, 
      mets.map = metabolites.map, 
      reac.to.gene = reaction.to.genes.df.reac, 
      reac.map = reactions.map
    )
  )
}

## only works for those mapping dfs from the ddbb
.mets_to_HMDB <- function(
    list.network  
) {
  metab.map <- list.network[[2]]
  list.network[[1]] <- list.network[[1]] %>% mutate(
    source = case_when(
      grepl("Metab__", source) ~ case_when(
        !is.na(
          metab.map[gsub("Metab__", replacement = "", x = source), "metHMDBID"]
        ) ~ paste0(
          "Metab__", 
          metab.map[gsub("Metab__", replacement = "", x = source), "metHMDBID"],
          "_", str_sub(source, start = nchar(source))
        ),
        !is.na(
          metab.map[gsub("Metab__", replacement = "", x = source), "metKEGGID"]
        ) ~ paste0(
          "Metab__", 
          metab.map[gsub("Metab__", replacement = "", x = source), "metKEGGID"],
          "_", str_sub(source, start = nchar(source))
        ), TRUE ~ source
      ), TRUE ~ source 
    ),
    target = case_when(
      grepl("Metab__", target) ~ case_when(
        !is.na(
          metab.map[gsub("Metab__", replacement = "", x = target), "metHMDBID"]
        ) ~ paste0(
          "Metab__", 
          metab.map[gsub("Metab__", replacement = "", x = target), "metHMDBID"],
          "_", str_sub(target, start = nchar(target))
        ),
        !is.na(
          metab.map[gsub("Metab__", replacement = "", x = target), "metKEGGID"]
        ) ~ paste0(
          "Metab__", 
          metab.map[gsub("Metab__", replacement = "", x = target), "metKEGGID"],
          "_", str_sub(target, start = nchar(target))
        ), TRUE ~ target
      ), TRUE ~ target 
    )
  )
  
  return(
    list(
      GSMM.PKN = list.network[[1]], 
      mets.map = metab.map, 
      reac.to.gene = list.network[[3]], 
      reac.map = list.network[[4]]
    )
  )
}


.genes_to_symbol <- function(
    list.network, 
    organism, 
    mapping.biomart = NULL,
    ont.from = "ensembl_gene_id", 
    ont.to = "external_gene_name"
) {
  dataset.biomart <- switch(
    as.character(organism), 
    "9606" = "hsapiens_gene_ensembl",
    "10090" = "mmusculus_gene_ensembl",
    "10116" = "rnorvegicus_gene_ensembl",
    "7955" = "drerio_gene_ensembl",
    "7227" = "dmelanogaster_gene_ensembl",
    "6239" = "celegants_gene_ensembl"
  )
  if (is.null(dataset.biomart)) 
    stop("Chosen organism is not recognizable")
  
  regex <- "(Gene\\d+__)|(_reverse)"
  regex.reverse <- "_reverse"
  ## getting biomart info
  genes.GSMM <- grep("Gene\\d+__", unlist(list.network[[1]]), value = TRUE) %>% 
    gsub(regex, "", .) %>% 
    ifelse(grepl("_", .), sapply(strsplit(., split = "_"), \(x) x), .) %>% 
    unlist()
  if (is.null(mapping.biomart)) {
    ensembl.link <- useEnsembl(biomart = "ensembl", dataset = dataset.biomart)
    ensembl.df <- getBM(
      filters = ont.from, 
      attributes = c('ensembl_gene_id', 'external_gene_name'),
      values = genes.GSMM,
      mart = ensembl.link
    ) %>% unique()
    rownames(ensembl.df) <- ensembl.df$ensembl_gene_id  
  } else {
    ensembl.df <- mapping.biomart %>% select(-ensembl_peptide_id) %>% 
      unique() %>% filter(
        !is.na(.data[[ont.from]]), !is.na(.data[[ont.to]]),
      ) 
    rownames(ensembl.df) <- ensembl.df[[ont.from]]
  }
  ## translating genes when possible (genes not found are not removed but kept 
  # with ENSEMBL IDs)
  genes.logical <- data.frame(
    source = grepl("Gene\\d+__", list.network[[1]]$source),
    target = grepl("Gene\\d+__", list.network[[1]]$target)
  ) %>% as.matrix()
  ## modifying the matrix in place seems to be faster than dplyr code / check 
  # checking_functions_COSMOS.R file in case I need to come back
  list.network[[1]] <- as.matrix(list.network[[1]])
  for (i in seq(nrow(list.network[[1]]))) {
    gene <- list.network[[1]][i, ][genes.logical[i, ]]
    new.character <- paste0(
      "Gene", 
      gsub("\\D", "", sapply(strsplit(x = gene, split = "__"), \(x) x[1])), 
      "__", 
      gsub(regex, "", gene) %>% 
        strsplit(x = ., split = "_") %>% 
        sapply(
          \(iter) {
            ifelse(
              test = is.na(ensembl.df[iter, ont.to]), 
              yes = iter, 
              no = ensembl.df[iter, ont.to]
            ) %>% paste(collapse = "_")
          }
        ),
      ifelse(grepl(regex.reverse, gene), yes = "_reverse", no = "")
    )
    list.network[[1]][i, ][genes.logical[i, ]] <- new.character
  }
  
  list.network[[1]] <- as.data.frame(list.network[[1]])
  
  list.network[[3]] <- list.network[[3]] %>% mutate(
    Gene = case_when(
      !is.na(
        ensembl.df[gsub(regex, replacement = "", x = Gene), ont.to]
      ) ~ paste0(
        "Gene", 
        gsub("\\D", "", sapply(strsplit(x = Gene, split = "__"), \(x) x[1])), 
        "__", 
        ensembl.df[gsub(regex, replacement = "", x = Gene), ont.to]
      ), 
      TRUE ~ Gene
    )
  )
  
  return(
    list(
      GSMM.PKN = list.network[[1]], 
      mets.map = list.network[[2]], 
      reac.to.gene = list.network[[3]], 
      reac.map = list.network[[4]]
    )
  )
}

.format_GSMM_COSMOS <- function(
    list.network, 
    verbose = TRUE
) {
  reaction.network <- list.network[[1]]
  enzyme_reacs <- unique(c(reaction.network$source, reaction.network$target))
  enzyme_reacs <- enzyme_reacs[grepl("^Gene", enzyme_reacs)]
  enzyme_reacs_reverse <- enzyme_reacs[grepl("_reverse",enzyme_reacs)]
  enzyme_reacs <- enzyme_reacs[!grepl("_reverse",enzyme_reacs)]
  
  if (verbose) message("\t>>> Step 1: Defining transporters")
  
  new_df_list <- sapply(
    X = enzyme_reacs, 
    FUN = function(enzyme_reac, reaction.network) {
      df <- reaction.network[which(
        reaction.network$source == enzyme_reac | 
          reaction.network$target == enzyme_reac
      ),]
      if (dim(df)[1] < 2) {
        return(NA)
      } else {
        if (dim(df)[1] < 3) {
          return(df)
        } else {
          for (i in 1:dim(df)[1]) {
            if (grepl("Metab__", df[i, 1])) {
              counterpart <- which(
                gsub("_[a-z]$", "", df[, 2]) == gsub("_[a-z]$", "", df[i, 1])
              )
              if (length(counterpart) > 0) {
                df[i, 2] <- paste0(df[i, 2], paste0("_TRANSPORTER", i))
                df[counterpart, 1] <- paste0(
                  df[counterpart, 1], paste0("_TRANSPORTER", i)
                )
              }
            }
          }
          return(df)
        }
      }
    }, 
    reaction.network = reaction.network
  ) 
  new_df <- as.data.frame(do.call(rbind, new_df_list))
  
  if (verbose) message("\t>>> Step 2: Defining reverse reactions")
  
  new_df_reverse <- sapply(
    X = enzyme_reacs_reverse, 
    FUN = function(enzyme_reac_reverse, reaction.network) {
      df <- reaction.network[which(
        reaction.network$source == enzyme_reac_reverse | 
          reaction.network$target == enzyme_reac_reverse
      ),]
      if (dim(df)[1] < 2) {
        return(NA)
      } else {
        if (dim(df)[1] < 3) {
          return(df)
        } else {
          for (i in 1:dim(df)[1]) {
            if(grepl("Metab__",df[i,1])) {
              counterpart <- which(
                gsub("_[a-z]$","",df[,2]) == gsub("_[a-z]$","",df[i,1])
              )
              if(length(counterpart) > 0) {
                transporter <- gsub("_reverse", "", df[i, 2])
                transporter <- paste0(
                  transporter, paste0(paste0("_TRANSPORTER", i), "_reverse")
                )
                df[i, 2] <- transporter
                df[counterpart, 1] <- transporter
              }
            }
          }
          return(df)
        }
      }
    }, reaction.network = reaction.network
  ) 
  new_df_reverse <- as.data.frame(do.call(rbind, new_df_list))
  reaction.network.new <- as.data.frame(rbind(new_df, new_df_reverse))
  reaction.network.new <- reaction.network.new[complete.cases(reaction.network.new),]
  ## filter metabolites in mapping mets
  metabs <- c(
    grep("Metab__", reaction.network.new[[1]], value = TRUE),
    grep("Metab__", reaction.network.new[[2]], value = TRUE)
  ) %>% unique() %>% gsub("(Metab__)|(_[a-z])", "", .)
  list.network[[2]] <- list.network[[2]] %>% 
    filter(metHMDBID %in% metabs | metBiGGID %in% metabs | mets %in% metabs)
  
  if (verbose) message("\n\tDONE")
  
  return(
    list(
      GSMM.PKN = reaction.network.new, 
      mets.map = list.network[[2]], 
      reac.to.gene = list.network[[3]], 
      reac.map = list.network[[4]]
    )
  )
}

.connecting_GSMM_omnipath <- function(
  GSMM.PKN, 
  omnipath.PKN, 
  verbose = TRUE  
) {
  elements <- unique(as.character(unlist(GSMM.PKN)))
  elements <- elements[grepl("^Gene\\d+__", elements)]
  elements <- gsub("(.*__)|(_TRANSPORTER[0-9]+)|(_reverse$)", "", elements)
  ## this function can be vectorized
  connectors.df <- sapply(
    X = elements, 
    FUN = function(ele) {
      if (grepl("_", ele)) {
        genes.sep <- str_split(string = ele, pattern = "_")[[1]]
        if(length(genes.sep) < 10) {
          genes_connector_list <- sapply(
            X = genes.sep,
            FUN = function(gene) {
              return(c(gene, ele))
            }
          )
          return(t(genes_connector_list))
        } 
      } else {
        return(c(ele, ele))
      }
    }
  ) %>% do.call(rbind, .) %>% as.data.frame()
  names(connectors.df) <- c("source", "target")
  connectors.df <- connectors.df[which(
    connectors.df$source %in% omnipath.PKN$source | 
      connectors.df$source %in% omnipath.PKN$target
  ),]
  network.df.new <- as.data.frame(rbind(GSMM.PKN, connectors.df))
  
  return(network.df.new)
}

.formatSTITCH <- function(
  actions.file, 
  links.file, 
  organism, 
  omnipath.PKN,
  mapping.biomart = NULL,
  threshold = 700, 
  verbose = TRUE
) {
  dataset.biomart <- switch(
    as.character(organism), 
    "9606" = "hsapiens_gene_ensembl",
    "10090" = "mmusculus_gene_ensembl",
    "10116" = "rnorvegicus_gene_ensembl",
    "7955" = "drerio_gene_ensembl",
    "7227" = "dmelanogaster_gene_ensembl",
    "6239" = "celegants_gene_ensembl"
  )
  if (is.null(dataset.biomart)) 
    stop("Chosen organism is not recognizable")
  if (verbose) {
    message("\t>>> Reading provided STITCH files\n")
  }
  links.detail <- as.data.frame(
    read_delim(
      links.file, "\t", escape_double = FALSE, 
      trim_ws = TRUE, progress = verbose
    )
  ) %>% filter(
    combined_score >= threshold, experimental >= threshold | database >= threshold
  ) %>% mutate(
    ID = paste(chemical, protein, sep = "_"),
    ID_reverse = paste(protein, chemical, sep = "_")
  )
  STITCH <- as.data.frame(
    read_delim(
      actions.file, "\t", escape_double = FALSE, 
      trim_ws = TRUE, progress = verbose
    )
  ) %>% filter(mode == "activation" | mode == "inhibition", a_is_acting) %>% 
    mutate(ID = paste(item_id_a, item_id_b, sep = "_")) %>% 
    filter(ID %in% links.detail$ID | ID %in% links.detail$ID_reverse)
  STITCH <- STITCH[,-7]
  ## df of proteins in STICH
  prots <- unique(c(STITCH$item_id_a, STITCH$item_id_b))
  prots <- prots[grepl(paste0(organism, "\\."), prots)]
  prots <- as.data.frame(cbind(prots, gsub(paste0(organism, "\\."), "", prots)))
  colnames(prots) <- c("original", "ensembl_prots")
  ## getting info from Biomart
  if (verbose) message("\n\n\t>>> Using info from BiomaRt")
  
  if (is.null(mapping.biomart)) {
    ensembl.link <- useEnsembl(biomart = "ensembl", dataset = dataset.biomart)
    ensembl.df <- getBM(
      filters = "ensembl_peptide_id", 
      attributes = c(
        "ensembl_peptide_id",'ensembl_gene_id', 'external_gene_name'# , 'entrezgene_id', "description"
      ),
      values = prots[[2]],
      mart = ensembl.link
    )
    colnames(ensembl.df)[1] <- "ensembl_prots"  
  } else {
    ensembl.df <- mapping.biomart %>% filter(
      ensembl_peptide_id %in% prots[[2]]
    )
    colnames(ensembl.df)[1] <- "ensembl_prots"  
  }
  
  prots <- merge(prots, ensembl.df, by = "ensembl_prots")
  ## external_gene_name for mouse, Idk in other cases, check this
  prots <- prots[prots$external_gene_name != "",]
  prots.vec <- prots$external_gene_name
  names(prots.vec) <- prots$original

  if (verbose) message("\n\t>>> Generating PKN network")
  
  STITCH <- STITCH %>% mutate(
    item_id_a = case_when(
      grepl("\\.", item_id_a) & (item_id_a %in% names(prots.vec)) ~ 
        prots.vec[item_id_a],
      grepl("^CID", item_id_a) ~ gsub("CID[a-z]0*", "Metab__", item_id_a),
      TRUE ~ item_id_a
    ),
    item_id_b = case_when(
      grepl("\\.", item_id_b) & (item_id_b %in% names(prots.vec)) ~ 
        prots.vec[item_id_b],
      grepl("^CID", item_id_b) ~ gsub("CID[a-z]0*", "Metab__", item_id_b),
      TRUE ~ item_id_b
    ),
    sign = case_when(action == "inhibition" ~ -1, TRUE ~ 1)
  ) %>% dplyr::select(1, 2, 7)
  colnames(STITCH) <- c("source", "target", "sign") 
  CIDs <- unique(as.character(unlist(STITCH[,c(1,3)])))
  CIDs <- CIDs[grepl("Metab__", CIDs)] %>% gsub("Metab__", "", .)
  
  ## Convert CID to HMDB Id when available
  metabolitesMapping.mod <- metabolitesMapping %>% 
    filter(CID %in% CIDs, !is.na(HMDB)) %>% 
    mutate(HMDB = paste0("Metab__", HMDB))
  metabolitesMapping.vec <- metabolitesMapping.mod$HMDB
  names(metabolitesMapping.vec) <- paste0("Metab__",metabolitesMapping.mod$CID)
  
  ## metabolites with no HMDB are kept
  STITCH <- STITCH %>% mutate(
    source = case_when(
      grepl("Metab__", source) & source %in% names(metabolitesMapping.vec) ~ 
        metabolitesMapping.vec[source],
      TRUE ~ source
    ),
    target = case_when(
      grepl("Metab__", target) & target %in% names(metabolitesMapping.vec) ~ 
        metabolitesMapping.vec[target],
      TRUE ~ target
    )
  )
  # TODO: at this point, STITCH contains metabolites in both columns of the dataframe
  ## this should be checked
  omn.prots <- unique(as.character(unlist(omnipath.PKN[,c(1, 2)])))
  STITCH <- unique(STITCH[which(STITCH$target %in% omn.prots),])
  
  STITCH$source <- paste(STITCH$source, "_c", sep = "") 
  
  if (verbose) message("\n\tDONE")
  
  return(STITCH)
}

.mixing_resources <- function(GSMM.network, omnipath.PKN, stitch.PKN) {
  ## connecting Omnipath and GSMM
  GSMM.network <- .connecting_GSMM_omnipath(
    GSMM.PKN = GSMM.network, omnipath.PKN = omnipath.PKN
  )
  GSMM.network$sign <- 1
  meta.PKN <- as.data.frame(
    rbind(omnipath.PKN, stitch.PKN, GSMM.network)
  ) %>% unique()
  meta.network <- meta.PKN[,c(1, 3, 2)]
  names(meta.network) <- c("source", "interaction", "target")
  # TODO: manual correction: difficult to generalize for different organisms, shall I remove it? 
  #probably erroneous interaction (WHY?? this only works for human / mouse)
  # meta.network <- meta.network[-which(
  #   meta.network$source == "Prkca" & meta.network$target == "Src"
  # ),]
  # meta.PKN <- meta.PKN[-which(
  #   meta.PKN$source == "Prkca" & meta.PKN$target == "Src"
  # ),]
  #probably erroneous interaction
  # meta_network <- meta_network[-which(meta.network$source == "Ltc45"),] 
  #I don't know where this interaction comes from, the sources are wrong (https://www.nature.com/articles/onc2008228)
  # meta_network <- meta.network[!(grepl("Cad_reverse", meta.network$source) | grepl("Cad_reverse", meta.network$target)) ,] 
  #redHuman confirms that the reaction is actually not reversible: NOT FOUND IN MOUSE EITHER
  
  return(list(meta.PKN, meta.network))
}
