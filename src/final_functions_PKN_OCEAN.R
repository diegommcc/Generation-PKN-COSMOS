################################################################################
## Description: set of functions to generate OCEAN PKN for different organisms 
## using GSMM information from <https://www.pnas.org/doi/10.1073/pnas.2102344118>
## Author: Diego Mañanes
## Date: 23/11/27
################################################################################

## set of dependencies
# suppressMessages(library("biomaRt"))
# suppressMessages(library("OmnipathR"))
# suppressMessages(library("readr"))
# suppressMessages(library("stringr"))
# suppressMessages(library("metaboliteIDmapping"))
# suppressMessages(library("R.matlab"))
# suppressMessages(library("dplyr"))

## SIF files

## exported
create_PKN_OCEAN <- function(
  GSMM.matlab.path,
  GSMM.reactions.map,
  GSMM.metabolites.map,
  KEGG.compounds = NULL,
  GSMM.reactions.map.col = "rxns", 
  GSMM.metabolites.map.col = "mets",
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
  translate.genes = FALSE,
  organism = NULL,
  genes.from = "ensembl_gene_id",
  genes.to = "external_gene_name",
  verbose = TRUE
) {
  
  ## checking if files exist
  if (!file.exists(GSMM.matlab.path)) {
    stop("GSMM.matlab.path file does not exist")
  }  else if (is.null(KEGG.compounds) & verbose) {
    message(">>> Using internal list of cofactors determined using KEGG")
  }
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
  ## for one of the organisms, tihs part will be needed, but not implemented here
  if (translate.genes){
    if (!organism %in% c(9606, 10090, 10116, 7955, 7227, 6239)) {
      stop(
        paste("organism is not recognizable. Supported organisms are", 
              paste(c(9606, 10090, 10116, 7955, 7227, 6239), collapse = ", "))
      )
    }
    gsmm.PKN.list <- .genes_to_symbol(
      gsmm.PKN.list, organism = organism, 
      ont.from = genes.from, ont.to = genes.to
    )
  }
  gsmm.PKN.list <- .remove_cofactors(
    list.network = gsmm.PKN.list, 
    KEGG.compounds = KEGG.compounds,
    verbose = verbose
  )
  if (verbose) message(">>> Compressing transporters...")
  gsmm.PKN.list <- .compress_transporters(list.network = gsmm.PKN.list)
  ## check this with Aurelien: removing during cofactor identification
  gsmm.PKN.list <- .split_transaminases(list.network = gsmm.PKN.list)
  if (verbose) message(">>> Translating metabolites...")
  gsmm.PKN.list <- .mets_to_HMDB(list.network = gsmm.PKN.list)
  
  if (verbose) message("\nDONE")
  
  return(gsmm.PKN.list)
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
  )
  reaction.ids <- unlist(
    mat.object[[which(attribs.mat == list.params.GSMM$reaction.ID.name)]]
  )
  ## reaction to genes df: 
  ##TODO: if length == 0, then no gene name, that's why sometimes there are genes with a number
  ## ask about this
  reaction.to.genes.df <- lapply(
    seq_along(reaction.list),
    \(idx) {
      genes.reac <- unlist(reaction.list[[idx]], recursive = FALSE)
      if (length(genes.reac) != 0) {
        genes <- unique(
          gsub(
            " and ", "_", 
            gsub(
              "[()]","", 
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
    if (as.vector(lbs[reac.idx, 4] == "forward")) {
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

.vec_cofactors <- c(
  'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00007', 'C00008', 
  'C00009', 'C00010', 'C00011', 'C00012', 'C00013', 'C00014', 'C00015', 'C00016', 
  'C00017', 'C00018', 'C00019', 'C00020', 'C00021', 'C00022', 'C00023', 'C00024', 
  'C00025', 'C00026', 'C00027', 'C00028', 'C00029', 'C00030', 'C00031', 'C00032', 
  'C00033', 'C00034', 'C00035', 'C00036', 'C00037', 'C00038', 'C00039', 'C00040', 
  'C00041', 'C00042', 'C00043', 'C00044', 'C00045', 'C00046', 'C00047', 'C00048', 
  'C00049', 'C00050', 'C00051', 'C00052', 'C00053', 'C00054', 'C00055', 'C00058', 
  'C00059', 'C00060', 'C00061', 'C00062', 'C00063', 'C00064', 'C00065', 'C00066', 
  'C00067', 'C00068', 'C00069', 'C00070', 'C00071', 'C00072', 'C00073', 'C00074', 
  'C00075', 'C00076', 'C00077', 'C00078', 'C00079', 'C00080', 'C00081', 'C00082', 
  'C00083', 'C00084', 'C00085', 'C00086', 'C00087', 'C00088', 'C00089', 'C00090', 
  'C00091', 'C00092', 'C00093', 'C00094', 'C00095', 'C00096', 'C00097', 'C00098', 
  'C00099', 'C00100', 'C00101', 'C00102', 'C00103', 'C00104', 'C00105', 'C00106', 
  'C00107', 'C00108', 'C00109', 'C00110', 'C00111', 'C00112', 'C00113', 'C00114', 
  'C00116', 'C00117', 'C00118', 'C00119', 'C00120', 'C00121', 'C00122', 'C00123', 
  'C00124', 'C00125', 'C00126', 'C00127', 'C00128', 'C00129', 'C00130', 'C00131', 
  'C00132', 'C00133', 'C00134', 'C00135', 'C00136', 'C00137', 'C00138', 'C00139', 
  'C00140', 'C00141', 'C00143', 'C00144', 'C00145', 'C00146', 'C00147', 'C00148', 
  'C00149', 'C00150', 'C00151', 'C00152', 'C00153', 'C00154', 'C00155', 'C00156', 
  'C00157', 'C00158', 'C00159', 'C00160', 'C00161', 'C00162', 'C00163', 'C00164', 
  'C00165', 'C00166', 'C00167', 'C00168', 'C00169', 'C00170', 'C00171', 'C00172', 
  'C00173', 'C00174', 'C00175', 'C00176', 'C00177', 'C00178', 'C00179', 'C00180', 
  'C00181', 'C00182', 'C00183', 'C00184', 'C00185', 'C00186', 'C00187', 'C00188', 
  'C00189', 'C00190', 'C00191', 'C00192', 'C00193', 'C00194', 'C00195', 'C00196', 
  'C00197', 'C00198', 'C00199', 'C00200', 'C00201', 'C00202', 'C00203', 'C00204', 
  'C00205', 'C00206', 'C00207', 'C00208', 'C00209', 'C00210', 'C00211', 'C00212', 
  'C00213', 'C00214', 'C00215', 'C00216', 'C00217', 'C00218', 'C00219', 'C00220', 
  'C00221', 'C00222', 'C00223', 'C00224', 'C00225', 'C00226', 'C00227', 'C00229', 
  'C00230', 'C00231', 'C00232', 'C00233', 'C00234', 'C00235', 'C00236', 'C00237', 
  'C00238', 'C00239', 'C00240', 'C00241', 'C00242', 'C00243', 'C00244', 'C00245',
  'C00246', 'C00247', 'C00248', 'C00249', 'C00250', 'C00251', 'C00252', 'C00253', 
  'C00254', 'C00255', 'C00256', 'C00257', 'C00258', 'C00259', 'C00261', 'C00262', 
  'C00263', 'C00264', 'C00265', 'C00266', 'C00267', 'C00268', 'C00269', 'C00270', 
  'C00272', 'C00273', 'C00275', 'C00279', 'C00280', 'C00282', 'C00283'
)

## this function removes cofactors based on BRITE information. 
# I have to improve it: check inputs, etc. 
# maybe creating a class for this could be a good idea... 
# transaminases are removed here. Ask Aurelien
.remove_cofactors <- function(
    list.network, 
    KEGG.compounds = NULL,
    mets.map.kegg.col = "metKEGGID",
    verbose = TRUE
) {
  kegg.ids.GSMM <- unique(
    na.omit(list.network[[2]][[mets.map.kegg.col]])
  )
  ## check if KEGG.compounds is not null: use internal lisst of cofactors
  if (!is.null(KEGG.compounds)) {
    ## check how many metabolites are shared in both objects
    valid.mets <- intersect(names(KEGG.compounds), kegg.ids.GSMM)
    if (verbose) {
      message(
        paste(
          ">>> Checking BRITE info of", 
          length(valid.mets), "metabolites"
        )
      )
    }
    ## getting bad mets
    bad_kegg_compounds <- Filter(
      f = \(compound) {
        brite <- compound$BRITE
        compound$NAME <- gsub(";", "", compound$NAME)
        if (!is.null(compound$ATOM)) {
          brite <- gsub("[ ]+", "", brite)
          if(("Cofactors" %in% brite | 
              "Nucleotides" %in% brite | 
              "CO2" %in% compound$NAME | 
              as.numeric(compound$ATOM[1]) <= 3 | 
              "ITP" %in% compound$NAME | 
              "IDP" %in% compound$NAME | 
              "NADH" %in% compound$NAME | 
              "NADPH" %in% compound$NAME | 
              "FADH2" %in% compound$NAME) &
             !"SAM" %in% compound$NAME)
          {
            return(TRUE)
          }
        }  
      },
      x = KEGG.compounds[valid.mets]
    )
    bad_entries <- sapply(bad_kegg_compounds, \(x) x$ENTRY[1])
    ## including diphosphate
    bad_entries <- unique(c(bad_entries, "C00013", "C00009"))
  } else {
    bad_entries <- intersect(.vec_cofactors, kegg.ids.GSMM)
  }
  if (verbose) {
    message(
      paste(
        "\n>>>", length(bad_entries), 
        "metabolites have been identified as cofactors"
      )
    )
  }
  ## removing cofactors from network
  mets.map.bad <- list.network[[2]] %>% 
    filter(.data[[mets.map.kegg.col]] %in% bad_entries)
  network.filt <- list.network[[1]] %>% mutate(
    source = case_when(
      grepl("Metab__", source) & 
        gsub("Metab__", "", source) %in% rownames(mets.map.bad) ~
        NA,
      TRUE ~ source
    ),
    target = case_when(
      grepl("Metab__", target) & 
        gsub("Metab__", "", target) %in% rownames(mets.map.bad) ~
        NA,
      TRUE ~ target
    )
  ) %>% filter(!is.na(source), !is.na(target))
  # reaction.network.no.cofact <- list.network[[1]][
  #   !(reaction.network$source %in% bad_entries) & 
  #     !(reaction.network$target %in% bad_entries),
  # ]
  mets.map.f <- list.network[[2]][!list.network[[2]][[mets.map.kegg.col]] %in% bad_entries, ]
  nodes <- unique(c(network.filt$source, network.filt$target))
  node.attributes <- data.frame(
    Node = nodes,
    Type = ifelse(grepl("Gene\\d+__", nodes), "reaction_enzyme", "metabolite")
  )
  return(
    list(
      GSMM.PKN = network.filt, 
      mets.map = mets.map.f, 
      reac.to.gene = list.network[[3]], 
      reac.map = list.network[[4]],
      attributes = node.attributes
    )
  )
}

.compress_transporters <- function(list.network) {
  test_1 <- paste(
    with(
      list.network[[1]], ifelse(
        grepl(pattern = "Metab__", source), 
        str_sub(source, start = 1L, end = 15L), 
        source
      )
    ),
    with(
      list.network[[1]], ifelse(
        grepl(pattern = "Metab__", target), 
        str_sub(target, start = 1L, end = 15L), 
        target
      )
    ),
    sep = "_"
  )
  test_2 <- paste(
    with(
      list.network[[1]], ifelse(
        grepl(pattern = "Metab__", target), 
        str_sub(target, start = 1L, end = 15L), 
        target
      )
    ),
    with(
      list.network[[1]], ifelse(
        grepl(pattern = "Metab__", source), 
        str_sub(source, start = 1L, end = 15L), 
        source
      )
    ),
    sep = "_"
  )
  transporters <- unique(
    c(list.network[[1]][test_1 %in% test_2, 1], 
      list.network[[1]][test_1 %in% test_2, 2])
  )
  network.transpor <- list.network[[1]] %>% mutate(
    source = case_when(
      source %in% transporters ~ "transporter", 
      TRUE ~ source
    ),
    target = case_when(
      target %in% transporters ~ "transporter", 
      TRUE ~ target
    ), 
    edgeId = paste(source, target, sep = "_")
  ) %>% filter(!duplicated(edgeId))
  edge_transporter <- grep("transporter", network.transpor$edgeId, value = TRUE)
  
  groups <- rep(0, length(edge_transporter))
  metab_memory <- c()
  for (i in 1:length(edge_transporter)) {
    if (grepl("Metab__", edge_transporter[i])) {
      metab <- edge_transporter[i]
      if (!metab %in% metab_memory) {
        metab_memory <- c(metab_memory, metab)
        idx <- grep(metab, edge_transporter)
        if (length(idx) != 0) groups[idx] <- i
      }
    }
  }
  names(groups) <- edge_transporter
  network <- network.transpor %>% mutate(
    source = case_when(
      source == "transporter" ~ paste(source, groups[edgeId], sep = ">"),
      TRUE ~ source
    ),
    target = case_when(
      target == "transporter" ~ paste(target, groups[edgeId], sep = ">"),
      TRUE ~ target
    )
  ) %>% filter(  ## only valid if symbol, not always true
    !grepl("^Slc", edgeId, ignore.case = TRUE) & 
      !grepl("_Slc", edgeId, ignore.case = TRUE)
  ) 
  network.transpor.def <- network[, -3]
  
  node.attributes <- list.network[[5]] %>% filter(
    Node %in% network.transpor.def[[1]] | Node %in% network.transpor.def[[2]]
  )
  new_transporters <- grep(
    "transporter", 
    unique(c(network.transpor.def[, 1], network.transpor.def[, 2])), 
    value = TRUE
  )
  new_transporters <- data.frame(
    Node = new_transporters, Type = "transporter"
  )
  node.attributes <- as.data.frame(rbind(node.attributes, new_transporters))
  
  return(
    list(
      GSMM.PKN = network.transpor.def, 
      mets.map = list.network[[2]], 
      reac.to.gene = list.network[[3]], 
      reac.map = list.network[[4]],
      attributes = node.attributes
    )
  )
}

## transaminases are removed during cofactor identification 
# therefore, I don't know if this makes sense...
.split_transaminases <- function(
    list.network
){
  ## check how many transaminases we have
  sub_net <- list.network[[1]]
  mets.trans <- list.network[[2]] %>% filter(
    grepl("C0002[56]", metKEGGID) | grepl("C00064", metKEGGID)
  ) 
  equiv <- list(
    C00025 = mets.trans %>% filter(metKEGGID == "C00025") %>% pull(mets) %>% 
      paste0("Metab__", .),
    C00026 = mets.trans %>% filter(metKEGGID == "C00026") %>% pull(mets) %>% 
      paste0("Metab__", .),
    C00064 = mets.trans %>% filter(metKEGGID == "C00064") %>% pull(mets) %>% 
      paste0("Metab__", .)
  )
  if (identical(mets.trans, character(0))) {
    message(">>> There are no transaminases in the network")
    return(list.network)
  }
  transaminases_net <- list.network[[1]] %>% filter(
    source %in% paste0("Metab__", mets.trans[["mets"]]) | 
      target %in% paste0("Metab__", mets.trans[["mets"]])
  )
  no_transaminases_net <- list.network[[1]] %>% filter(
    !source %in% paste0("Metab__", mets.trans[["mets"]]) | 
      !target %in% paste0("Metab__", mets.trans[["mets"]])
  )
  
  transaminases_potential <- grep(
    "Metab__", unique(c(transaminases_net$source, transaminases_net$target)), 
    value = TRUE, invert = TRUE
  )
  transaminases_net <- list.network[[1]] %>% filter(
    source %in% transaminases_potential | target %in% transaminases_potential
  )
  
  splitted_transaminase <- sapply(
    X = transaminases_potential, 
    FUN = function(transaminase, transaminases_net) {
      sub_net_transaminase <- transaminases_net[transaminase == transaminases_net$source |
                                                  transaminase == transaminases_net$target,]
      columns <- c(1,2)
      elements <- unique(c(sub_net_transaminase[,1],sub_net_transaminase[,2]))
      if (sum(elements %in% equiv[["C00025"]]) & sum(elements %in% equiv[["C00026"]])) {
        for(i in 1:length(sub_net_transaminase[, 1])) {
          for(j in 1:2) {
            if (sub_net_transaminase[i,j] %in% unlist(equiv[c("C00025", "C00026")])) {
              sub_net_transaminase[i, columns[-j]] <- paste(
                sub_net_transaminase[i, columns[-j]], "_gluakg", sep = ""
              )
            }
          }
        }
      }
      if (sum(elements %in% equiv[["C00025"]]) & sum(elements %in% equiv[["C00064"]])) {
        for(i in 1:length(sub_net_transaminase[,1])) {
          for(j in 1:2) {
            if(sub_net_transaminase[i,j] %in% unlist(equiv[c("C00025", "C00064")])) {
              sub_net_transaminase[i,columns[-j]] <- paste(
                sub_net_transaminase[i, columns[-j]], "_glugln", sep = ""
              )
            }
          }
        }
      }
      return(sub_net_transaminase)
    }, 
    transaminases_net = transaminases_net, USE.NAMES = TRUE, simplify = F
  ) %>% do.call(rbind, .)
  
  reaction.net <- unique(rbind(no_transaminases_net, splitted_transaminase))
  ifelse(
    grepl("Metab__", list.network[[1]][,1]), "metabolite", "reaction_enzyme"
  )
  node.attributes <- data.frame(
    Node = unique(c(reaction.net[, 1], reaction.net[, 2]))
  ) %>% mutate(
    case_when(
      grepl("Metab__", Node) ~ "metabolite", 
      TRUE ~ "reaction_enzyme"
    )
  )
  
  return(
    list(
      GSMM.PKN = reaction.net, 
      mets.map = list.network[[2]], 
      reac.to.gene = list.network[[3]], 
      reac.map = list.network[[4]],
      attributes = node.attributes
    )
  )
}

## only works for those mapping dfs from the ddbb
.mets_to_HMDB <- function(list.network) {
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
    ont.from = "ensembl_gene_id", ## check if correect?? funciton is not for users...
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
  ## getting biomart info
  genes.GSMM <- grep("Gene\\d+__", unlist(list.network[[1]]), value = TRUE) %>% 
    gsub(regex, "", .) %>% 
    ifelse(grepl("_", .), sapply(strsplit(., split = "_"), \(x) x), .) %>% 
    unlist()
  ensembl.link <- useEnsembl(biomart = "ensembl", dataset = dataset.biomart)
  ensembl.df <- getBM(
    filters = ont.from, 
    attributes = c('ensembl_gene_id', 'external_gene_name'),
    values = genes.GSMM,
    mart = ensembl.link
  ) %>% unique()
  rownames(ensembl.df) <- ensembl.df$ensembl_gene_id
  ## translating genes when possible (not found genes are not removed)
  ## when complexes are present (several genes concatenated), this code does not work
  list.network[[1]] <- list.network[[1]] %>% mutate(
    source = case_when(
      ## cases with a single gene
      grepl("Gene\\d+__", source) ~ case_when(
        !is.na(
          ensembl.df[gsub(regex, replacement = "", x = source), ont.to]
        ) ~ paste0(
          "Gene", 
          gsub("\\D", "", sapply(strsplit(x = source, split = "__"), \(x) x[1])), 
          "__", 
          ensembl.df[gsub(regex, replacement = "", x = source), ont.to]
        ), 
        ## cases with complexes: more than 1 gene
        grepl("[0-9]_[E]", source) ~ 
          paste0(
            "Gene", 
            gsub("\\D", "", sapply(strsplit(x = target, split = "__"), \(x) x[1])), 
            "__",
            unlist(
              strsplit(
                gsub(
                  pattern = "reverse", replacement = "", 
                  grep("[0-9]_[E]", source, value = T)[1]
                ), 
                split = "_"
              )
            )[-c(1:2)] %>% ensembl.df[., ont.to] %>% paste(collapse = "_")
          ),
        TRUE ~ source
      ), TRUE ~ source
    ),
    target = case_when(
      ## cases with a single gene
      grepl("Gene\\d+__", target) ~ case_when(
        !is.na(
          ensembl.df[gsub(regex, replacement = "", x = target), ont.to]
        ) ~ paste0(
          "Gene", 
          gsub("\\D", "", sapply(strsplit(x = target, split = "__"), \(x) x[1])), 
          "__", 
          ensembl.df[gsub(regex, replacement = "", x = target), ont.to]
        ), 
        ## cases with complexes: more than 1 gene
        grepl("[0-9]_[E]", target) ~ 
          paste0(
            "Gene", 
            gsub("\\D", "", sapply(strsplit(x = target, split = "__"), \(x) x[1])), 
            "__",
            unlist(
              strsplit(
                gsub(
                  pattern = "reverse", replacement = "", 
                  grep("[0-9]_[E]", target, value = T)[1]
                ), 
                split = "_"
              )
            )[-c(1:2)] %>% ensembl.df[., ont.to] %>% paste(collapse = "_")
          ),
        TRUE ~ target
      ), TRUE ~ target
    )
  ) 
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


#### this function is not used, it removes those links with no KEGG metabolites
## this function returns the reaction network including only metabolites
# with KEGG ID + included in KEGG compounds. In addition, the 
# mets IDs are transformed into KEGG instead of keeping MetAtlas
.metsKEGG <- function(
    reacction.network, 
    mets.map,
    KEGG.compounds,
    mets.map.id.col = "mets",
    mets.map.kegg.col = "metKEGGID"
) {
  ## check that KEGG.compounds is correctly built
  
  ## 
  compounds <- unique(c(reacction.network$source, reacction.network$target))
  compounds <- compounds[grepl("Metab__", x = compounds)] %>% 
    gsub(pattern = "Metab__", replace = "", x = .) %>% as.data.frame()
  kegg.in.map <- mets.map[[mets.map.kegg.col]][!is.na(mets.map[[mets.map.kegg.col]])]
  
  ## filtering reaction network acconrding to KEGG in mapping df
  which.met <- which(mets.map[[mets.map.kegg.col]] %in% names(KEGG.compounds))
  mets.map.kegg <- mets.map[which.met, ]
  reacction.network.kegg <- reacction.network[
    reacction.network$source %in% 
      paste0("Metab__", mets.map.kegg[[mets.map.id.col]]) | 
      reacction.network$target %in% 
      paste0("Metab__", mets.map.kegg[[mets.map.id.col]]),
  ] 
  ## transforming mets IDs into kegg IDs 
  mm.source <- match(
    reacction.network.kegg$source, 
    paste0("Metab__", mets.map.kegg[[mets.map.id.col]])
  )
  mm.source <- mm.source[!is.na(mm.source)]
  mm.target <- match(
    reacction.network.kegg$target, 
    paste0("Metab__", mets.map.kegg[[mets.map.id.col]])
  )
  mm.target <- mm.target[!is.na(mm.target)]
  
  reacction.network.kegg.v <- reacction.network.kegg
  reacction.network.kegg.v[
    grepl("Metab__", reacction.network.kegg.v$source), "source"
  ] <- mets.map.kegg[mm.source, mets.map.kegg.col]
  reacction.network.kegg.v[
    grepl("Metab__", reacction.network.kegg.v$target), "target"
  ] <- mets.map.kegg[mm.target, mets.map.kegg.col]
  
  ## making unique reactions: there are some cases where diff MetAtlas IDs points
  # toward the same KEGG ID
  reacction.network.kegg.v <- unique(reacction.network.kegg.v)
  
  return(list(reacction.network.kegg.v, mets.map.kegg))
}
