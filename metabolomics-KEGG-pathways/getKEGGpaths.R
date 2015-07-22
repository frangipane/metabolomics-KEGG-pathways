library(KEGGREST)

## given a dataset with metabolite KEGG codes as row names, order of calls:

# row.names(dataset) = metab
# info.metab = getKEGGpaths(metab)
# metabolite.count = metabolitesPerPath(info.metab)
# keep.paths = filterPathways(metabolite.count, thresh=5, species="mmu")
# path.attrib = manualfilterPathways(keep.paths)
# metabolite_set = trimMetaboliteSets(path.attrib, metab)

##===================================================================================
## get paths per metabolite

## INPUT: metab: a vector of metabolites as KEGG codes
## RETURN: info.metab: list of KEGG info per metabolite
getKEGGpaths = function(metab) {
  if(!is.vector(metab)) stop("Metabolites must be provided as a vector.")
  
  ## check for duplicated compounds.  only query KEGG for unique list of compounds
  if(sum(duplicated(metab))) metab = unique(metab)
  
    ## keggGet only allows up to 10 queries at a time, so fetch info in batches of 10
  batch.idxs = c(seq(from=1, to=length(metab), by=10), length(metab)+1)
  info.metab = vector(mode="list", length=0)
  for (i in seq_along(batch.idxs[1:length(batch.idxs)-1])) {
    batch = keggGet(metab[batch.idxs[i]:(batch.idxs[i+1]-1)])
    #info.metab[batch.idxs[i]:(batch.idxs[i+1]-1)] = batch
    info.metab = c(info.metab, batch)
  }
  return(info.metab)
}

##===================================================================================
## get modules per metabolite
## consider analyzing modules instead of pathways?
# vec.modules = vector(mode="list", length=24)
# for (i in seq_along(fullquery)) {
#   vec.modules[[i]] = fullquery[[i]]$MODULE
# }
# sum(sapply(vec.modules, is.null))
## returns 6
## a significant number of compounds lack associated modules, but all exist in at least one pathway,
## example: no module contains ASN

##===================================================================================
## count number of metabolites per path

## INPUT: info.metab: a list containing KEGG data per metabolite (including paths involved)
## RETURN: number of metabolites per pathway as a table object
metabolitesPerPath = function(info.metab) {
  
  ## get names of pathways "map-----" for which a metabolite is present; apply to all metabolites in list
  pathnames = vector(mode="list", length=length(info.metab))
  for (i in seq_along(info.metab)) {
    pathnames[[i]] = info.metab[[i]]$PATHWAY
  }
  ## get a count of the number of metabolites in the dataset present in each pathway
  metabolite.count = table(names(unlist(pathnames)))
  
  ## show distribution of metabolite count
  print(hist(metabolite.count))

  return(metabolite.count)
}

##===================================================================================
## apply threshold for keeping pathways based on number of metabolites measured, and
## if they exist for for the species specified

## INPUT: -metabolite.count: count of metabolites per pathway as table object
##        -thresh: an integer threshold for minimum number of measured metabolites in order to keep pathway
## RETURN: -keep.paths: a list of KEGG paths with associated data, thresholded and
##          checked for existence in user-specified species
filterPathways = function(metabolite.count, thresh, species="mmu") {
  
  ## how many pathways contain at least thresh measured metabolites?
  npath = sum(metabolite.count >= thresh)
  print(paste0(npath, " pathways contain at least ", thresh, " metabolites."))
  
  ## get map names of the pathways containing at least n=thresh measured metabolites
  map.thresh = names(metabolite.count)[metabolite.count >= thresh]
  
  ## substitute "mmu" for "map" in pathway names
  mmu.thresh = gsub("map", species, map.thresh)
  
  ## check that all pathways exist specifically for the specified species
  path.exists = sapply(mmu.thresh, function(mmupath) tryCatch(keggGet(mmupath), error=function(e) {NA}))
  
  ## only keep paths that exist in "mmu"
  keep.paths = path.exists[!is.na(path.exists)]
  
  return(keep.paths)
}

##===================================================================================
## manually filter additional paths that are global or irrelevant

## INPUT: a list of KEGG paths with associated data, already thresholded and
## checked for existence in species
## RETURN: list of manually filtered paths containing {entry = KEGG pathway identifier,
## pathname = english name for pathway, compounds = metabolites present in pathway}
manualfilterPathways = function(keep.paths) {
  ## get names and compounds per path
  path.attrib = lapply(keep.paths,
                       function(x) return(list(entry = x$ENTRY,
                                               pathname = x$PATHWAY_MAP, 
                                               compounds = x$COMPOUND)))
  ## remove global paths (3 paths).
  ## global paths do not have compounds list, so their compounds property will be NULL
  ## mmu01100 = "metabolic pathways"
  ## mmu01230 = "biosynthesis of amino acids"
  ## mmu01200 = "carbon metabolism"
  global.idxs = sapply(path.attrib, function(x) is.null(x$compounds))
  path.attrib = path.attrib[!global.idxs]
  
  ## manually remove some other irrelevant pathways
  ## mmu04974 = protein digestion and absorption
  ## mmu04978 = mineral absorption
  ## mmu02010 = ABC transporters
  ## mmu00970 = Aminoacyl-tRNA biosynthesis
  irrel.idxs = sapply(path.attrib, function(x) x$entry %in% c("mmu04974", "mmu04978", "mmu02010","mmu00970"))
  path.attrib = path.attrib[!irrel.idxs]

  print(paste0(length(path.attrib), " paths remaining after thresholding and manual filtering."))
  return(path.attrib)
}

##===================================================================================
## return metabolite sets

## INPUT: -list of manually filtered paths containing {entry = KEGG pathway identifier,
##        pathname = english name for pathway, compounds = metabolites present in pathway} returned
##        by manualfilterPathways()
##        -metab: a vector of metabolites as KEGG codes measured in the dataset
## RETURN: -metabolite_set: a list of metabolites per pathway, trimmed to only include
##          metabolites measured in the dataset
trimMetaboliteSets = function(path.attrib, metab) {
  ## compounds per path
  cpds_per_path = sapply(path.attrib, function(x) x$compounds)
  
  ## remove compounds in each path that are not present in the dataset
  metabolite_set = sapply(cpds_per_path, 
                          function(x) {
                            idxs = names(x) %in% metab
                            return(x[idxs])
                          })
  return(metabolite_set)
}