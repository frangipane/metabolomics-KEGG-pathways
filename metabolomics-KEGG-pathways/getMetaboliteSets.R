## getMetaboliteSets.R
## metaboliteSets() gets KEGG metabolic pathways associated with a user-provided vector of 
## metabolites provided as KEGG codes (e.g. Serine's KEGG code is "C00065").  The list of pathways
## retrieved from KEGG are filtered according to a threshold minimum number of measurements per pathway,
## the existence of the pathway in organism of interest, and some manual filtering.  Metabolite sets
## are constructed from the pathways and trimmed to only contain the measured metabolites.
##
##
## INPUT:
##  -metabolites: a vector of metabolites as KEGG codes
##    example: metabolites = c("C00065","C00222","C00788",...)
##  -threshold: an integer lower threshold for the number of metabolites
##    that must be measured to keep a pathway for analysis
##  -species: character (usually 3-letter) KEGG code for a particular species
##    (e.g. "hsa" = human, "mmu" = mouse, see http://www.genome.jp/kegg/catalog/org_list.html for more)
## OUTPUT:
##  -a histogram of metabolite.count, the number of metabolites per pathway
## RETURNS:
##  -metabolite_set: a list of metabolites per pathway, trimmed to only include
##    metabolites measured in the dataset
##  -path.attrib: list of manually filtered metabolic paths corresponding to the
##    returned metabolite_set, containing {entry = KEGG pathway identifier,
##    pathname = english name for pathway, compounds = metabolites present in pathway}
###########################################################################################################

source("getKEGGpaths.R")

getMetaboliteSets = function(metabolites, threshold=5, species="hsa") {
  
  ## get paths per metabolite (and other metadata)
  info.metab = getKEGGpaths(metabolites)
  
  ## count number of metabolites per path
  metabolite.count = metabolitesPerPath(info.metab)
  
  ## apply threshold for keeping pathways based on number of metabolites measured, and
  ## if they exist for for the species specified
  keep.paths = filterPathways(metabolite.count, threshold, species="hsa")
  
  ## manually filter additional paths that are global or irrelevant
  path.attrib = manualfilterPathways(keep.paths)
  
  ## return metabolite sets, trimmed only to contain metabolites that were measured in the
  ## dataset
  metabolite_set = trimMetaboliteSets(path.attrib, metabolites)
  
  return(list(path.attrib = path.attrib, metabolite_set = metabolite_set))
}
