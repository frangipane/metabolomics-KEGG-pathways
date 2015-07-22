## viewpathways.R
## viewpathways() plots KEGG pathways associated with a dataset provided by the user using the
## "pathview" library. The KEGG pathways are selected by getMetaboliteSets().
##
## INPUT:
##  -filename, a dataset in .csv format organized in the following way:
##    -Metabolites are rows, labeled by KEGG codes (e.g. Serine's KEGG code is "C00065")
##    -Samples are columns
##    -First row contains header names ("row.names", "vehicle1", "vehicle2", "treatment1", etc.)
##    -Actually, only the first column containing KEGG codes is necessary 
##      (since nodes are colored according to a dummy sample column)
##  -threshold: an integer lower threshold for the number of metabolites
##    that must be measured to keep a pathway for analysis
##  -species: character (usually 3-letter) KEGG code for a particular species
##    (e.g. "hsa" = human, "mmu" = mouse, see http://www.genome.jp/kegg/catalog/org_list.html for more)
## OUTPUT:
##  -one image file (.png) per pathway, saved to working directory

###########################################################################################################

library("pathview")
source("getMetaboliteSets.R")

viewpathways = function(filename, threshold=5, species="hsa") {
  
  #filename = "vector_kegg_codes.csv"
  dataset = read.table(filename, sep=",", header=T, row.names=1)
  
  ## user-specified properties for selecting metabolic pathways for the dataset
  metabolites = row.names(dataset)
  
  ## get names of metabolite pathways e.g. "hsa04974" (for protein digestion)
  infoset = getMetaboliteSets(metabolites, threshold=threshold, species=species)
  path.names = names(infoset$metabolite_set)
  
  ## for calling pathview function, remove species prefix from path.names
  path.names = gsub("[A-Za-z]+", "", path.names)
  
  ## create dummy data frame for coloring in nodes
  dummy.df = data.frame(dummyvar = rep(1, nrow(dataset)), row.names=row.names(dataset))
  
  dir.create("pathview_graphs")
  
  ## Generate a kegg.native graph per pathway. Nodes are colored 
  ## using dummy dataframe.
  for (path_i in path.names) {
    pv.out = pathview(cpd.data=dummy.df, 
                      cpd.idtype = "kegg",
                      pathway.id = path_i,
                      species= species, 
                      out.suffix="mydata",
                      kegg.native=T,
                      multi.state=F,
                      limit=list(gene=1,cpd=1),
                      kegg.dir="./pathview_graphs/",
                      plot.col.key = F)
  }
}