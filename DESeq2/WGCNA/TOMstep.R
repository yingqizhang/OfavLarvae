library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
library(flashClust)

lnames = load(file="SamplesAndTraits_OutliersRemoved.RData") 

softPower=4
adjacency2019 = adjacency(datExprOut2019, power=softPower,type="signed")
# translate the adjacency into topological overlap matrix (TOM) and calculate the corresponding dissimilarity
TOM2019= TOMsimilarity(adjacency2019,TOMType = "signed")
dissTOM2019= 1-TOM2019
# Generate a clustered gene tree
geneTree2019= flashClust(as.dist(dissTOM2019), method="average")

adjacency2021 = adjacency(datExprOut2021, power=softPower,type="signed")
TOM2021= TOMsimilarity(adjacency2021,TOMType = "signed")
dissTOM2021= 1-TOM2021
geneTree2021= flashClust(as.dist(dissTOM2021), method="average")

save(dissTOM2019, dissTOM2021, geneTree2019, geneTree2021, file="TOM_Output.RData")

