#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library(dplyr)
library(gplots)
# - - - - - - - - - - - - - - -  Open and format files - - - - - - - - - - - - - -:
h<-function(DF, x, y) {return(DF[c(x:y),c(x:y)])}

d=read.table("myTIMrun/allNodesCounts.txt", header=T, sep="\t")

d2=d %>% group_by(taxa) %>% mutate(corrPval=p.adjust(Pval, method = args[1], n = length(Pval))) %>% as.data.frame

write.table(d2, sep="\t", file="correctedPvalues.txt", row.names=F, quote=F)

# test itol:
d3=d2 %>% filter(corrPval<0.05 & nbLeaves>1 & (nbUniqPolBperTaxa/nbLeaves>0.5 | dist_average<0.02)) %>% select(node, taxa)
#d3=d2 %>% filter(corrPval<0.05 & nbLeaves>1) %>% select(node, taxa)
#d3=d2 %>% filter(corrPval<0.05 & nbLeaves>1 & (nbUniqPolBperTaxa/nbLeaves>0.5 | dist_farthest<0.025)) %>% select(node, taxa)

print("number of node with enrichement:")
print(nrow(d3))

taxaNames=unique(d3$taxa)
taxaColors=rainbow(length(taxaNames))
names(taxaColors)<-taxaNames

for (row in 1:nrow(d3)) {
	node=d3[row, "node"]
	taxa=as.character(d3[row, "taxa"])
	#color=col2hex(as.character(taxaColors[taxa]))
	color=as.character(taxaColors[taxa])
	system(sprintf("cat myTIMrun/NODES/node_%s.leaves.txt | sed 's/$/ %s '%s'/g' | grep '^[0-9]' > itol_node_%s_%s.txt", node, color, taxa, node, taxa))
	system(sprintf("cat scripts/headerItol.txt | sed 's/EukGrp/%s_%s/g' > tmpHeaderItol.txt", node, taxa))
	system(sprintf("cat tmpHeaderItol.txt itol_node_%s_%s.txt > plotItol_node_%s_%s.txt", node, taxa, node, taxa))
}


