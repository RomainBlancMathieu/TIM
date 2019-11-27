#!/Users/romain/anaconda_ete/bin/python
from __future__ import division
import os
import sys
import random

os.system("ls -ltrh plotItol_node_* | sed 's/plot/@/g' | cut -d '@' -f 2 > inputForPieChart.txt")
f=open("inputForPieChart.txt", "r")
seenNodes=[]
allGrp=set()
dico={}
for l in f:
	t=l.strip().split(".")[0].split("_")
	node=int(t[2])
	eukGrp=t[3]
	if node in dico.keys():
		dico[node].append(eukGrp)
		print "node %d seen more than once\n" % node
	else:
		dico[node]=[eukGrp]
	allGrp.add(eukGrp)
f.close()
print allGrp

# assign random color to euk grp
dicoCol={}
for euk in allGrp:
	color = "#%06x" % random.randint(0, 0xFFFFFF)
	dicoCol[euk]=color


fo=open("PIECHART_ITOL.txt", "w")
fo.write("DATASET_PIECHART\n")
fo.write("SEPARATOR COMMA\n")
fo.write("DATASET_LABEL,TIM\n")

fo.write("FIELD_COLORS,")
tmp=""
for k in dicoCol:
        tmp+="%s," % dicoCol[k]
fo.write("%s\n" % tmp[:-1])

fo.write("FIELD_LABELS,")
tmp=""
for k in dicoCol:
        tmp+="%s," % k
fo.write("%s\n" % tmp[:-1])

fo.write("MAXIMUM_SIZE,10\n")
fo.write("DATA\n")

# - - - - Write DATA - --
for node in dico:
	vector=''
	nbEukGrp=len(dico[node])
        for k in dicoCol.keys():
		if k in dico[node]:
			vector+="%f," % (1/nbEukGrp)
		else:
			vector+="0,"    

	fo.write("INT_%d,0.5,2,%s\n" % (node,vector[:-1]))

print "finito"
