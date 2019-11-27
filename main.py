#!/Users/romain/anaconda_ete/bin/python
from __future__ import division
import os
import sys
from ete3 import Tree
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import subprocess
import scipy.stats as stats

if len(sys.argv) <= 1:
 exit('enter name of a newick tree (without internal nodes ID)')
t = Tree(sys.argv[1])

if len(sys.argv) <= 2:
 exit('enter name of file containing connections between query and subject')
C = "%s" % sys.argv[2]

if len(sys.argv) <= 3:
 exit("What connections to analyze?\nfor positives enter 'POS'\nfor negatives enter 'NEG'\nfor both enter 'ALL'")
else:
 DIRECTION=sys.argv[3]

if DIRECTION=="POS":
 os.system("awk '$3>0' %s > tmpConnections.txt" % C)
elif DIRECTION=="NEG":
 os.system("awk '$3<0' %s > tmpConnections.txt" % C)
elif DIRECTION=="ALL":
 os.system("cp %s tmpConnections.txt" % C)
else:
 exit("Second argument should be 'POS', 'NEG' or 'ALL'")

USED_FILE="tmpConnections.txt"


def getNameGivenTAXIDandRank(TAXID, RANK):
        li=ncbi.get_lineage(TAXID)
        dicoRank=ncbi.get_rank(li)
        #print dicoRank
        if RANK in dicoRank.values():
                taxid_for_targetedRank=dicoRank.keys()[dicoRank.values().index(RANK)]
                name_for_targetedRank=ncbi.get_taxid_translator([taxid_for_targetedRank]).values()[0]
        else:
                name_for_targetedRank="NA"

        return name_for_targetedRank

def getAllNames():
        allMyNCBItaxid=ncbi.get_descendant_taxa('root', intermediate_nodes=True)
        allMyNCBInames=ncbi.translate_to_names(allMyNCBItaxid)
        SET=set()
        for i in allMyNCBInames:
                SET.add(i)
        return SET

# - - function to get average diatance between leaves and a node:
def avLeavesToNodeDist(NODE):
        sumDist=0.0
        cpt=0
        for leaf in NODE.get_leaves():
                dist = node.get_distance(leaf)
                sumDist+=dist
                cpt+=1
        avDist=sumDist/cpt
        return avDist


# - - function to read connections file:
def getCoPerEukGrp(FILE):
        # - - READ THE CONNECTION FILE FOR A NODE - -
        f=open(FILE, "r")
        level="order"; level2="class"; level3="phylum"
        listLevel=[]
        for  l in f:
                #OM-RGC.v1.000065429    4202    0.5330  0.0099  Eukaryota|...|Cryptosporidium+fragile   core Apicomplexa
                t=l.strip().split("\t")
                polb=t[0]; otu=t[1]; lineage=t[4]
                LML=lineage.split("|")[-1] #LML mean lowest meaningfull level of tax anno
                if "_" in LML: LML=LML.split("_")[0]
                if "+" in LML: LML=LML.split("+")[0]

                if LML in ncbiNames:
                        LML_taxid=ncbi.get_name_translator([LML]).values()[0][0]
                        taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level)
                        #if taxoGrp=="NA":
			#	taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level1)
                        if taxoGrp=="NA":
                                taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level2)
                        if taxoGrp=="NA":
                                taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level3)
                        listLevel.append(taxoGrp)
                elif LML not in seenLML:
                        fwarning.write("%s" % l)

                seenLML.append(LML)
	f.close()
	return listLevel	

def getNbUniqPolBperTaxa(FILE): # IT IS NEARLY THE SAME FUN AS ABOVE BUT HERE I COUNT HOW MANY polb are connected to given taxa 
	f=open(FILE, "r")
	level="order"; level2="class"; level3="phylum"
	dico={}
	for  l in f:
              	#OM-RGC.v1.000065429    4202    0.5330  0.0099  Eukaryota|...|Cryptosporidium+fragile   core Apicomplexa
       	        t=l.strip().split("\t")
                polb=t[0]; lineage=t[4]
               	LML=lineage.split("|")[-1] #LML mean lowest meaningfull level of tax anno
       	        if "_" in LML: LML=LML.split("_")[0]
                if "+" in LML: LML=LML.split("+")[0]

                if LML in ncbiNames:
                        LML_taxid=ncbi.get_name_translator([LML]).values()[0][0]
                        taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level)
                        #if taxoGrp=="NA":
                        #       taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level1)
                        if taxoGrp=="NA":
                                taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level2)
                        if taxoGrp=="NA":
                                taxoGrp=getNameGivenTAXIDandRank(LML_taxid, level3)
                        
			if taxoGrp in dico.keys():
				dico[taxoGrp]+=1
			else:
				dico[taxoGrp]=1
			#print dico
	return dico

# - - load all taxa names in ncbi tree - - 
ncbiNames=getAllNames()

# - -  TRAVERSE THE TREE - - 
fwarning=open("taxaNotInNCBI.txt", "w") #function getCoPerEukGrp will write in this file
seenLML=[] # function getCoPerEukGrp will fill in this list
nodeNb=0
for node in t.iter_descendants("preorder"):
	nodeNb+=1
	# - - Do some analysis on node
	fout=open("node_tmp.leaves.txt", "w")
        tmp=node.get_leaves()
	leaves=0
	for feuille in tmp:
		fout.write("%s" % feuille)
		leaves+=1

	tmp1=len(node)
	if tmp1!=leaves:
		exit("\terror: difference in the way to compte leaves under a node!\n")
	fout.close()
	
	# info on leaves distance
	farthest, dist_farthest = node.get_farthest_leaf()
	averageDist=avLeavesToNodeDist(node)
	
	os.system("sed 's/--//g' node_tmp.leaves.txt | sed '/^$/d' > node_%d.leaves.txt" % nodeNb)
	os.system("grep -w -f node_%d.leaves.txt %s > node_%d.connections.txt" % (nodeNb, USED_FILE, nodeNb))
        os.system("grep -w -f node_%d.leaves.txt %s -v > others.connections.txt" % (nodeNb, USED_FILE))	
	
	IC=getCoPerEukGrp("node_%d.connections.txt" % nodeNb) #IC=in clade
	#print IC
	NIC=getCoPerEukGrp("others.connections.txt") #NIC=not in clade
	dico_nbUniqPolBperTaxa=getNbUniqPolBperTaxa("node_%d.connections.txt" % nodeNb)

	# - - - find out which Euk group receives the greatest number of association
	setLevel=set(IC)
	#print setLevel
	modeName=""
	mode=0
	fout0=open("node_%d.counts.txt" % nodeNb, "w")
	for lev in setLevel:
		# number of uniq viral sequence connected to this taxa in this clade
		nbUniqPolBperTaxa=dico_nbUniqPolBperTaxa[lev]
		# Fisher extact test:
		levIC=IC.count(lev)
		levNIC=NIC.count(lev)
		otherIC=len(IC)-levIC
		otherNIC=len(NIC)-levNIC
		oddsratio, pvalue = stats.fisher_exact([[levIC, otherIC], [levNIC, otherNIC]])
		#get rank
		if lev=="NA":
			rank="no rank"
		else:
			taxid=ncbi.get_name_translator([lev]).values()[0]
			rank=ncbi.get_rank(taxid).values()[0]
		fout0.write("%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%f\t%f\n" % (nodeNb, leaves, lev, rank, levIC, otherIC, levNIC, otherNIC, pvalue, nbUniqPolBperTaxa, dist_farthest, averageDist))	
		
		if levIC>mode:
			modeName=lev	
			mode=levIC
	#print "mode name is %s" % modeName
	fout0.close()
fwarning.close()


# CLEANING
os.system("echo 'node\tnbLeaves\ttaxa\trank\ttaxaIC\totherIC\ttaxaNIC\totherNIC\tPval\tnbUniqPolBperTaxa\tdist_farthest\tdist_average' > tmp")
os.system("cat tmp node_*counts.txt > allNodesCounts.txt")
os.system("mkdir NODES")
os.system("mv node_*.txt NODES/")
os.system("mkdir myTIMrun")
os.system("mv allNodesCounts.txt taxaNotInNCBI.txt NODES/  myTIMrun/.")
os.system("rm tmp tmpConnections.txt others.connections.txt")

print "finito"




