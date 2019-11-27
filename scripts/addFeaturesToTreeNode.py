#!/Users/romain/anaconda_ete/bin/python
from __future__ import division
import os
import sys
from ete3 import Tree
import subprocess
import scipy.stats as stats

treeName=sys.argv[1]
t = Tree("%s" % treeName)

# - -  TRAVERSE THE TREE - - 
nodeNb=0
for node in t.iter_descendants("preorder"):
        #print node.support
        nodeNb+=1
        #print nodeNb
        if not node.is_leaf():
                node.name = "INT_%d" % nodeNb #If I do this only I loose node support..
                #node.add_feature("label", "INT_%d" % nodeNb)
                #node.add_feature("confidence", "%f" % node.support)

t.write(format=1, outfile="treeWithNodeID_forItolPlot.nwk")
#t.write(format=9, outfile="treeWithNodeFeatures.nwk",features=["label", "confidence"]) 


print "finito"
