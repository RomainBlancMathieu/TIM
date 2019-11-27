# TIM - Taxon Interaction Mapper
The Taxon Interaction Mapper (TIM) maps taxon interactions onto a phylogenetic tree of a target group of organism using the organism’s associations predicted from a species co-occurence-based network (such as one generated by FlashWeave or SparCC)

It assumes that evolutionarily related organisms or viruses (refer to as ```query```) infect evolutionary related organisms (```subject```)

## Prerequisites

Listed softwares are required by TIM.
Please install these softwares and add their directories to the ```PATH``` environment variable or specify them in ```Pipeline.sh```.
* [ETE3 v3](http://etetoolkit.org/download/) 
* [Python 2.7.11 or greater](https://www.python.org/downloads/release/python-2711/)
* [R](https://www.r-project.org/)

### Input files

The TIM pipeline needs two inputs: <br />
* A phylogenetic tree for the ```query``` <br /> 
* A table file containing connection between ```query``` and ```subject``` with follwing tab separated columns:<br />
&nbsp;&nbsp;Query ID (same as leaves's names in the tree) <br />
&nbsp;&nbsp;Subject ID <br />
&nbsp;&nbsp;Direction of connections (typically the weight, positive or neagtive, in a co-occurence-based network) <br />
&nbsp;&nbsp;P-value (facultative can be NA if you pre-filtred the connections) <br />
&nbsp;&nbsp;Genealogy of the subject ID (ideally in NCBI taxonomy terms, if a taxomic name is not found the tool will report it) <br />


## Usage
TIM runs in two main steps scripts: <br />
```main.py tree.nwk connections.tsv [ALL, POS, NEG]``` will ... <br />
```downstream.py``` ... <br />

## Output files
taxaNotInNCBI.txt : list taxa's name in your connection file that were not found in the NCBI taxonomy <br />

## List of reference and script files
### Scripts
| Filename | Description |
| ---- | :--- |
|```main.py```|report number of connection for a taxonomic group wiht leaves in the tree|
|```taxo_functions.py```|function for taxonomy|
|```downstream.py```|filters out the output of ```main.py``` and prepare ITOL formated files for visualization|
|```filt_form.R```|used by ```downstream.py``` for filtering and files formating|
|```addFeaturesToTreeNode.py```|add node ID to the tree used by assignHost. Can be use to visualized the results on ITOL|

## Authors

* **Romain Blanc-Mathieu**  - romain.blancmathieu@gmail.com

## References
Blanc-Mathieu R, Kaneko H, Endo H, Chaffron S, Hernández-Velázquez R, Nguyen CH, Mamitsuka H, Henry N, Vargas C de, Sullivan MB, et al. 2019. Viruses of the eukaryotic plankton are predicted to increase carbon export efficiency in the global sunlit ocean. bioRxiv:710228.
https://www.biorxiv.org/content/10.1101/710228v1.full

## Acknowledgments

* ...

