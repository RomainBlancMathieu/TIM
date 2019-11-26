# TIM
Taxon Interaction Mapper


The Taxon Interaction Mapper (TIM) is a set of scripts to...

## Prerequisites

Listed softwares are required by TIM.
Please install these softwares and add their directories to the ```PATH``` environment variable or specify them in ```Pipeline.sh```.
* [ETE3 v3](http://etetoolkit.org/download/) 
* [Python 2.7.11 or greater](https://www.python.org/downloads/release/python-2711/)
* [R](https://www.r-project.org/)

## Usage
TIM runs in two main steps corresponding to scripts ```run.sh``` and ```run.py```.
* ```-maindir MAIN_DIR``` : Path of main directory, which should at least contain ```Scripts``` directory and ```References``` directory<br>


## List of reference and script files
### Scripts
| Filename | Description |
| ---- | :--- |
|```run.sh```|main script|
|```assignHost.py```|report number of connection for a taxonomic group wiht leaves of a given node in the tree  |
|```rrrr.R```|test for connetion enrichement for a given taxonomic group at a given node in the tree|
|```addFeaturesToTreeNode.py```|add node ID to the tree used by assignHost. Can be use to visualized the results on ITOL|

## Authors

* **Romain Blanc-Mathieu**  - romain.blancmathieu@gmail.com

## References
Blanc-Mathieu R, Kaneko H, Endo H, Chaffron S, Hernández-Velázquez R, Nguyen CH, Mamitsuka H, Henry N, Vargas C de, Sullivan MB, et al. 2019. Viruses of the eukaryotic plankton are predicted to increase carbon export efficiency in the global sunlit ocean. bioRxiv:710228.
https://www.biorxiv.org/content/10.1101/710228v1.full

## Acknowledgments

* ...

