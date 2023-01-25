Data and code used in Diniz and Aguiar (2022). Scientific Reports (under review). 

"Beyond a 'perfect fit': the interplay between morphology and spatiotemporal overlap as determinants of microstructure in a bat-flower network"

Version: 25.01.2023


Authors: Ugo M. Diniz (Technische Universität München), Ludmilla M. S. Aguiar (Universidade de Brasília).

E-mail for contact: ugo.diniz@tum.de


DOI: 



-> Folders
	- Main folder (Analysis): contains the main script, network and .txt and .xls files for most analysis, except the files used for the likelihood analysis and the felis used in Gephi, both of which have their own folders.
		
		1. morph_bats.xlsx > morphological data from bats separated by modules to build desity plots, and for the ANOVA analysis.

		2. morph_plants.xlsx > morphological data from plants separated by modules to build desity plots, and for the ANOVA analysis.

		3. network.txt > main network matrix containing the interactions between bats and plants. 

		4. samp.bats.xlsx > bat abundance data per sampling event for the construction of rarefaction curves.

		5. sampplantschiro.xlsx > plant (chiropterophilous only) abundance data per sampling event for the construction of rarefaction curves.

		6. sampplantsfull.xlsx > plant (all species) abundance data per sampling event for the construction of rarefaction curves.

		7. script.R > main R code file to reproduce step-by-step the analyses in the paper.

		8. spatial.xlsx > data on bat relative abundance per habitat type.
	 
		9. temporal.xlsx > data on bat relative abundance per month of the year.


	-likelihood models: contains three separate folders, which need to be accessed at different points of the code.
		
		1. full > contains all matrices (observed and probabilities) used in the likelihood analysis in .txt and .xlsx formats, with all species from the network.

		2. nectar > contains all matrices (observed and probabilities) used in the likelihood analysis in .txt and .xlsx formats, with nectarivores and their interactions only.

		3. other > contains all matrices (observed and probabilities) used in the likelihood analysis in .txt and .xlsx formats, with bats from other guilds and their interactions only.

 
	-Gephi files: contains the node and edge data to make Gephi graph.

		1. gephi_edges.xlsx

		2. gephi_nodes.xlsx		 		


R version used: 4.1.0 (2021-05-18) "Camp Potanezen".


*OBS* Original sources: the code uses functions and code sections created and kindly made public by other authors, which are acknowledged within the script itself. These are:  


1. The functions PosteriorProb and RestNullModel and the corresponding script sections to calculate WNODA and its variants. Source: Pinheiro, R. B., Felix, G. M., & Lewinsohn, T. M. (2022). Hierarchical compound topology
	uncovers complex structure of species interaction networks. Journal of Animal Ecology, 91(11), 2248-2260. Data and code: https://doi.org/10.5281/zenodo.7007952   

2. The function mlik. Source: P. Vázquez, Diego; Chacoff, Natacha P.; Cagnolo, Luciano (2016): Evaluating multiple determinants of the structure of plant-animal mutualistic networks. Wiley. Collection.
	Data and code: https://doi.org/10.6084/m9.figshare.c.3301145.v1 


## Reference

Diniz, U. M., Aguiar, L. M. S. (2022) "Beyond a 'perfect fit': the interplay between morphology and spatiotemporal overlap as determinants of microstructure in a bat-flower network". Scientific Reports (Under review)
