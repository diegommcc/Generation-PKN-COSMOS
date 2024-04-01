# Generation of prior knowledge network for COSMOS

This repository contains a set of functions to generate organism-specific prior knowledge networks (PKNs) for the [cosmosR](https://www.bioconductor.org/packages/release/bioc/html/cosmosR.html) R package. It is structured as follows: 

* `data`: it contains the starting data needed to generate PKNs for every organism. These data come from [this repository](https://github.com/SysBioChalmers/), a platform of genome-scale metabolic models (GEM) for different organisms. For more details about how these GEMs were generated, see Wang et al., 2021. In this repository, only GEMs for human (_Homo sapiens_), mouse (_Mus musculus_) and rat (_Ratus norvegicus_) are available. In addition, information from KEGG to determine cofactors is also already downloaded and stored as an RDS file. Importantly, COSMOS also uses information of chemical-protein interactions obtained from the STITCH database (Kuhn et al., 2008). However, these files are not available in this repository, and must be downloaded from [here](http://stitch.embl.de/). 
* `notebooks`: Rmd with a brief analysis of COSMOS results comparing the old (available at the [cosmosR](https://www.bioconductor.org/packages/release/bioc/html/cosmosR.html) R package) and new PKN for human data. 
* `output`: final PKNs as RDS files and COSMOS results. 
* `reports`: HTML and plots (png) generated after rendering Rmds in `notebooks`. 
* `src`: 
    * `final_functions_PKN_COSMOS.R`: set of functions in charge of generating the PKNs. 
    * `generation_networks_COSMOS.R`: script which loads `final_functions_PKN_COSMOS.R` functions and write the final PKNs as RDS files in `output`. 
    * `run_COSMOS_*`: scripts used to run COSMOS with data from <https://github.com/saezlab/Sciacovelli_Dugourd_2021_paper> and the old and new PKNs. The suffix determines the PKN and `maximum_network_depth` parameter used. 

These functions will be available on the [OmnipathR](https://bioconductor.org/packages/release/bioc/html/OmnipathR.html) R package.  

# References

* Kuhn M, von Mering C, Campillos M, Jensen LJ, Bork P. STITCH: interaction networks of chemicals and proteins. Nucleic Acids Research. 2008 Jan 1;36(suppl_1):D684–8.
* Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley PE, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. Proceedings of the National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.
