# Synechococcus-metabolic-modelling  
This archive contains code files and datasets utilized in the publication: "A hybrid multi-omic modeling and machine learning pipeline to identify adaptation mechanisms of cyanobacteria".

The model.xml model file for the Synechococcus PCC 7002 was previously published in https://www.sciencedirect.com/science/article/abs/pii/S0960852416302747
and converted into .mat format for modelling - this is saved as SynechococcusPCC7002.mat

The folder transcriptomic_data contains all RNA sequencing data downloaded from Cyanomics, a database based on the results of Synechococcus sp. PCC 7002 omics studies:
https://lag.ihb.ac.cn/cyanomics (link no longer active) 
https://academic.oup.com/database/article/doi/10.1093/database/bau127/2433127

 - The initial .xls datasets containing RPKM values for each gene/locus (Dataset1split.xls and Dataset2split.xls) were imported into     Matlab as matrices (Dataset1RPKM and Dataset2RPKM) that were converted into fold change values centred around 1 by dividing each condition by the mean of three standard controls (Dataset1newFC and Dataset2newFC). Combining both of these matrices gives the single matrix DatasetsnewFC.
 - All other filenames ending in "...newFC" are separate vectors for each growth condition converted into expression profiles that are called by evaluate_objective_minNorm.m when running RUN_all.m.
 
The simulation begins by running the RUN_all script, where regularized flux balance analysis is conducted for three different pairs of flux objectives: Biomass - ATP maintenance, Biomass - Photosystem I and Biomass - Photosystem II. All simulations are run using the Cobra Toolbox in MATLAB R2019b with the Gurobi Optimizer 9.0 as a solver 
https://opencobra.github.io/cobratoolbox/stable/

- All outputs are converted into absolute values and flux values < 10^-4 are set to zero to account for solver error.
- Prior to combining transcriptomic and fluxomic data in a common matrix, fold change is performed on fluxes by dividing all conditions by the standard control flux.

PCA in conducted in R using the script PCA_script.R

k-means clustering is run using the statistics_on_genes.m script, which also calls mdscale_robust.m, a script that applies multidimensional scaling to avoid co-location of data points during clustering.

The folder lasso contains the script lasso_script.m for running LASSO regularization in Matlab with subsets of transcript/flux data serving as predictor data (x) and growth rates measured across 12 growth conditions as responses (y).

Flux maps in Fig 5 of the main text were generated using Escher http://escher.github.io/.
The model file SynechococcusPCC7002.mat is converted into SynPCC7002_model.json using cobrapy https://opencobra.github.io/cobrapy/.
The JSON model and map are saved as SynPCC7002_model.json and SynPCC7002_map.json.
Reaction data were loaded using the python script flux_comparison_json.py, which produces output files for various growth conditions.








