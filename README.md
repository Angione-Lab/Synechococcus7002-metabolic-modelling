# Synechococcus7002-metabolic-modelling  
This archive contains code files and data utilized in the publication: "A hybrid multi-omic modeling and machine learning pipeline to elucidate metabolic response of cyanobacteria to different growth conditions".

The `modelXML.xml` model file for the Synechococcus sp. PCC 7002 was previously published in https://www.sciencedirect.com/science/article/abs/pii/S0960852416302747
and converted into .mat format for modelling - this is saved as `SynechococcusPCC7002.mat`

The folder transcriptomic_data contains all RNA sequencing data downloaded from Cyanomics, a database based on the results of Synechococcus sp. PCC 7002 omics studies:
https://lag.ihb.ac.cn/cyanomics (link no longer active).
https://academic.oup.com/database/article/doi/10.1093/database/bau127/2433127.

Currently, the relevant files are available at the NCBI Sequence Read Archive:
- https://www.ncbi.nlm.nih.gov/sra?term=SRP007372 
- https://www.ncbi.nlm.nih.gov/sra?term=SRP013965 
- https://www.ncbi.nlm.nih.gov/sra/?term=SRP066851

 - The initial .xls datasets containing RPKM values for each gene/locus (`Dataset1split.xls` and `Dataset2split.xls`) were imported into     Matlab as matrices (`Dataset1RPKM.mat` and `Dataset2RPKM.mat`) that were converted into fold change values centred around 1 by dividing each condition by the mean of three standard controls (`Dataset1newFC.mat` and `Dataset2newFC.mat`). Combining both of these matrices gives the single matrix `DatasetsnewFC.mat`.
 - All other filenames ending in "...newFC.mat" are separate vectors for each growth condition converted into expression profiles that are called by `evaluate_objective_minNorm.m` when running `RUN_all.m`.
 
The simulation is initialized by running the `RUN_all.m` script, where regularized flux balance analysis is conducted for three different pairs of flux objectives: Biomass - ATP maintenance, Biomass - Photosystem I and Biomass - Photosystem II. Directions in the `RUN_all.m` script must be carefully followed to manually change the pair of objectives optimized for FBA in each of these three cases.
All simulations were run using the Cobra Toolbox in MATLAB R2019b with the Gurobi Optimizer 9.0 as a solver:
https://opencobra.github.io/cobratoolbox/stable/. Loading `bounds.mat` ensured adjustment of specific upper and lower bounds according to growth media and other requirements specific to each condition (mainly nutrient and photon uptake).

- All outputs were converted into absolute values and flux values < 10^-4 are set to zero to account for solver error.
- Prior to combining transcriptomic and fluxomic data in a common matrix, fold change was performed on fluxes by dividing all conditions by the standard control flux.

PCA was conducted in R using the script `PCA_script.R`, which uses the FactoMineR and factoextra packages for analysis.

k-means clustering was run using the `statistics_on_genes.m` script, which also calls `mdscale_robust.m`, a script that applies multidimensional scaling to avoid co-location of data points during clustering: https://github.com/jooh/matlab-plotting/blob/master/mdscale_robust.m.

The folder lasso contains the script `lasso_script.m` for running LASSO regularization in Matlab with subsets of transcript/flux data serving as predictor data (x) and growth rates measured across 12 growth conditions as responses (y).

Flux maps in Fig 5(a) & 5(b) of the main text were generated using Escher: https://escher.github.io/.
The model file `SynechococcusPCC7002.mat` was converted into `SynPCC7002_model.json` using cobrapy https://opencobra.github.io/cobrapy/.
The JSON map was saved as `SynPCC7002_map.json`.
Reaction data were loaded using the python script `flux_comparison_json.py`, which produces output files for various growth conditions.

The Pearson correlation coefficients are calculated using the script `corrcoef_tf_gr.m`, along with their respective p-values, and the lower and upper bounds of the 95% confidence interval.
`sort_subsys.m` is a script used to sort fluxes by their unique subsystem names for plotting the mean Pearson correlation coefficient (PCC)  according to model subsystems in Fig 5(c) of the main text.

The contributions by Supreeta Vijayakumar (S.V.) and Claudio Angione (C.A.) to the scripts listed above were as follows:
- `associate_genes_reactions.m`,`compute_reaction_expression.m`,`evaluate_objective_minNorm.m`,`flux_balance_minNorm.m`,`statistics_on_genes.m` - C.A.
- `corrcoef_tf_gr.m`,`flux_comparison_json.py`, `lasso_script`, `PCA_script.R`, `sort_subsys.m` - S.V.
- `RUN_all.m` - C.A and S.V
