# Synechococcus-metabolic-modelling  
This archive contains code files and datasets utilized in the publication: "A hybrid multi-omic modeling and machine learning pipeline to identify adaptation mechanisms of cyanobacteria".

The model.xml model file for the Synechococcus PCC 7002 was previously published in https://www.sciencedirect.com/science/article/abs/pii/S0960852416302747
and converted into .mat format for modelling - this is saved as SynechococcusPCC7002.mat

The folder transcriptomic_data contains all RNA sequencing data downloaded from Cyanomics:
 - The initial .xls datasets contain RPKM values for each gene/locus (Dataset1split.xls and Dataset2split.xls) were imported into     Matlab as matrices (Dataset1RPKM and Dataset2RPKM) that were converted into fold change values centred around 1 by dividing each condition by the mean standard control (Dataset1newFC and Dataset2newFC). Combining both of the matrices gives the single matrix DatasetsnewFC.
 - Filenames ending in "...newFC" are matlab structures to be directly loaded for FBA



