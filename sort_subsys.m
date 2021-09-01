% Load model and flux data
 %load('SynechococcusPCC7002.mat');
 %load('protocol_10.03.21.mat');
% 
% Use absolute flux values
 %all_atp_flux = abs(all_atp_flux);
 %all_p1_flux = abs(all_p1_flux);
 %all_p2_flux = abs(all_p2_flux);

% all_atp_flux(all_atp_flux<=0.0001)=0;
% all_p1_flux(all_p1_flux<=0.0001)=0;
% all_p2_flux(all_p2_flux<=0.0001)=0;

% Make a query list of unique subsystems
% unique_subsystems = {'Amino Acid Biosynthesis','Amino Acid Metabolism','Aminoacyl-tRNA biosynthesis','Biomass Synthesis','Biosynthesis of other secondary metabolites','Biotin Biosynthesis','Carbohydrate metabolism','Carotenoid Biosynthesis','Cell wall','Coenzymes & prosthetic groups','Cyanophycin Metabolism','Energy metabolism','Exchange','Fatty Acid Synthesis','Folate Metabolism','Glycan biosynthesis & metabolism','Glycerophospholipid metabolism','Hydrogen Metabolism','Lipid metabolism','Membrane bioenergetics','Metabolism of cofactors & vitamins','Metabolism of other amino acids','Metabolism of terpenoids & polyketides','Modelling','Nucleotide Metabolism','Nucleotides & nucleic acids','Peptidoglycan biosynthesis','Phenylmercury acetate degradation','Proline Biosynthesis','Purine metabolism','Pyridine Metabolism','Quinone Biosynthesis','Riboflavin Metabolism','Salt tolerance','Thiamine metabolism','Transport','Transport, Extracellular','Unassigned','Vitamin E biosynthesis'}';
% sort fluxes according to subsystems
% list_subsystems = unique([fbamodel.subSystems{:}])'; % list all subsystems
% list_subsystems = unique([new_subsystems{:}])';
%ixs_subsystems = cell(numel(list_subsystems),1); % create cell variables for pathway indices
%ixs_subsystems = cell(length(list_subsystems),1); % create cell variables
cardinality_subsystems = zeros(numel(list_subsystems),1);
sum_contrib_subsystems = zeros(numel(list_subsystems),1);
avg_contrib_subsystems = zeros(numel(list_subsystems),1);

% for i=1:numel(list_subsystems)
%     ixs_subsystems{i} = find(strcmp(list_subsystems{i},fbamodel.subSystems));
%     sum_flux_subsystems(i) = sum(abs(all_atp_flux(ixs_subsystems{i}))); % select each group of fluxes individually i.e. all_atp_flux, all_p1_flux or all_p2_flux
%     cardinality_subsystems(i) = numel(ixs_subsystems{i});
% end

 for i=1:numel(list_subsystems)
%for i=1:length(list_subsystems)
    %ixs_subsystems{i} = find(strcmpi(list_subsystems{i},fbamodel.subSystems));
    sum_contrib_subsystems(i) = sum(Dim_1_and_2(ixs_subsystems{i},1)); % column 1 for Dim 1, column 2 for Dim 2
    cardinality_subsystems(i) = numel(ixs_subsystems{i});
    avg_contrib_subsystems(i) = mean(sum_contrib_subsystems(i)./cardinality_subsystems(i));
end