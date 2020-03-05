% Load model and flux data
% load('SynechococcusPCC7002.mat');
% load('matlab_10.02.20.mat');
% 
% Use absolute flux values
% all_atp_flux = abs(all_atp_flux);
% all_p1_flux = abs(all_p1_flux);
% all_p2_flux = abs(all_p2_flux);

all_atp_flux(all_atp_flux<=0.0001)=0;
all_p1_flux(all_p1_flux<=0.0001)=0;
all_p2_flux(all_p2_flux<=0.0001)=0;

% sort fluxes according to subsystems
list_subsystems = unique(fbamodel.subSystems); % list all subsystems
ixs_subsystems = cell(numel(list_subsystems),1); % create cell variables
cardinality_subsystems = zeros(numel(list_subsystems),1);
sum_flux_subsys= zeros(numel(list_subsystems),1);

for i=1:numel(list_subsystems)
    ixs_subsystems{i} = find(strcmp(list_subsystems{i},fbamodel.subSystems));
    sum_flux_subsys(i) = sum(abs(all_atp_flux(ixs_subsystems{i}))); % select each group of fluxes individually i.e. all_atp_flux, all_p1_flux or all_p2_flux
    cardinality_subsystems(i) = numel(ixs_subsystems{i});
end