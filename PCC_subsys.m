%% Calculate average Pearson correlation coefficient for metabolic pathways
% Load table of PCC values for flux data (ATP/P1/P2) 
% corr_ATP = array2table(corr_ATP,'VariableNames',{'PCC'});

% Select column containing PCC values and convert into array
% ATP_PCC = table2array(corr_ATP_table(:,3));

% Convert NaNs into 0
ATP_PCC(isnan(ATP_PCC)) = 0;
P1_PCC(isnan(P1_PCC)) = 0;
P2_PCC(isnan(P2_PCC)) = 0;

% Convert all coefficients into absolute values
ATP_PCC_abs = abs(ATP_PCC);
P1_PCC_abs = abs(P1_PCC);
P2_PCC_abs = abs(P2_PCC);

% Load reaction indices for subsystems 
% load('ixs_subsystems.mat');
% Load number of reactions in subsystems
% load('cardinality_subsystems.mat');
cardinality_old_subsystems = zeros(numel(ixs_subsystems_uniq_old),1);
cardinality_subsystems = zeros(numel(ixs_subsystems),1);

% Create vectors to store sum and average PCC values
%ATP_PCC_sum = zeros(numel(ixs_subsystems),1);
ATP_PCC_mean = zeros(numel(ixs_subsystems),1);
P1_PCC_mean = zeros(numel(ixs_subsystems),1);
P2_PCC_mean = zeros(numel(ixs_subsystems),1);

% Calculate mean PCC for subsystems by calculating PCC sums then dividing by no. of reactions (cardinality)
for c = 1:numel(ixs_subsystems)
%ATP_PCC_sum(c) = sum(ATP_PCC_abs(ixs_subsystems{c},1));
cardinality_subsystems(c) = numel(ixs_subsystems{c},1);
ATP_PCC_mean(c) = mean(ATP_PCC_abs(ixs_subsystems{c},1));
P1_PCC_mean(c) = mean(P1_PCC_abs(ixs_subsystems{c},1));
P2_PCC_mean(c) = mean(P2_PCC_abs(ixs_subsystems{c},1));
%ATP_PCC_mean_old(c) = ATP_PCC_sum(c)./cardinality_old_subsystems(c);
end

%% Combine all PCC values
PCC_all = horzcat(ATP_PCC_mean,P1_PCC_mean,P2_PCC_mean);

% Plot mean PCC values against subsystems
X = categorical(list_subsystems);
bar(X,PCC_all);
xlabel('Subsystems'); 
ylabel('Mean PCC');
ax.FontSize = 8;
hold on
set(gca, 'XTickLabelRotation',45);
legend({'Biomass - ATP maintenance','Biomass - Photosystem I', 'Biomass - Photosystem II'},'Location','northwest');
