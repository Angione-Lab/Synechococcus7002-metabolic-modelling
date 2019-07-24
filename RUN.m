%addpath(genpath('C:\Program Files\MATLAB\R2015b\toolbox\cobra'));
%addpath(genpath('C:\Users\Supreeta\OneDrive\Synechococcus code'));
%initCobraToolbox

%load('geni_names.mat');
load('reaction_expression.mat');
load('SynechococcusPCC7002.mat'); %fbamodel
load('pos_genes_in_react_expr.mat');
load('ixs_geni_sorted_by_length.mat');
load('Syn7002_IDs.mat'); % list of gene IDs extracted from transcriptomic reads file

%load conditions i.e. RPKM fold change values normalised 0-1
load('AmmonianewFC.mat')
load('DarkanoxicnewFC.mat')
load('DarkoxicnewFC.mat')
load('FelimnewFC.mat')
load('HeatshocknewFC.mat')
load('HighlightnewFC.mat')
load('HighsaltnewFC.mat')
load('lowCO2newFC.mat')
load('lowO2newFC.mat')
load('LowsaltnewFC.mat')
load('MixotrophicnewFC.mat')
load('NitratenewFC.mat')
load('NlimnewFC.mat')
load('OD04newFC.mat')
load('OD10newFC.mat')
load('OD30newFC.mat')
load('OD50newFC.mat')
load('OxstressnewFC.mat')
load('PlimnewFC.mat')
load('SlimnewFC.mat')
load('T22newFC.mat')
load('T30newFC.mat')
load('UreanewFC.mat')

genes = fbamodel.genes; 
genes_in_dataset = Syn7002_IDs;

M = 2; %number of objectives
V = numel(genes); %number of variables

ix_f = find(fbamodel.f==1); %check current primary objective
ix_g = find(fbamodel.g==1); %check current secondary objective
ix_new_f = 735; % set new main objective = standard biomass or 73=c-lim, 74=nlim, 75=llim?

%Set new secondary objective g
ix_new_g = find(ismember(fbamodel.rxnNames,'ATP maintenance requirment')==1);
%ix_new_g = find(ismember(fbamodel.rxnNames,'Photosystem I Reaction (cytochrome c6)')==1);
%ix_new_g = find(ismember(fbamodel.rxnNames,'photosystem II reaction')==1);
%ix_new_g = find(ismember(fbamodel.rxns,'DGDGST_C181C161_SYN')==1); %digalactosyldiglycerol (glycolipid) synthase rxn no. 662
%ix_new_g = find(ismember(fbamodel.rxns,'pgly_syn')==1); % fatty acid synthesis rxn no. 640
% ix_new_g = find(ismember(fbamodel.rxns,'SQDG_syn')==1); % fatty acid synthesis rxn no. 646
% ix_new_g = find(ismember(fbamodel.rxns,'MGDG_syn')==1); % fatty acid synthesis rxn no. 657
% ix_new_g = find(ismember(fbamodel.rxns,'DGDG_syn')==1); % fatty acid synthesis rxn no. 665
% ix_new_g = find(ismember(fbamodel.rxns,'LIPSYN_SYN')==1); % fatty acid synthesis rxn no. 666


fbamodel.f(ix_f) = 0;
fbamodel.f(ix_new_f) = 1;
fbamodel.g(ix_g) = 0;
fbamodel.g(ix_new_g) = 1;

%% Flux distribution control
%%
x = ones(numel(genes),1);  %we start from the all-one configuration
[v1_control, f_out_control] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length);   

%% Flux distribution in condition
%%
expr_profile = DarkoxicnewFC(:,1); %choose condition
pos_genes_in_dataset = zeros(numel(genes),1);

expression = '[.]\d';
replace = '';
genes_truncated = regexprep(genes,expression,replace);  %we remove the last two characters (e.g. '.1') because in the model we have the transcripts indicated with '.1', while these are not in the dataset
%genes_in_dataset = strtrim(cellstr(num2str(Boston_CD4_Entrez_IDs))); %converts numerical array into cell array of strings

for i=1:numel(genes)
    position = find(strcmp(genes_truncated{i},genes_in_dataset)); 
    if ~isempty(position)
        pos_genes_in_dataset(i) = position;    
        x(i) = expr_profile(pos_genes_in_dataset(i));
    end
end

V = numel(genes);
[v1, f_out] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length);   
