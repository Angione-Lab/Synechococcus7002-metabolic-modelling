%RUN_all script runs norm-2 regularized bi-level optimization, integrating RNA-seq transcriptomic data from studies conducted by Ludwig and Bryant deposited in Cyanomics.
%Add code repositories to MATLAB path
addpath(genpath('C:/Users/****/Synechococcus-metabolic-modelling-master'));
addpath(genpath('C:/Users/****/Synechococcus-metabolic-modelling-master/transcriptomic_data'));
addpath(genpath('C:/Users/****/Synechococcus-metabolic-modelling-master/lasso'));
savepath();

load('reaction_expression.mat');
load('SynechococcusPCC7002.mat'); %fbamodel
load('pos_genes_in_react_expr.mat');
load('ixs_geni_sorted_by_length.mat');
load('Syn7002_IDs.mat'); % list of gene IDs extracted from transcriptomic reads file

% load conditions i.e. RPKM fold change values normalised 0-1
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

% Check current objective functions
ix_f = find(fbamodel.f==1); %check current primary objective
ix_g = find(fbamodel.g==1); %check current secondary objective

% Set new primary objective f
ix_new_f = 735; % set new main objective = standard biomass reaction (735) 

% In each simulation, a secondary objective must be specified *manually* by uncommenting the selected reaction and commenting the other two reactions when not in use.
% Therefore, in line 665, the name of the outputs should be be renamed to specify which secondary objective was in use i.e. all_atp_flux, all_p1_flux or all_p2_flux.
% Set new secondary objective g
ix_new_g = find(ismember(fbamodel.rxnNames,'ATP maintenance requirment')==1); % select this objective to produce all_atp_flux
% ix_new_g = find(ismember(fbamodel.rxnNames,'Photosystem I Reaction (cytochrome c6)')==1); % select this objective to produce all_p1_flux
% ix_new_g = find(ismember(fbamodel.rxnNames,'photosystem II reaction')==1); % select this objective to produce all_p2_flux

% Select new objective functions for simulation

fbamodel.f(ix_f) = 0;
fbamodel.f(ix_new_f) = 1;
fbamodel.g(ix_g) = 0;
fbamodel.g(ix_new_g) = 1;

%% Model constraints
%% Boundary constraints to simulate growth medium and record experimentally feasible growth rates
%% Load new flux bounds 
%% Load list of variables including rxn names, indices and new values for lower and upper bounds in the model
load('bounds.mat')

%% Solver
%% Set Gurobi as the solver for linear and quadratic problems
changeCobraSolver('gurobi','LP');
changeCobraSolver('gurobi','QP');
changeCobraSolverParams('QP', 'method', 1); %avoid solver feasibility error

%% Set new bounds for standard control condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,1);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,1);
% fbamodel_control=fbamodel;

%% Flux distribution control
%
x = ones(numel(genes),1);  %we start from the all-one configuration
[v1_control, f_out_control] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length);   

%% Set new bounds for dark oxic condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,2);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,2);
% fbamodel_do=fbamodel;
%% Flux distribution in dark oxic condition
% %
expr_profile = DarkoxicnewFC(:,1); %choose growth condition
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
[v1_do, f_out_do] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length);   
%% Set bounds for dark anoxic condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,3);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,3);
% fbamodel_da=fbamodel;
% 
%% Flux distribution in dark anoxic condition
% %
expr_profile = DarkanoxicnewFC(:,1); %choose growth condition
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
[v1_da, f_out_da] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length);   
%% Set bounds for high light condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,4);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,4);
% fbamodel_hl=fbamodel;
% 
%% Flux distribution in high light condition
% %%
expr_profile = HighlightnewFC(:,1); %choose growth condition
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
[v1_hl, f_out_hl] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for OD 0.4 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,5);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,5);
% fbamodel_od04=fbamodel;
% 
%% Flux distribution in OD 0.4 condition
%%
expr_profile = OD04newFC(:,1); %choose growth condition
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
[v1_od04, f_out_od04] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for OD 1.0 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,6);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,6);
% fbamodel_od10=fbamodel;

%% Flux distribution in OD 1.0 condition
% %%
expr_profile = OD10newFC(:,1); %choose growth condition
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
[v1_od10, f_out_od10] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for OD 3.0 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,7);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,7);
% fbamodel_od30=fbamodel;

%% Flux distribution in OD 3.0 condition
%%
expr_profile = OD30newFC(:,1); %choose growth condition
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
[v1_od30, f_out_od30] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for OD 5.0 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,8);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,8);
% fbamodel_od50=fbamodel;
% % 
%% Flux distribution in OD 5.0 condition
%%
expr_profile = OD50newFC(:,1); %choose growth condition
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
[v1_od50, f_out_od50] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for low O2 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,9);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,9);
% fbamodel_lo2=fbamodel;
% 
%% Flux distribution in low O2 condition
%%
expr_profile = lowO2newFC(:,1); %choose growth condition
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
[v1_lo2, f_out_lo2] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for low CO2 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,10);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,10);
% fbamodel_lco2=fbamodel;
% 
%% Flux distribution in low CO2 condition
%%
expr_profile = lowCO2newFC(:,1); %choose growth condition
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
[v1_lco2, f_out_lco2] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for N-limited condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,11);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,11);
% fbamodel_nlim=fbamodel;
% 
%% Flux distribution in N-limited condition
%
expr_profile = NlimnewFC(:,1); %choose growth condition
pos_genes_in_dataset = zeros(numel(genes),1);

expression = '[.]\d';
replace = '';
genes_truncated = regexprep(genes,expression,replace);  %we remove the last two characters (e.g. '.1') because in the model we have the transcripts indicated with '.1', while these are not in the dataset
% genes_in_dataset = strtrim(cellstr(num2str(Boston_CD4_Entrez_IDs))); %converts numerical array into cell array of strings

for i=1:numel(genes)
    position = find(strcmp(genes_truncated{i},genes_in_dataset)); 
    if ~isempty(position)
        pos_genes_in_dataset(i) = position;    
        x(i) = expr_profile(pos_genes_in_dataset(i));
    end
end

V = numel(genes);
[v1_nlim, f_out_nlim] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for S-limited condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,12);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,12);
% fbamodel_slim=fbamodel;

%% Flux distribution in S-limited condition
% %%
expr_profile = SlimnewFC(:,1); %choose growth condition
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
[v1_slim, f_out_slim] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for P-limited condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,13);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,13);
% fbamodel_plim=fbamodel;
% 
%% Flux distribution in P-limited condition
%%
expr_profile = PlimnewFC(:,1); %choose growth condition
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
[v1_plim, f_out_plim] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for Fe-limited condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,14);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,14);
% fbamodel_felim=fbamodel;
% 
%% Flux distribution in Fe-limited condition
% %%
expr_profile = FelimnewFC(:,1); %choose growth condition
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
[v1_felim, f_out_felim] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for NO3- condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,15);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,15);
% fbamodel_no3=fbamodel;
% 
%% Flux distribution in NO3- condition
% % %%
expr_profile = NitratenewFC(:,1); %choose growth condition
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
[v1_no3, f_out_no3] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for NH3 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,16);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,16);
% fbamodel_nh3=fbamodel;

%% Flux distribution in NH3 condition
%
expr_profile = AmmonianewFC(:,1); %choose growth condition
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
[v1_nh3, f_out_nh3] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for CO(NH2)2 condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,17);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,17);
% fbamodel_urea=fbamodel;
% % 
%% Flux distribution in CO(NH2)2 condition
% %%
expr_profile = UreanewFC(:,1); %choose growth condition
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
[v1_urea, f_out_urea] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for heat shock condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,18);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,18);
% fbamodel_heat=fbamodel;
% 
%% Flux distribution in heat shock condition
%%
expr_profile = HeatshocknewFC(:,1); %choose growth condition
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
[v1_heat, f_out_heat] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for 22C condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,19);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,19);
% fbamodel_22c=fbamodel;
% 
%% Flux distribution in 22C condition
%%
expr_profile = T22newFC(:,1); %choose growth condition
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
[v1_22c, f_out_22c] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for 30C condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,20);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,20);
% fbamodel_30c=fbamodel;
% 
%% Flux distribution in 30C condition
%%
expr_profile = T30newFC(:,1); %choose growth condition
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
[v1_30c, f_out_30c] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for oxidative stress condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,21);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,21);
% fbamodel_oxs=fbamodel;
% 
%% Flux distribution in oxidative stress condition
%
expr_profile = OxstressnewFC(:,1); %choose growth condition
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
[v1_oxs, f_out_oxs] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for mixotrophic condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,22);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,22);
% fbamodel_mix=fbamodel;
% 
%% Flux distribution in mixotrophic condition
% % 
expr_profile = MixotrophicnewFC(:,1); %choose growth condition
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
[v1_mix, f_out_mix] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for low salinity condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,23);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,23);
% fbamodel_ls=fbamodel;
% % 
%% Flux distribution in low salinity condition
% % % %
expr_profile = LowsaltnewFC(:,1); %choose growth condition
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
[v1_ls, f_out_ls] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
%% Set bounds for high salinity condition
fbamodel.lb(new_lb_ixs)=new_lb_val(1:15,24);
fbamodel.ub(new_ub_ixs)=new_ub_val(1:2,24);
% fbamodel_hs=fbamodel;
% 
%% Flux distribution in high salinity condition
% %%
expr_profile = HighsaltnewFC(:,1); %choose growth condition
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
[v1_hs, f_out_hs] = evaluate_objective_minNorm(x,M,V,fbamodel,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length); 
% % % 
% % %% Concatenate flux vectors for all growth conditions
all_atp_flux = [v1_control,v1_do,v1_da,v1_hl,v1_od04,v1_od10,v1_od30,v1_od50,v1_lo2,v1_lco2,v1_nlim,v1_slim,v1_plim,v1_felim,v1_no3,v1_nh3,v1_urea,v1_heat,v1_22c,v1_30c,v1_oxs,v1_mix,v1_ls,v1_hs];
% Change all_atp_flux to all_p1_flux or all_p2_flux when changing ix_new_g (secondary flux objective) in lines 55-57.
