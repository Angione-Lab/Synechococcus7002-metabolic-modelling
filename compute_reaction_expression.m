%load genes_reactions.mat
%load geo_s_genes_reactions.mat


%load('geo_s_react.mat');
load('SynechococcusPCC7002.mat');


%the grules are already written nicely and ordered by reaction, so there
%is no need to assign each geneset to a reaction, so there is no need to
%use the following code. We only need to convert "AND" into "and", "OR"
%into "or" because of how associate_genes_reactions.m is written

%% WARNING: we can remove brackets only in Recon because brackets are only used for gene names and not for associative rules between AND and OR. Strangely there is no associativity in the grRules!
% %However, in general, to avoid problems I am substituting all the '((' and
% %'))' (that are the useful brackets because the first '(' is part of the associative rules, while the second '(' is part of the gene name) with '[' and ']' to save them
% % then I am removing all the '(' and ')', 
% %and finally I am substituting back '[' and ']' with '(' and ')' respectvely. 
% % This way, all the useless brackets part of the gene names are
% %removed, while the useful brackets are kept
% grules = regexprep(grules,'((','['); 
% grules = regexprep(grules,'))',']'); 
% grules = regexprep(grules,'(',''); 
% grules = regexprep(grules,')',''); 
% grules = regexprep(grules,'[','('); 
% grules = regexprep(grules,']',')'); 

geni = fbamodel.genes;
grules = fbamodel.grRules;
for i = 1:numel(geni)
    gene = geni{i};
    grules = strrep(grules,['(' gene ')'],gene);
end

%%
grules = regexprep(grules,' AND ',' and '); 
grules = regexprep(grules,' OR ',' or '); 
%grules = upper(grules);  %upper('str') converts any lowercase characters in the string str to the corresponding uppercase characters and leaves all other characters unchanged. We need it because some genes in the grules are YorXXX while in the fbamodel are YORXXX

% %Boolean_grules = hacky_genes_reactions_conversion('geo_s_react');
%Boolean_grules = fbamodel.grRules;  %fbamodel.grRules contains the grules (also fbamodel.rules, but it contains only gene numbers and not their names)

% Boolean_grules = strcat('(',Boolean_grules,')');
% Boolean_grules = strcat('(',Boolean_grules,')');
% Boolean_grules = regexprep(Boolean_grules,', ',' and '); 
% Boolean_grules = regexprep(Boolean_grules,'_AND_',' and '); 
% %Boolean_grules = regexprep(Boolean_grules,'_OR_',' or '); 
% Boolean_grules = regexprep(Boolean_grules,'_OR_',') or (');  %brackets are necessary to force associative rules when these are not included in the original Boolean grules
% 
% Boolean_grules = regexprep(Boolean_grules,'  ',' '); 
% Boolean_grules = regexprep(Boolean_grules,'  ',' '); 
% Boolean_grules = regexprep(Boolean_grules,'  ',' '); 
% 
% 

% full_G = full(fbamodel.G); %geneset-to-reaction matrix
% fbamodel.pts; %grules
% 
% grules = cell(size(full_G,2),1);   %we need to map the grules to the reactions (all the grules are associated to a reaction, but non every reaction in the model has a corresponding geneset (e.g. the exchange or sink reactions do not have grules)
% grules(:) = {''};
% for j = 1 : size(full_G,2)
%     i = find(full_G(:,j) == 1);  %find the geneset j responsible for reaction i (which means G(j,i)=1)
%     if ~isempty(i)
%         grules(j) = Boolean_grules(i);
%     end
% end


reaction_expression = cell(length(grules),1);  %initializes an empty cell array, where each cell will be a string (it's the only way to create an array of strings in matlab)
reaction_expression(:) = {''};

%parfor_progress(length(grules)); %initialise
parfor i = 1:length(grules)
    %i
    aux = associate_genes_reactions(grules{i});
    reaction_expression{i} = aux;  
    %parfor_progress;
end

reaction_expression=strrep(reaction_expression,' ',''); %removes white spaces
reaction_expression




for i=1:size(geni)
    lung(i)=length(geni{i});
end
[lung_sorted, ixs_geni_sorted_by_length] = sort(lung,'descend');

reaction_expression_aux = reaction_expression;

for i = 1:numel(ixs_geni_sorted_by_length)
    j = ixs_geni_sorted_by_length(i);
    matches = strfind(reaction_expression_aux,geni{j});     %this and the following instruction find the locations of the gene 'bXXXX' in the array reaction_expression
    pos_genes_in_react_expr{j} = find(~cellfun('isempty', matches));
    reaction_expression_aux(pos_genes_in_react_expr{j}) = strrep(reaction_expression_aux(pos_genes_in_react_expr{j}),geni{j},'');   %replaced with empty char so it's not found again later (this avoid for instance looking for 'HGNC:111' and replacing also partially another gene named for instance 'HGNC:1113'
end

%% everything follows is useless because it's all in fbamodel.genes
% we also prepare the array of genes
%geni = gpr2genes(fbamodel.grRules);
%geni = unique(geni); %removes duplicates, we need unique array of genes, without duplicates
%reaction_expression = compute_geneset_expression('geo_s_react',hacky_genes_reactions_conversion('geo_s_react'));

