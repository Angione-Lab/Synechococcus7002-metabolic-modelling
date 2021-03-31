% Create cell arrays to store subsystems from each column of the cell array fbamodel.subSystems or new_subsystems
first_subsystems = cell(numel(list_subsystems),1);
second_subsystems = cell(numel(list_subsystems),1);
third_subsystems = cell(numel(list_subsystems),1);
fourth_subsystems = cell(numel(list_subsystems),1);
fifth_subsystems = cell(numel(list_subsystems),1);

% Retrieve subsystem names for each subsystem ordinal
for k = 1 : length(fbamodel.subSystems)
  thisCellContents = fbamodel.subSystems{k};
  % Get the first subsystem for all reactions
  first_subsystems{k} = thisCellContents{1};
  if length(thisCellContents) > 1
  % Get the second subsystem if present
  second_subsystems{k} = thisCellContents{2};
  else
  % If there is only one subsystem for the reaction, assign the second a blank []
  second_subsystems{k} = [];
  end
  if length(thisCellContents) > 2
  % Get the third subsystem if present
  third_subsystems{k} = thisCellContents{3};
  else
  % If there are no more than two subsystems for the reaction, assign the third a blank []
  third_subsystems{k} = [];
  end
  if length(thisCellContents) > 3
  % Get the fourth subsystem if present
  fourth_subsystems{k} = thisCellContents{4};
  else
  % If there are no more than three subsystems for the reaction, assign the fourth a blank []
  fourth_subsystems{k} = [];
  end
  if length(thisCellContents) > 4
  % Get the fifth subsystem if present
  fifth_subsystems{k} = thisCellContents{5};
  else
  % If there are no more than four subsystems for the reaction, assign the fifth a blank []
  fifth_subsystems{k} = [];
  end
end

%% Find indices of reactions within each subsystem
% Specify the number of unique subsystems
%N = size(list_subsystems,1);, %or
N = length(list_subsystems);
% Find the number of reactions that have a certain number of subsystems by indexing all non-blank cells
% All reactions have atleast one subsystem, but the number of reactions decreases for each successive ordinal 
% R_1 = size(find(~cellfun(@isempty,first_subsystems)),1); % first/second/third/fourth/fifth

% Create empty cell arrays (with length of list_subsystems) to store reaction indices of each number of subsystems
ix_first = cell(N,1);
ix_second = cell(N,1);
ix_third = cell(N,1);
ix_fourth = cell(N,1);
ix_fifth = cell(N,1);

% Retrieve reaction indices for each unique subsystem that matches each set of subsystems
for s = 1:N
    ix_first{s} = find(strcmpi(list_subsystems{s},first_subsystems));
    ix_second{s} = find(strcmpi(list_subsystems{s},second_subsystems));
    ix_third{s} = find(strcmpi(list_subsystems{s},third_subsystems));
    ix_fourth{s} = find(strcmpi(list_subsystems{s},fourth_subsystems));
    ix_fifth{s} = find(strcmpi(list_subsystems{s},fifth_subsystems));
end
% or an alternative loop:

% for s = 1:N
%     ix_first{s} = find(ismember(first_subsystems,list_subsystems{s}));
% end

% Concatenate all five columns 
ix_all = horzcat(ix_first,ix_second,ix_third,ix_fourth,ix_fifth);
% Merge columns to get total list of indices for each subsystem
ixs_subsystems = cell(length(ix_all),1);
for a=1:length(ixs_subsystems)
   ixs_subsystems{a}=vertcat(ix_all{a,:});
end

%% Find number of reactions in each subsystem then calculate sum and average of contributions for principal components 1 and 2
cardinality_subsystems = zeros(numel(list_subsystems),1);
sum_contrib_subsystems_PC1 = zeros(numel(list_subsystems),1);
sum_contrib_subsystems_PC2 = zeros(numel(list_subsystems),1);
avg_contrib_subsystems_PC1 = zeros(numel(list_subsystems),1);
avg_contrib_subsystems_PC2 = zeros(numel(list_subsystems),1);

for i=1:numel(ixs_subsystems)
    sum_contrib_subsystems_PC1(i) = sum(Dim_1_and_2(ixs_subsystems{i},1)); 
    sum_contrib_subsystems_PC2(i) = sum(Dim_1_and_2(ixs_subsystems{i},2)); 
    cardinality_subsystems(i) = numel(ixs_subsystems{i});
    avg_contrib_subsystems_PC1(i) = sum_contrib_subsystems_PC1(i)./cardinality_subsystems(i);
    avg_contrib_subsystems_PC2(i) = sum_contrib_subsystems_PC2(i)./cardinality_subsystems(i);
end