% FLUX_BALANCE Flux-balance analysis of FBA model
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL) performs a basic flux-balance
%    analysis of the FBA model FBAMODEL and returns a biomass-maximizing
%    flux distribution in the vector V.  The maximimum and minimum 
%    synthetic objective possible in a biomass-maximizing flux distribution
%    are given in FMAX and FMIN, respectively.
%
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL, QUIET) performs the analysis
%    and supresses screen output if QUIET is set to true.

function [FBAsolution, vmax, fmax, vbiomass] = flux_balance_minNorm(fbamodel, quiet)

% if number of arguments
if nargin < 2
    quiet = false;
end

param.tmlim  = -1;
param.msglev = 1;
param.save   = 0;
minNorm = 1e-6;

nrxn   = numel(fbamodel.rxns);
nmetab = numel(fbamodel.mets);

yt = ones(nrxn,1); %the old yt = fbamodel.present, meaning that all the reactions must be considered as active in the model

A = [ fbamodel.S; 
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ]; %stoichiometric matrix
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1) ]; % steady state so b all 0
ctype = char('S' * ones(1, nmetab + nnz(~yt))); % 'S' Fixed Variable (A(i,:)*x = b(i)).
vartype = char('C' * ones(1, nrxn)); % C = continuous variable
% v is fluxes for all reactions, vbiomass is the optimal flux for the
% objective function
[v, vbiomass] = glpk(fbamodel.f, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1); 
%%
A = [ fbamodel.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      fbamodel.f' ]; % add first objective as a metabolite constraint 
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass ]; % change the first objective concentration target to be its optimum max value
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 1)); %S = equality for first objective
%[vmin, fmin] = glpk(fbamodel.g, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, 1);
[vmax, fmax] = glpk(fbamodel.g, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1); %vmax = fluxes, fmax = optimal flux for 2nd objective

    
if ~quiet
    fprintf('Flux balance Biomass flux:    %f\n', vbiomass)
    fprintf('Flux balance secondary obj flux:  %f\n', fmax)
end

%% minNorm
% Set up the LP problem struct
[nMets,nRxns] = size(A);
LPproblem.A = A;
LPproblem.c = fbamodel.g; % second objective
LPproblem.b = b;
LPproblem.lb = fbamodel.lb;
LPproblem.ub = fbamodel.ub;
LPproblem.csense(1:nMets,1) = 'E'; % set constraints to equalities (S = b)
LPproblem.osense = 1; %minimise the fluxes other than the first and second objectives 



%if length(minNorm)> 1 || minNorm > 0
    if nnz(LPproblem.c)>1 %nnz number of non-zero matrix elements
        error('Code assumes only one non-negative coefficient in linear part of objective');
    end
    % quadratic minimization of the norm.
    % set previous optimums as constraint.
    LPproblem.A = [LPproblem.A;
        (LPproblem.c'~=0 + 0)];%new constraint must be a row with a single unit entry
    LPproblem.csense(end+1) = 'E';
    
    LPproblem.b = [LPproblem.b; fmax];
    LPproblem.c = zeros(size(LPproblem.c)); % no need for c anymore.
    %Minimise Euclidean norm using quadratic programming
    if length(minNorm)==1
        minNorm=ones(nRxns,1)*minNorm;
    end
    LPproblem.F = spdiags(minNorm,0,nRxns,nRxns);
    
    %quadratic optimization will get rid of the loops unless you are maximizing a flux which is
    %part of a loop. By definition, exchange reactions are not part of these loops, more
    %properly called stoichiometrically balanced cycles.
    solution = solveCobraQP(LPproblem);  
    
    FBAsolution = solution.full(1:nRxns); % fluxes values
    
    %objective = solution.obj;
    %FBAsolution.f = model.c'*solution.full(1:nRxns); %objective from original optimization problem.
    %if abs(FBAsolution.f - objective) > .01  
     %   error('optimizeCbModel.m: minimizing Euclidean norm did not work')
    %end
