function solution = simulateChemostat(model,gRate,pos,minProt)
% Function that sets chemostat constraints in the ecModel and runs an FBA
% first by minimizing glucose uptake rate and then minimizing the protein
% usage
%
%   model       (struct) MATLAB structure
%   pos         (vector) index for main carbon source uptake reaction and 
%               growth rxn
%   gRate       Dilution rate for the desired chemostat condition
%   minProt     Minimize protein usage after after fixing a minimal carbon
%               source uptake (Default = False)
%
% Usage:    solution = simulateChemostat(model,gRate,pos,minProt)
%
% Function adapted from GECKO:https://github.com/SysBioChalmers/GECKO
%
% Last modified.    Carl Malina 2020-08-30

if nargin < 4
    minProt = false;
end
cSource = pos(1);
gPos    = pos(2);
protPos = find(contains(model.rxnNames,'prot_pool'));

% Set growth rate and incorporate the growth rate in constraints
model.lb(gPos) = gRate;
model.ub(gPos) = 1.01*gRate;
modelWithMuConstraints = incorporateMuInConstraints2(model,gRate);

% Set minimization of carbon source uptake as objective
modelWithMuConstraints.c(:) = 0;
modelWithMuConstraints.c(cSource) = 1;

% Run Optimization
solution = optimizeCbModel(modelWithMuConstraints,'min','one');
%solution = optimizeCbModel(modelWithMuConstraints,'min');
solution = solution.x;
if ~isempty(solution)
    if minProt
        modelWithMuConstraints.lb(cSource) = 0.99999*solution(cSource);
        modelWithMuConstraints.ub(cSource) = 1.00001*solution(cSource);
        % Set protein usage rxn as objective
        modelWithMuConstraints.c(:)       = 0;
        modelWithMuConstraints.c(protPos) = 1;
        try
            solution_tmp = optimizeCbModel(modelWithMuConstraints,'min','one');
        catch
            solution_tmp.x = [];
        end
        %solution_tmp = optimizeCbModel(modelWithMuConstraints,'min','one');
        if ~isempty(solution_tmp.x)
            solution = solution_tmp.x;
        else
            disp('Protein minimization not feasible')
        end
    end
else
    disp('Chemostat conditions too stringent')
end

end