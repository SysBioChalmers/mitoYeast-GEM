%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printRxnFormulaRAVEN(model,rxns,metNameFlag, printFlag)
%
% Inputs:
%   model: a RAVEN compatible model structure (1x1 struct)
%   rxns:  Cell array containing reaction IDs
%
% Optional inputs:          
%      metNameFlag: print full metNames instead of metIDs 
%                   (default = false)
%      printFlag: print formula or just return them (default = true)
%      grRuleFlag: print grRule or just return them (default = true)
%
% Output:
%    formulas: Cell array containing formulas of specified reactions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function formulas = printRxnFormulaRAVEN(model,rxns,metNameFlag,printFlag,grRuleFlag)

if nargin < 3
    metNameFlag = false;
    printFlag   = true;
    grRuleFlag  = false;
end

if nargin < 4
    printFlag   = true;
    grRuleFlag  = false;
end
if nargin < 5
    grRuleFlag  = false;
end

formulas = cell(length(rxns),1);
if grRuleFlag
    grRules = cell(length(rxns),1);
end

for i = 1: length(rxns)
    rxnID = rxns{i};
    rxnIndex = strcmp(model.rxns,rxnID);
   
    rxnS = full(model.S(:,rxnIndex));
   
    % Substrates
    subInd = rxnS < 0;
    subS = rxnS(subInd);
    SmetCompIndices = model.metComps(subInd);
    if metNameFlag
        Smets       = model.metNames(subInd);
        for j = 1:length(Smets)
            Smets{j} = [Smets{j} ' [' model.compNames{SmetCompIndices(j)} ']'];
        end
    else
        Smets       = model.mets(subInd);
        for j = 1:length(Smets)
            Smets{j} = [Smets{j} ' [' model.comps{SmetCompIndices(j)} ']'];
        end
    end
    
    % Products
    prodInd       = rxnS > 0;
    prodS         = rxnS(prodInd);
    PmetCompIndices = model.metComps(prodInd);
    if metNameFlag
        Pmets         = model.metNames(prodInd);
        for j = 1:length(Pmets)
            Pmets{j} = [Pmets{j} ' [' model.compNames{PmetCompIndices(j)} ']'];
        end
    else
        Pmets     = model.mets(prodInd);
        for j = 1:length(Pmets)
            Pmets{j} = [Pmets{j} ' [' model.comps{PmetCompIndices(j)} ']'];
        end
    end
    
    formulaStr = '';
    for j = 1:length(Smets)
        if (j > 1)
            formulaStr = sprintf('%s+ ', formulaStr);
        end
        if (abs(subS(j)) ~= 1)
            formulaStr = sprintf('%s%g %s ', formulaStr, abs(subS(j)), Smets{j});
        else
            formulaStr = sprintf('%s%s ', formulaStr, Smets{j});
        end
    end
    
    if model.lb(rxnIndex) < 0
        formulaStr = sprintf('%s <=> ',formulaStr);
    else
        formulaStr = sprintf('%s -> ',formulaStr);
    end
    
    for j = 1:length(Pmets)
        if (j > 1)
            formulaStr = sprintf('%s+ ', formulaStr);
        end
        if (abs(prodS(j)) ~= 1)
            formulaStr = sprintf('%s%g %s ', formulaStr, abs(prodS(j)), Pmets{j});
        else
            formulaStr = sprintf('%s%s ', formulaStr, Pmets{j});
        end
    end
    
    if grRuleFlag
        grRules{i} = model.grRules{rxnIndex};
    end
    
    formulas{i} = formulaStr;
end

if grRuleFlag
    for i = 1:length(formulas)
        formulas{i} = [formulas{i} '  ' grRules{i}];
    end
end

if printFlag
    for i = 1:length(formulas)
        disp([rxns{i} '  ' formulas{i}]);
    end
    %else
        %for i = 1:length(formulas)
            %disp([rxns{i} '  ' formulas{i}]);
        %end
    %end
    
end

end