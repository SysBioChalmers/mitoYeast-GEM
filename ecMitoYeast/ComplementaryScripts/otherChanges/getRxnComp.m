%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rxnComps = getRxnComp(model,rxnList)
%
%
%
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rxnComps = getRxnComp(model,rxnList)

% get the compartments of all metabolites involved in the reaction(s)
if iscell(rxnList)
    [tmp,rxnInd] = ismember(rxnList,model.rxns);
    rxnComps = cell(length(rxnList),1);
    for i = 1:length(rxnInd)
        rxnID        = rxnInd(i);
        metPos       = model.S(:,rxnID) ~= 0;
        metCompInd   = model.metComps(metPos);
        metComps     = unique(model.comps(metCompInd));
        rxnComps{i}  = metComps;
    end
else
    rxnInd = find(strcmp(model.rxns,rxnList));
    if isempty(rxnInd)
        rxnInd = 0;
    else
        metPos     = model.S(:,rxnInd) ~= 0;
        metCompInd = model.metComps(metPos);
        metComps   = unique(model.comps(metCompInd));
        rxnComps   = metComps;
    end
end

end