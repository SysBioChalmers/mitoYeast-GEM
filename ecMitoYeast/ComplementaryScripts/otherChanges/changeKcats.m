function model = changeKcats(model,protIDs,newValues,isImportProtein)
% changeKcats
%   function that changes the stoichiometric coefficients of enzymes based
%   on a kcats value
if ~iscell(protIDs)
    protIDs = {protIDs};
end

if nargin < 4
    isImportProtein = false;
end

% Udate kcats
for i = 1:length(protIDs)
    protPos  = contains(model.mets,protIDs{i});
    rxnPos   = find(model.S(protPos,:)~=0);
    % Get only metabolic reactions or translocation rxns
    if isImportProtein
        uniprotID = strsplit(protIDs{i},'_');
        uniprotID = uniprotID{3};
        rxns      = model.rxns(rxnPos);
        rxnPos    = rxnPos(~contains(rxns,['prot_' uniprotID]));
    else
        rxnPos   = contains(model.rxns(rxnPos),'No');
    end
    newValue = newValues(i);
    model.S(protPos,rxnPos) = -1/(newValue*3600); % 1/s -> 1/h
end

