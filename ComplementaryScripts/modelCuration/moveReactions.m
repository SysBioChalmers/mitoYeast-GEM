%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = moveLipidReactions(model)
%
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = moveReactions(model)
model = creategrRulesField(model);
% Read file with rxn information
cd ../../ComplementaryData/modelCuration/
fid  = fopen('rxnsToMove.tsv');
info = textscan(fid,'%s%s%s%s','HeaderLines',1,'Delimiter','\t');
fclose(fid);

rxnInfo.rxnIDs   = info{1};
rxnInfo.rxnName  = info{2};
rxnInfo.newComps = info{4};

% Add mitochondrial outer membrane compartment
model.comps{end+1}      = 'om';
model.compNames{end+1}  = 'mitochondrial outer membrane';

% go through rxn list and add all necessary metabolites
for i = 1:length(rxnInfo.rxnIDs)
    newComp     = rxnInfo.newComps{i};
    rxnID       = rxnInfo.rxnIDs{i}; 
    rxnPos      = strcmp(model.rxns,rxnID);
    metPos      = model.S(:,rxnPos) ~= 0;
    metIDs      = model.mets(metPos);
    metNames    = model.metNames(metPos);
    metFormulas = model.metFormulas(metPos);
    metChEBIID  = model.metChEBIID(metPos);
    metKEGGID   = model.metKEGGID(metPos);
    metCharge   = model.metCharges(metPos);
    cd ../../ComplementaryScripts/otherChanges/
    for j = 1:length(metIDs)
        newID           = getNewIndex(model.mets);
        [baseMetName,~] = extractBaseMetNames(metNames(j));
        if strcmp(newComp,'mitochondrial outer membrane')
            newMetName      = [baseMetName{1} '[mitochondrial outer membrane]'];
            ID              = ['s_' num2str(newID) '[om]'];
        elseif strcmp(newComp,'mitochondrial membrane')
            newMetName      = [baseMetName{1} '[mitochondrial membrane]'];
            ID              = ['s_' num2str(newID) '[mm]'];
        end
        if ~strcmp(model.metNames,newMetName)
            model = addMetabolite(model,ID,newMetName, ...
                                  metFormulas{j},metChEBIID{j}, ...
                                  metKEGGID{j},'','',metCharge(j));
        end    
    end
end

% Add reactions
metsToRemove = {};
for i = 1:length(rxnInfo.rxnIDs)
    newComp      = rxnInfo.newComps{i};
    rxnID        = rxnInfo.rxnIDs{i};
    rxnPos       = strcmp(model.rxns,rxnID);
    metPos       = model.S(:,rxnPos) ~= 0;
    metNames     = model.metNames(metPos);
    metIDs       = model.mets(metPos);
    metsToRemove = [metsToRemove;metIDs];
    stoichCoeffs = zeros(size(metNames));
    baseMetNames = extractBaseMetNames(metNames);
    newMetNames  = cell(size(metNames));
    newMetIDs    = cell(size(metIDs));
    for j = 1:length(metNames)
        stoichCoeffs(j) = model.S(strcmp(model.mets,metIDs{j}),rxnPos);
        newMetName      = [baseMetNames{j} '[' newComp ']'];
        newMetNames{j}  = newMetName;
        newMetIDs{j}    = model.mets{strcmp(model.metNames,newMetName)};
    end
    rxnName      = model.rxnNames{rxnPos};
    if model.lb(rxnPos) < 0 && model.ub > 0
        rxnRev = true;
    else
        rxnRev = false;
    end
    rxnLB        = model.lb(rxnPos);
    rxnUB        = model.ub(rxnPos);
    rxnSubSystem = model.subSystems{rxnPos};
    rxngrRule    = model.grRules{rxnPos};
    model        = addReaction(model, rxnID, ...
                            'reactionName',rxnName, ...
                            'metaboliteList',newMetIDs, ...
                            'stoichCoeffList',stoichCoeffs, ...
                            'reversible',rxnRev, ...
                            'lowerBound',rxnLB, ...
                            'upperBound',rxnUB, ...
                            'subSystem',rxnSubSystem, ...
                            'geneRule',rxngrRule, ...
                            'checkDuplicate',1);
end
metsToRemove = unique(metsToRemove);
% remove unused metabolites
for i = 1:length(metsToRemove)
    metID  = metsToRemove{i};
    metPos = strcmp(model.mets,metID);
    if sum(model.S(metPos,:)) == 0
        model = removeMetabolites(model,metID);
    end
end
cd ../modelCuration/
end