%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addPMFcostToRxns(model)
%
% Modifies the the transport reactions affecting the proton motive force
% (PMF) to include co-transport of a pseudometabolite, PMF, representing
% the effect on the PMF. Also couples the respiratory chain and ATP
% synthase to the proton motive force.
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,energyResults,redoxResults,growth_rate] = addPMFcostToRxns(model)

% open file with changes to proton motive force (PMF)- effecting rxns

cd ../../ComplementaryData/otherChanges/

fid = fopen('PMFrxns.tsv');
PMFrxnInfo  = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
info.rxnIDs   = PMFrxnInfo{1};
info.metCoeff = cellfun(@str2num, PMFrxnInfo{2});
info.metIDs   = PMFrxnInfo{3};
info.mettype  = PMFrxnInfo{4};
info.metComp  = PMFrxnInfo{5};
fclose(fid);

% change coefficient for reactants
for i = 1:length(info.rxnIDs)
    if strcmp(info.mettype{i},'reactant')
        info.metCoeff(i) = -info.metCoeff(i);
    end
end

% Define compartment
compData = cat(2,model.compNames,model.comps);
leftBracket = ' [';
lBracket = '[';
rightBracket = ']';
space = ' ';
[row, col] = size(compData);
for i = 1:row
    compName = char(compData(i,1));
    for j = 1:length(info.rxnIDs)
        metComp = info.metComp{j};
        if strcmp(metComp,compName)
            info.newComps(j,1) = compData(i,2);
        end
    end
end

for i = 1:length(info.rxnIDs)
    info.metNames{i,1} = [info.metIDs{i} leftBracket info.metComp{i} rightBracket];
    info.newComps{i,1}    = [lBracket info.newComps{i'} rightBracket];
end

% map metabolites to model.mets and add new metabolites
cd ../../ComplementaryScripts/otherChanges/
for i = 1:length(info.metNames)
    [~,metIndex] = ismember(info.metNames(i),model.metNames);
    if metIndex ~= 0
        info.mets{i} = model.mets{metIndex};
    elseif metIndex == 0
        newID        = getNewIndex(model.mets);
        info.mets{i} = ['s_' newID info.newComps{i}];
        model        = addMetabolite(model,info.mets{i},'metName', ...
                                     info.metNames{i});
    end
end

% Add metabolite information
cd ../../ComplementaryData/otherChanges/
fid = fopen('PMFmetAnnotation.tsv');
metAnnotation = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t', ...
                         'HeaderLines',1);
metInfo.metNames    = metAnnotation{1};
metInfo.metFormulas = metAnnotation{2};
metInfo.metCharges  = cellfun(@str2num,metAnnotation{3});
metInfo.metNotes    = metAnnotation{7};
fclose(fid);

for i = 1:length(metInfo.metNames)
    [~,metInd] = ismember(metInfo.metNames(i),model.metNames);
    if metInd ~= 0
        model.metFormulas{metInd} = metInfo.metFormulas{i};
        model.metCharges(metInd)  = metInfo.metCharges(i);
        model.metNotes{metInd}    = metInfo.metNotes{i};
    end
end

% Open file with reaction properties
fid = fopen('PMFrxnProp.tsv');
rxnProp = textscan(fid,'%s %s %s %s %s %s %s', 'Delimiter','\t','HeaderLines',1);
rxnInfo.ID           = rxnProp{1};
rxnInfo.rev          = cellfun(@str2num,rxnProp{2});
rxnInfo.GPR          = rxnProp{3};
rxnInfo.rxnNames     = rxnProp{4};
rxnInfo.rxnKEGGid    = rxnProp{5};
rxnInfo.notes        = rxnProp{6};
rxnInfo.references   = rxnProp{7};
fclose(fid);
    
% Modify reactions in model
cd ../../ComplementaryScripts/otherChanges/
energyResults = cell(size(rxnInfo.ID));
redoxResults = cell(size(rxnInfo.ID));
growth_rate = cell(size(rxnInfo.ID));
for i = 1:length(rxnInfo.ID)
    rxnID = rxnInfo.ID{i};
    if sum(strcmp(model.rxns,rxnID)) > 0
        ID = rxnID;
    elseif sum(strcmp(model.rxnNames,rxnID)) > 0
        ID = model.rxns{strcmp(model.rxnNames,rxnID)};
    else
        newID = getNewIndex(model.rxns);
        ID = ['r_' newID];
    end
    rxnInd = find(strcmp(info.rxnIDs,rxnInfo.ID{i}));
    mets = info.mets(rxnInd);
    coeff = transpose(info.metCoeff(rxnInd));
    model = addReaction(model, ID, ...
                        'reactionName',rxnInfo.rxnNames{i},...
                        'metaboliteList',mets, ...
                        'stoichCoeffList',coeff, ...
                        'reversible',rxnInfo.rev(i,1), ...
                        'geneRule',rxnInfo.GPR{i}, ...
                        'checkDuplicate',1);
    % Check if reaction was successfully added to model
    if ~strcmp(model.rxns,ID)
        continue
    end
    if strcmp(rxnID,'r_1128_REV')
        model.rxns{end} = 'r_1128_REV';
        ID = 'r_1128_REV';
    end
    [energyRes,redoxRes,growth] = checkATPandRedox(model,ID);
    energyResults{i} = energyRes;
    redoxResults{i}  = redoxRes;
    growth_rate{i}   = growth;
end

% Add standard gene names for new genes
cd ../../ComplementaryData/databases/
fid = fopen('SGDgeneNames.tsv');
SGDgeneAnnotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

geneIndex = zeros(1,1);
for i = 1:length(model.genes)
    geneIndex = strcmp(SGDgeneAnnotation{1},model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(SGDgeneAnnotation{2}{geneIndex})
        model.geneNames{i} = SGDgeneAnnotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = ['COBRAprotein' num2str(i)];
end

% Add reaction annotation
for i = 1:length(rxnInfo.ID)
    [~,rxnID] = ismember(rxnInfo.ID,model.rxnNames);
    if rxnID ~= 0
        model.rxnNames{rxnID}  = rxnInfo.rxnNames{i};
        model.notes{i}         = rxnInfo.notes{i};
        model.rxnReferences{i} = rxnInfo.references{i};
    end
end
cd ../../ComplementaryScripts/otherChanges/
end
