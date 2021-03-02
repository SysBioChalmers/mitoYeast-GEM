%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script adds new mitochondrial reactions and curations of already
% existing reactions
% Input: model, newRxnMatrix.tsv, newRxnMetAnnotation.tsv, newRxnProp.tsv 
% and SGDgeneNames.txt
%
% Function adapted from yeastGEM: https://github.com/SysBioChalmers/yeast-GEM
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = addMitochondrialRxns(model)

cd ../../ComplementaryData/modelCuration/

% open file with new reactions
fid = fopen('newRxnMatrix.tsv');
rxnInfo = textscan(fid,'%s %s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1);
info.rxnIDs = rxnInfo{1};
info.metCoeff = cellfun(@str2num, rxnInfo{2});
info.metIDs = rxnInfo{3};
info.mettype = rxnInfo{4};
info.metcompartments = rxnInfo{5};
fclose(fid);

% Change coefficient of reactants
for i = 1:length(info.rxnIDs)
    if strcmp(info.mettype(i),'reactant')
        info.metCoeff(i) = info.metCoeff(i)*-1;
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
    compName = compData(i,1);
    compName = char(compName);
    for j = 1:length(info.rxnIDs)
        metComp = info.metcompartments(j,1);
        metComp = char(metComp);
        if strcmp(metComp,compName)
            info.newComps(j,1) = compData(i,2);
        end
    end
end
for i = 1:length(info.rxnIDs)
    info.metnames(i) = strcat(info.metIDs(i),leftBracket,info.metcompartments(i),rightBracket);
    info.newComps(i) = strcat(lBracket,info.newComps(i),rightBracket);
end

% Define metabolite s_ index for new metabolites
cd ../../ComplementaryScripts/otherChanges/
for i = 1:length(info.metnames)
    [~,metindex] = ismember(info.metnames(i),model.metNames);
    if metindex ~= 0
        info.mets(i) = model.mets(metindex);
    elseif metindex == 0
        newID = getNewIndex(model.mets);
        info.mets(i) = strcat('s_',newID, info.newComps(i));
        model = addMetabolite(model,char(info.mets(i)), ...
            'metName',info.metnames(i));
    end
end


% Add metabolite information 
cd ../../ComplementaryData/modelCuration/
fid = fopen('newMetAnnotation.tsv');
metAnnotation = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
metInfo.metNames = metAnnotation{1};
metInfo.metFormulas = metAnnotation{2};
metInfo.metCharges = cellfun(@str2num,metAnnotation{3});
metInfo.metKEGGID = metAnnotation{5};
metInfo.metChEBIID = metAnnotation{6};
metInfo.metNotes = metAnnotation{7};
fclose(fid);

for i = 1:length(metInfo.metNames)
    [~, metID] = ismember(metInfo.metNames(i),model.metNames);
    if metID ~= 0
        model.metFormulas{metID} = metInfo.metFormulas{i};
        model.metCharges(metID)  = metInfo.metCharges(i);
        model.metKEGGID{metID}   = metInfo.metKEGGID{i};
        model.metChEBIID{metID}  = metInfo.metChEBIID{i};
        model.metNotes{metID}    = metInfo.metNotes{i};
    end
end

% Load file with reaction properties
fid = fopen('newRxnProp.tsv');
rxnProp = textscan(fid,'%s %s %s %s %s %s %s', 'Delimiter','\t','HeaderLines',1);
newRxn.ID           = rxnProp{1};
newRxn.rev          = cellfun(@str2num,rxnProp{2});
newRxn.GPR          = rxnProp{3};
newRxn.rxnNames     = rxnProp{4};
newRxn.rxnECNumbers = rxnProp{5};
newRxn.subSystems   = rxnProp{6};
newRxn.rxnKEGGID    = rxnProp{7};
fclose(fid);

% Add new reactions to model
cd ../../ComplementaryScripts/otherChanges/
for i = 1:length(newRxn.ID)
    if sum(strcmp(model.rxns,newRxn.ID{i})) > 0
        ID = newRxn.ID{i};
    else
        newID = getNewIndex(model.rxns);
        ID = ['r_' newID];
    end
    rxnInd = find(strcmp(info.rxnIDs,newRxn.ID{i}));
    met    = info.mets(rxnInd);
    Coeff  = transpose(info.metCoeff(rxnInd));
    [model,rxnID]  = addReaction(model, ...
                         ID, ...
                         'reactionName', newRxn.ID{i}, ...
                         'metaboliteList', met, ...
                         'stoichCoeffList', Coeff, ...
                         'reversible', newRxn.rev(i,1), ...
                         'geneRule', newRxn.GPR{i}, ...
                         'checkDuplicate',1);
end

% Add standard gene names for new genes
cd ../../ComplementaryData/databases/
fid = fopen('SGDgeneNames.tsv');
SGDgeneAnnotation = textscan(fid,'%s%s','Delimiter','\t','HeaderLines',1);
fclose(fid);

geneIndex = zeros(1,1);
for i = 1:length(model.genes)
    geneIndex = strcmp(SGDgeneAnnotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(SGDgeneAnnotation{2}{geneIndex})
        model.geneNames{i} = SGDgeneAnnotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAprotein',num2str(i));
end

% Add reaction annotation
for i = 1:length(newRxn.ID)
    [~,rxnID] = ismember(newRxn.ID(i),model.rxnNames);
    if rxnID ~= 0
        model.rxnNames{rxnID}        = newRxn.rxnNames{i};
        model.rxnECNumbers(rxnID)    = newRxn.rxnECNumbers(i);
        model.rxnKEGGID(rxnID)       = newRxn.rxnKEGGID(i);
    end
end
cd ../../ComplementaryScripts/
end