%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addIMScomp(model)
%
% This function creates a compartment for the mitochondrial intermembrane space
% and adds reactions occuring in the compartment
% Input: model, newIMSrxnMatrix.tsv, newIMSetAnnotation.tsv, newIMSrxnProp.tsv 
% SGDgeneNames.txt
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = addIMScomp(model)

% Add compartments name and abbreviation for intermembrane space
model.comps{end+1}     = 'ims';
model.compNames{end+1} = 'mitochondrial intermembrane space';

% Load file containing rxn data
cd ../../ComplementaryData/modelCuration/
% remove previously deined occurences of rxnInfo
clear rxnInfo info
fid           = fopen('newIMSrxnMatrix.tsv');
rxnInfo       = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
info.rxnIDs    = rxnInfo{1};
info.metCoeff = cellfun(@str2num, rxnInfo{2});
info.metIDs   = rxnInfo{3};
info.mettype  = rxnInfo{4};
info.metCompartments = rxnInfo{5};
fclose(fid);

% Change the sign of reactant coefficients
for i = 1:length(info.rxnIDs)
    if strcmp(info.mettype(i),'reactant')
        info.metCoeff(i) = info.metCoeff(i)*-1;
    end
end

% Define the compartment
compData = cat(2,model.compNames,model.comps);
leftBrack   = ' [';
lBrack      = '[';
rightBrack  = ']';
space       = ' ';
[row,col]   = size(compData);
for i = 1:row
    compName = compData(i,1);
    compName = char(compName);
    for j = 1:length(info.rxnIDs)
        metComp = info.metCompartments(j,1);
        metComp = char(metComp);
        if strcmp(metComp,compName)
            info.newComps(j,1) = compData(i,2);
        end
    end
end

for i = 1:length(info.rxnIDs)
    info.metNames(i) = strcat(info.metIDs(i),leftBrack,info.metCompartments(i),rightBrack);
    info.newComps(i) = strcat(lBrack,info.newComps(i),rightBrack);
end

% Define new metabolite s_ index for new metabolites
cd ../../ComplementaryScripts/otherChanges/
for i = 1:length(info.metNames)
    [~,metIndex] = ismember(info.metNames(i), model.metNames);
    if metIndex  ~= 0
        info.mets(i) = model.mets(metIndex);
    elseif metIndex == 0
        newIndex = getNewIndex(model.mets);
        info.mets(i) = strcat('s_',newIndex,info.newComps(i));
        model = addMetabolite(model,char(info.mets(i)), ...
            'metName',info.metNames(i));
    end
end

% Add metabolite information
cd ../../ComplementaryData/modelCuration/
clear metAnnotation metInfo
fid = fopen('newIMSmetAnnotation.tsv');
metAnnotation = textscan(fid,'%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
metInfo.metNames = metAnnotation{1};
metInfo.metFormulas = metAnnotation{2};
metInfo.metCharges = cellfun(@str2num,metAnnotation{3});
metInfo.metKEGGID = metAnnotation{5};
metInfo.metChEBIID = metAnnotation{6};
metInfo.metNotes = metAnnotation{7};
fclose(fid);

for i = 1:length(metInfo.metNames)
    [~,metID] = ismember(metInfo.metNames,model.metNames);
    if metID ~= 0
        model.metFormulas{metID} = metInfo.metFormulas{i};
        model.metKEGGID{metID} = metInfo.KEGGID{i};
        model.metChEBIID{metID} = metInfo.metChEBIID{i};
        model.metNotes{metID} = metInfo.metNotes{i};
    end
end

% Load file with reaction properties
clear rxnProp newRxn
fid = fopen('newIMSrxnProp.tsv');
rxnProp = textscan(fid,'%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
newRxn.ID = rxnProp{1};
newRxn.rev = cellfun(@str2num,rxnProp{2});
newRxn.GPR = rxnProp{3};
newRxn.rxnNames = rxnProp{4};
newRxn.rxnECnumbers = rxnProp{5};
newRxn.subSystems = rxnProp{6};
newRxn.rxnKEGGID = rxnProp{7};
fclose(fid);

% Add new reactions to model
cd ../../ComplementaryScripts/otherChanges/
for i = 1:length(newRxn.ID)
    newID = getNewIndex(model.rxns);
    rxnInd = find(strcmp(info.rxnIDs,newRxn.ID{i}));
    met = info.mets(rxnInd);
    Coeff = transpose(info.metCoeff(rxnInd));
    model = addReaction(model, ...
                        ['r_' newID], ...
                        'reactionName', newRxn.ID{i}, ...
                        'metaboliteList',met, ...
                        'stoichCoeffList', Coeff, ...
                        'reversible',newRxn.rev(i,1), ...
                        'geneRule',newRxn.GPR{i}, ...
                        'checkDuplicate',1);
end

% Add reaction annotation
for i = 1:length(newRxn.ID)
    [~,rxnID] = ismember(newRxn.ID(i),model.rxnNames); % CHECK if it should
    % really be newRxn.ID(i) or if it should be newRxn.rxnNames
    if rxnID ~= 0
        model.rxnNames{rxnID} = newRxn.rxnNames{i};
        model.rxnECNumbers{rxnID} = newRxn.rxnECnumbers{i};
        model.rxnKEGGID{rxnID} = newRxn.rxnKEGGID{i};
    end
end

cd ../
end


