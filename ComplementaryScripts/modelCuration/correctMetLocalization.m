%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = correctMetLocalization(model)
%
% This function corrects the metabolite localization for reactions occuring
% over the inner mitochondrial membrane. It requires that
% changeMitochondrialCompartments.m, and addIMSreactions have been run
% prior to running this script.
% Input: model, rxnsOccurringOverMembrane.tsv
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = correctMetLocalization(model)

% Load and read file containing information on changes
cd ../../ComplementaryData/modelCuration/
fid = fopen('rxnMetChanges.tsv');
% use delimiter q to deal with " introduced when reading tsv file
changes = textscan(fid,'%s%q%q','Delimiter','\t','HeaderLines',1);
rxnChanges.rxn = changes{1};
rxnChanges.mets2change = changes{2};
rxnChanges.newMets = changes{3};
fclose(fid);

% Change the metabolites used in reactions
for i = 1:length(rxnChanges.rxn)
    rxn = rxnChanges.rxn{i};
    mets2change_str = rxnChanges.mets2change{i};
    mets2change = transpose(strsplit(mets2change_str,';'));
    newMets_str = rxnChanges.newMets{i};
    newMets = transpose(strsplit(newMets_str,';'));
    % find the metabolite IDs based on the metabolite name
    newMetIDs = cell(size(newMets));
    for j = 1:length(newMets)
        newMet = newMets{j};
        oldMet = mets2change{j};
        if sum(strcmp(model.metNames,newMet)) ~= 1
            disp(strcat('could not find metabolite: ',{' '},newMet))
        elseif sum(strcmp(model.mets,oldMet)) ~= 1
            disp(strcat('could not find metabolite: ',{' '},oldMet))   
        else
            newMetID = char(model.mets(strcmp(model.metNames,newMet)));
            newMetIDs{j} = newMetID;
        end
    end
    [model,modifiedRxns] = changeRxnMets(model,mets2change,newMetIDs,rxn);
    disp(modifiedRxns)
end
cd ../../ComplementaryScripts/
end