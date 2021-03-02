%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMitochondrialGPR(model)
% Change the GPR according to the available localization evidence of 
% enzymes based on the literature, SGD and UniProt
%
% NOTE: Changes are listed in updateMitoGPR.tsv
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = updateMitoGPR(model)
% Load file with updated GPRs
cd ../../ComplementaryData/modelCuration/
fid             = fopen('updateMitoGPRs.tsv');
updateGPR       = textscan(fid, '%s%s%s%s', 'Delimiter', '\t', 'HeaderLines',1);
newGPR.rxnID    = updateGPR{1};
newGPR.rxnName  = updateGPR{2};
newGPR.oldGPR   = updateGPR{3};
newGPR.GPR      = updateGPR{4};
fclose(fid);

% update GPR
for i = 1:length(newGPR.rxnID)
    rxnID  = newGPR.rxnID{i};
    rxnInd = strcmp(model.rxns,rxnID);
    if sum(rxnInd) == 1 && strcmp(model.rxnNames{rxnInd},newGPR.rxnName{i})
        model  = changeGeneAssociation(model, model.rxns{rxnInd}, newGPR.GPR{i});
    else
        disp(['RxnID ' rxnID ' does not match with the rxnName specified']);
    end
end

% Delete any unused genes
model = removeUnusedGenes(model);

% Add standard gene names for any new genes
cd ../databases/
fid            = fopen('SGDgeneNames.tsv');
SGD_annotation = textscan(fid, '%s%s', 'Delimiter', '\t', 'HeaderLines', 1);
model.geneNames = cell(size(model.genes));
for i = 1:length(model.genes)
    geneInd = strcmp(SGD_annotation{1},model.genes{i});
    if sum(geneInd) == 1 && ~isempty(SGD_annotation{2}{geneInd})
        model.geneNames{i} = SGD_annotation{2}{geneInd};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein names for any new genes
for i = 1:length(model.genes)
   model.proteins{i} = strcat('COBRAProtein',num2str(i)); 
end

cd ../../ComplementaryScripts/
end
