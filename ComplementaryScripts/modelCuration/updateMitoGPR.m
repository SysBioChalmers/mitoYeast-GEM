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
updateGPR       = textscan(fid, '%s%s%s', 'Delimiter', '\t', 'HeaderLines',1);
newGPR.rxnID    = updateGPR{1};
newGPR.oldGPR   = updateGPR{2};
newGPR.GPR      = updateGPR{3};
fclose(fid);

% update GPR
for i = 1:length(newGPR.rxnID)
   rxnInd = find(strcmp(model.rxns, newGPR.rxnID(i)));
   model  = changeGeneAssociation(model, model.rxns{rxnInd}, newGPR.GPR{i});
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
