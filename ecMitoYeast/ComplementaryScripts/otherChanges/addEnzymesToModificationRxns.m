%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addEnzymesToModificationRxns(model)
%
%
%
%
% Carl Malina       Last edited. 2020-05-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addEnzymesToModificationRxns(model)
%Read kcat data:
cd ../../ComplementaryData/cofactors/
fID           = fopen('kcatData_cofactorModifications.txt');
data          = textscan(fID,'%s %s %s %s %s %f','delimiter','\t', ...
                         'HeaderLines',1);
modifiedProt  = data{1};
structure     = data{3};
protGenes     = data{5};
kcats         = data{6}.*3600;
fclose(fID);

% Load data on mitochondrial enzymes
cd ../proteinImport/
fID = fopen('enzymeLocalization.tsv');
data = textscan(fID,'%s%s%s%s%s%s%s%s%s%s','Delimiter','\t',...
                     'HeaderLines',1);
enzymeLoc.uniprots          = data{1};
enzymeLoc.genes             = data{2};
enzymeLoc.geneShort         = data{3};
enzymeLoc.localization      = data{4};
enzymeLoc.pathway           = data{5};
enzymeLoc.MPPcleavageSite   = cellfun(@str2num,data{7},'UniformOutput',false);
enzymeLoc.OCT1cleavageSite  = cellfun(@str2num,data{8},'UniformOutput',false);
enzymeLoc.ICP55cleavageSite = cellfun(@str2num,data{9},'UniformOutput',false);
fclose(fID);

% Read uniprot and kegg data
cd ../../../../GECKO/databases/
data          = load('ProtDatabase.mat');
swissprot     = data.swissprot;
kegg          = data.kegg;

cd ../../mitoYeast-GEM/ecMitoYeast/ComplementaryScripts/customGECKO/change_model/
current = pwd;

% collect uniprot IDs from data
for i = 1:length(modifiedProt)
    cd(current)
    prot   = modifiedProt{i};
    rxn    = ['prot_' prot '_modification'];
    rxnPos = contains(model.rxns,rxn);
    kcat   = kcats(i);
    % get proteins required for modification
    uniprots = strsplit(structure{i},' + ');
    
    for j = 1:length(uniprots)
        uniprot = uniprots{j};
        % check if enzyme exists in model
        enzInModel = strcmp(model.enzymes,uniprot);
        % check if enzyme is mitochondrial
        isMitochondrial = strcmp(enzymeLoc.uniprots,uniprot);
        if sum(enzInModel) == 0
            model = addEnzymesToModel(model,{uniprot},kcat);
            % check if enzyme is mitochondrial and if so, add translocation
            % rxn
            if sum(isMitochondrial) ~= 0
                cd ../../otherChanges/
                model = addTranslocationRxns(model,enzymeLoc);
            end 
        end
        % get enzyme position
        if sum(isMitochondrial) == 1
            comp = lower(enzymeLoc.localization{isMitochondrial});
            if strcmp(comp,'im')
                comp = 'mm';
            end
            enzID = ['holo_prot_' uniprot '_' comp];
        else
            enzID = ['holo_prot_' uniprot];
        end
        % add enzymes to reaction
        enzPos = strcmp(model.mets,enzID);
        coeff           = -1/kcat;
        model.S(enzPos,rxnPos) = coeff;    
    end
    cd ../../miscellaneous/
    printRxnFormulaRAVEN(model,model.rxns(rxnPos),true,true,false);
end
% remove P48445 (BPL1) from acetyl-CoA carboxylase, reaction (No1) (r_0109No1)
genePos                = strcmp(model.geneShortNames,'BPL1');
enzPos                 = strcmp(model.mets,'holo_prot_P48445');
rxnPos                 = strcmp(model.rxns,'r_0109No1');
model.S(enzPos,rxnPos) = 0;
model.grRules{rxnPos}  = 'YNR016C'; % Only leave ACC1 in grRule
model.rxnGeneMat(rxnPos,genePos) = 0;

% Update rxnGeneMat
model = buildRxnGeneMat(model);

cd(current);

end