%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addTranslocationComplexes to model
%
%
%
% Carl Malina   Last edited. 2020-06-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addTranslocationComplexesToModel(model)
tic
% Read info for translocation complexes
current = pwd;
cd ../../ComplementaryData/proteinImport/
fid                          = fopen('importMachineryComponents.tsv');
fileData                     = textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
importMachinery.uniprots     = fileData{1};
importMachinery.genes        = fileData{2};
importMachinery.geneShort    = fileData{3};
importMachinery.component    = fileData{5};
importMachinery.localization = fileData{6};
fclose(fid);

% Add genes that are not present in model already
% Check which genes are already present in model
inModel = ismember(importMachinery.genes,model.genes);
genesToAdd.genes = importMachinery.genes(~inModel);
genesToAdd.geneShortNames = importMachinery.geneShort(~inModel);

model = addGenesRaven(model,genesToAdd);

% read kcat info for complexes
fid = fopen('kcatsTranslocationComplexes.txt');
data = textscan(fid,'%s%s%s','Delimiter','\t','HeaderLines',1);
importComplexInfo.complex = data{1};
importComplexInfo.kcats   = cellfun(@str2double,data{2});
fclose(fid);

% read file with localization info for import complex subunits
fid = fopen('enzymeLocalization.tsv');
fileData = textscan(fid,'%s%s%s%s%s%s%s%s%s%s','Delimiter','\t',...
                     'HeaderLines',1);
enzymeLoc.uniprots         = fileData{1};
enzymeLoc.genes             = fileData{2};
enzymeLoc.geneShort         = fileData{3};
enzymeLoc.localization      = fileData{4};
enzymeLoc.pathway           = fileData{5};
enzymeLoc.MPPcleavageSite   = cellfun(@str2num,fileData{7},'UniformOutput',false);
enzymeLoc.OCT1cleavageSite  = cellfun(@str2num,fileData{8},'UniformOutput',false);
enzymeLoc.ICP55cleavageSite = cellfun(@str2num,fileData{9},'UniformOutput',false);
fclose(fid);
cd ../../ComplementaryScripts/customGECKO/change_model/
for i = 1:length(importComplexInfo.complex)
    complex = importComplexInfo.complex{i};
    kcat = importComplexInfo.kcats(i);
    if kcat ~= 0
        % Get list of complex subunits
        componentPos = strcmp(importMachinery.component,complex);
        uniprotIDs = importMachinery.uniprots(componentPos);
        kcats = kcat*ones(size(uniprotIDs));
        % add proteins to model
        model = addEnzymesToModel(model,uniprotIDs,kcats);
    end    
end

% Add translocation reaction for import machinery complexes
cd(current)
model = addTranslocationRxns(model,enzymeLoc);

% Add usage of subunits of import machinery to translocation reactions
cd ../miscellaneous/
for i = 1:length(importComplexInfo.complex)
    complex = importComplexInfo.complex{i};
    % Enable matching complex name to pathways in enzymeLoc
    if strcmpi(complex,'Disulfide relay')
        complexName = 'MIA40-ERV1';
    elseif strcmpi(complex,'MPP')
        complexName = 'TIM23';
    else
        complexName = complex;
    end
    kcat = importComplexInfo.kcats(i);
    if kcat ~= 0
        coeff = -1/kcat;
        componentPos =strcmp(importMachinery.component,complex);
        genes = importMachinery.genes(componentPos);
        uniprotIDs = importMachinery.uniprots(componentPos);
        % Loop through each subunit of complex
        for j = 1:length(genes)
            gene = genes{j};
            if strcmp(gene,'YJR045C')
                disp('YJR045C')
            end
            genePos = strcmp(enzymeLoc.genes,gene);
            uniprot = uniprotIDs{j};
            comp    = lower(enzymeLoc.localization{genePos});
            if strcmpi(comp,'im')
                comp = 'mm';
            end
            metID = ['holo_prot_' uniprot '_' comp];
            metInd = strcmp(model.mets,metID);
            translocationRxnInd = find(contains(model.rxns,' translocation'));
            y = 0;
            % Identify translocation reaction requiring complex subunit
            for k = 1:length(translocationRxnInd)
                rxnInd = translocationRxnInd(k);
                rxnID = model.rxns{rxnInd};
                enzyme = strsplit(rxnID,'_');
                enzyme = strsplit(enzyme{3},' ');
                enzyme = enzyme{1};
                if strcmp(enzyme,'P0CS90') && strcmp(gene,'YJR045C')
                    disp('SSC1')
                end
                enzymePos = strcmp(enzymeLoc.uniprots,enzyme);
                pathway   = enzymeLoc.pathway{enzymePos};
                % Add protein usage to reactions
                if contains(model.grRules{rxnInd},gene) && contains(pathway,complexName) %&& ~contains(rxnID,uniprot)
                    if strcmp(enzyme,uniprot)
                        continue
                    end
                    model.S(metInd,rxnInd) = coeff;
                    printRxnFormulaRAVEN(model,{rxnID},true,true,false);
%                     sol = optimizeCbModel(model);
%                     %disp(num2str(sol.obj))
%                     if sol.obj < 0.001
%                         disp(num2str(0))
%                         warning(['Including protein ' uniprot ' renders model unable to grow'])
%                     else
%                         disp(num2str(sol.obj))
%                     end
                elseif strcmp(complex,'TOM') && ~contains(pathway,'OM') && ~contains(pathway,'MIM1')
                    if strcmp(enzyme,uniprot)
                        continue
                    end
                    model.S(metInd,rxnInd) = coeff;
                    printRxnFormulaRAVEN(model,{rxnID},true,true,false);
%                     sol = optimizeCbModel(model);
%                     if sol.obj < 0.001
%                         disp(num2str(0))
%                         warning(['Including protein ' uniprot ' renders model unable to grow'])
%                     else
%                         disp(num2str(sol.obj))
%                     end
                elseif contains(model.grRules{rxnInd},gene) && strcmp(gene,'YKL134C') && ~isempty(enzymeLoc.OCT1cleavageSite{enzymePos})
                    y = y+1;
                    if strcmp(enzyme,uniprot)
                        continue
                    end
                    model.S(metInd,rxnInd) = coeff;
                    printRxnFormulaRAVEN(model,{rxnID},true,true,false);
%                     sol = optimizeCbModel(model);
%                     if sol.obj < 0.001
%                         disp(num2str(0))
%                         warning(['Including protein ' uniprot ' renders model unable to grow'])
%                     else
%                         disp(num2str(sol.obj))
%                     end
                elseif contains(model.grRules{rxnInd},gene) && strcmp(gene,'YER078C') && ~isempty(enzymeLoc.ICP55cleavageSite{enzymePos})
                    if strcmp(enzyme,uniprot)
                        continue
                    end
                    model.S(metInd,rxnInd) = coeff;
                    printRxnFormulaRAVEN(model,{rxnID},true,true,false);
%                     sol = optimizeCbModel(model);
%                     if sol.obj < 0.001
%                         disp(num2str(0))
%                         warning(['Including protein ' uniprot ' renders model unable to grow'])
%                     else
%                         disp(num2str(sol.obj))
%                     end
                end
            end
        end            
    end
end
% update rxnGeneMat based on grRules field
if isfield(model,'rules')
    model = rmfield(model,'rules');
end
model = buildRxnGeneMat(model);
% remove unused genes from model
genesToRemoveInd = find(sum(model.rxnGeneMat) == 0);
model.rxnGeneMat(:,genesToRemoveInd) = [];
model.genes(genesToRemoveInd) = [];
model.geneShortNames(genesToRemoveInd) = [];
cd(current)
toc
end