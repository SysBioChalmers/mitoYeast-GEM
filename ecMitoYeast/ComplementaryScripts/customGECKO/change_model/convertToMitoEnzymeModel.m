%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = convertToMitoEnzymeModel(irrevModel,Genes,uniprots,kcats)
% Converts standard GEM to GEM accounting for enzymes as pseudo
% metabolites, with -(1/kcat) as the corresponding stoich. coeffs. Also
% adds reaction for coupling enzymes to cofactors; iron-sulfur clusters,
% lipoic acid various hemes. Also adds iport reactions for mitochondrial
% proteins
%
% INPUT:
% model             The GEM structure (1x1 struct)
% Genes             Genes that were matched to uniprot codes for each Rxn
% uniprots          uniprot codes of all metabolic reactions
% kcats             kcats for each enzyme/reaction
%
% OUTPUT:
% eModel            Modified GEM structure (1x1 struct)
%
% Function adapted from the GECKO toolbox:
% https://github.com/SysBioChalmers/GECKO
% 
% Carl Malina.  Last edited: 2019-12-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = convertToMitoEnzymeModel(irrevModel,Genes,uniprots,kcats)

eModel = irrevModel;
current = pwd;
% get location on change_model in GECKO
cd ../../../../../GECKO/geckomat/change_model/
change_model_GECKO = pwd;
cd(current)
% Read files with submitochondrial localization info for mitochondrial
% enzymes
cd ../../../ComplementaryData/proteinImport/
fid = fopen('enzymeLocalization.tsv');
fileData = textscan(fid,'%s%s%s%s%s%s%s%s%s%s','Delimiter','\t',...
                    'HeaderLines',1);
fclose(fid);

% extract file data
enzymeLoc.uniprots         = fileData{1};
enzymeLoc.genes             = fileData{2};
enzymeLoc.geneShort         = fileData{3};
enzymeLoc.localization      = fileData{4};
enzymeLoc.pathway           = fileData{5};
enzymeLoc.MPPcleavageSite   = cellfun(@str2num,fileData{7},'UniformOutput',false);
enzymeLoc.OCT1cleavageSite  = cellfun(@str2num,fileData{8},'UniformOutput',false);
enzymeLoc.ICP55cleavageSite = cellfun(@str2num,fileData{9},'UniformOutput',false);

% Update kcat information with values for mitochondrial carriers
fid      = fopen('mitoCarriers.txt');
fileData = textscan(fid,'%s%s%s%s%s%s%s','Delimiter','\t',...
                    'HeaderLines',1);
fclose(fid);
% Extract information
mitoCarriers.uniprots = fileData{3};
fileData{6}           = strrep(fileData{6},',','.');
mitoCarriers.kcats    = str2double(fileData{6})*3600; % s-1 -> h-1
[m,n] = size(uniprots);
% Update array with kcat information
for i = 1:length(mitoCarriers.uniprots)
    uniprotID = mitoCarriers.uniprots{i};
    for j = 1:n
        % skip looping over columns lacking data
        column = uniprots(:,j);
        % replace empty ([]) cells
        column(cellfun(@isempty,column)) = {'empty'};
        carrierInd = contains(string(column),uniprotID);
        if sum(carrierInd) == 0
            continue
        else
            kcats(carrierInd,j) = mitoCarriers.kcats(i);
        end
    end
end

cd(current)

% Add enzymes to model. All enzymes are added to cytoplasm, with the
% exception of mitochondrially synthesized enzymes. Enzymes are added as
% holo-proteins
eModel = addEnzymesToModel(eModel,uniprots,kcats);

% Add cofactor modification reactions
cd ../../otherChanges/
eModel = coFactorModifications(eModel,enzymeLoc);

% Add translocation reactions for mitochondrial enzymes
% Match model enzymes to list of mitochondrial enzymes
eModel = addTranslocationRxns(eModel,enzymeLoc);

[m,n]   = size(uniprots);

for i = 1:m
   rxnID = irrevModel.rxns{i};
   
   x = 0;
   for j = 1:n
       % Count the number of isozymes (x):
       if ~isempty(uniprots{i,j}) && kcats(i,j) > 0
           uniprots{i,j} = strsplit(uniprots{i,j},' ');
           x = x+1;
       end
   end
   
   if x > 0
       %>1 enzyme: Will include an "arm reaction" for controlling the total
       %flux thorugh system of parallel rxns.
       cd(change_model_GECKO)
       if x > 1
           eModel = addArmReaction(eModel,rxnID);
       end
       x = 0;
       % For each enzyme adds a new parallel rxn with that enzyme as pseudo
       % metabolite
       cd(current)
       newIDs = {};
       for j = 1:n
           if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 % if kcat=0 no mets will be added
               x = x + 1;
               newID   = [rxnID 'No' num2str(x)];
               newIDs = [newIDs;newID];
               newName = [irrevModel.rxnNames{i} ' (No' num2str(x) ')'];
               kvalues = ones(size(uniprots{i,j}))*kcats(i,j);
               newMets = cell(size(uniprots{i,j}));
               comp = cell(size(uniprots{i,j}));
               % check compartment of initial reaction
               cd ../../otherChanges/
               rxnComps = getRxnComp(eModel,rxnID);
               for k = 1:length(newMets)
                   newMets{k} = ['holo_prot_' uniprots{i,j}{k}];
                   isMitoEnzyme = strcmp(enzymeLoc.uniprots,uniprots{i,j}{k});
                   if length(rxnComps) == 1 && strcmp(rxnComps,'c')
                       comp{k} = 'c';
                   elseif sum(isMitoEnzyme) == 1
                       comp{k} = lower(enzymeLoc.localization{isMitoEnzyme});
                   else
                       comp{k} = 'c';
                   end
               end
               grRule = strrep(Genes{i,j},' ',' and ');
               cd ../customGECKO/change_model/
               eModel = addEnzymesToRxn(eModel,kvalues,rxnID,newMets,comp,{newID,newName},grRule);
               %printRxnFormulaRAVEN(eModel,{rxnID},true,true,true);
               %printRxnFormulaRAVEN(eModel,{newID},true,true,true);
           end
       end
       eModel = removeReactions(eModel,{rxnID}); %Remove the original rxn
   end
   if rem(i,100) == 0 || i == m
       disp(['Adding enzymes to rxns: Ready with rxn ' num2str(i)])
   end
end

end