%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = coFactorModifications(eModel, mitoEnzymes)
%
%
%
% CarlMalina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = coFactorModifications(eModel,mitoEnzymes)
current = pwd;
% change LB of biotin transport to mitochondrion
eModel.lb(strcmp(eModel.rxnNames,'biotin transport')) = -1000;

% Read file with co-factor content info
cd ../../../mitoYeast/ComplementaryData/proteinInfo/
fid    = fopen('coFactorInfo.tsv');
cfInfo = textscan(fid,'%s%s%s%s%s%s%s%s','HeaderLines',1,'Delimiter','\t');
fclose(fid);
coFactorInfo.genes   = cfInfo{1};
coFactorInfo.type    = cfInfo{4};
coFactorInfo.stoich  = cellfun(@str2num,cfInfo{5});
coFactorInfo.pool    = cfInfo{6};
coFactorInfo.grRules = cfInfo{7};

% Read file with template reactions for cofactor modification
cd(current)
filename = '../../ComplementaryData/cofactors/coFactorModificationTemplateRxns.txt';
templateRxns = readTemplateRxnFile(filename);

% Read file with gene information, based on SGD, to add new genes
cd ../../../mitoYeast/ComplementaryData/databases/
fid = fopen('SGDgeneNames.tsv');
SGDgeneAnnotation = textscan(fid,'%s%s','Delimiter','\t','HeaderLines',1);
SGDgeneName = SGDgeneAnnotation{1};
SGDgeneShortName = SGDgeneAnnotation{2};
fclose(fid);
cd(current)

% Find cofactor containing enzymes
enzGenes = eModel.enzGenes;
% get list of all unique rxns in model
for i = 1:length(enzGenes)
    gene    = enzGenes{i};
    genePos = ismember(coFactorInfo.genes,gene);
    if sum(genePos) > 0
       disp(num2str(i))
       % get all unique rxns catalyzed by enzyme
       rxns      = eModel.rxns(~cellfun(@isempty,strfind(eModel.grRules,gene)));
       rxns      = strrep(rxns,'_REV','');
       rxns      = unique(rxns); 
       % get (unique) compartment of rxns
       rxnComps = cell(size(rxns));
       for k = 1:length(rxns)
           rxn  = rxns{k};
           if ~contains(rxn,'holo_prot_')
               rxnComp = getRxnComp(eModel,rxn);
               % check if rxn comp is already in rxnComps
               compareComps = cellfun(@(x) isequal(x,rxnComp),rxnComps,...
                                     'UniformOutput',false);
               isPresent    = cellfun(@(x) any(x),compareComps);
               if sum(isPresent) == 0
                   rxnComps{k}  = rxnComp;
               end
           end
       end
       rxnComps = rxnComps(~cellfun(@isempty,rxnComps));
       % Check which compartment(s) modifications should be performed in
       modificationComps = cell(size(rxnComps));
       mitoComps = {'m','mm','ims','om'};
       for l = 1:length(rxnComps)
           comps = rxnComps{l};
           if length(comps) == 1
               comp = comps{1};
               if sum(strcmp(mitoComps,comp)) > 0
                   % Check the subcompartment localization
                   geneInd = strcmp(mitoEnzymes.genes,gene);
                   loc = lower(mitoEnzymes.localization{geneInd});
                   if strcmp(loc,'im')
                       % change compartment symbol to match the model
                       loc = 'mm';
                   end
               elseif strcmp(comp,'c')
                   % check if enzyme is localized to outer membrane
                   mitoGene = strcmp(mitoEnzymes.genes,gene);
                   if sum(mitoGene) > 0 && strcmpi(mitoEnzymes.localization{mitoGene},'om')
                       loc = 'om';
                   else
                       loc = comp;
                   end
               else
                   % All enzymes that are not mitochondrial are added to 'c'
                   loc = 'c'; 
               end
           else
               if sum(ismember(mitoComps,comps)) > 0
                   % Check mitochondrial subcompartment localization
                   geneInd = strcmp(mitoEnzymes.genes,gene);
                   loc = lower(mitoEnzymes.localization{geneInd});
                   if strcmp(loc,'im')
                       % change compartment symbol to match the model
                       loc = 'mm';
                   end
               else
                   % All non-mitochondrial enzymes are added to 'c'
                   loc = 'c';
               end
           end
           if ~sum(strcmp(modificationComps,loc)) > 0
               modificationComps{l} = loc;
           end
       end
       modificationComps = modificationComps(~cellfun(@isempty,modificationComps));
       % Add co-factor modifications
       uniprotID = eModel.enzymes{i};
       holoProt = ['holo_prot_' uniprotID];
       apoProt  = ['apo_prot_' uniprotID];
       oldRxnName = ['holo_prot_' uniprotID '_exchange'];
       newRxnName = ['apo_prot_' uniprotID '_exchange'];
       rxnPos = strcmp(eModel.rxnNames,oldRxnName);
       for m = 1:length(modificationComps)
           comp = modificationComps{m};
           % Add metabolites apo_prot_X to model
           if sum(strcmp(mitoComps,comp))
               if sum(strcmp(modificationComps,'c')) == 0
                   % change cytosolic exchange rxns to provide apo_protein_X
                   metPos   = strcmp(eModel.mets,holoProt);
                   % replace metID and metName with apo_prot_x
                   eModel.mets{metPos}     = apoProt; 
                   eModel.metNames{metPos} = apoProt;
                   % Change ID and name of exchange rxn
                   eModel.rxns{rxnPos} = newRxnName;
               end
               % Add metabolites apo_protein_X and holo_protein_X to
               % mitochondrial compartments
               apoProt_mito           = [apoProt '_' comp];
               holoProt_mito          = [holoProt '_' comp];
               metsToAdd.mets         = {apoProt_mito,holoProt_mito};
               metsToAdd.metNames     = {apoProt,holoProt};
               metsToAdd.compartments = comp;
               eModel = addMets(eModel,metsToAdd);
           else
               % Add apo_prot_X met to cytoplasm
               metsToAdd.mets         = {apoProt};
               metsToAdd.metNames     = {apoProt};
               metsToAdd.compartments = 'c';
               eModel = addMets(eModel,metsToAdd);
               % Change exchange rxn to provide the apoenzyme
               eModel.S(:,rxnPos)   = 0;
               eModel.S(end,rxnPos) = 1;
               % Change rxnID and name
               eModel.rxns{rxnPos}     = newRxnName;
               eModel.rxnNames{rxnPos} = newRxnName;
           end
       
           cofactors = coFactorInfo.type(genePos);
           cofactorStoich = coFactorInfo.stoich(genePos);
           for j = 1:length(cofactors)
               cofactor = cofactors{j};
               % convert type to corresponding metName in model
               if strcmp(cofactor,'fe2s2')
                   cofactors{j} = '[2Fe-2S] iron-sulfur cluster';
               elseif strcmp(cofactor,'fe4s4')
                   cofactors{j} = '[4Fe-4S] iron-sulfur cluster';
               elseif strcmp(cofactor,'fe3s4')
                   cofactors{j} = '[3Fe-4S] iron-sulfur cluster';
               elseif strcmp(cofactor,'heme_a')
                   cofactors{j} = 'heme a';
               elseif strcmp(cofactor,'ferroheme_b')
                   cofactors{j} = 'ferroheme b';
               end
           end
           compPos = find(strcmp(eModel.comps,comp));
           clear subs prods
           lb = 0;
           ub = 1000;
           grRule = coFactorInfo.grRules(genePos);
           grRule = grRule(~cellfun(@isempty,grRule));
           if sum(strcmp(cofactors,'lipoate')) > 0
               lipoylationRxnPos = strcmp(templateRxns.rxnAbbreviation,'prot_X_lipoylation (lumped)');
               templateRxnEq = templateRxns.rxnEquation{lipoylationRxnPos};
               grRules_modification = templateRxns.grRules{lipoylationRxnPos};
               
               % parse template reaction equation and extract info
               [metList,stoichCoeffList,compList,~] = parseTemplateRxn(templateRxnEq);
               metIDs = cell(length(metList),1);
               % get metIDs from metList
               for j = 1:length(metList)
                   met = metList{j};
                   compartment = compList{j};
                   if contains(met,'prot_')
                       % update to specific protein
                       metIDs{j} = strrep(metList{j},'X',[uniprotID '_m']);
                   else
                       % find metabolite in specific compartment
                       metID = findCompSpecificMetIDs(eModel,met,compartment);
                       metIDs{j} = metID;
                   end
                end

               rxnToAdd.mets         = metIDs;
               rxnToAdd.stoichCoeffs = str2double(stoichCoeffList);
               rxnToAdd.rxns         = cellstr(['prot_' uniprotID '_modification_' comp]);
               rxnToAdd.rxnNames      = cellstr(['prot_' uniprotID '_modification']);
               rxnToAdd.lb           = lb;
               rxnToAdd.ub           = ub;
               rxnToAdd.grRules      = cellstr(grRules_modification);
            
            elseif sum(strcmp(cofactors,'biotin')) > 0
                biotinPool = coFactorInfo.pool{genePos};
                if strcmp(biotinPool,'mitochondrial')
                    templateRxnPos = strcmp(templateRxns.rxnAbbreviation,'prot_X_biotinylation (mitochondrial)');
                else
                    templateRxnPos = strcmp(templateRxns.rxnAbbreviation,'prot_X_biotinylation (cytosolic)');
                end
                templateRxnEq = templateRxns.rxnEquation{templateRxnPos};
                grRules_modification = templateRxns.grRules{templateRxnPos};
                % parse template rxn and extract info
                [metList,stoichCoeffList,compList,~] = parseTemplateRxn(templateRxnEq);
                metIDs = cell(length(metList),1);
                % get metIDs from metList
                for j = 1:length(metList)
                    met = metList{j};
                    compartment = compList{j};
                    if contains(met,'prot_') && ~sum(strcmp(mitoComps,compartment)) > 0
                        % update to specific protein
                        metIDs{j} = strrep(metList{j},'X',uniprotID);
                    elseif contains(met,'prot_') && sum(strcmp(mitoComps,compartment)) > 0 
                        metIDs{j} = strrep(metList{j},'X',[uniprotID '_m']);
                    else
                        % find metabolite in specific compartment
                        metID = findCompSpecificMetIDs(eModel,met,compartment);
                        metIDs{j} = metID;
                    end
                end
                 
                % Define rxn ID
                if sum(strcmp(mitoComps,comp)) > 0
                    rxn = ['prot_' uniprotID '_modification_' comp];
                else
                    rxn = ['prot_' uniprotID '_modification'];
                end

               rxnToAdd.rxns         = cellstr(rxn);
               rxnToAdd.rxnNames     = cellstr(rxn);
               rxnToAdd.mets         = metIDs;
               rxnToAdd.stoichCoeffs = str2double(stoichCoeffList);
               rxnToAdd.grRules      = cellstr(grRules_modification);
               rxnToAdd.lb           = lb;
               rxnToAdd.ub           = ub;
           
           else
               cofactorPool = unique(coFactorInfo.pool(genePos));
               if strcmp(cofactorPool,'mitochondrial') && ~strcmp(comp,'c')
                   templateRxnPos = strcmp(templateRxns.rxnAbbreviation,'prot_X_cofactor_modification (mitochondrial)');
               elseif strcmp(cofactorPool,'non-mitochondrial') && strcmp(comp,'m')
                   templateRxnPos = strcmp(templateRxns.rxnAbbreviation,'prot_X_cofactor_modification (mitochondrial)');
               else
                   templateRxnPos = strcmp(templateRxns.rxnAbbreviation,'prot_X_cofactor_modification (cytosolic)');
               end
               templateRxnEq = templateRxns.rxnEquation{templateRxnPos};
               % parse template reaction equation and extract info
               [metList,stoichCoeffList,compList,~] = parseTemplateRxn(templateRxnEq);
               metIDs = cell(length(metList),1);
               % get metIDs from metList
               for j = 1:length(metList)
                   met = metList{j};
                   compartment = compList{j};
                   if contains(met,'prot_')
                       % update to specific protein
                       if sum(strcmp(mitoComps,comp)) > 0
                           metIDs{j} = strrep(metList{j},'X',[uniprotID '_' comp]);
                           rxn = ['prot_' uniprotID '_modification_' comp];
                       else
                           metIDs{j} = strrep(metList{j},'X',uniprotID);
                           rxn = ['prot_' uniprotID '_modification'];
                       end
                   else
                       % find metabolite in specific compartment
                       metID = findCompSpecificMetIDs(eModel,met,compartment);
                       metIDs{j} = metID;
                   end
               end
               
               % Update stoichiometric coefficients %%%%%%%%%%%%%%%%%%%%%%%
               stoichCoeffList = str2double(stoichCoeffList);
               for j = 1:length(metList)
                   metabolite = metList{j};
                   metabolitePos = strcmp(cofactors,metabolite);
                   if contains(metabolite,'prot_')
                       continue
                   elseif sum(metabolitePos) == 1
                       stoichCoeffList(j) = -cofactorStoich(metabolitePos);
                   else
                       stoichCoeffList(j) = 0;
                   end
               end
               

               rxnToAdd.rxns         = cellstr(rxn);
               rxnToAdd.rxnNames     = cellstr(rxn);
               rxnToAdd.mets         = metIDs;
               rxnToAdd.stoichCoeffs = stoichCoeffList;
               rxnToAdd.grRules      = grRule;
               rxnToAdd.lb           = lb;
               rxnToAdd.ub           = ub;
               % Add any gene in grRule that is not present in model
               if ~isempty(grRule)
                   geneNames = strsplit(grRule{1}, ' and ');
                   genesToAdd.genes = {};
                   genesToAdd.geneShortNames = {};
                   for j = 1:length(geneNames)
                       gene_name = geneNames{j};
                       if sum(strcmp(eModel.genes,gene_name)) == 0
                           genesToAdd.genes = [genesToAdd.genes;gene_name];
                           geneShort = SGDgeneShortName{strcmp(SGDgeneName,gene_name)};
                           genesToAdd.geneShortNames = [genesToAdd.geneShortNames;geneShort];
                       end
                   end
                   if ~isempty(genesToAdd.genes)
                       eModel = addGenesRaven(eModel,genesToAdd);
                   end
               else
                   rxnToAdd.grRules = {''};
               end
           end
           eModel = addRxns(eModel,rxnToAdd,1);
           % check that rxn can carry flux
           cd ../miscellaneous/
           printRxnFormulaRAVEN(eModel,eModel.rxns(end),true,true,true);
           cd(current)
           testModificationRxn(eModel,eModel.enzymes{strcmp(enzGenes,gene)},rxn);
       end
    end
end

end
