%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = addTranslocationRxns(eModel, enzymeLoc)
%
%
%
%
%
% Carl Malina 2019-10-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = addTranslocationRxns(eModel,enzymeLoc)
current = pwd;
% read file with info about integral inner membrane proteins
cd ../../ComplementaryData/proteinImport/
fid = fopen('integralMembraneProteins.tsv');
fileData = textscan(fid,'%s%s%s%s%s%s%s%s%s%s','Delimiter','\t',...
                    'HeaderLines',1);
fclose(fid);

% extract data for integral inner membrane proteins
IMproteins.uniprotID             = fileData{1};
IMproteins.totalImportedSequence = cellfun(@str2double,fileData{7});

% check if translocation reactions already added to model
if sum(contains(eModel.rxns,' translocation')) == 0
    % load information about import machinery components
    fid = fopen('importMachineryComponents.tsv');
    fileData = textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fid);

    % Extract gene information
    importComponents.genes = fileData{2};
    importComponents.geneShortNames = fileData{3};

    % Check which genes are already present in model
    inModel = ismember(importComponents.genes,eModel.genes);

    % update importComponents to contain only new genes
    importComponents.genes = importComponents.genes(~inModel);
    importComponents.geneShortNames = importComponents.geneShortNames(~inModel);

    % Add genes to model
    eModel = addGenesRaven(eModel,importComponents);

    % Add gene information to newly added genes, based on SGD
    cd ../../../mitoYeast/ComplementaryData/databases/
    fid = fopen('SGDgeneNames.tsv');
    SGDgeneAnnotation = textscan(fid,'%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fid);

    for i = 1:length(eModel.genes)
        geneIndex = strcmp(SGDgeneAnnotation{1}, eModel.genes{i});
        if sum(geneIndex) == 1 && ~isempty(SGDgeneAnnotation{2}{geneIndex})
            eModel.geneNames{i} = SGDgeneAnnotation{2}{geneIndex};
        else
            eModel.geneNames{i} = eModel.genes{i};
        end
    end
end
cd(current)
% read file with template reactions
filename = '../../ComplementaryData/proteinImport/templateReactionsTranslocation.tsv';
templateRxns = readTemplateRxnFile(filename);

% All presequence related import pathways produce a presequence
% add pseudo-metabolite prot_presequence to model to represent the
% presequnce being cleaved off upon import
preSeq_metID = 'prot_presequence';
if sum(strcmp(eModel.mets,preSeq_metID)) == 0
    metsToAdd.mets 	  = {preSeq_metID};
    metsToAdd.metNames = {preSeq_metID};
    metsToAdd.compartments = 'm';
    eModel = addMets(eModel,metsToAdd);
end

% Add protein import
for i = 1:length(eModel.enzymes)
    enzyme = eModel.enzymes{i};
    % Check if enzyme is mitochondrial and if translocation reaction
    % already present in model
    enzPos = strcmp(enzymeLoc.uniprots,enzyme);
    if sum(enzPos) == 1 && sum(contains(eModel.rxns,['prot_' enzyme ' translocation'])) == 0
        comp = lower(enzymeLoc.localization{enzPos});
        if strcmpi(comp,'im')
            comp = 'mm';
        end
        pathway = enzymeLoc.pathway{enzPos};
        enzGene = eModel.enzGenes{strcmp(eModel.enzymes,enzyme)};
        % initialize/reset preSeq_pathway (for checking whether enzyme 
        % depends on the presequence pathway)
        preSeq_pathway = false;
        
        %%%%%%% Check if enzyme is localized to inner membrane here %%%%%%%%
        if strcmpi(pathway,'tim23-im') || strcmpi(pathway,'tim23-oxa')
            IMenzPos = strcmpi(IMproteins.uniprotID,enzyme);
            totalImportedSequence = IMproteins.totalImportedSequence(IMenzPos);
            % assign enzyme to pathway based on "proline rule"
            if ~isnan(totalImportedSequence)%~isnan(IMproteins.totalImportedSequence(strcmpi(IMproteins.uniprotID,enzyme)))
                pathway = 'TIM23-OXA';
                %totalImportedSequence = IMproteins.totalImportedSequence(strcmpi(IMproteins.uniprotID,enzyme));
            else
                pathway = 'TIM23-IM';
            end
        end
       %%%%% ------ %%%%%%%%%%

        sequence = eModel.sequences{i};
        % Check which import pathway protein goes through
        pathwayPos = strcmp(templateRxns.rxnAbbreviation,pathway);
        templateRxnEq = templateRxns.rxnEquation{pathwayPos};
        grRules_translocation = templateRxns.grRules{pathwayPos};
        % parse template reaction equation and extract info
        [metList,stoichCoeffList,compList,~] = parseTemplateRxn(templateRxnEq);
        metIDs = cell(length(metList),1);
        % get metIDs from metList
        for j = 1:length(metList)
            met = metList{j};
            compartment = compList{j};
            if contains(met,'prot_') && ~contains(met,'presequence')
                % protein compartment will be updated separately
                metIDs{j} = met;
            else
                % find metabolite in specific compartment
                metID = findCompSpecificMetIDs(eModel,met,compartment);
                metIDs{j} = metID;
            end
        end

       	% check if protein is apoprotein or holoprotein and collect
       	% metabolite info
       	if sum(strcmp(eModel.mets,['apo_prot_' enzyme])) >0
           	% enzyme requires cofactor modification, imported as apoprotein
           	% enzyme already added to mitochondrial comp when adding
           	% cofactor modifications
           	metID = ['apo_prot_' enzyme];
           	metID_mito = [metID '_' comp];
       	else
           	% enzyme imported as holoprotein, need to add enzyme to
           	% mitochondrial compartment
           	metID = ['holo_prot_' enzyme];
           	metID_mito = [metID '_' comp];
           	metsToAdd.mets = {metID_mito};
           	metsToAdd.metNames = {metID};
           	metsToAdd.compartments = comp;
           	eModel = addMets(eModel,metsToAdd);
       	end
       	% Add specific protein IDs to template reaction
       	proteinPos = strcmp(metIDs,'prot_X');
       	metIDs(proteinPos) = {metID;metID_mito};

       	if contains(pathway,'TIM23')
            % import depends on preSeq_pathway
            preSeq_pathway = true;
            % reset variables for indicating cleavage site
            OCT1_cleavage = false;
            ICP55_cleavage = false;
            % calculate proton motive force (PMF) cost for presequence import
            % get sequence length from protease cleavage site
            if ~isempty(enzymeLoc.ICP55cleavageSite{enzPos})
                preSeqLength = enzymeLoc.ICP55cleavageSite{enzPos};
                ICP55_cleavage = true;
            elseif ~isempty(enzymeLoc.OCT1cleavageSite{enzPos})
                preSeqLength = enzymeLoc.OCT1cleavageSite{enzPos};
                OCT1_cleavage = true;
            else
                preSeqLength = enzymeLoc.MPPcleavageSite{enzPos};
            end

            % Update grRules for translocation based on peptidase
            % cleavage site OCT1 (YKL134C) and/or ICP55 (YER078C)
            
            if OCT1_cleavage && ICP55_cleavage
            	grRules_translocation = [grRules_translocation ' and YKL134C and YER078C'];
            elseif OCT1_cleavage
            	grRules_translocation = [grRules_translocation ' and YKL134C'];
            elseif ICP55_cleavage 
            	grRules_translocation = [grRules_translocation ' and YER078C'];
            end
            
            % Add enzyme gene to grRule for translocation
            grRules_translocation = [grRules_translocation ' and ' enzGene];
            
            % Calculate the proton motive force (PMF) cost of import
            if strcmpi(pathway,'tim23-pam')
                % entire protein imported, calculate cost based on overall
                % net charge of the protein
                netCharge = calculateNetCharge(sequence);
                % presequence required for ATP cost calculations
                preSequence = sequence(1:preSeqLength);
            else
                %Only presequence imported, calculate cost based on net
                % charge of presequence
                preSequence = sequence(1:preSeqLength);
                netCharge = calculateNetCharge(preSequence);
            end
            
            % membrane poteintial contributes 90 % of total PMF
            PMFcost = netCharge * 0.9;
            
            % Check the compartment/pathway
            if strcmpi(pathway,'tim23-oxa')
                % calculate the ATP cost for import of transmembrane segments
                % 1 binding site for mtHsp70 on average every 25 amino acids
                % therefore assumed ro require 1 ATP/25 amino acids
                ATPcost = totalImportedSequence/25;
          % elseif strcmpi(pathway,'tim23-ims')
            elseif strcmpi(pathway,'tim23-pam')
                % calculate ATP cost of import
                totalImportedSequence = length(sequence(preSeqLength+1:end));
                ATPcost = totalImportedSequence/25;
            elseif strcmpi(pathway,'tim23-bcs1')
                % calculate ATP cost of import
                totalImportedSequence = length(sequence(preSeqLength+1:end));
                ATPcost = totalImportedSequence/25;
            end
            % update stoichiometric coefficients for ATP, H2O, ADP and Pi
            % and PMF pseudo-metabolite
            for j = 1:length(stoichCoeffList)
            	coeff = stoichCoeffList{j};
            	if ischar(coeff)
            		if contains(coeff,'x')
            			coeff = strrep(coeff,'x',num2str(ATPcost));
            		elseif contains(coeff,'y')
                        if PMFcost < 0
                            coeff = strrep(coeff,'y',num2str(-PMFcost));
                        else
            			    coeff = strrep(coeff,'y',num2str(PMFcost));
                        end
            		end
            		stoichCoeffList{j} = str2double(coeff);
            	end
            end
            % convert stoichCoeffList to double
            stoichCoeffList = cell2mat(stoichCoeffList);
        else
            % convert stoichCeoffList to double
            stoichCoeffList = cellfun(@str2double,stoichCoeffList);
        end

        % create structure for adding import rxns
        translocationRxn 	= [metID ' translocation'];
        %if preSeq_pathway
        rxnsToAdd.rxns 		   = {translocationRxn};%,preSeq_degradation_rxn};
        rxnsToAdd.mets 		   = metIDs;%,metIDs_preSeq_Deg};
        rxnsToAdd.stoichCoeffs = {stoichCoeffList};%,stoichCoeffs_preSeq};
        rxnsToAdd.rxnNames 	   = {translocationRxn};%,preSeq_degradation_rxn};
        rxnsToAdd.lb 		   = zeros(length(rxnsToAdd.rxns),1);
        rxnsToAdd.ub 		   = 1000*ones(length(rxnsToAdd.rxns),1);
        rxnsToAdd.grRules 	   = {grRules_translocation};%,grRules_preSeq_Deg};
       
        % Add translocation reaction and pre-sequence degradation reaction 
        %(if applicable) to model
        eModel = addRxns(eModel,rxnsToAdd,1);
        % print reactions added
        cd ../miscellaneous/
        if preSeq_pathway
            printRxnFormulaRAVEN(eModel,eModel.rxns(end-1:end),true);
        else
            printRxnFormulaRAVEN(eModel,eModel.rxns(end),true);
        end
    end
    disp(['Successfully went through ' num2str(i) ' enzymes']);
    cd(current);
end

% Add a sink reaction for prot_presequence representing presequence
% degradation
preSeq_degradation_rxn   = 'prot_presequence_degradation';
if sum(strcmp(eModel.rxns,preSeq_degradation_rxn)) == 0
    preSeq_degRxnPos         = strcmpi(templateRxns.rxnAbbreviation,preSeq_degradation_rxn);
    rxnsToAdd.rxns           = {preSeq_degradation_rxn};
    rxnsToAdd.rxnNames       = {preSeq_degradation_rxn};
    preSeq_degradation_rxnEq = templateRxns.rxnEquation{preSeq_degRxnPos};
    grRules_preSeq_Deg       = templateRxns.grRules{preSeq_degRxnPos};
    [metList_preSeq,stoichCoeffs_preSeq,~,~] = parseTemplateRxn(preSeq_degradation_rxnEq);
    rxnsToAdd.mets 		   = metList_preSeq; %,metIDs_preSeq_Deg};
    rxnsToAdd.stoichCoeffs = cellfun(@str2num,stoichCoeffs_preSeq);
    rxnsToAdd.lb 		   = 0;
    rxnsToAdd.ub 		   = 1000;
    rxnsToAdd.grRules 	   = {grRules_preSeq_Deg};

    eModel = addRxns(eModel,rxnsToAdd,1);
end

end