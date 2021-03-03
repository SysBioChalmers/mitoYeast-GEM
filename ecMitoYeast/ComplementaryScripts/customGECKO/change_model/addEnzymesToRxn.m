%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,comp,newRxnName,protGenes)
% Adds new metabolite to the left side of a selected reaction in the model.
% If the reaction does not exist it will create a new one.
%
% INPUT:
% model        the GEM structure (1x1 struct)
% kvalues      kcat values of the enzyme/complex
% rxn          the reaction original ID   (string)
% newMets      name of the new pseudo-metabolite (enzyme)
% comp         compartment of the new pseudo-metabolite (enzyme)
% newRxnName   {ID,name} of the new reaction
% protGenes    Gene ID associated with the provided enzyme
% 
% OUTPUTS:
% model             Modified GEM structure (1x1 struct)
% 
% Function adapted from the GECKO toolbox:
% https://github.com/SysBioChalmers/GECKO
%
% Carl Malina. Last edited: 2020-01-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addEnzymesToRxn(model,kvalues,rxn,newMets,comp,newRxnName,protGenes)

if nargin < 7
    protGenes = '';
end

% Define all parts necessary for the new (or changed) rxn
rxnIndex = strcmp(model.rxns,rxn);
metS     = model.mets(model.S(:,rxnIndex) < 0)';
metP     = model.mets(model.S(:,rxnIndex) > 0)';
coeffsS  = model.S(model.S(:,rxnIndex)<0,rxnIndex)';
coeffsP  = model.S(model.S(:,rxnIndex)>0,rxnIndex)';

% Find enzyme in the correct compartment, if enzyme is mitochondrial, and
% update metabolite (enzyme) ID
mitoEncodedGenes = {'Q0105','Q0045','Q0250','Q0275','Q0085',...
    'Q0080','Q0130'};
for i = 1:length(comp)
    compartment = comp{i};
    mitoGene = false;
    if strcmpi(compartment,'om')
        newMets{i} = [newMets{i} '_om'];
    elseif strcmpi(compartment,'im') || strcmpi(compartment,'mm')
        newMets{i} = [newMets{i} '_mm'];
    elseif strcmpi(compartment,'ims')
        newMets{i} = [newMets{i} '_ims'];
    elseif strcmpi(compartment,'m') %&& sum(strcmp(mitoEncodedGenes,protGenes{i})) == 0
        % mitchondrially synthesized enzymes lack _m in metID
        for j = 1:length(mitoEncodedGenes)
            gene = mitoEncodedGenes{j};
            if contains(protGenes,gene)
                mitoGene = true;
            end
        end
        % Update metID if enzyme not mitochondrially synthesized
        if ~mitoGene  
            newMets{i} = [newMets{i} '_m'];
        end
    end
end

% Include enzymes in reactions
rxnToAdd.mets = [metS,newMets,metP];
rxnToAdd.stoichCoeffs = [coeffsS,-kvalues.^-1,coeffsP];
if ismember(newRxnName{1},model.rxns)
    model = changeRxns(model,newRxnName(1),rxnToAdd);
else
    rxnToAdd.rxns     = newRxnName(1);
    rxnToAdd.rxnNames = newRxnName(2);
    rxnToAdd.lb       = model.lb(rxnIndex);
    rxnToAdd.ub       = model.ub(rxnIndex);
    rxnToAdd.obj      = model.c(rxnIndex);
    if ~isempty(protGenes)
        rxnToAdd.grRules = {protGenes};
    end
    if isfield(model,'subSystems')
        rxnToAdd.subSystems = model.subSystems(rxnIndex);
    end
    model = addRxns(model,rxnToAdd,1);
end
    
end