%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveMitoYeastModel(model)
%
% Saves model as a .xml, .txt and .mat file.
%
% Function adapted from yeast-GEM:https://github.com/SysBioChalmers/yeast-GEM
%
% Carl Malina 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveMitoYeastModel(model,version)
model.S = full(model.S);
% change model ID, description and name
model.description = 'mitoYeastGEM.xml';
model.modelName   = 'yeastGEM with expanded mitochondrial metabolism';
model.modelID     = ['mitoYeastGEM_v' version];
model             = creategrRulesField(model); 
grRules = model.grRules; % used when creating txt file for model

% Remove grRules field from model
model = rmfield(model,'grRules');

% Update SBO terms in model
cd ../../../yeast-GEM/ComplementaryScripts/missingFields/
model = addSBOterms(model);
cd ../../../mitoYeast-GEM/mitoYeast/ComplementaryScripts/

% check if model is a valid SBML structure
writeCbModel(model,'sbml','tempModel.xml');
[~,errors] = TranslateSBML('tempModel.xml');
if ~isempty(errors)
    delete('tempModel.xml');
    error('Model should be a valid SBML structure. Please fix all errors before saving.')
end

% Update .xml, .txt and .mat model files:
% .xml
copyfile('tempModel.xml','../ModelFiles/xml/mitoYeastGEM.xml');
delete('tempModel.xml');
% .txt
cd ../ModelFiles/txt/
rxnFormulas = printRxnFormula(model,model.rxns,false,true,true);
subSystems = cell(size(model.rxns));
for i = 1:length(model.rxns)
    subSystem = model.subSystems{i};
    subSystems{i} = strjoin(subSystem,';');
end 
colNames    = {'rxnID','rxnName','Formula','ECnumber', ...
               'GeneAssociation','LB','UB','subSystem','Notes', ...
               'Reference','ConfidenceScore'};
modelInfo = table(model.rxns,model.rxnNames,rxnFormulas,model.rxnECNumbers, ...
                  grRules,model.lb,model.ub,subSystems,model.rxnNotes, ...
                  model.rxnReferences, model.rxnConfidenceScores, ...
                  'VariableNames',colNames);
writetable(modelInfo,'mitoYeastGEM.txt','FileType','text', ...
           'Delimiter','\t');
cd ../mat/
% .mat
save('mitoYeastGEM.mat','model');
cd ../../ComplementaryScripts/

end