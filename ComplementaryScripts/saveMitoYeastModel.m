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
model.modelName   = 'yeastGEM with expanded mitochondrial compartments';
model.modelID     = ['mitoYeastGEM_v' version];

if isfield(model,'grRules')
    model = rmfield(model,'grRules');
end

% Update SBO terms in model
cd ../../yeast-GEM/ComplementaryScripts/missingFields/
model = addSBOterms(model);
cd ../../../mitoYeast-GEM/ComplementaryScripts/

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
writeCbModel(model,'text','../ModelFiles/txt/mitoYeastGEM.txt');
% .mat
save('mitoYeastGEM.mat','model');

end