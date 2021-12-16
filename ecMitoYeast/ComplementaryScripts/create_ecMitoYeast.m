%% Create ecMitoYeast

% 1. Import the metabolic model
cd ../../mitoYeast/ModelFiles/mat
initCobraToolbox
changeCobraSolver('ibm_cplex');
model = load('mitoYeastGEM.mat');
model = model.model;
model = creategrRulesField(model);
modelID = strsplit(model.modelID,'_');
modelName = modelID{1};
modelVersion = modelID{2};
model_raven = ravenCobraWrapper(model);
cd ../../../../GECKO/geckomat/change_model/
[model,~,~] = preprocessModel(model_raven,modelName,modelVersion); 

% Collect kcats and assign to reactions, using GECKO toolbox
cd ../../../mitoYeast-GEM/ecMitoYeast/ComplementaryScripts/customGECKO/get_enzyme_data/
[model_data,kcats] = getKcatsFromGecko(model);

% Integrate enzymes in the model:
cd customGECKO/change_model/
ecModel = readKcatDataGECKO(model_data,kcats);
[ecModel,modifications] = manualModifications_modified(ecModel);

% Add enzymes to cofactor modification reactions
cd ../../otherChanges/
ecModel = addEnzymesToModificationRxns(ecModel);

% Adjust co-factor content after accounting for co-factor modfications in
% model
cd ../../otherChanges/
ecModel = takeOutCofactorsAccountedForInModel(ecModel,ecModel.enzGenes);

% Perform additional manual modifications to iron-sulfur cluster biogenesis
cd ../otherChanges/
ecModel = modificationsISCbiosynthesis(ecModel);
% Add translocation complexes to model
ecModel = addTranslocationComplexesToModel(ecModel);
% Perform additional manual modifications protein import
ecModel = modificationsProteinImport(ecModel);
% limit metabolite secretion
cd ../simulation/
ecModel = limitMetaboliteSecretion(ecModel);
% Correct respiratory chain efficiency
cd ../otherChanges/
ecModel_corr_resp_efficiency = changeRespChainEfficiency(ecModel);
% Get constrained model
cd ../customGECKO/limit_proteins/
[ecMitoYeastGEM_batch_1,optSigma] = getConstrainedModel_modified(ecModel_corr_resp_efficiency,modifications,'ecMitoYeast_1');
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
