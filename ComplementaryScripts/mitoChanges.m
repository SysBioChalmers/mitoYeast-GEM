%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes modification to the mitochondrial compartment of the
% yeast8 model, including adding reactions, an additional compartment 
% representing the intermembrane space (IMS) and updating the GPR rules for
% mitochondrial reactions
%
% It is dependent on functions from the yeast-GEM project and SLIMEr:
% https://github.com/SysBioChalmers/yeast-GEM
% https://github.com/SysBioChalmers/SLIMEr
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize COBRA toolbox and load model
initCobraToolbox()

cd ../../yeast-GEM/ModelFiles/xml/
model = readCbModel('yeastGEM.xml');
model = buildRxnGeneMat(model);

% Change the lipid pseudoreaction to include additional lipid species
cd ../../../mitoYeast-GEM/ComplementaryScripts/otherChanges/
model = changeLipidPseudoRxns(model);

% Add new mitochondrial reactions
cd modelCuration/
model = addMitochondrialRxns(model);

% Add biotin and ubiquinol to cofactor pseudoreaction
cd otherChanges/
model = addCofactorsToBiomass(model);

% Add mitochondrial aa-tRNAs and prosthetic groups to protein
% pseudoreaction
cd otherChanges/
model = changeProteinPseudoReaction(model);

% correct GPR for mitochondrial reactions
cd modelCuration/
model = updateMitoGPR(model);

% Add mitochondrial intermembrane space (IMS) compartment
cd modelCuration/
model = addIMScomp(model);

% update compartment of metabolites in IMS reactions
cd modelCuration/
model = correctMetLocalization(model);

cd modelCuration/
% Remove incorrect reactions
[model,notFound,growthRate] = removeIncorrectRxns(model);

% Move reactions to correct compartment
% This includes adding the mitochondrial outer membrane compartment
model = moveReactions(model);

cd ../otherChanges/
% Represent the proton motive force (PMF) in model
[model,energy,redox,growthRate] = addPMFcostToRxns(model);

%%
% Save model
saveMitoYeastModel(model,'1.0.0');

