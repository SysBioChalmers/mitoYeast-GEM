%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addCofactorsToBiomass(model)
%
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addCofactorsToBiomass(model)

% Load data on biomass components
cd ../../ComplementaryData/Physiology/
fid = fopen('biomassComponents.tsv');
bd  = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
biomassData.mets   = bd{1};
biomassData.MWs    = cellfun(@str2num,bd{4});
biomassData.groups = bd{5};
fclose(fid);

% calculate total mitochondrial protein amount (used for estimating
% ubiquinone (CoQ) content
cd ../../ComplementaryScripts/otherChanges/
mitoProtFrac  = calculateMitoProteinFraction(model); % [g mitochondrial protein/gDW]
CoQ_amount    = 5.4e-3; % [mmol/g mitochondrial protein] taken from DOI:10.1016/B978-0-12-510150-9.50010-9
CoQ_abundance = CoQ_amount * mitoProtFrac; %[mmol/gDW]
CoQpos        = strcmp(model.metNames,'ubiquinol-6 [mitochondrion]');
CoQ_MW        = 592.905; % [g/mol]


% add ubiquinone to cofactors pseudoreaction

cofactorRxnPos = strcmp(model.rxnNames,'cofactor pseudoreaction');
model.S(CoQpos,cofactorRxnPos)    = -CoQ_abundance;
% Remove heme a from cofactor pseudoreaction, added as prosthetic group
% instead
hemeApos = strcmp(model.metNames,'heme a [cytoplasm]');
model.S(hemeApos,cofactorRxnPos) = 0;

% balance out biomass after adding additional cofactors

% add info for CoQ and biotin to biomassData
biomassData.mets   = [biomassData.mets;model.mets{CoQpos}];
biomassData.MWs    = [biomassData.MWs;CoQ_MW;];
biomassData.groups = [biomassData.groups;{'cofactor'}];
%cd ../../../yeast-GEM/ComplementaryScripts/modelCuration/
[X,~,C,~,~,~,~,F] = sumBioMass(model,biomassData);

%cd ../../../mitoYeast-GEM/ComplementaryScripts/otherChanges/

delta = X - 1;
fC = (C -delta)/C;
model = rescalePseudoReaction(model,'carbohydrate',fC);
%cd ../../../yeast-GEM/ComplementaryScripts/modelCuration/
sumBioMass(model,biomassData);

cd ../

end