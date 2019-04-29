%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeLipidPseudoRxns(model)
%
% Carl Malina 2019-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = changeLipidPseudoRxns(model)
% save a copy of original model
model_old = model;

% read data from Ejsing et al. 2009
cd ../../../SLIMEr/data
data = readEjsingData(1);
data = convertEjsingData(data,model,true);
lipidData = data.lipidData;
chainData = data.chainData;

% Find the IDs in model for all lipid backbones
backboneIDs = cell(size(lipidData.metNames));
for i = 1:length(backboneIDs)
   lipid = lipidData.metNames{i};
   if ~strcmp(lipid,'ergosterol')
       backboneName   = [lipid ' backbone [cytoplasm]'];
       backbonePos    = strcmp(model.metNames,backboneName);
       backboneIDs(i) = model.mets(backbonePos);
   else
       backboneName   = [lipid ' [cytoplasm]'];
       backbonePos    = strcmp(model.metNames,backboneName);
       backboneIDs(i) = model.mets(backbonePos);
   end
end

% Get the fields for the lipid backbone pseudoreaction 
rxnName     = 'lipid backbone pseudoreaction';
rxnID       = model.rxns{strcmp(model.rxnNames,rxnName)};
backbonePos = strcmp(model.metNames,'lipid backbone [cytoplasm]');
backboneID  = model.mets(backbonePos);

% Modify the existing backbone pseudoreaction according to the data
model = addReaction(model, rxnID, ...
                    'reactionName', rxnName, ...
                    'metaboliteList', [backboneIDs;backboneID], ...
                    'stoichCoeffList', [-lipidData.abundance;1], ...
                    'reversible', false, ...
                    'lowerBound', 0, ...
                    'upperBound', 1000);
                     
printRxnFormula(model,rxnID,true,true,true);
model.rxnConfidenceScores(strcmp(model.rxns,rxnID)) = 1;
% remove grRules field automatically created by COBRA function
% addReaction()
if isfield(model,'grRules')
    model = rmfield(model,'grRules');
end

% Find IDs in model for all lipid chains
chainIDs = cell(size(chainData.metNames));
for i = 1:length(chainIDs)
   chainName   = [chainData.metNames{i} ' [cytoplasm]'];
   chainPos    = strcmp(model.metNames,chainName);
   chainIDs(i) = model.mets(chainPos);
end

% Get the fields for the lipid chain pseudoreaction
rxnName  = 'lipid chain pseudoreaction';
rxnID    = model.rxns{strcmp(model.rxnNames,rxnName)};
chainPos = strcmp(model.metNames,'lipid chain [cytoplasm]');
chainID  = model.mets(chainPos);

% Modify the lipid chain pseudoreaction according to the data
model = addReaction(model, rxnID, ...
                    'reactionName', rxnName, ...
                    'metaboliteList', [chainIDs;chainID], ...
                    'stoichCoeffList', [-chainData.abundance;1], ...
                    'reversible', false, ...
                    'lowerBound', 0, ...
                    'upperBound', 1000);
           
printRxnFormula(model,rxnID,true,true,true);
model.rxnConfidenceScores(strcmp(model.rxns,rxnID)) = 1;

% remove grRules field automatically created by COBRA function
% addReaction()
if isfield(model,'grRules')
    model = rmfield(model,'grRules');
end

% Read file with biomass composition data
cd ../../mitoYeast-GEM/ComplementaryData/Physiology/

fid = fopen('biomassComponents.tsv');
bd  = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
biomassData.mets   = bd{1};
biomassData.MWs    = cellfun(@str2num,bd{4});
biomassData.groups = bd{5};
fclose(fid);

% calculate the current composition
cd ../../../yeast-GEM/ComplementaryScripts/modelCuration/

[X,P,C,R,~,~,~,~] = sumBioMass(model,biomassData);

% Correct biomass composition with data
otherData.mets       = data.otherData.metIDs;
otherData.abundances = data.otherData.abundance;
for i = 1:length(otherData.mets)
   met = otherData.mets{i};
   abundance = otherData.abundances(i);
   if strcmp(met,'protein')
       fP = abundance/P;
       model = rescalePseudoReaction(model,'protein',fP);
   elseif strcmp(met,'RNA')
       fR = abundance/R;
       model = rescalePseudoReaction(model,'RNA',fR);
   end
end

[X,P,C,R,~,~,~,~] = sumBioMass(model,biomassData);

% Balance out mass with carbohydrate content
delta = X - 1;
fC = (C - delta)/C;
model = rescalePseudoReaction(model,'carbohydrate',fC);
sumBioMass(model,biomassData);

% Unblock exchange of lipid backbones and chains
model.ub(strcmp(model.rxnNames,'lipid chain exchange')) = 1000;
model.ub(strcmp(model.rxnNames,'lipid backbone exchange')) = 1000;

% Scale abundances of backbones and chains to be consistent
cd ../../../SLIMEr/models
[model,k] = scaleAbundancesInModel(model,data,'tails');

% correct the other components of the biomass again, accounting for the new
% lipid content
cd ../../yeast-GEM/ComplementaryScripts/modelCuration/
[X,~,C,~,~,L,~,~] = sumBioMass(model,biomassData);
% balance out content with carbohydrates
delta = X - 1;
fC = (C - delta)/C;
model = rescalePseudoReaction(model,'carbohydrate',fC);
sumBioMass(model,biomassData);

%Block exchange of tails and backbones:
 posT  = strcmp(model.rxnNames,'lipid chain exchange');
 posB  = strcmp(model.rxnNames,'lipid backbone exchange');
 model = changeRxnBounds(model,model.rxns(posT),0,'b');
 model = changeRxnBounds(model,model.rxns(posB),0,'b');

% Compare model simulations after modifying the lipid pseudoreactions
sol_1        = optimizeCbModel(model_old);
growthRate_1 = sol_1.x(strcmp(model_old.rxnNames,'growth'));
qs_1         = sol_1.x(strcmp(model_old.rxnNames,'D-glucose exchange'))/1000*180;
disp(['Biomass yield for original model: ' num2str(-growthRate_1/qs_1) ' gDW/g(glucose)'])
sol_2        = optimizeCbModel(model);
growthRate_2 = sol_2.x(strcmp(model.rxnNames,'growth'));
qs_2         = sol_2.x(strcmp(model.rxnNames,'D-glucose exchange'))/1000*180;
disp(['Biomass yield for model with altered lipid pseudoreactions: ' num2str(-growthRate_2/qs_2) ' gDW/g(glucose)'])

cd ../../../mitoYeast-GEM/ComplementaryScripts/
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function adapted from yeastGEM: https://github.com/SysBioChalmers/yeast-GEM

function model = rescalePseudoReaction(model,metName,f)

rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
for i = 1:length(model.mets)
    S_ir   = model.S(i,rxnPos);
    isProd = strcmp(model.metNames{i},[metName ' [cytoplasm]']);
    if S_ir ~= 0 && ~isProd
        model.S(i,rxnPos) = f*S_ir;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

