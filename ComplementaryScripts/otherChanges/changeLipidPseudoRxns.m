%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeLipidPseudoRxns(model)
%
% Carl Malina. Last updated: 2020-08-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = changeLipidPseudoRxns(model)
% save a copy of original model
model_old = model;

% read data from Ejsing et al. 2009
cd ../../../SLIMEr/data
data = readEjsingData(1);
cd ../../mitoYeast-GEM/ComplementaryScripts/otherChanges/
data = convertEjsingData(data,model,true);
lipidData = data.lipidData;
chainData = data.chainData;
%read data from Lahtvee et al. 2017 - to be used for ergosterol ester (SE) 
%and free fatty acids (FFA)
cd ../../../SLIMEr/data/
LahtveeData = readLahtveeData(1);
lipidDataLahtvee = LahtveeData.lipidData;
chainDataLahtvee = LahtveeData.chainData;
% Add abundance data for SE and FFA to lipidData
SE_pos = strcmp(lipidDataLahtvee.metAbbrev,'SE');
FFA_pos = strcmp(lipidDataLahtvee.metAbbrev,'FFA');
% % Calculate the percentage of total lipids constituted by SE and FFA, to be
% % used for calculating the acyl chain abundance in these two species
fractionSE            = lipidDataLahtvee.abundance(SE_pos)/sum(lipidDataLahtvee.abundance);
fractionFFA           = lipidDataLahtvee.abundance(FFA_pos)/sum(lipidDataLahtvee.abundance);
% Combine data
combinedPos           = [find(SE_pos);find(FFA_pos)];
lipidData.metNames    = [lipidData.metNames;lipidDataLahtvee.metNames(combinedPos)];
% % Calculate a scaling factor to account for the 8% total lipid content from
% % Lahtvee study
scalingFactor1 = (0.08-0.08*(fractionSE+fractionFFA))/sum(lipidData.abundance);
scalingFactor2 = 0.08/sum(lipidDataLahtvee.abundance);
lipidData.abundance   = double([lipidData.abundance*scalingFactor1;lipidDataLahtvee.abundance(combinedPos)*scalingFactor2]);
lipidData.std         = [lipidData.std;0;0];
acylChainAbundanceFFA = (chainDataLahtvee.abundance+chainDataLahtvee.std)*fractionFFA;
acylChainAbundanceSE  = zeros(size(chainDataLahtvee.abundance));
palmitoleate_pos      = strcmp(chainDataLahtvee.metNames,'C16:1 chain');
oleate_pos            = strcmp(chainDataLahtvee.metNames,'C18:1 chain');
acylChainAbundanceSE(palmitoleate_pos) = (chainDataLahtvee.abundance(palmitoleate_pos)+chainDataLahtvee.std(palmitoleate_pos))*fractionSE;
acylChainAbundanceSE(oleate_pos)       = (chainDataLahtvee.abundance(oleate_pos)+chainDataLahtvee.std(oleate_pos))*fractionSE;
combinedFraction = fractionSE + fractionFFA;
chainData.abundance = double(chainData.abundance + acylChainAbundanceFFA + acylChainAbundanceSE);%(chainDataLahtvee.abundance+chainDataLahtvee.std)*combinedFraction);

% incorporate changes into data
data.lipidData = lipidData;
data.chainData.abundance = chainData.abundance;

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
cd ../../ComplementaryScripts/otherChanges/

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

% Check if model can grow given the current lipid composition
sol = optimizeCbModel(model);
% if model doesn't predict growth, try scaling acyl chain pseudoreaction
if sol.f == 0
    scaleAll = false;
    lipidChainPseudoRxnPos = strcmp(model.rxnNames,'lipid chain pseudoreaction');
    acylChainPos = find(model.S(:,lipidChainPseudoRxnPos)<0);
    increment = 1.5:-0.01:0.5;
    scalingModel = model;
    optScalingFactor = 1;
    for j = 1:length(increment)
        % scale lipid chain pseudoreaction
        scalingFactor = increment(j);
        scalingModel.S(acylChainPos,lipidChainPseudoRxnPos) = model.S(acylChainPos,lipidChainPseudoRxnPos)*scalingFactor;
        %printRxnFormula(testModel,'r_4063',true,true,true);
        sol = optimizeCbModel(scalingModel);
        if sol.f > 0.01
            disp(['Scaling acyl chains by factor ' num2str(scalingFactor) ' restores growth to: ' num2str(sol.f)]);
            try
                scalingModel = scaleAbundancesInModel(scalingModel,data,'tails');
                if abs(1-scalingFactor) < optScalingFactor
                    optScalingFactor = scalingFactor;
                end
            catch
                warning('Scaling does not work')
                continue
            end
        end
    end
    % If scaling acyl chains restores growth, scale acyl chains. Otherwise, try
    % scaling individual acyl chains
    if optScalingFactor ~= 1
        scaleAll = true;
    else
        increment = 0.95:0.01:1.05;
        optScalingFactor = 1;
        for j = 1:length(increment)
            % scale individual lipid chains in lipid chain pseudoreaction
            scalingFactor = increment(j);
            for z = 1:length(acylChainPos)
                scalingModel = model;
                acylChain_pos = acylChainPos(z);
                scalingModel.S(acylChain_pos,lipidChainPseudoRxnPos) = model.S(acylChain_pos,lipidChainPseudoRxnPos)*scalingFactor;
                %printRxnFormula(testModel,'r_4063',true,true,true);
                sol = optimizeCbModel(scalingModel);
                if sol.f > 0.01
                    %disp(['Scaling acyl chains by factor ' num2str(scalingFactor) ' restores growth to: ' num2str(sol.f)]);
                    disp(['Scaling acyl chain ' scalingModel.metNames{acylChain_pos} ' by factor: ' num2str(scalingFactor) ' restores growth to ' num2str(sol.f)])
                    try
                        scalingModel = scaleAbundancesInModel(scalingModel,data,'tails');
                        if optScalingFactor == 1
                            optScalingFactor = scalingFactor;
                        elseif abs(1-scalingFactor) < abs(1-optScalingFactor)
                            optScalingFactor = scalingFactor;
                            acylChainToScalePos = acylChain_pos;
                            solution = optimizeCbModel(scalingModel);
                            disp(['After scaling, model grows at: ' num2str(solution.f)])
                        end
                        %continue
                    catch
                        warning('Scaling does not work')
                        continue
                    end
                end
            end
        end
    end
end
if scaleAll
    model.S(acylChainPos,lipidChainPseudoRxn) = model.S(acylChainPos,lipidChainPseudoRxn)*optScalingFactor;
else
    model.S(acylChainToScalePos,lipidChainPseudoRxnPos) = model.S(acylChainToScalePos,lipidChainPseudoRxnPos)*optScalingFactor;
end

% Scale abundances of backbones and chains to be consistent
[model,k] = scaleAbundancesInModel(model,data,'tails');

% correct the other components of the biomass again, accounting for the new
% lipid content
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

