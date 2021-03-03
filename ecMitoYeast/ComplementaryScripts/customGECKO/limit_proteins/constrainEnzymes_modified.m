%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainEnzymes(model,Ptot,sigma,f,GAM,pIDs,data,gRate,GlucUptake)
% 
% Function adapted from the GECKO toolbox:
% https://github.com/SysBioChalmers/GECKO
%
% Carl Malina. Last edited: 2020-09-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,enzUsages,modifications,GAM,massCoverage] = constrainEnzymes_modified(model,f,GAM,Ptot,pIDs,data,gRate,c_UptakeExp)
cd ..
parameters = getModelParameters_modified;
sigma      = parameters.sigma;
c_source   = parameters.c_source;
cd ../../../../GECKO/geckomat/limit_proteins/
%Compute f if not provided:
if nargin < 2
    [f,~] = measureAbundance(model.enzymes);
else
    if isempty(f)
        [f,~] = measureAbundance(model.enzymes);
    end
end
%Leave GAM empty if not provided (will be fitted later):
if nargin < 3
    GAM = [];
end
% Load Ptot if not provided
if nargin < 4
    Ptot = parameters.Ptot;
end
% Scale biomass composition to match Ptot
cd ../../../mitoYeast-GEM/ecMitoYeast/ComplementaryScripts/otherChanges/
ATPpos        = strcmp(model.mets,findCompSpecificMetIDs(model,'ATP','c'));
biomassRxnPos = strcmpi(model.rxnNames,'biomass pseudoreaction');
cd ../miscellaneous/
[~,Pbase,Cbase,R,D,~,~,~,~] = sumBioMass(model);
if isfield(parameters,'pol_cost')
    cost   = parameters.pol_cost;
    GAMpol = Pbase*cost(1) + Cbase*cost(2) + R*cost(3) + D*cost(4);
end
currentGAM    = -model.S(ATPpos,biomassRxnPos) - GAMpol;
cd ../customGECKO/limit_proteins/
[model,~]     = scaleBioMass_modified(model,Ptot,currentGAM);

%No UB will be changed if no data is available -> pool = all enzymes(FBAwMC)
if nargin < 5
    pIDs = cell(0,1);
    data = zeros(0,1);
end
%Remove zeros or negative values
data = cleanDataset(data);
%Assign concentrations as UBs [mmol/gDW]:
model.concs = nan(size(model.enzymes));      %OBS: min value is zero!!
disp('Matching data to enzymes in model...')
for i = 1:length(model.enzymes)
    match = false;
    for j = 1:length(pIDs)
        if strcmpi(pIDs{j},model.enzymes{i}) && ~match
            model.concs(i) = data(j)*model.MWs(i); %g/gDW
            rxn_name       = ['prot_' model.enzymes{i} '_exchange'];
            pos            = ~cellfun(@isempty,strfind(model.rxns,rxn_name));
            model.ub(pos)  = data(j);
            match          = true;
        end
    end
end
%Count mass of non-measured enzymes:
measured       = ~isnan(model.concs);
concs_measured = model.concs(measured);
Pmeasured      = sum(concs_measured);
%Get protein content in biomass pseudoreaction:
cd ../../miscellaneous/
Pbase = sumProtein(model);
cd ../../../../GECKO/geckomat/limit_proteins/
if Pmeasured > 0
    %Calculate fraction of non measured proteins in model out of remaining mass:
    [fn,~] = measureAbundance(model.enzymes(~measured));
    fm     = Pmeasured/Ptot;
    f      = fn/(1-fm);
    %Discount measured mass from global constrain:
    fs = (Ptot - Pmeasured)/Pbase*f*sigma;
else
    fs = f*sigma;
end
cd ../../../mitoYeast-GEM/ecMitoYeast/ComplementaryScripts/customGECKO/limit_proteins/
%Constrain the rest of enzymes with the pool assumption:
if sum(strcmp(model.rxns,'prot_pool_exchange')) == 0
    model = constrainPool(model,~measured,full(fs*Pbase));
end
if sum(data)==0
   %Modify protein/carb content and GAM:
   model = scaleBioMass_modified(model,Ptot,GAM);
end
%Display some metrics:
disp(['Total protein amount measured = '     num2str(Pmeasured)              ' g/gDW'])
disp(['Total enzymes measured = '            num2str(sum(measured))          ' enzymes'])
disp(['Enzymes in model with 0 g/gDW = '     num2str(sum(concs_measured==0)) ' enzymes'])
disp(['Total protein amount not measured = ' num2str(Ptot - Pmeasured)       ' g/gDW'])
disp(['Total enzymes not measured = '        num2str(sum(~measured))         ' enzymes'])
disp(['Total protein in model = '            num2str(Ptot)                   ' g/gDW'])
enzUsages = [];
if nargin > 7
    [tempModel,enzUsages,modifications] = flexibilizeProteins(model,gRate,c_UptakeExp,c_source);
    Pmeasured = sum(tempModel.concs(~isnan(tempModel.concs)));
    model     = updateProtPool(tempModel,Ptot,f*sigma);
end
massCoverage  = Pmeasured/Ptot;
if isempty(enzUsages)
    enzUsages     = table({},zeros(0,1),'VariableNames',{'prot_IDs' 'usage'});
    modifications = table({},zeros(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'protein_IDs' 'previous_values' 'modified_values' 'flex_mass'});
else
    plotHistogram(enzUsages,'Enzyme usage [-]',[0,1],'Enzyme usages','usages')
end
%Plot histogram (if there are measurements):
plotHistogram(concs_measured,'Protein amount [mg/gDW]',[1e-3,1e3],'Modelled Protein abundances','abundances')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotHistogram(variable,xlabelStr,xlimits,titleStr,option)
if iscell(variable)
    cell2mat(variable);
end
if sum(variable) > 0
    figure
    if strcmpi(option,'abundances')
        hist(variable*1e3,10.^(-3:0.5:3))
        set(gca,'xscale','log')
    else
        hist(variable,(0:0.05:1))
    end
    xlim(xlimits)
    xlabel(xlabelStr)
    ylabel('Frequency');
    title(titleStr)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = cleanDataset(data)
for i=1:length(data)
    if data(i)<=0
        data(i) = NaN;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%