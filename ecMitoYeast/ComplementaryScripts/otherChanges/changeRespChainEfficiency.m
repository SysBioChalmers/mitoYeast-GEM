function model = changeRespChainEfficiency(model)
% changeRespChainEfficiency
% Function that changes the efficiency of respiratory chain to the
% theoretical P/O of 1.5
%
% Usage: model = changeRespEfficiency(model)
%
% Carl Malina   2020-09-18
%

% Get position of relevant metabolites and reactions
PMFpos_ims            = strcmp(model.mets,'s_4324');
PMFpos_m              = strcmp(model.mets,'s_4323');
complexIIIpos         = contains(model.rxns,'r_0439No');
arm_rxn_complexIIIpos = contains(model.rxns,'arm_r_0439');
complexIVpos          = contains(model.rxns,'r_0438No');
arm_rxn_complexIVpos  = contains(model.rxns,'arm_r_0438');
% Change PMF stoichiometry of reactions
scalingFactor = 1.5/1.266;
model.S(PMFpos_m,arm_rxn_complexIIIpos) = model.S(PMFpos_m,arm_rxn_complexIIIpos)*scalingFactor; % 0.1% electron leakge to O2 included
model.S(PMFpos_ims,complexIIIpos)       = model.S(PMFpos_ims,complexIIIpos)*scalingFactor;
model.S(PMFpos_m,arm_rxn_complexIVpos)  = model.S(PMFpos_m,arm_rxn_complexIVpos)*scalingFactor;
model.S(PMFpos_ims,complexIVpos)        = model.S(PMFpos_ims,complexIVpos)*scalingFactor;

% Add pseudometabolite used for coupling efficiency
% get new met ID
newID                  = getNewInd(model.mets,'met');
metID                  = ['s_' newID];
metsToAdd.mets         = {metID};
metsToAdd.metNames     = {'coupling'};
metsToAdd.compartments = {'m'};
model                  = addMets(model,metsToAdd,false);
% Couple r_0438 (complex IV) to production of coupling pseudometabolite
metPos                       = strcmp(model.mets,metID);
model.S(metPos,complexIVpos) = +1;
% Add reaction for empirical coupling efficiency to adjust P/O ratio
% Flux through r_0438 is 2x flux through r_0439 -> stoichiometry 2/3 in
% coupling reaction leads to same flux as through respiratory chain to
% balance coupling pseudometabolite
newID                  = getNewInd(model.rxns,'rxn');
rxnID                  = ['r_' newID];
rxnsToAdd.rxns         = {rxnID};
rxnsToAdd.mets         = {metID,'s_4324','s_4323'};
rxnsToAdd.stoichCoeffs = [-(2/3);-0.12;0.12];
rxnsToAdd.rxnNames     = {'Empirical coupling efficiency'};
rxnsToAdd.lb           = 0;
rxnsToAdd.ub           = 1000;
model                  = addRxns(model,rxnsToAdd,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newID = getNewInd(IDs,type)
% get index type
if strcmp(type,'met')
    % extract metabolites that are not proteins
    metIDs = IDs(contains(IDs,'s_'));
    metIDs = regexprep(metIDs,'[^(\d*)]','');
    metIDs = str2double(metIDs);
    newID  = max(metIDs)+1;
    
elseif strcmp(type,'rxn')
    % extract reactions that are not related to protein translocation or
    % modification
    rxnIDs = IDs(startsWith(IDs,'r_'));
    rxnIDs = regexprep(rxnIDs,'No(\d{1})','');
    rxnIDs = rxnIDs(~endsWith(rxnIDs,'REV'));
    rxnIDs = regexprep(rxnIDs,'[^(\d*)]','');
    rxnIDs = str2double(rxnIDs);
    newID  = max(rxnIDs)+1;
else
    error('type must be met or rxn') 
end

newID = num2str(newID);

end
