%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainPool(model,non_measured,UB)
%
% Function adapted from the GECKO toolbox:
% https://github.com/SysBioChalmers/GECKO
%
% Carl Malina. Last edited: 2020-01-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = constrainPool(model,non_measured,UB)

%Find default compartment:
cytIndex = strcmpi(model.compNames,'cytoplasm');
if sum(cytIndex) == 1
    comp = model.comps{cytIndex};	%Protein pool in cytosol
else
    comp = model.comps{1};
end

%Define new rxns: For each enzyme, add a new rxn that draws enzyme from the
%enzyme pool (a new metabolite), and remove previous exchange rxn. The new
%rxns have the following stoichiometry (T is the enzyme pool):
% MW[i]*P[T] -> P[i]
for i = 1:length(model.enzymes)
    if non_measured(i)
        % Find the name of the exchange reaction for enzyme i and enzyme name 
        exc_rxn_pos = ~cellfun(@isempty,strfind(model.rxns,['prot_' model.enzymes{i} '_exchange']));
        enz_name    = model.mets{find(model.S(:,exc_rxn_pos))};
        rxnToAdd.rxns         = {['draw_prot_' model.enzymes{i}]};
        rxnToAdd.rxnNames     = rxnToAdd.rxns;
        rxnToAdd.mets         = {'prot_pool' enz_name};
        rxnToAdd.stoichCoeffs = [-model.MWs(i) 1];
        rxnToAdd.lb           = 0; % ub is taken from model default, otherwise inf
        rxnToAdd.grRules      = model.enzGenes(i);
        model = addRxns(model,rxnToAdd,1,comp,true);
        model = removeReactions(model,model.rxns(exc_rxn_pos));
    end
end

%Finally, constraint enzyme pool by fixed value:
rxnToAdd.rxns         = {'prot_pool_exchange'};
rxnToAdd.rxnNames     = rxnToAdd.rxns;
rxnToAdd.mets         = {'prot_pool'};
rxnToAdd.stoichCoeffs = 1;
rxnToAdd.lb           = 0;
rxnToAdd.ub           = UB;
rxnToAdd.grRules      = {''};
model = addRxns(model,rxnToAdd);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%