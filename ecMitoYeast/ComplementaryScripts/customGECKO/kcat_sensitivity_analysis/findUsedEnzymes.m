%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usedEnzymeIndxs = findUsedEnzymes(model,base_sol,fluxCarriers)
%
% INPUT:  model        - enzyme constrained model (1x1 struct)
%         base_sol     - solution from FBA (1x1 struct)
%         fluxCarriers - original rxn indexes of the flux carrier enzymes 
%                       (enzyme usages)
%
% OUTPUT: usedEnzymeIndexes - list of metabolite indexes for used enzymes
%
% Carl Malina       Last edited 2020-04-24.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function usedEnzymesIndxs = findUsedEnzymes(model,base_sol,fluxCarriers)
usedEnzymesIndxs = [];
for i = 1:length(fluxCarriers)
    enzUsageRxn = model.rxns{fluxCarriers(i)};
    rxnSplit = strsplit(enzUsageRxn,'_');
    uniprot = rxnSplit{3};
    enzMets = model.mets(~cellfun(@isempty,strfind(model.mets,['holo_prot_' uniprot])));
    for j = 1:length(enzMets)
        met = enzMets{j};
        metIndx = strcmp(model.mets,met);
        enzRxnsIndxs = find(model.S(metIndx,:));
        enzRxns = model.rxns(enzRxnsIndxs);
        enzUsed = false;
        for k = 1:length(enzRxns)
            rxnIndx = enzRxnsIndxs(k);
            rxn = enzRxns{k};
            if contains(rxn,'r_') && base_sol.x(rxnIndx) > 0
                enzUsed = true;
            end 
        end
        if enzUsed
            usedEnzymesIndxs = [usedEnzymesIndxs;find(metIndx)];
        end
    end
end
end
