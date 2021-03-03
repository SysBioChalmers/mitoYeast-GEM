%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = incorporateMuInConstraints(model,mu)
% This functions incorporates the growth rate in the stoichiometric
% coefficients of enzyme usages (in metabolic rxns) and the enzyme exchange
% or draw_prot_X rxns in the case of no protein abundance data provided
%
% INPUT:
% model        an enzyme-constrained model (1x1 struct)
% mu           growth rate
%
% OUTPUT:
% model        a model where the growth rate is accounted for in
%              constraints
%
% Carl Malina. Last edited: 2020-01-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = incorporateMuInConstraints(model,mu)

% update all enzyme-involving reactions
for i = 1:length(model.rxns)
    rxn = model.rxns{i};
    subsPos = find(model.S(:,i) < 0);
    
    % check type of reaction and update coefficient accordingly
    if contains(rxn,'draw_prot_')
        % MW prot_pool --> prot_X. Only one substrate
        model.S(subsPos,i) = model.S(subsPos,i) / mu;
    elseif contains(rxn,'prot_') && contains(rxn,'exchange') && ~contains(rxn,'_pool_')
        % protein exchange rxn in case of proteomics data provided
        model.ub(i) = model.ub(i)/mu;
    elseif startsWith(rxn,'r_') && contains(rxn,'No')
        % metabolic reaction, update coefficient of enzymes
        subsMets = model.mets(subsPos);
        for j = 1:length(subsMets)
            met = subsMets{j};
            if contains(met,'prot_')
                model.S(subsPos(j),i) = model.S(subsPos(j),i) * mu;
            end
        end
    else
        continue
    end
    
%    enzyme  = model.enzymes{i};
%    enzGene = model.enzGenes{i};
%    genePos = strcmp(model.genes,enzGene);
%    % Get list of rxns involving the enzyme
%    rxnPos  = find(model.rxnGeneMat(:,genePos));
%    rxnList = model.rxns(rxnPos);
%    if sum(strcmp(model.mets,'prot_pool')) == 1
%       prot_pool_pos = strcmp(model.mets,'prot_pool'); 
%    end
   % go through list of reactions and update rxns involving enzymes
%    for j = 1:length(rxnList)
%       rxn     = rxnList{j};
%       rxn_pos = rxnPos(j);
%       if contains(rxn,'draw_prot')
%           % update stoichiometric coefficient
%           model.S(prot_pool_pos,rxn_pos) = model.S(prot_pool_pos,rxn_pos) / mu;
%       elseif contains(rxn,'exchange')
%           % update upper bound of reaction
%           model.ub(rxn_pos) = model.ub(rxn_pos)/mu;
%       elseif startsWith(rxn,'r_')
%          subs_pos = find(model.S(:,rxn_pos) < 0);
%          subs     = model.mets(subs_pos);
%          for k = 1:length(subs)
%              sub = subs{k};
%              if contains(sub,enzyme)
%                 subPos = subs_pos(k);
%                 % Update stoichiometric coefficient of enzyme in rxn
%                 model.S(subPos,rxn_pos) = model.S(subPos,rxn_pos) * mu;
%              end
%          end  
%       end
%    end
   % Check if model is feasible
   %sol = optimizeCbModel(model);
   %disp(num2str(sol.obj));
   %if ~strcmpi(sol.origStat,'optimal')
    %   disp(['Problem after constraining enzyme ' enzyme ' ' model.enzNames{i}])
     %  error('Constraining enzyme renders model infeasible')
   %end
end

end