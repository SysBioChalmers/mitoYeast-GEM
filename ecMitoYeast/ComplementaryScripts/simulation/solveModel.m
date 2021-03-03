%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solveModel
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu,sol] = solveModel(model,mu_low,mu_high)

% unconstrain glucose uptake rate and set minimizing glucose uptake rate as objective
   glcUptakePos = ~cellfun(@isempty,strfind(model.rxnNames,'D-glucose exchange (reversible)'));
   model.ub(glcUptakePos) = 1000;
   model.c = zeros(size(model.rxns));
   model.c(glcUptakePos) = 1;
   growth_pos = strcmpi(model.rxnNames,'growth');

while mu_high - mu_low > 0.0001
   mu_mid = (mu_high + mu_low)/2;
   % set growth rate to mu_mid
   model.lb(growth_pos) = mu_mid;
   model.ub(growth_pos) = mu_mid;
   
   % incorporate growth rate in constraints
   modelWithMuConstraint = incorporateMuInConstraints(model,mu_mid);
   
   % try solving the model with minimizing glucose uptake as objective
   sol = optimizeCbModel(modelWithMuConstraint,'min','one');
  
   if sol.stat == 1
       mu_low     = mu_mid;
       %fluxes_tmp = sol.x;
       sol_tmp = sol;
   else
       mu_high = mu_mid;
   end
   
end

mu = mu_low;
%fluxes = fluxes_tmp;
if exist('sol_tmp','var')
    sol = sol_tmp;
end

% get index of protein pool exchange
protPoolPos = strcmp(model.rxns,'prot_pool_exchange');
% fix glucose uptake rate and minimize enzyme usage
modelWithMuConstraint.ub(glcUptakePos) = sol.f*1.001;
modelWithMuConstraint.lb(glcUptakePos) = sol.f*1.001;
modelWithMuConstraint.c                = zeros(size(model.rxns));
modelWithMuConstraint.c(protPoolPos)   = 1;
   
sol = optimizeCbModel(modelWithMuConstraint,'min','one');
if sol.stat ~= 1
    sol = sol_tmp;
end
end