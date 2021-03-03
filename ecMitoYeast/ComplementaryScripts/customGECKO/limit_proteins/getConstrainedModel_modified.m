%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel_batch,OptSigma] = getConstrainedModel_modified(ecModel,c_source,sigma,Ptot,gRate,modifications,name)
%
% Function that gets a GEM with kinetic data and returns an enzyme 
% constrained model, either with individual enzyme levels or with the total
% measured protein content.
% 
% The model parameters are automatically curated (querying the BRENDA
% files) for reaching the experimental maximal growth rate on glucose
% minimal media (if it was overconstrained). Then, the average saturation 
% factor for enzymes is fitted for the same growth conditions. Finally the 
% fitted model is simulated on the same conditions and the top ten used 
% enzymes (mass-wise) are saved in a file as an output, this can be used for 
% suggesting further parameters curation (enzyme usages > 10% of the total 
% proteome).
%
% INPUT:
%   ecModel         An enzyme constrained model.
%	c_source        The name of the exchange reaction that supplies the model
%                   with carbon.
%   sigma           Average saturation factor.
%   Ptot            Total protein content [g/gDW].
%   gRate           Minimum growth rate the model should grow at [1/h]. For
%                   finding the growth reaction, GECKO will choose the
%                   non-zero coeff in the objective function.
%   modifications	List of manually changed kcats (will be skipped in the
%                   kcat flexibilization).
%
% OUTPUT:
%   ecModel_batch	The enzyme constrained model under batch conditions.
%   OptSigma        Optimized saturation factor.
%
% Function adapted from GECKO toolbox: https://github.com/SysBioChalmers/GECKO
%
% Ivan Domenzain        2018-09-27
% Benjamin J. Sanchez   2018-12-11
% Carl Malina           2020-09-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ecModel_batch,OptSigma] = getConstrainedModel_modified(ecModel,modifications,name)
	
    %Get f (estimated mass fraction of enzymes in model)
    cd ..
    parameters = getModelParameters_modified;
    c_source   = parameters.c_source;
    sigma      = parameters.sigma;
    Ptot       = parameters.Ptot;
    gRate      = parameters.gR_exp;
    % Get f (estimated mass fraction of enzymes in model)
    cd ../../../../GECKO/geckomat/limit_proteins/
    [f,~] = measureAbundance(ecModel.enzymes);
    %Change media to batch conditions:
    cd ../../../mitoYeast-GEM/ecMitoYeast/ComplementaryScripts/customGECKO/kcat_sensitivity_analysis/
    ecModel = changeMedia_batch(ecModel,c_source);
    cd ../limit_proteins/
    %Get a preliminary enzyme constrained model for performing the Kcats
    %sensitivity analysis
    [ecModel_batch,~,~] = constrainEnzymes_modified(ecModel,f);
	%solution            = solveLP(ecModel_batch,1);
    cd ../../simulation/
    [mu,sol] = solveModel(ecModel_batch,0,gRate);
    if ~isempty(sol.f)
        %Set the media according to the experimental conditions
        cd ../customGECKO/kcat_sensitivity_analysis/
        %ObjIndex = find(ecModel_batch.c);
        % If the model is overconstrained
        if (gRate-mu)>0 
            fprintf('\n')
            disp('***************************************************************')
            disp('                The ECmodel is overconstrained                 ')
            %Perform a sensitivity analysis on the objective function with 
            %respect to the individual Kcat coefficients, the algorithm will 
            %iterate replacing the top limiting value according to the maximum 
            %value available in BRENDA for the same EC number until the objective
            %is no longer underpredicted 
            ecModel_batch = modifyKcats(ecModel_batch,gRate,modifications,name);
        else
            fprintf('\n')
            disp('***************************************************************')
            disp('              The ECmodel is not overconstrained               ')
        end    
        %The sigma factor is reffited for the specified conditions (constraints in the model)
        fprintf('\n')
        disp('***************************************************************')
        disp('        Fitting the average enzymes saturation factor          ')
        OptSigma          = sigmaFitter(ecModel_batch,Ptot,gRate,f);
        enzymePos         = strcmp(ecModel_batch.rxns,'prot_pool_exchange');
        currentEnzymeUB   = ecModel_batch.ub(enzymePos);
        ecModel_batch     = setParam(ecModel_batch,'ub','prot_pool_exchange', ...
                                     currentEnzymeUB*OptSigma/sigma);
        
        %Simulate growth on minimal media and export to the output folder:
        % 1) the exchange fluxes to the file "exchangeFluxes.txt"
        % 2) the top ten used enzymes to the file "topUsedEnzymes.txt"
        cd ../../simulation/
        [~,solution] = solveModel(ecModel_batch,0,gRate);
        if ~isempty(solution.x)
            disp('Saving simulation results files...')
            fluxFileName = ['../../../ModelFiles/' name '_exchangeFluxes.txt'];
            printFluxes(ecModel_batch,solution.x,true,10^-6,fluxFileName);
            cd ../customGECKO/kcat_sensitivity_analysis/
            topUsedEnzymes(solution.x,ecModel_batch,{'Min_glucose'},name);
            cd ../limit_proteins/
        end
    else
        disp('ecModel with enzymes pool constraint is not feasible')
    end
end
