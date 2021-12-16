%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ecModelBatch = modifyKcats(ecModelBatch,expVal,modifiedKcats,name)
%
% Function that gets the limiting Kcat values in an EC model (according to
% a sensitivity analysis), then it modifies each of those values according to 
% the maximum available values in the BRENDA files (Kcats and SA*Mw) when a
% manual curated option is not specified.
% 
% The algorithm iterates until the simulated value for the objective 
% function agrees with the provided by the user 
% (batch growth on glucose minimal media recommended).
%
% INPUTS
%   - ecModelBatch:  Enzyme-constrained GEM with the total protein pool
%                    global constraint.
%   - expVal:        Experimentally measured value for the objective
%                    function.
%   - modifiedKcats: Cell array containing IDs for the previously manually 
%                    modified kcats ('UniprotCode_rxnIndex')
%   - name:          String containing the name for the model files
%                    provided by the user.
% OUTPUTS
%   - ecModel:       Enzyme-constrained GEM with the automatically curated
%                    Kinetic parameters.
%
% Ivan Domenzain    Last edited. 2018-10-09
% Carl Malina       Last edited. 2020-04-24 - Adapt to ME-like concept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ecModelBatch = modifyKcats(ecModelBatch,expVal,modifiedKcats,name)

changes = [];
error   = -100; 
current = pwd;
%Load BRENDA data:
cd ../../../../../GECKO/geckomat/get_enzyme_data/
[BRENDA,SA_cell] = loadBRENDAdata;
cd(current)
%Iterates while the objective value is being underpredicted
disp('******************* Limiting Kcats curation *******************')
i=1; 
% Tolerance of 10% underprediction for allowing a sigma factor readjustment
while error<=-6%10
    cd(current)
    %Get the top limiting enzyme (uniprot code basis)
    [limKcat,breakFlag] = findTopLimitations(ecModelBatch,modifiedKcats,0,expVal);   
    if breakFlag == false
        disp(['*Iteration #' num2str(i)])
        [ecModelBatch,data] = changeKcat(ecModelBatch,limKcat,expVal,BRENDA,SA_cell);
        %Saves the parameter modification information
        changes = [changes; data];
        if ~isempty(data{1,9})
            error = data{1,9};
        end
        %Add a string with the uniprot code and the rxn number in order to 
        %keep track of the modified coefficients
        str           = {horzcat(data{1},'_',num2str(limKcat{3}))};
        modifiedKcats = [modifiedKcats; str];
        disp(horzcat('  Protein:',data{1},' Rxn#:',num2str(limKcat{1,3}),' name: ',limKcat{6}{1}))
        disp(['  prev_Kcat:' num2str(data{1,7}) ' new_Kcat:' num2str(data{1,8}) ...
              ' CC:' num2str(limKcat{1,5}) ' Err:' num2str(error) '%'])
        i = i+1;
    else
        break
    end
    fprintf('\n')
end
%Create a .txt file with all the modifications that were done on the
%individual Kcat coefficients
cd(current)
if ~isempty(changes)
    varNamesTable = {'Unicode','enz_pos','rxn_pos','Organism','Modified',...
                    'Parameter','oldValue','newValue','error','ControlCoeff'};  
    changes = cell2table(changes,'VariableNames',varNamesTable);
    
    changes = truncateValues(changes,[7:10]);
    %Write results in a .txt for further exploration.
    writetable(changes,['../../../ModelFiles/' name '_kcatModifications.txt']);
else
    %If the model is not feasible then the analysis is performed in all the 
    %Kcats matched either to: option 1 -> each of the enzymatic rxns, 
    %option 2 -> each of the individual enzymes
    [limRxns,~] = findTopLimitations(ecModelBatch,modifiedKcats,1,expVal);
    [limEnz, ~] = findTopLimitations(ecModelBatch,modifiedKcats,2,expVal);
    
    if ~isempty(limRxns)
        varNamesTable = {'rxnNames','rxnPos','ControlCoeff'};
        changes = cell2table(limRxns,'VariableNames',varNamesTable);
        changes = truncateValues(changes,3);
        writetable(changes,['../../../ModelFiles/' name '_limitingRxns.txt']);
    end
    if ~isempty(limEnz)
        varNamesTable = {'EnzNames','EnzPos','ControlCoeff'};
        changes = cell2table(limEnz,'VariableNames',varNamesTable);
        changes = truncateValues(changes,3);
        writetable(changes,['../../../ModelFiles/' name '_limitingEnzymes.txt']);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,output] = changeKcat(model,limKcats,objVal,BRENDA,SA_cell)
current_dir = pwd;
cd ../customGECKO/kcat_sensitivity_analysis/
% Gets the Unicode
prot = limKcats{1}{1};
prot_split = strsplit(prot,'_');
UniCode = prot_split{3};
%UniCode = limKcats{1}{1}(strfind(limKcats{1}{1},'_')+1:end);
% Map the UNIPROT code (kcat)
[ECnumber, ~] = findECnumber(UniCode);
enzIndx = limKcats{2}(1);
rxnIndx = limKcats{3}(1);
output  = {UniCode,[],[],[],[],[],[],[],[],[]};
error   = 0;
if ~isempty(ECnumber)
    flag           = false;
    previous_value = -1/(3600*model.S(enzIndx,rxnIndx)); %[1/s]
    disp([' Automatic search // ' 'EC#: ' ECnumber])
    %Looks for the maximal value available for the respective EC
    %number (Kcats and SA*Mw if available)
    [Kcat,org, match] = findMaxValue(ECnumber,BRENDA,SA_cell);
    coeff             = -1/(Kcat);
    %If a higher value was found then change the kinetic coefficient
    if coeff > model.S(enzIndx,rxnIndx)
        flag = true;
        model.S(enzIndx,rxnIndx) = coeff;
    end
    new_value = -1/(3600*model.S(enzIndx,rxnIndx));
    %After changing the i-th kcat limiting value a simulation is
    %performed and the objective value and absolute error are saved
    model_sim               = model;
    growthPos               = strcmpi(model.rxnNames,'growth');
    %model_sim.c          = zeros(size(model_sim.c));
    %model_sim.c(objPos)  = 1;
    %solution             = solveLP(model_sim);
    cd(current_dir)
    [~,solution]           = solveModel(model_sim,0,objVal);
    %model_sim.lb(growthPos) = 0.999*solution.x(growthPos);
    %model_sim.ub(growthPos) = solution.x(growthPos);
    %solution             = solveLP(model_sim);
    %[~,solution]            = solveModel(model_sim,0,solution.x(growthPos));
    %Calculate relative error
    error  = ((solution.x(growthPos)-objVal)/objVal)*100;
    output = {UniCode,limKcats{2}(1),limKcats{3}(1),org,flag,match,...
              previous_value,new_value,error,limKcats{5}(1)};
end
%cd(current)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ECnumber, Mw] = findECnumber(Unicode)
current = pwd;
load('../../../../../GECKO/databases/ProtDatabase.mat')
DB1{1} = swissprot(:,1);DB1{2} = swissprot(:,4);DB1{3} = swissprot(:,5);
DB2{1} = kegg(:,1);     DB2{2} = kegg(:,4);     DB2{3} = kegg(:,5);
ECnumber = {};
Mw       = {};
%First search for the UNIPROT ID in the swissprot DB structure
matching = find(strcmpi(DB1{1},Unicode));
if ~isempty(matching)
    ECnumber = DB1{2}{matching};
    Mw       = DB1{3}{matching};
end
%If nothing comes up then look into the KEGG DB structure
if isempty(ECnumber)
    matching = find(strcmpi(DB2{1},Unicode));
    if ~isempty(matching)
        ECnumber = DB2{2}{matching};
        Mw       = DB2{3}{matching};
    end
end
cd (current)
end
