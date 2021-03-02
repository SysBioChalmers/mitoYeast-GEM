%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [energy, redox] = checkATPandRedox(model,rxn,energyResults,redoxResults)
%
% This function checks what effect the added/modified reaction has on the 
% ATP production and redox balance
%
% Function modified from yeastGEM: https://github.com/SysBioChalmers/yeast-GEM
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energy,redox,growth] = checkATPandRedox(model,rxn)

[~,rxnInd] = ismember(rxn,model.rxns);

if rxnInd ~= 0
    model_temp = model;
    
    % change to minimal media
    exchangeRxns = findExcRxns(model_temp);

    desiredExchange = {'D-glucose';
                   'ammonium';
                   'oxygen';
                   'biotin';
                   'iron(2+)';
                   'H+';
                   'water';
                   'chloride';
                   'Mn(2+)';
                   'Zn(2+)';
                   'Mg(2+)';
                   'sodium';
                   'Cu2(+)';
                   'Ca(2+)';
                   'potassium';
                   'sulphate';
                   'phosphate'};
               
    exchangesToBlock = {'bicarbonate';
                        'lipid backbone';
                        'lipid chain'};
               
    model_temp.lb(exchangeRxns) = 0;
    model_temp.ub(exchangeRxns) = 1000;
    
    for i = 1:length(desiredExchange)
        excRxn = [desiredExchange{i} ' exchange'];
        excRxnInd = strcmp(model_temp.rxnNames,excRxn);
        
        if strcmp(excRxn,'D-glucose exchange')
            model_temp.lb(excRxnInd) = -1;
        else
            model_temp.lb(excRxnInd) = -1000;
        end
    end
     
    for i = 1:length(exchangesToBlock)
       excRxn = [exchangesToBlock{i} ' exchange'];
       excRxnInd = ismember(model_temp.rxnNames,excRxn);
       model_temp.lb(excRxnInd) = 0;
       model_temp.ub(excRxnInd) = 0;
    end
    count = 0;
    if count == 0
        FBAsol = optimizeCbModel(model_temp);
        gr = FBAsol.obj;
        disp(['After changing exc rxn bounds, model grows at: ' num2str(gr)])
    end
    % Check growth
    sol = optimizeCbModel(model_temp);
    growth = [model.rxns{rxnInd} ' ' num2str(sol.obj)];
    
    % Add reaction hydrolyzing ATP

    testATP_mets   = {'ATP [cytoplasm]';'H2O [cytoplasm]';'ADP [cytoplasm]'; ...
                      'H+ [cytoplasm]';'phosphate [cytoplasm]'};
    [~,metInd]         = ismember(testATP_mets,model_temp.metNames);
    testATP_metIDs = model_temp.mets(metInd);
    coeffs         = [-1,-1,1,1,1];
    model_temp = addReaction(model_temp,{'ATP hydrolysis','leakTest'}, ...
                         testATP_metIDs,coeffs,false,0,1000);
    model_temp = changeObjective(model_temp,'ATP hydrolysis',1);
    sol = optimizeCbModel(model_temp);
    
    if sol.obj <= 100 && sol.obj > 0
        energy = [model.rxns{rxnInd} ' pass ' num2str(sol.obj)];
    elseif sol.obj > 100
        energy = [model.rxns{rxnInd} ' fail ' num2str(sol.obj)];
    else
        energy = [model.rxns{rxnInd} ' error ' 'error'];
    end
    
    % Add reaction oxidizing NADH
    testNADH_mets   = {'NADH [cytoplasm]';'H+ [cytoplasm]';'NAD [cytoplasm]'};
    [~,metInd]      = ismember(testNADH_mets,model_temp.metNames);
    testNADH_metIDs = model_temp.mets(metInd);
    coeffs = [-1,-1,1];
    model_temp = addReaction(model_temp,{'NADH oxidation','leakTest2'}, ...
                             testNADH_metIDs,coeffs,false,0,1000);
    model_temp = changeObjective(model_temp,'NADH oxidation',1);
    sol = optimizeCbModel(model_temp);
    
    if sol.obj <= 100 && sol.obj > 0
        redox = [model.rxns{rxnInd} ' pass ' num2str(sol.obj)];
    elseif sol.obj > 100
        redox = [model.rxns{rxnInd} ' fail ' num2str(sol.obj)];
    elseif sol.obj <= 0
        redox = [model.rxns{rxnInd} ' error ' ' error '];
    end
else
    redox = ['Already exists' ' skip ' ' skip '];
end

end