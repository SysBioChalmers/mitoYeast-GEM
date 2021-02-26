%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = adjustModel(model,k,block,scaling)
%
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Carl Malina. Last update: 2020-07-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = adjustModel(model,k,block,scaling)

%Block exchange of tails and backbones:
if block
    posT  = strcmp(model.rxnNames,'lipid chain exchange');
    posB  = strcmp(model.rxnNames,'lipid backbone exchange');
    model = changeRxnBounds(model,model.rxns(posT),0,'b');
    model = changeRxnBounds(model,model.rxns(posB),0,'b');
end

%Switch what to rescale depending on flag:
switch scaling
    case 'backbones'
        rxnName = 'lipid backbone pseudoreaction';
    case 'tails'
        rxnName = 'lipid chain pseudoreaction';
end

%Find positions:
scaleRxn  = strcmp(model.rxnNames,rxnName);
scaleMets = model.S(:,scaleRxn) < 0;

%Change stoich coeffs:
model.S(scaleMets,scaleRxn) = k*model.S(scaleMets,scaleRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%