%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L,I,F] = sumBioMass(model,data)
% Calculates breakdown of biomass
%
% model     metabolic model in COBRA format
% data      structure with at least the following 2 fields:
%   mets    Cell array with metabolite ids
%   MWs     Numeric array with molecular weights for each metabolite
%
% X         Total biomass fraction [gDW/gDW]
% P         Protein fraction [g/gDW]
% C         Carbohydrate fraction [g/gDW]
% R         RNA fraction [g/gDW]
% D         DNA fraction [g/gDW]
% L         Lipid fraction [g/gDW]
% F         cofactor [g/gDW]
% I         ion [g/gDW]
% 
%
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Benjamin Sanchez. Update: 2018-09-04
% Feiran Li. Last update: 2018-09-24 -Update the sumBiomass to include ions
% and cofactors
% Carl Malina. 2020-03-06. Update to include prosthetic groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L,I,F,G] = sumBioMass(model,data)
if nargin < 2
    %Load data for biomass composition:
    fid = fopen('biomassComponents_ecModel.txt');
    d = textscan(fid,'%s %s %f32 %f32 %s','Delimiter','\t','HeaderLines',1);
    data.mets       = d{1};
    data.abundances = double(d{3});
    data.MWs        = double(d{4});
    data.groups     = d{5};
    fclose(fid);
end
%Get main fractions:
[P,X] = getFraction(model,data,'P',0);
[C,X] = getFraction(model,data,'C',X);
[R,X] = getFraction(model,data,'R',X);
[D,X] = getFraction(model,data,'D',X);
[L,X] = getFraction(model,data,'L',X);
[I,X] = getFraction(model,data,'I',X);
[F,X] = getFraction(model,data,'F',X);
[G,X] = getFraction(model,data,'G',X);

%disp(['X -> ' num2str(X) ' gDW/gDW'])

% Simulate growth:
%sol = optimizeCbModel(model);
%disp(['Growth = ' num2str(sol.f) ' 1/h'])
%disp(' ')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,data,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'N','biomass');
rxnName = strrep(rxnName,'L','lipid backbone');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'I','ion');
rxnName = strrep(rxnName,'F','cofactor');
rxnName = strrep(rxnName,'G','prosthetic groups');


%Add up fraction:
rxnPos = strcmp(model.rxnNames,rxnName);
if ~all(rxnPos==0)
    isSub   = model.S(:,rxnPos) < 0;        %substrates in pseudo-rxn
    if strcmp(compType,'L')
        F = -sum(model.S(isSub,rxnPos));   %g/gDW
    else
        F = 0;
        %Add up all components:
        for i = 1:length(model.mets)
            pos = strcmp(data.mets,model.mets{i});
            if isSub(i) && sum(pos) == 1
                if strcmp(compType,'I') || strcmp(compType,'F') || strcmp(compType,'G')
                    MW = data.MWs(pos);
                else
                    MW = data.MWs(pos)-18;
                end
                abundance = -model.S(i,rxnPos)*MW/1000;
                F         = F + abundance;
            end
        end
    end
    X = X + F;

    %disp([compType ' -> ' num2str(F) ' g/gDW'])
else
    disp([compType ' do not exist '])
        F = 0;
    X = X + F;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%