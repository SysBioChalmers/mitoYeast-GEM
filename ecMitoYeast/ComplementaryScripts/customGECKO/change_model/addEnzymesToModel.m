%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = addEnzymesToModel(model,uniprots,kcats)
% Adds enzymes and exchange rxns for enzymes to model
% 
% Function adapted from: https://github.com/SysBioChalmers/GECKO/
%
% Carl Malina. Last edited: 2020-01-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eModel = addEnzymesToModel(model,uniprots,kcats)
current = pwd;
% Load databases
cd ../../../../../GECKO/databases/
data      = load('ProtDatabase.mat');
swissprot = data.swissprot;
kegg      = data.kegg;

eModel   = model;
enzymes  = cell(5000,1);
[m,n]    = size(uniprots);
y        = 0;

if isfield(model,'enzymes')
    enzModel = true;
else
    enzModel = false;
end
%Create a vector for enzymes
for i = 1:m
    for j = 1:n
       if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 % if kcat = 0 no enzyme will be added
           uniprots{i,j} = strsplit(uniprots{i,j},' ');
           for k = 1:length(uniprots{i,j})
               y = y + 1;
               enzymes{y} = uniprots{i,j}{k};
           end
       end
    end
end

% Eliminate duplicated uniprots from enzymes vector
enzymes(y+1:end) = [];
enzymes          = unique(enzymes);

% Add additional fields to model
if ~enzModel
    eModel.enzymes   = cell(0,1);
    eModel.enzGenes  = cell(0,1);
    eModel.enzNames  = cell(0,1);
    eModel.MWs       = zeros(0,1);
    eModel.sequences = cell(0,1);
    eModel.pathways  = cell(0,1);
end
cd(current)
for i = 1:length(enzymes)
    enzyme = enzymes{i};
    if enzModel
        % check if enzyme already in model before adding enzymes to model
        enzInModel = strcmp(model.enzymes,enzyme);
        if sum(enzInModel) == 0
            eModel = addProteinToModel(eModel,enzyme,kegg,swissprot);
        end
    else
        eModel = addProteinToModel(eModel,enzyme,kegg,swissprot);
    end
end

end
