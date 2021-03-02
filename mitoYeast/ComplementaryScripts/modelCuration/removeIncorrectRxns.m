%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,notFound] = removeIncorrectRxns(model)
%
%
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model, notFound, growthRate] = removeIncorrectRxns(model)
% Load file with information on reactions to remove

cd ../../ComplementaryData/modelCuration/

fid     = fopen('rxnsToRemove.tsv');
rxnfile = textscan(fid,'%s%s%s%s%s%s%s%s','HeaderLines',1,'Delimiter','\t');
fclose(fid);

rxnInfo.genes    = rxnfile{1};
rxnInfo.rxnIDs   = rxnfile{3};
rxnInfo.rxnNames = rxnfile{4};

% remove reactions
notFound = {};
growthRate = zeros(length(rxnInfo.rxnIDs)+1,1);
FBAsol_orig = optimizeCbModel(model);
growthRate(1) = FBAsol_orig.obj;
for i = 1:length(rxnInfo.rxnIDs)
   rxnID   = rxnInfo.rxnIDs{i};
   rxnName = rxnInfo.rxnNames{i};
   rxnPos  = strcmp(model.rxns,rxnID);
   if sum(rxnPos) == 1 && sum(strcmp(model.rxnNames,rxnName)) > 0 
       model = removeRxns(model,{rxnID});
       FBAsol_deletion = optimizeCbModel(model);
       growthRate(i+1) = FBAsol_deletion.obj;
   % Handle exception for r_4249    
   elseif strcmp(rxnID,'r_4249') && contains(rxnName,'O3-acetyl-L-serine')
       model = removeRxns(model,{rxnID});
       FBAsol_deletion = optimizeCbModel(model);
       growthRate(i+1) = FBAsol_deletion.obj;
   else
       notFound = [notFound;rxnID];
       disp(['Rxn ID ' rxnID ' does not match with the specified rxn name'])
   end
end

model = removeUnusedGenes(model);
cd ../../ComplementaryScripts/modelCuration/
end
        