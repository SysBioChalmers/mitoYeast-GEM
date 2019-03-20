%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = rescalePseudoReaction(model,metName,f)
%
% Function adapted from yeastGEM: https://github.com/SysBioChalmers/yeast-GEM
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = rescalePseudoReaction(model,metName,f)

rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
for i = 1:length(model.mets)
   S_ir = model.S(i,rxnPos);
   isProd = strcmp(model.metNames{i},[metName ' [cytoplasm]']);
   if S_ir ~= 0 && ~isProd
      model.S(i,rxnPos) = f*S_ir; 
   end
end

end