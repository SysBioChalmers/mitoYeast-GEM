%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = readKcatDataGECKO(model_data,kcats)
% Reads the output of the getBRENDAdata module, and creates the modified
% enzyme model, with additional metabolites (the enzymes) and reactions
% (for all isoenzymes and also the enzyme exchange reactions).
%
% INPUT:
% model_data        model and EC numbers and substrates/products from each
%                   reaction (output from "getECnumbers.m")
% kcats             kcats for each reaction/enzyme (output from
%                   "matchKcats.m")
%
% OUTPUT:
% eModel            modified model accounting for enzymes
%
% Function adapted from the GECKO toolbox:
% https://github.com/SysBioChalmers/GECKO
% 
% Carl Malina. Last edited: 2020-01-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = readKcatDataGECKO(model_data,kcats)

%Get kcat value for both directions:
Fkcat = kcats.forw.kcats;
Bkcat = kcats.back.kcats;
rev   = logical(model_data.model.rev);
kcats = [Fkcat;Bkcat(rev,:)];

%Update uniprots with both directions:
uniprots = [model_data.uniprots; model_data.uniprots(rev,:)];

%Update matched genes with both directions:
matchedGenes = [model_data.matchedGenes; model_data.matchedGenes(rev,:)];

%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
model = convertToIrrev(model_data.model);

%Convert original model to enzyme model according to uniprots and kcats:
% Uses a function adapted to add co-factor modifications and mitochondrial
% protein import reactions
eModel = convertToMitoEnzymeModel(model,matchedGenes,uniprots,kcats);

%Leave all UB = +Inf:
eModel.ub(eModel.ub == 1000) = +Inf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%