%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_data,kcats] = getKcatsfromGECKO(model)
%
% This function collect the kcats for the enzymes in the model. This
% function depends on functions in GECKO: https://github.com/SysBioChalmers/GECKO
%
% INPUT:
% model: A RAVEN-compatible GEM structure (1x1 struct)
%
% OUTPUT:
%           model_data, which contains:
%           *model:      The standardized GEM
%           *substrates: Substrates associated for each rxn
%           *products:   Products associated, when rxn is reversible
%           *uniprots:   All possible uniprot codes, for each rxn
%           *EC_numbers: All possible EC numbers, for each uniprot
%           *count(1):   #rxns with data from Swissprot
%           *count(2):   #rxns with data from KEGG
%           *count(3):   #exchange/transport rxns with no GPRs
%           *count(4):   #other rxns
%           
%           kcats, which contains:
%           *forw.kcats:   kcat values for the forward reactions (mxn)
%           *forw.org_s:   Number of matches for organism - substrate in
%                          forward reaction (mxn)
%           *forw.rest_s:  Number of matches for any organism - substrate
%                          in forward reaction (mxn)
%           *forw.org_ns:  Number of matches for organism - any substrate
%                          in forward reaction (mxn)
%           *forw.rest_ns: Number of matches for any organism - any
%                          substrate in forward reaction (mxn)
%           *forw.org_sa:  Number of matches for organism - using s.a.
%                          in forward reaction (mxn)
%           *forw.rest_sa: Number of matches for any organism - using s.a.
%                          in forward reaction (mxn)
%           *back.kcats:   kcat values for the backward reactions (mxn)
%           *back.org_s:   Number of matches for organism - substrate in
%                          backwards reaction (mxn)
%           *back.rest_s:  Number of matches for any organism - substrate
%                          in backwards reaction (mxn)
%           *back.org_ns:  Number of matches for organism - any substrate
%                          in backwards reaction (mxn)
%           *back.rest_ns: Number of matches for any organism - any
%                          substrate in backwards reaction (mxn)
%           *back.org_sa:  Number of matches for organism - using s.a.
%                          in backwards reaction (mxn)
%           *back.rest_sa: Number of matches for any organism - using s.a.
%                          in backwards reaction (mxn)
%           *tot.queries:  The total amount of ECs matched (1x1)
%           *tot.org_s:    The amount of ECs matched for the organism & the
%                          substrate (1x1)
%           *tot.rest_s:   The amount of ECs matched for any organism & the
%                          substrate (1x1)
%           *tot.org_ns:   The amount of ECs matched for the organism & any
%                          substrate (1x1)
%           *tot.rest_ns:  The amount of ECs matched for any organism & any
%                          substrate (1x1)
%           *tot.org_sa:   The amount of ECs matched for the organism & 
%                          using s.a. (1x1)
%           *tot.rest_sa:  The amount of ECs matched for any organism & 
%                          using s.a. (1x1)
%
% Carl Malina. Last edited: 2020-01-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model_data, kcats] = getKcatsFromGecko(model)
% Define the organism name
org_name = 'saccharomyces cerevisiae';
format short e

% retrieve kcats for each reaction in the model
cd ../../../../../GECKO/geckomat/get_enzyme_data/
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);
%kcats      = kcats.forw.kcats(:,1);
%rxnID      = model.rxns;
%grRules    = model.grRules;
%M_kcats    = [rxnID grRules num2cell(kcats)];

cd ../../../mitoYeast-GEM/ecMitoYeast/ComplementaryScripts/

end

