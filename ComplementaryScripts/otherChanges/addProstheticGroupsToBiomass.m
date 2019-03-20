%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addProstheticGroupsToBiomass(model)
%
% This function creates a pseudoreaction for protein prosthetic groups,
% resulting in the formation of the pseudometabolite 'prosthetic groups
% [cytoplasm]'
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addProstheticGroupsToBiomass(model)

% Calculate the content of prosthetic groups
cofactorContent = calculateCofactorContent(model);

coFactors = fieldnames(cofactorContent);

% collect info in an array
cofactorData = {};
for i = 1:length(coFactors)
   cofactor = coFactors{i};
   abundance = cofactorContent.(cofactor);
   cofactorAbundance = {cofactor,abundance};
   cofactorData = [cofactorData;cofactorAbundance];
end

% merge values for heme a [m] and [c] and ferroheme b [m] and [c], since 
% both taken from same compartment in model
heme_a      = 0;
ferroheme_b = 0;
biotin      = 0;
for i = 1:length(cofactorData)
   cofactor = cofactorData{i,1};
   if contains(cofactor,'heme_a')
       heme_a = heme_a + cofactorData{i,2};
   elseif contains(cofactor,'ferroheme_b')
      ferroheme_b = ferroheme_b + cofactorData{i,2};
   elseif contains(cofactor,'biotin')
       biotin = biotin + cofactorData{i,2};
   end
end

% get position of mets in model
fe2s2_m_ID     = model.mets{strcmp(model.metNames,'[2Fe-2S] iron-sulfur cluster [mitochondrion]')};
fe4s4_m_ID     = model.mets{strcmp(model.metNames,'[4Fe-4S] iron-sulfur cluster [mitochondrion]')};
fe3s4_m_ID     = model.mets{strcmp(model.metNames,'[3Fe-4S] iron-sulfur cluster [mitochondrion]')};
fe2s2_c_ID     = model.mets{strcmp(model.metNames,'[2Fe-2S] iron-sulfur cluster [cytoplasm]')};
fe4s4_c_ID     = model.mets{strcmp(model.metNames,'[4Fe-4S] iron-sulfur cluster [cytoplasm]')};
lipoate_m_ID   = model.mets{strcmp(model.metNames,'lipoate (protein bound) [mitochondrion]')};
hemeA_ID       = model.mets{strcmp(model.metNames,'heme a [cytoplasm]')};
ferroheme_b_ID = model.mets{strcmp(model.metNames,'ferroheme b [mitochondrion]')};
siroheme_c_ID  = model.mets{strcmp(model.metNames,'siroheme [cytoplasm]')};
biotin_ID      = model.mets{strcmp(model.metNames,'biotin (protein bound) [cytoplasm]')};

% Add pseudoreaction for prosthetic groups
metID   = ['s_' getNewIndex(model.mets) '[c]'];
metName = 'prosthetic groups [cytoplasm]';
model = addMetabolite(model,metID,'metName',metName);

metIDs       = {fe2s2_m_ID;fe4s4_m_ID;fe3s4_m_ID;fe2s2_c_ID;fe4s4_c_ID;
                lipoate_m_ID;hemeA_ID;ferroheme_b_ID;siroheme_c_ID;biotin_ID;
                metID};
stoichCoeffs = -[cofactorData{strcmp(cofactorData(:,1),'fe2s2_m'),2};
                cofactorData{strcmp(cofactorData(:,1),'fe4s4_m'),2};
                cofactorData{strcmp(cofactorData(:,1),'fe3s4_m'),2};
                cofactorData{strcmp(cofactorData(:,1),'fe2s2_c'),2};
                cofactorData{strcmp(cofactorData(:,1),'fe4s4_c'),2};
                cofactorData{strcmp(cofactorData(:,1),'lipoate_m'),2};
                heme_a;ferroheme_b;
                cofactorData{strcmp(cofactorData(:,1),'siroheme_c'),2};
                biotin];
            
rxnID = ['r_' getNewIndex(model.rxns)];
rxnName = 'prosthetic groups pseudoreaction';

model = addReaction(model,rxnID,...
                    'reactionName',rxnName,...
                    'metaboliteList',metIDs,...
                    'stoichCoeffList',[stoichCoeffs;1],...
                    'reversible',false,...
                    'checkDuplicate',true);

% remove heme a from cofactor pseudoreaction
cofactorPseudoRxnPos = strcmp(model.metNames,'cofactor pseudoreaction');
hemeApos = strcmp(model.metNames,'heme a [cytoplasm]');
model.S(hemeApos,cofactorPseudoRxnPos) = 0;

% Allow biotin uptake, since present in standard minimal media
biotinExchange           = strcmp(model.rxnNames,'biotin exchange');
model.lb(biotinExchange) = -1000;
biotinUptake             = strcmp(model.rxnNames,'biotin uptake');
model.lb(biotinUptake)   = -1000;

% Add prosthetic groups [cytoplasm] to protein pseudoreaction
prostheticGroupsPos = strcmp(model.metNames,'prosthetic groups [cytoplasm]');
proteinRxnPos       = strcmp(model.rxnNames,'protein pseudoreaction');
model.S(prostheticGroupsPos,proteinRxnPos) = -1;

end