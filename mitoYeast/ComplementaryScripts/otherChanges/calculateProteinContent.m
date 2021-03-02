%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ptot = calculateProteinContent(model,rxnName,aa_abbr,aa_MWs)
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ptot = calculateProteinContent(model,rxnName,aa_abbr,aa_MWs)
F = 0;
proteinRxn = strcmp(model.rxnNames,rxnName);
isSub = model.S(:,proteinRxn) < 0;
substrate_IDs = model.mets(isSub);
substrate_names = model.metNames(isSub);
coeffs = model.S(isSub,proteinRxn);
for j = 1:length(substrate_names)
    aatRNA = substrate_names{j};
    abbr = aatRNA(1:3);
    MW = aa_MWs(strcmp(aa_abbr,abbr));
    abundance = -coeffs(j)*MW/1000;
    F = F + abundance;
end
Ptot = F;
end