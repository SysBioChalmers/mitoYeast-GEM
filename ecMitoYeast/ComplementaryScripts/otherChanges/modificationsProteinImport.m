%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modificationsProteinImport(model)
%
%
%
% Carl Malina. Last updated: 2020-12-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modificationsProteinImport(model)
% Update kcats for complexes involved in mitochondrial protein import

% Translocatse of the outer membrane (TOM). No kcat found in literature,
% but roughly 1/3 of all mitochondrial proteins can be co-translationally
% imported PMID: 31030976. Therefore, the translocation rate (9.5 aa/s) and
% the median protein length (358 aa) was used to calculate the rate of
% translocation. TOM complexes have been shown to dimerize and form two
% protein translocation channels.
newValue = (10/358)*[1,1,1,2,1,1,2];%(9.5/358)*[1,1,1,2,1,1,2];
model = changeKcats(model,{'holo_prot_P80967_om','holo_prot_P33448_om',...
    'holo_prot_P53507_om','holo_prot_P35180_om','holo_prot_P49334_om',...
    'holo_prot_P23644_om','holo_prot_P07213_om'},newValue,true);

% Mitochondrial inner membrane insertase OXA1 (P39952). Involved in
% insertion of protein into the mitochondrial inner membrane.
% Co-translational insertion of mitochondrially synthesized proteins.
% Rate set to translation rate (9.5 aa/s) over median protein length (358).
newValue = 10/358;%9.5/358;
model = changeKcats(model,'holo_prot_P39952_mm',newValue,true);

% BCS1 (P32839). Involved in membrane insertion of complex III subunit
% Rip1. Assumed to operate at the same rate as OXA1
newValue = 10/358;%9.5/358;
model = changeKcats(model,'holo_prot_P32839_m',newValue,true);

% Octapeptidyl aminopeptidase (P35999-EC3.4.24.59). Only values reported in
% BRENDA are from Human. Use the median turnover number (0.059 s-1).
model = changeKcats(model,'holo_prot_P35999_m',0.059,true);

% Intermediate cleaving peptidase (P40051-EC3.4.11.26). Turnover number
% for S. cerevisiae reported in BRENDA is orders of magnitude higher than
% values reported for octapeptidyl aminopeptidase. Assume similar rate as
% that enzyme (0.059 s-1).
model = changeKcats(model,'holo_prot_P40051_m',0.059,true);

% Mitochondrial processing peptidase (MPP) (P10507/P11914-EC3.4.24.64). Use
% median turnover number from Human available in BRENDA
model = changeKcats(model,{'holo_prot_P10507_m','holo_prot_P11914_m'},[0.059,0.059],true);

% Translocase of the inner membrane TIM23. Value based on ATPase activity
% of mtHsp70 (SSC1, 1 ATP/min), average length between binding sites within 
% protein (25 aa), median protein length (386) and stimulation of ATPase 
% activity by co-chaperones and target protein. 
newValue = 9/386;%(1/60/358)*25*11;
% Presequence translocase-associated motor (PAM). P38523 functions as
% a dimer
model    = changeKcats(model,{'holo_prot_P0CS90_m','holo_prot_Q01852_m',...
    'holo_prot_P42949_m','holo_prot_Q07914_m','holo_prot_P38523_m'},...
    newValue*[1,1,1,1,0.5],true);
% TIM23
model = changeKcats(model,{'holo_prot_P39515_mm','holo_prot_P53220_mm', ...
                           'holo_prot_P32897_mm','holo_prot_Q02776_mm'},...
                           newValue*ones(4,1),true);
end