%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modificationsTranslocationISCbiogenesis(model)
%
%
%
% Carl Malina. Last updated: 2020-11-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modificationsISCbiosynthesis(model)
% Make curations to iron-sulfur cluster biogenesis
current = pwd;
% Add Ssq1 (Q05931) and Jac1 (P53193) to model
% Load file with information about enzyme localization (used to add
% enzymes and translocation of the enzymes)
cd ../../ComplementaryData/proteinImport/
fid = fopen('enzymeLocalization.tsv');
fileData = textscan(fid,'%s%s%s%s%s%s%s%s%s%s','Delimiter','\t',...
                     'HeaderLines',1);
enzymeLoc.uniprots         = fileData{1};
enzymeLoc.genes             = fileData{2};
enzymeLoc.geneShort         = fileData{3};
enzymeLoc.localization      = fileData{4};
enzymeLoc.pathway           = fileData{5};
enzymeLoc.MPPcleavageSite   = cellfun(@str2num,fileData{7},'UniformOutput',false);
enzymeLoc.OCT1cleavageSite  = cellfun(@str2num,fileData{8},'UniformOutput',false);
enzymeLoc.ICP55cleavageSite = cellfun(@str2num,fileData{9},'UniformOutput',false);
fclose(fid);
% Add enzymes to model
cd ../../ComplementaryScripts/customGECKO/change_model/
uniprots = enzymeLoc.uniprots(ismember(enzymeLoc.geneShort,{'JAC1','SSQ1'}));
kcats = [1;1]; % kcats will be updated in a subsequent step
model = addEnzymesToModel(model,uniprots,kcats);
% Add translocation reactions
cd(current)
model = addTranslocationRxns(model,enzymeLoc);

% Add enzymes to [2Fe-2S] iron-sulfur cluster transfer. No EC number
% identified. kcat for HSP70 (SSQ1,PMID: 32397253) used for reaction. 
% Transfer involves proteins Grx5 (Q02784), Ssq1 (Q05931), Jac1 (P53193) 
% and Mge1 (P38523).
coeff = -1/(0.035*3600);
model.S(ismember(model.mets,{'holo_prot_P53193_m','holo_prot_Q05931_m'}),...
    strcmp(model.rxnNames,'[2Fe-2S] iron-sulfur cluster transfer')) = coeff;
% Account for stoichiometry of GRX5 (2/[2Fe-2S] cluster) and MGE1 (functional form is dimeric)
model.S(ismember(model.mets,{'holo_prot_Q02784_m','holo_prot_P38523_m'}),strcmp(model.rxnNames,'[2Fe-2S] iron-sulfur cluster transfer')) = coeff*2;
% Update rxn ID and name to reflect inclusion of proteins
rxnPos = strcmp(model.rxnNames,'[2Fe-2S] iron-sulfur cluster transfer');
model.rxnNames{rxnPos} = ['[2Fe-2S] iron-sulfur cluster transfer' ' (No1)'];
model.rxns{rxnPos} = [model.rxns{rxnPos} 'No1'];

% Block reaction r_4217, which otherwise causes production of oxygen
model.ub(strcmp(model.rxns,'r_4217_REVNo1')) = 0;

% Cysteine desulfurase (EC2.8.1.7). According to Complex Portal, the
% cysteine desulfurase complex is heterotrimeric. kcat used is for the
% persulfide transferl to scaffold for cluster synthesis PMID: 31395877
coeff      = -1/((2.4)*60); % min-1 -> h-1
model.S(ismember(model.mets,{'holo_prot_P25374_m','holo_prot_P32463_m','holo_prot_Q6Q560_m'}),strcmp(model.rxns,'r_4173No1')) = coeff;

% Update [2Fe-2S] iron-sulfur cluster synthesis rxns
% remove use of cysteine desulfurase complex as these are represented in a
% separate reaction
model.S(ismember(model.mets,{'holo_prot_P25374_m','holo_prot_P32463_m','holo_prot_Q6Q560_m'}),contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster synthesis')) = 0;
% Add requirement of yeast ferredoxin (YAH1, Q12184) and update kcat of
% reaction [2Fe-2S] iron-sulfur cluster synthesis (No1)/(No2). kcat taken
% from PMID: 31395877
coeff = -1/((2.5/2)*60); % min-1 -> h-1
model.S(ismember(model.mets,{'holo_prot_Q03020_m','holo_prot_Q12184_m','holo_prot_Q07540_m'}),contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster synthesis (No1)')) = coeff;
model.S(ismember(model.mets,{'holo_prot_Q12056_m','holo_prot_Q12184_m','holo_prot_Q07540_m'}),contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster synthesis (No2)')) = coeff;
% Update kcat for Yah1 in [2Fe-2S] iron-sulfur cluster biosynthesis
coeff = -1/((4.5/2)*60); % min-1 -> h-1
model.S(ismember(model.mets,{'holo_prot_Q12184_m'}),contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster synthesis (No1)')) = coeff;
model.S(ismember(model.mets,{'holo_prot_Q12184_m',}),contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster synthesis (No2)')) = coeff;

% Change kcat of adrenodoxin reductase homolog rxn (r_4630No1), assuming
% that it operates at the same rate as reduction of the persulfide
% reduction step of [2Fe-2S] cluster synthesis. Value from PMID: 31395877
coeff = -1/(4.5*60); % min-1 -> h-1
model.S(ismember(model.mets,'holo_prot_Q12184_m'),contains(model.rxnNames,'adrenodoxin:NADP+ oxidoreductase (No1)')) = 0; % remove use of Yah1 since used in subsequent reactions
model.S(ismember(model.mets,'holo_prot_P48360_m'),contains(model.rxnNames,'adrenodoxin:NADP+ oxidoreductase (No1)')) = coeff;

% Update kcat for [4Fe-4S] cluster reduction. kcat assumed to be the same
% as for [2Fe-2S] cluster assembly PMID: 31395877.
coeff = -1/(4.5*60); % min-1 -> h-1
model.S(ismember(model.mets,{'holo_prot_Q02784_m','holo_prot_Q07821_m','holo_prot_Q12425_m',...
    'holo_prot_P47158_m','holo_prot_Q12184_m'}),contains(model.rxnNames,'[4Fe-4S] iron-sulfur cluster reduction (No1)')) = coeff;
% Account for stoichiometry of Yah1 (2 ferredoxin required) and Grx5
model.S(ismember(model.mets,{'holo_prot_Q12184_m','holo_prot_Q02784_m'}),...
    contains(model.rxnNames,'[4Fe-4S] iron-sulfur cluster reduction (No1)')) = coeff*2;

% Change kcat of [2Fe-2S] iron-sulfur cluster transfer 2 -
% co-chaperone independent, therefore assumed to be slower. Use kcat
% without co-chaperone stimulation from PMID: 32397253.
coeff = -1/(0.0087*3600);
model.S(ismember(model.mets,{'holo_prot_Q02784_m','holo_prot_Q07821_m','holo_prot_Q12425_m',...
    'holo_prot_P47158_m'}),contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster transfer 2 (No1)')) = coeff;
model.S(contains(model.mets,'holo_prot_Q02784_m'),...
    contains(model.rxnNames,'[2Fe-2S] iron-sulfur cluster transfer 2 (No1)')) = coeff*2;

% Remove usage of proteins in r_4642No1 [4Fe-4S] iron-sulfur cluster release
% since protein usage included in preceeding rxns
model.S(ismember(model.mets,{'holo_prot_Q02784_m','holo_prot_Q07821_m','holo_prot_Q12425_m',...
    'holo_prot_P47158_m'}),contains(model.rxnNames,'[4Fe-4S] iron-sulfur cluster release (No1)')) = 0;
end