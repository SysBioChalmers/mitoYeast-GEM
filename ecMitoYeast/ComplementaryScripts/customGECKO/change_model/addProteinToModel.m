%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addProteinToModel(model,P,kegg,swissprot)
% Adds an exchange reaction for protein P and updates model.enzymes,
% model.MWs and model.pathways to account for P.
%
% Function adapted from GECKO toolbox: 
%
% INPUT:
% model             Model with enzymes
% p                 Uniprot code of the protein
% kegg              KEGG database
% swissprot         Swissprot database
%
% OUTPUTS:
% model             Model with the added protein
% 
% Function adapted from the GECKO toolbox:
% https://github.com/SysBioChalmers/GECKO
% 
% Carl Malina. Last edited: 2020-01-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addProteinToModel(model,P,kegg,swissprot)

prot_name      = ['holo_prot_' P ];
model.enzymes = [model.enzymes;P];
pos_e          = strcmp(model.enzymes,P);
mitoProteins = {'P00163';'P00401';'P00410';'P00420';'P00854';'P00856'; ...
                'P61829'};

% Update model.MWs & model.sequences vectors:
match_gene     = false;
match_enzName  = false;
match_MW       = false;
match_seq      = false;
gene           = [];

for i = 1:length(swissprot)
    % Gene name
    if strcmp(P,swissprot{i,1}) && ~isempty(swissprot{i,3}) && ~match_enzName
        match_enzName     = true;
        geneNames         = swissprot{i,3};
        geneIDs           = strsplit(geneNames,' ');
        [geneSwissprot,~] = intersect(geneIDs,model.genes);
        if ~isempty(geneSwissprot)
            match_gene = true;
            gene       = geneSwissprot{1};
        end
        model.enzNames{pos_e,1} = geneIDs{1};
    end
    % Molecular weight
    if strcmp(P,swissprot{i,1}) && swissprot{i,5} ~= 0 && ~match_MW
        match_MW       = true;
        model.MWs(pos_e,1) = swissprot{i,5}/1000; % g/mmol
    end
    % Sequence
    if strcmp(P,swissprot{i,1}) && ~isempty(swissprot{i,6}) && ~match_seq
        match_seq                 = true;
        model.sequences{pos_e,1} = swissprot{i,6};
    end
end

if ~match_enzName
    model.enzName{pos_e,1} = '-';
end
if ~match_MW
    model.MWs(pos_e,1) = mean(cell2mat(swissprot(:,5)))/1000; % average MW [g/mmol]
end
if ~match_seq
    model.sequences{pos_e,1} = '-';
end

% Update model.genes & model.pathways
match_path = false;
for i = 1:length(kegg)
    if strcmp(P,kegg{i,1})
        % Gene
        if ~isempty(kegg{i,3}) && ~match_gene
            match_gene = true;
            gene       = kegg{i,3};
        end
        % Pathway
        if ~isempty(kegg{i,6}) && ~match_path
            match_path = true;
            model.pathways{pos_e,1} = kegg{i,6};
        end
        % Molecular way (if nothing found in uniprot)
        if kegg{i,5} > 0 && ~match_MW
            match_MW            = true;
            model.MWs(pos_e,1) = kegg{i,5}/1000; % g/mmol
        end
        % Sequence (if nothing found in uniprot)
        if ~isempty(kegg{i,7}) && ~match_seq
            match_seq                 = true;
            model.sequences(pos_e,1) = kegg(i,7);
        end
    end
end

if ~match_gene
    unknowns = ~cellfun(@isempty,strfind(model.genes,'unknown_'));
    if sum(unknowns) == 0
        idx = 0;
    else
        unknowns  = model.genes(unknowns);
        pos_final = strfind(unknowns{end},'_')+1;
        idx       = str2double(unknowns{end}(pos_final:end));
    end
    gene = ['unknown_' num2str(idx+1)];
end
if ~match_path
    model.pathways{pos_e,1} = '-';
end

model.enzGenes{pos_e,1} = gene;
% Add gene to list if non-existing previously
if ~ismember(gene,model.genes)
    geneToAdd.genes         = {gene};
    geneToAdd.geneShortName = model.enzNames{pos_e,1};
    model                  = addGenesRaven(model,geneToAdd);
end

% Add exchange reaction of protein: -> P
rxnID = [prot_name '_exchange'];
% all enzymes are added to the cytoplasm, except the ones synthesized
% in mitochondria
if sum(strcmp(mitoProteins,P)) > 0
    mitoIndex = strcmpi(model.compNames,'mitochondrion');
    comp      = model.comps{mitoIndex};
else
    cytIndex = strcmpi(model.compNames,'cytoplasm');
    comp    = model.comps{cytIndex};
end

%Add exchange reaction of protein: -> P
rxnToAdd.rxns         = {rxnID};
rxnToAdd.rxnNames     = {rxnID};
rxnToAdd.mets         = {prot_name};
rxnToAdd.stoichCoeffs = 1;
rxnToAdd.lb           = 0;
rxnToAdd.grRules      = {gene};
model                = addRxns(model,rxnToAdd,1,comp,true);

end