%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addCofactorsToBiomass(model)
%
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coFactorContent = calculateCofactorContent(model)
cd ../../ComplementaryData/databases/
% read PaxDb abundance data
fid = fopen('proteinAbundancePaxDb.txt');
pdb = textscan(fid,'%s %s %s','Delimiter','\t','HeaderLines',13);
paxDb.genes = strrep(pdb{2},'4932.','');
paxDb.abundances = cellfun(@str2num,pdb{3})*1e-6;
fclose(fid);

% read uniprot protein information data
fid = fopen('uniprot.tsv');
up  = textscan(fid,'%s %s %s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
uniprot.gene_ids = up{3};
uniprot.MWs       = cellfun(@str2num, up{5})*1000; % kDa --> g/mol
fclose(fid);

% read file with cofactor content information
cd ../proteinInfo/
fid    = fopen('coFactorInfo.tsv');
coFact = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t', ...
                  'HeaderLines',1);
coFactorInfo.genes    = coFact{1};
coFactorInfo.complex  = coFact{3};
coFactorInfo.type     = coFact{4};
coFactorInfo.stoich   = cellfun(@str2num,coFact{5});
coFactorInfo.group     = coFact{6};
fclose(fid);

% read file with stoichiometry info for mitochondrial complexes
fid     = fopen('complexStoichiometry.tsv');
compl   = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
complexInfo.genes    = compl{1};
complexInfo.complex  = compl{3};
complexInfo.stoich   = cellfun(@str2num,compl{4});
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Estimate the cofactor content based on protein abundance %%%%%%%%%%

% calculate the protein content in the model
cd ../Physiology/
clear aa_data aa_usage
fid = fopen('mitoAminoAcidUsage.tsv');
aa_usage          = textscan(fid,'%s %s %f %f %f %f','Delimiter','\t','HeaderLines',1);
aa_data.name      = aa_usage{1}; % amino acid name
aa_data.abbr      = aa_usage{2}; % 3 letter abbreviations
aa_data.MWs       = aa_usage{3}; % molecular weight in polymer
fclose(fid);

cd ../../ComplementaryScripts/otherChanges/
rxnName = 'protein pseudoreaction';
Ptot = calculateProteinContent(model,rxnName,aa_data.abbr,aa_data.MWs);

% calculate the abundance of cofactor-containing enzymes (mmol/gDW)
coFactorInfo.abundances = zeros(size(coFactorInfo.genes));
for i = 1:length(coFactorInfo.genes)
    protein = coFactorInfo.genes{i};
    complex = coFactorInfo.complex{i};
    if contains(complex,'complex')
       % find genes in complex
       protein_stoich = complexInfo.stoich(strcmp(complexInfo.genes,protein)); % Stoichiometry of actual subunit
       genePos        = strcmp(complexInfo.complex,complex);
       genes          = complexInfo.genes(genePos);
       subunit_stoich = complexInfo.stoich(genePos);
       abundance = 0;
       for j = 1:length(genes)
           gene            = genes{j};
           stoichCoeff     = subunit_stoich(j);
           pos             = strcmp(paxDb.genes,gene);
           paxDb_abundance = paxDb.abundances(pos)/stoichCoeff;
           if paxDb_abundance > abundance
               abundance = paxDb_abundance;
           end
       end
       abundance = abundance*protein_stoich;
    else
        pos = strcmp(paxDb.genes,protein);
        abundance = paxDb.abundances(pos);
    end

    % calculate cofactor content (mmol/gDW)
    for k = 1:length(uniprot.gene_ids)
        gene_id = uniprot.gene_ids{k};
        if contains(gene_id,protein)
            MW = uniprot.MWs(k);
            coFactorContent = Ptot*abundance*1000/MW; % mmol/gDW
            coFactorInfo.abundances(i) = coFactorContent;
        end
    end
end

% collect compartments specific total cofactor content
coFactorContent = struct();
for i = 1:length(coFactorInfo.type)
    
   type = coFactorInfo.type{i};
   group = coFactorInfo.group{i};
   content = coFactorInfo.abundances(i);
   % create structure field for cofactor
   switch group
       case 'non-mitochondrial'
           comp = 'c';
       case 'mitochondrial'
           comp = 'm';
   end
   
   cofactor = [type '_' comp];
   
   if ~isfield(coFactorContent,cofactor)
       coFactorContent.(cofactor) = content;
   else
      coFactorContent.(cofactor) = coFactorContent.(cofactor) + content;
   end
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end