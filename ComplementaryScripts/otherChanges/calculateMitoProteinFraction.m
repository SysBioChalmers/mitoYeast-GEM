%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculateMitoProteinContent
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mitoProteinFraction = calculateMitoProteinFraction(model)
% Read file with PaxDb protein abundance
cd ../../ComplementaryData/databases
fid = fopen('proteinAbundancePaxDb.txt');
pdb = textscan(fid,'%s %s %s','Delimiter','\t','HeaderLines',12);
paxDb.genes = strrep(pdb{2},'4932.','');
paxDb.abundances = cellfun(@str2num,pdb{3})*1e-6;
fclose(fid);

% read file with stoichiometry info for mitochondrial complexes
cd ../
fid     = fopen('complexStoichiometry.tsv');
compl   = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
complexInfo.genes    = compl{1};
complexInfo.complex  = compl{3};
complexInfo.stoich   = cellfun(@str2num,compl{4});
complexInfo.integralMembrane = compl{5};
fclose(fid);

% read file with info on mitochondrial proteins
fid = fopen('mitoProteins.tsv');
mgenes = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s', ...
                  'Delimiter','\t','HeaderLines',1);
mitoGenes = mgenes{1};
fclose(fid);

% read file with amino acid info
fid               = fopen('mitoAminoAcidUsage.tsv');
aa_usage          = textscan(fid,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
aa_data.name      = aa_usage{1}; % amino acid name
aa_data.abbr      = aa_usage{2}; % 3 letter abbreviations
aa_data.MWs       = cellfun(@str2num,aa_usage{3}); % molecular weight in polymer
fclose(fid);

% get the total protein content in model
cd ../../ComplementaryScripts/otherChanges/
rxnName = 'protein pseudoreaction';
Ptot = calculateProteinContent(model,rxnName,aa_data.abbr,aa_data.MWs);


% calculate total mitochondrial fraction of the proteome
totalMitoProtFrac = 0;
for i = 1:length(mitoGenes)
    gene      = mitoGenes{i};
    abundance = 0;
    inComplex = strcmp(complexInfo.genes,gene);
    if sum(inComplex) == 1 && complexInfo.integralMembrane{inComplex} == 'Y'
        gene_stoich    = complexInfo.stoich(inComplex); % stoichiometry of subunit
        complex        = complexInfo.complex{inComplex};
        complexPos     = strcmp(complexInfo.complex,complex);
        complexGenes   = complexInfo.genes(complexPos);
        complex_stoich = complexInfo.stoich(complexPos); % stoichiometry of complex
        for j = 1:length(complexGenes)
           complexGene     = complexGenes{j};
           stoichiometry   = complex_stoich(j);
           posInPaxDb      = strcmp(paxDb.genes,complexGene);
           paxDb_abundance = paxDb.abundances(posInPaxDb)/stoichiometry;
           if paxDb_abundance > abundance
              abundance = paxDb_abundance; 
           end
        end
        abundance = abundance * gene_stoich;
    else
        pos = strcmp(paxDb.genes,gene);
        % Handle genes not found in paxDb data
        if sum(pos) == 1
            abundance = paxDb.abundances(pos);
        end
    end
    totalMitoProtFrac = totalMitoProtFrac + abundance;
    if isempty(totalMitoProtFrac)
        disp(num2str(i))
        break
    end
end

mitoProteinFraction = totalMitoProtFrac*Ptot;

end
