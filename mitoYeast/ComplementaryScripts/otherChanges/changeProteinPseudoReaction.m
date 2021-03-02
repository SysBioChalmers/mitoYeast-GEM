%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addMitoProteinRxn(model)
%
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = changeProteinPseudoReaction(model)
% Read data for mitochondrial amino acid usage
cd ../../ComplementaryData/Physiology/
fid               = fopen('mitoAminoAcidUsage.tsv');
aa_usage          = textscan(fid,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
aa_data.name      = aa_usage{1}; % amino acid name
aa_data.abbr      = aa_usage{2}; % 3 letter abbreviations
aa_data.MWs       = cellfun(@str2num,aa_usage{3}); % molecular weight in polymer
aa_data.fraction = cellfun(@str2num,aa_usage{6}); % g amino acid/g protein
fclose(fid);

% read file with data for mitochondrially synthesized proteins
cd ../proteinInfo/
fid = fopen('mitoEncodedProteins.tsv');
mp = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
mitoProt.genes   = mp{3};     % 
mitoProt.complex = mp{5};
fclose(fid);

% read file with stoichiometry info for RC complexes and ATP synthase
fid     = fopen('complexStoichiometry.tsv');
compl   = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
complexInfo.genes    = compl{1};
complexInfo.complex  = compl{3};
complexInfo.stoich   = cellfun(@str2num,compl{4});
fclose(fid);

cd ../databases/
% read PaxDb abundance data
fid = fopen('proteinAbundancePaxDb.txt');
pdb = textscan(fid,'%s %s %s','Delimiter','\t','HeaderLines',13);
paxDb.genes = strrep(pdb{2},'4932.','');
paxDb.abundances = cellfun(@str2num,pdb{3})*1e-6;
fclose(fid);

cd ../../ComplementaryScripts/otherChanges/
rxnName = 'protein pseudoreaction';
Ptot = calculateProteinContent(model,rxnName,aa_data.abbr,aa_data.MWs);

disp(['Protein content in model is ' num2str(Ptot) ' g/gDW']);
 
% calculate the abundance of mitochondrial proteins
mitoProt.abundance = zeros(size(mitoProt.genes));
for i = 1:length(mitoProt.genes)
    protein = mitoProt.genes{i};
    complex = mitoProt.complex{i};
    if contains(complex,'complex')
        % find genes in complex
        protein_stoich = complexInfo.stoich(strcmp(complexInfo.genes,protein)); % Stoichiometry of actual subunit
        genePos        = strcmp(complexInfo.complex,complex);
        genes          = complexInfo.genes(genePos);
        subunit_stoich = complexInfo.stoich(genePos);
        abundance      = 0;
        for j = 1:length(genes)
            gene = genes{j};
            stoichCoeff = subunit_stoich(j);
            pos = strcmp(paxDb.genes,gene);
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
    mitoProt.abundance(i) = abundance;
end

% calculate the content of mitochondrially synthesized proteins
tot_mitoProt = sum(mitoProt.abundance.*Ptot); % g/gDW

% get metabolite IDs for amino acids and aminoacyl-tRNAs
aa_data.metID_m  = cell(size(aa_data.name));
aa_data.metID_c  = cell(size(aa_data.name));
aa_data.tRNA_IDs = cell(size(aa_data.abbr));
for i = 1:length(aa_data.name)
    metName_m = [aa_data.name{i} ' [mitochondrion]'];
    metName_c = [aa_data.name{i} ' [cytoplasm]'];
    metID_m = model.mets{strcmp(model.metNames,metName_m)};
    metID_c = model.mets{strcmp(model.metNames,metName_c)};
    aa_data.metID_m{i,1} = metID_m;
    aa_data.metID_c{i,1} = metID_c;
    aa_abbr = aa_data.abbr{i};
    for j = 1:length(model.metNames)
        metName = model.metNames{j};
        if startsWith(metName,[aa_abbr '-tRNA']) && endsWith(metName,'[cytoplasm]')
           aa_data.tRNA_IDs{i,1} = model.mets{j}; 
        end
    end
end

% calculate the abundance (mmol/gDW) of amino acids in mito synthesized
% proteins
aa_data.abundance = zeros(size(aa_data.fraction));
for i = 1:length(aa_data.fraction)
   fraction = aa_data.fraction(i);
   abundance = fraction * tot_mitoProt * 1000 / aa_data.MWs(i);
   aa_data.abundance(i) = abundance;
end

% create a pseudoreaction for protein prosthetic groups
model = addProstheticGroupsToBiomass(model);

% Modify protein pseudoreaction

% get the ID all mitochondrial aminoacyl-tRNAs (excluding fMet-tRNA(fMet),
% L-glytamyl-tRNA(Gln) and L-Aspartyl-tRNA(Asn))
aa_tRNAs_name = {};
tRNAs_name    = {};
for i = 1:length(model.metNames)
    met = model.metNames{i};
    if contains(met,'-tRNA') && contains(met,'[mitochondrion]') && ~contains(met,'fMet') && ~contains(met,'L-')
        aa_tRNAs_name = [aa_tRNAs_name; met];
    elseif contains(met,'tRNA(') && contains(met,'[mitochondrion]') && ~contains(met,'-tRNA')
        tRNAs_name = [tRNAs_name;met];
    end
end
aa_tRNAs_name = sort(aa_tRNAs_name);
tRNAs_name    = sort(tRNAs_name);

% subtract mitochondrial amino acid amount from protein pseudoreaction
%cd ../
proteinRxnPos = strcmp(model.rxnNames,'protein pseudoreaction');
printRxnFormula(model,model.rxns(proteinRxnPos),true,true,true);
for i = 1:length(aa_tRNAs_name)
    metName = aa_tRNAs_name(i);
    [baseMetName,~] = extractBaseMetNames(metName);
    baseMetName = char(baseMetName);
    aa_abbr = baseMetName(1:3);
    mitoAbundance = aa_data.abundance(strcmp(aa_data.abbr,aa_abbr));
    aa_tRNApos = strcmp(model.metNames,[baseMetName '[cytoplasm]']);
    oldStoichCoeff = abs(model.S(aa_tRNApos,proteinRxnPos));
    newStoichCoeff = abs(oldStoichCoeff)-mitoAbundance;
    metPos = abs(model.S(:,proteinRxnPos)) == oldStoichCoeff;
    metIndex = find(metPos);
    if sum(metPos) ~= 2
        error('did not find tRNA(aa) pair')
    else
        for j = 1:2
           posInS = metIndex(j);
           coeff = model.S(posInS,proteinRxnPos);
           if coeff < 0
               model.S(posInS,proteinRxnPos) = -newStoichCoeff;
           else
               model.S(posInS,proteinRxnPos) = newStoichCoeff;
           end
               
        end
    end
end

printRxnFormula(model,model.rxns(proteinRxnPos),true,true,true);

% change reaction bounds for r_1874 and r_1935 to allow for uptake of
% L-alanine and L-methionine in mitochondria
model.lb(strcmp(model.rxns,'r_1874')) = -1000;
model.lb(strcmp(model.rxns,'r_1935')) = -1000;

% add mitochondrial tRNA pairs to protein pseudoreaction
for i = 1:length(aa_tRNAs_name)
    aa_tRNA_name = aa_tRNAs_name(i);
    [baseMetName,~] = extractBaseMetNames(aa_tRNA_name);
    baseMetName = char(baseMetName);
    aa_abbr = baseMetName(1:3);
    tRNA = ['tRNA(' aa_abbr ') [mitochondrion]'];
    tRNAid       = model.mets{strcmp(model.metNames,tRNA)};
    tRNA_pos    = strcmp(model.mets,tRNAid);
    aa_tRNA     = model.mets{strcmp(model.metNames,aa_tRNA_name)};
    aa_tRNA_pos = strcmp(model.mets,aa_tRNA);
    abundance = aa_data.abundance(i);
    % change coefficient in S-matrix for tRNA pair
    model.S(aa_tRNA_pos,proteinRxnPos) = -abundance;
    model.S(tRNA_pos,proteinRxnPos)    = abundance;  
end

printRxnFormula(model,model.rxns(proteinRxnPos),true,true,true);

cd ../

end