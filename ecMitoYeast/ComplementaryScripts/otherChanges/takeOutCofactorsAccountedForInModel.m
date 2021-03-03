%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function model = takeOutCofactorsAccountedForInModel(model,geneList)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = takeOutCofactorsAccountedForInModel(model,geneList)

% read file with biomass composition
cd ../../ComplementaryData/physiology/
fid = fopen('biomassComponents_ecModel.txt');
d = textscan(fid,'%s %s %f32 %f32 %s','Delimiter','\t','HeaderLines',1);
data.mets       = d{1};
data.abundances = double(d{3});
data.MWs        = double(d{4});
data.groups     = d{5};
fclose(fid);

% calculate the amount of cofactors included in the model

% Load data from PaxDb used to calculate protein abundance
cd ../../../mitoYeast/ComplementaryData/databases/
fid = fopen('proteinAbundancePaxDb.txt');
pdb = textscan(fid,'%s %s %s','Delimiter','\t','HeaderLines',13);
paxDb.genes = strrep(pdb{2},'4932.','');
paxDb.abundances = cellfun(@str2num,pdb{3})*1e-6;
fclose(fid);

% read file with cofactor information
cd ../proteinInfo/
fid    = fopen('coFactorInfo.tsv');
coFact = textscan(fid,'%s %s %s %s %s %s %s %s','Delimiter','\t', ...
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

cd ../../../ecMitoYeast/ComplementaryScripts/miscellaneous/
current = pwd;
% calculate protein content in model
[~,P,~,~,~,~,~,~,~] = sumBioMass(model,data);

% Create a strcuture for storing cofactor content information
cofactorContent.cofactorIDs = model.mets(model.S(:,strcmp(model.rxnNames,'prosthetic groups pseudoreaction')) < 0);
cofactorContent.abundance = zeros(size(cofactorContent.cofactorIDs));

% Calculate the amount of cofactors accounted for by the proteins
% represented in the model

for i = 1:length(geneList)
    cd(current);
    gene = geneList{i};
    % Check if gene contains cofactors
    genePos = strcmp(coFactorInfo.genes,gene);
    if sum(genePos) > 0
        % Check if gene encoded complex subunit
        if length(find(genePos)) == 1
            inComplex = contains(coFactorInfo.complex(genePos),'complex');
        else
            inComplex = sum(strcmp(coFactorInfo.complex(genePos),'complex'));
        end

        if sum(inComplex) > 0
            % find genes in complex
            if inComplex > 1
                complex = unique(coFactorInfo.complex(genePos));
                complex = complex{1};
            else
                complex = coFactorInfo.complex{genePos};
            end
           
            complexPos    = strcmp(complexInfo.complex,complex);
            complexGenes  = complexInfo.genes(complexPos);
            complexStoich = complexInfo.stoich(complexPos);
            subunitStoich = complexInfo.stoich(strcmp(complexInfo.genes,gene));
           
            % Calculate the abundance using PaxDB dataset, using the most
            % abundant subunit as proxy for abundance of complex and
            % accounting for stoichiometry of subunits
            abundance = 0;
            for j = 1:length(complexGenes)
                complexGene     = complexGenes{j};
                stoichCoeff     = complexStoich(j);
                paxDbPos        = strcmp(paxDb.genes,complexGene);
                paxDb_abundance = paxDb.abundances(paxDbPos)/stoichCoeff;
                if paxDb_abundance > abundance
                    abundance = paxDb_abundance;
                end
            end
            abundance = abundance * subunitStoich;
        else
            paxDbPos  = strcmp(paxDb.genes,gene);
            abundance = paxDb.abundances(paxDbPos);
        end
        MW = model.MWs(i) * 1000; % kDa --> Da (g/mol)
        % Calculate cofactor content in mmol/gDW
        cofactorAbundance = P * abundance * 1000 / MW;
        % Collect info on cofactors in structure
        type = coFactorInfo.type(genePos);
        pool = coFactorInfo.group(genePos);
        for j = 1:length(type)
            cofactor = type{j};
            group    = pool{j};
            switch cofactor
                case 'fe2s2'
                    cofactorName = '[2Fe-2S] iron-sulfur cluster';
                case 'fe4s4'
                    cofactorName = '[4Fe-4S] iron-sulfur cluster';
                case 'fe3s4'
                    cofactorName = '[3Fe-4S] iron-sulfur cluster';
                case 'heme_a'
                    cofactorName = 'heme a';
                case 'lipoate'
                    cofactorName = 'lipoate (protein bound)';
                case 'ferroheme_b'
                    cofactorName = 'ferroheme b';
                case 'biotin'
                    cofactorName = 'biotin (protein bound)';
                case 'siroheme'
                    cofactorName = 'siroheme';
            end
       
            % Update compartment 
            if strcmp(group,'non-mitochondrial')
                if strcmp(cofactorName,'[3Fe-4S] iron-sulfur cluster') || strcmp(cofactorName,'lipoate (protein bound)') || strcmp(cofactorName,'ferroheme b')
                    comp = 'm';
                else
                    comp = 'c';
                end
            elseif strcmp(group,'mitochondrial')
                if strcmp(cofactorName,'heme a') || strcmp(cofactorName,'siroheme') || strcmp(cofactorName,'biotin (protein bound)')
                    comp = 'c';
                else
                    comp = 'm';
                end
            end
            cd ../otherChanges/
            metID = findCompSpecificMetIDs(model,cofactorName,comp);
       
            % Update abundance of cofactor in cofactor structure
            cofactorPos = strcmp(cofactorContent.cofactorIDs,metID);
            cofactorContent.abundance(cofactorPos) = cofactorContent.abundance(cofactorPos) + cofactorAbundance;
        end
    end
end

% Remove content covered by enzymes included in model
for i = 1:length(cofactorContent.cofactorIDs)
    cofactorID   = cofactorContent.cofactorIDs{i};
    cofactorInd  = strcmp(model.mets,cofactorID);
    pseudoRxnInd = strcmp(model.rxnNames,'prosthetic groups pseudoreaction');
    coeff        = model.S(cofactorInd,pseudoRxnInd);
    newCoeff     = coeff + cofactorContent.abundance(i);
    % prevent positive stoichiometric coefficients
    if newCoeff > 0
        newCoeff = 0;
    end
    model.S(cofactorInd,pseudoRxnInd) = newCoeff;
    % Update biomass composition structure with new abundances
    data.abundances(strcmp(data.mets,cofactorID)) = -newCoeff;
end

% sum biomass
cd(current);
[X,~,C,~,~,~,~,~,~] = sumBioMass(model,data);

delta = X - 1;
fC    = (C - delta)/C;
cd ../../../mitoYeast/ComplementaryScripts/otherChanges/
model = rescalePseudoReaction(model,'carbohydrate',fC);
cd(current);

sumBioMass(model,data);

end



