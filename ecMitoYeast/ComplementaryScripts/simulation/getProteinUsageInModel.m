function protUsage = getProteinUsageInModel(model,uniprotIDs,fluxes,scaledUsage,protDataIncorporated,writeFile,name)
% getProteinUsageInModel
%
%	Function that takes an ecModel, a list of uniprot IDs, and an FBA solution vector
%	and calculates the usage of the proteins specified (uniprotIDs) in a mass-wise way.
%	If desired, results are written to a .txt file
%
%	model 					ecModel structure used to generate the FBA solution
%	uniprotIDs 				Cell array with uniprot IDs for the enzymes of interest
%	fluxes 					FBA solution vector
%	protDataIncorporated 	Logical, whether proteomics data has been used to 
%							constrain model. (opt, default false)
%	writeFile 				Logical, whether results should be written to a file
% 							in location createMitoYeastGEM/Results (opt, default false)
%
%	Usage: protUsage = getProteinUsageInModel(model,uniprotIDs,fluxes,scaledUsage,protDataIncorporated,writeFile,name)
%
% 	Carl Malina 	Last edited 2020-10-07
%

current = pwd;

if ~iscell(uniprotIDs)
	uniprotIDs = {uniprotIDs};
end

if nargin < 4
    scaledUsage = false;
elseif nargin < 5
	writeFile = false;
	protDataIncorporated = false;
elseif nargin < 6
	protDataIncorporated = false;
end

if protDataIncorporated
	usages = find(contains(model.rxns,'prot_'));
	% remove prot_pool_exchange
	usages = usages(1:end-1);
	% remove translocation and modification
	usageRxns = model.rxns(usages);
	usages    = usages(~contains(usageRxns,'translocation') && ~contains(usageRxns,'modification'));
	usageRxns = usageRxns(usages);
else
	usages = find(contains(model.rxns,'draw_prot_'));
	usageRxns = model.rxns(usages);
end
gRate      = fluxes(strcmpi(model.rxnNames,'growth'));
enzUsages  = fluxes(usages);
cd ../miscellaneous/
totalUsage = sumProtein(model); % Get total protein content %sum((enzUsages.*model.MWs)./gRate);
cd(current)
if scaledUsage
    scaledUsageInd = find(contains(usageRxns,'_scaled'));
    enzUsages = enzUsages(scaledUsageInd);
    usageRxns = usageRxns(scaledUsageInd);
end
massFraction = zeros(size(uniprotIDs));
for i = 1:length(uniprotIDs)
	prot 	 		= uniprotIDs{i};
	protPos  		= strcmp(model.enzymes,prot);
	MW 		 		= model.MWs(protPos);
	usagePos 		= contains(usageRxns,prot);
	usage    		= enzUsages(usagePos);
	enzUsage 		= (usage*MW)/gRate;
    % If protein usage has been scaled, scale back
    if scaledUsage
        enzUsage  = enzUsage/1000;
    end
	enzUsage 		= enzUsage/totalUsage;
	massFraction(i) = enzUsage;
end
colNames 	 = {'uniprotID','mass_fraction_of_proteome'};
data 	     = [uniprotIDs,num2cell(massFraction)];
fileData 	 = cell2table(data,'VariableNames',colNames);
if writeFile
	writetable(fileData,['Results/' name '.txt'],'Delimiter','tab');
end
protUsage = fileData;
end