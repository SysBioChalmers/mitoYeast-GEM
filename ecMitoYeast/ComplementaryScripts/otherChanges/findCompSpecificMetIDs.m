%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: IDs = findCompSpecificMetIDs(metNames,comp)
%
% Input:
%	  model: 	  a model structure with fields metComps and comps
%     metNames:   cell array containing metabolite names
%	  compartment: compartment symbol matching model.comps. At the moment
%				   this has to be the same size as metNames, except in the
%				   case of finding a single metabolite in multiple comps
%
% Output:
%     metIDs: 		Cell array containing the metabolite IDs 
%
% Carl Malina 2019-11-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metIDs = findCompSpecificMetIDs(model,metNames,compartment)

if ~iscell(metNames)
	metNames = {metNames};
end

if ~iscell(compartment)
	compartment = {compartment};
end

if length(metNames) == 1
	% find metabolite ID in multiple comps
    metIDs = cell(length(compartment),1);
    metName = metNames{1};
    for i = 1:length(compartment)
		comp = compartment{i};
		compID = find(strcmp(model.comps,comp)); 
		metID = model.mets{strcmp(model.metNames,metName) & model.metComps == compID};
		metIDs{i} = metID;
    end
else
	% find multiple metIDs
    metIDs = cell(length(metNames),1);
    for i = 1:length(metNames)
		metName = metNames{i};
		if length(compartment) > 1
			comp = compartment{i};
		else
			comp = compartment{1};
		end
		compID = find(strcmp(model.comps,comp));
		% find compartment specific metID
		metID = model.mets{strcmp(model.metNames,metName) & model.metComps == compID};
		metIDs{i} = metID;
    end
end

% Output a string in case of single metabolite
if length(metIDs) == 1
    metIDs = cell2mat(metIDs);
end

end 