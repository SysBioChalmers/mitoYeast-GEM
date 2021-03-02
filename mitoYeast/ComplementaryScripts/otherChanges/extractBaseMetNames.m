%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [baseMetNames, compNames] = extractBaseMetNames(metNames)
%
% Carl Malina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [baseMetNames, compNames] = extractBaseMetNames(metNames)

for i = 1:length(metNames)
    % Separate metabolite ID and compartment symbol
    metName         = metNames{i};
    [tokens, tmp]   = regexp(metName,'(.+)\[(.+)\]','tokens','match');
    baseMetName     = tokens{1}{1};
    compName        = tokens{1}{2};
    baseMetNames{i} = baseMetName;
    compNames{i}  = compName;
end

baseMetNames = columnVector(baseMetNames);
compNames    = columnVector(compNames);

end