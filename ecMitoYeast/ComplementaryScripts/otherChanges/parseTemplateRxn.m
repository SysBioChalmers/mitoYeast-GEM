%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [metList,stochCoeffList,compList,revFlag] = parseTemplateRxn(rxn)
%
% INPUT:
%    formula:   Reaction formula, may contain symbols '+', '=>', '<=>' in
%               addition to stoichiometric coefficients and metabolite names
%               If no stoichiometric coefficient is provided, it is assumed
%               to be = 1.
%               Can also contain characters as stoichiometric coefficients
%
% % OUTPUTS:
%    metList:            Cell array with metabolite names
%    stoichCoeffList:    List of stoichiometric coefficients (cell array in
%                        case of non-numerical stoichCoeffs in rxn string 
%    compList:           Cell array of compartment symbols
%    revFlag:            Indicates whether the reaction is reversible (true) or not (false)
%
% Carl Malina 2019-09-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [metList,stoichCoeffList,compList,revFlag] = parseTemplateRxn(rxn)

characters = strsplit(rxn)';

% deal with metabolites containing a space in name (e.g ferricytochrome c)
% there's probably a better solution but for lack of another one, this
% works
expr = '^\w+\[\w{1,3}\]';
indices = find(~cellfun(@isempty,regexp(characters,expr)));
    
expr2 = '$\W+';
expr3 = '^\[\d{1,1}\w+';
expr4 = '^\d+';

x = 1;
newMets = cell(size(indices));
for b = 1:length(indices)
    nextMet = false;
    newMet = characters{indices(b)};
    if b ~= 1
        x = indices(b-1)+1;
    end
    for a = indices(b)-1:-1:x
        character = characters{a};
        if strcmp(character,'+') || strcmp(character,'=>')
            nextMet = true;
        elseif strcmp(character,'y') || strcmp(character,'x')
            continue
        end
        if isempty(regexp(character,expr2)) && isempty(regexp(character,expr4)) && ~nextMet 
            newMet = [character ' ' newMet];
            characters{a} = [];
        elseif ~isempty(regexp(characters{a},expr3)) && ~nextMet
            newMet = [character ' ' newMet];
            characters{a} = [];
        end
    end
    newMets{b} = newMet;
end
characters(indices) = newMets;
characters = characters(~cellfun(@isempty,characters));

% pattern = '^\w{1,1}\[\w+\]';
% 
% indices = find(~cellfun(@isempty, regexp(characters,pattern)));
% if sum(indices) > 0
%     for i = 1:length(indices)
%         index = indices(i);
%         part1 = characters{index-1};
%         part2 = characters{index};
%         characters{index-1} = [part1 ' ' part2];
%     end
% end
% characters(indices) = [];

% initialize variables
metList = {};
stoichCoeffList = {};
compList = {};
revFlag = true;

% Define the start of a new metabolite
newMetFlag = true;
% Disinguish substrates and products
productFlag = false;
for i = 1:length(characters)
    c = characters{i}; 
    % check if newMet   
    if strcmp(c,'+')
        newMetFlag = true;
    % check if irreversible   
    elseif strcmp(c,'=>')
        revFlag = false;
        newMetFlag = true;
        productFlag = true;
    % check if reversible  
    elseif strcmp(c,'<=>')
        revFlag = true;
        newMetFlag = true;
        productFlag = true;
    else
        % check stoichCoeff
        coeff = str2double(c);
        if (~isnan(coeff))
            if ~productFlag
                coeff = -coeff;
            end
            stoichCoeffList{end+1,1} = num2str(coeff);
            newMetFlag = false;
        elseif ischar(c) && length(c) == 1
            if ~productFlag
                coeff = ['-' c];
            else
                coeff = c;
            end
            stoichCoeffList{end+1,1} = num2str(coeff);
            newMetFlag = false;
        else
            % check metabolite name
            [baseMetName,compSymbol] = extractBaseMetNames(c);
            metList(end+1,1) = baseMetName;
            compList(end+1,1) = compSymbol;
            if newMetFlag
                if ~productFlag
                    stoichCoeffList{end+1,1} = '-1';
                else
                    stoichCoeffList{end+1,1} = '1';
                end
                %newMetFlag = true;
            end       
        end           
    end
end
% convert to double array in case of only numeric coeffs
if sum(cellfun(@ischar,stoichCoeffList)) == 0
    stoichCoeffList = cell2mat(stoichCoeffList);
end
    
end