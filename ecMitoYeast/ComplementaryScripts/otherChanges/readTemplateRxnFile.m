%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: templateRxns = readTemplateRxnFile(filename)
%
% Input:
%
%   filename    Name of a .tsv file containing info about template
%               reactions. The file should for columns in the following 
%               order: the rxnAbbreviation, reaction formula 
%               (e.g 'A + 2 B => C'), grRules and references
%
% Output: 
%   templateRxns    A structure containing the fields: rxnAbbreviation
%                                                      rxnEquation 
%                                                      grRules
%                                                      references
%
%
% Carl Malina 2019-09-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function templateRxns = readTemplateRxnFile(filename)

% read file
fid = fopen(filename);
fileData = textscan(fid,'%s%s%s%s','Delimiter','\t','HeaderLines',1);
fclose(fid);

% extract data
templateRxns.rxnAbbreviation = fileData{1};
templateRxns.rxnEquation = fileData{2};
templateRxns.grRules = fileData{3};
templateRxns.references = fileData{4};

end