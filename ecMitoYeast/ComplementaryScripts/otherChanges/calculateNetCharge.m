%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: netCharge = calculateNetCharge(sequence)
%
% Input:
%     sequence:   A string with the amino acid sequence
%
% Output:
%     netCharge: The net charge of the protein sequence  
%
% Carl Malina 2019-09-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function netCharge = calculateNetCharge(sequence)
% Define table with amino acid properties at physiological pH (7.4)
aas = {'L-alanine'          'A'     0      % A     Alanine
       'L-cysteine'         'C'     0      % C     Cysteine
       'L-aspartate'        'D'     -1     % D     Aspartic acid
       'L-glutamate'        'E'     -1     % E     Glutamic acid
       'L-phenylalanine'	'F'     0      % F     Phenylalanine
       'L-glycine'          'G'     0      % G     Glycine
       'L-histidine'        'H'     0      % H     Histidine
       'L-isoleucine'       'I'     0      % I     Isoleucine
       'L-lysine'           'K'     1      % K     Lysine
       'L-leucine'          'L'     0      % L     Leucine
       'L-methionine'       'M'     0      % M     Methionine
       'L-asparagine'       'N'     0      % N     Asparagine
       'L-proline'          'P'     0      % P     Proline
       'L-glutamine'        'Q'     0      % Q     Glutamine
       'L-arginine'         'R'     1      % R     Arginine
       'L-serine'       	'S'     0      % S     Serine
       'L-tryptophan'       'T'     0      % T     Threonine
       'L-valine'           'V'     0      % V     Valine
       'L-tryptophan'       'W'     0      % W     Tryptophan
       'L-tyrosine'         'Y'     0};    % Y     Tyrosine

sequence = upper(sequence);
netCharge = 0;

for i = 1:length(sequence)
    aa = sequence(i);
    aa_pos = strcmp(aas(:,2),aa);
    if ~sum(aa_pos) == 1
        error('Sequence contains invalid amino acid characters')
    else
        netCharge = netCharge + aas{aa_pos,3};
    end
end


end