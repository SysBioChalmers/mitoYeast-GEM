function model = limitMetaboliteSecretion(model)
%limitMetaboliteSecretion Function that limits the secretion of metabolites
%other then measured exchange fluxes

limitSecretion = {'r_1663'; ... % bicarbonate
                  'r_1549'; ... % 2,3-butanediol
                  'r_1631'; ... % acetaldehyde
                  'r_1810'; ... % glycine
                  'r_4526'; ... % O-phophonatoxy-D-serine
                  'r_1793'; ... % formate
                  'r_1815'; ... % glyoxylate
                  'r_1994'};    % palmitoleate
              
for i = 1:length(limitSecretion)
    excPos           = strcmp(model.rxns,limitSecretion{i});
    model.ub(excPos) = 1e-5;
end

% Also limit r_1048_REV in order to avoid glycolysis bypass
model.ub(strcmp(model.rxns,'arm_r_1048_REV')) = 0; % D-erythrose 4-phosphate + D-fructose 6-phosphate => glyceraldehyde 3-phosphate + sedoheptulose 7-phosphate

end

