%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModificationRxn(model,eznyme)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testModificationRxn(model,enzyme,modRxnName)
enzymePos = ~cellfun(@isempty,strfind(model.mets,enzyme));
enzymeMets = model.mets(enzymePos);
for i = 1:length(enzymeMets)
   enzymeMet = enzymeMets{i};
   mitoComps = {'m','mm','om','ims'};
   enzIDsplit = strsplit(enzymeMet,'_');
   enzComp = enzIDsplit{end};
   mitoProt = sum(strcmp(mitoComps,enzComp));
   apoProt = contains(enzymeMet,'apo_');
   if mitoProt > 0 && apoProt
       rxnID = ['apo_prot_' enzComp '_' enzyme '_exchange'];
       rxnToAdd.rxns         = {rxnID};
       rxnToAdd.rxnNames     = {rxnID};
       rxnToAdd.mets         = {enzymeMet};
       rxnToAdd.stoichCoeffs = 1;
       rxnToAdd.lb           = 0;
       model                 = addRxns(model,rxnToAdd,1);
   elseif ~apoProt
       model = addSinkReactions(model,{enzymeMet},0,1000);
   end     
end

% Maximize flux through modification rxn
rxnName = modRxnName;
rxnPos = strcmp(model.rxns,rxnName);
model.c = zeros(size(model.rxns));
model.c(rxnPos) = 1;
sol = optimizeCbModel(model);
disp(['Maximum flux = ' num2str(sol.obj)])

end