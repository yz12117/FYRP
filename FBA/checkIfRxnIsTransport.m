function isTrans=checkIfRxnIsTransport(substrates,products)
model = load ('iCHOv1.mat'); % load iCHO1766 model
model = model.iCHOv1;
compartments = ['c';'l';'m';'n';'r';'x'];
substrates_compartment = [];
products_compartment = [];
for i = 1:length(compartments)
    substrates_compartment = horzcat(substrates_compartment,append(append(substrates,'_'),compartments(i)));
    products_compartment = horzcat(products_compartment,append(append(products,'_'),compartments(i)));
end
[rxnList, rxnFormulaList] = findRxnsFromMets(model,[substrates_compartment;products_compartment]);
[transRxns, nonTransRxns, transRxnsBool]  = findTransRxns(model);
for i = 1:length(rxnList)
    if 
        
    end
end
end