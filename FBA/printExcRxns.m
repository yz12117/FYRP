function printExcRxns(model)
    excRxns = findExcRxns(model);
    rxnText = printRxnFormula(model,model.rxns(excRxns),false);
    lb = model.lb(excRxns);
    ub = model.ub(excRxns);
    for k = 1:length(rxnText)
        fprintf('%20s :\tlb = % 5d\tub = % 5d\n',rxnText{k},lb(k),ub(k))
    end
end