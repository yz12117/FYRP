function [model,upt] = PseudoRxns(model,drains,rxnRev,flag,uptake,ir2r)
% flag = 1 -> Reversible model, has model.rev
% flag = 2 -> Irreversible model, has model.matchRev

% rxns = model.rxns;
num_rxns = size(model.rxns,1);
j = 1;
if flag == 1 
    for i=1:num_rxns
        if drains(i) 
            n = model.rxns(i);
            if ismember(n,uptake)
                rxnNEW(i) = {sprintf('Uptake%i',i)};
                upt(j) = i;
                j = j + 1;
            else
                rxnNEW(i) = {sprintf('E%i',i)};
            end
        else
            rev = model.rev;
            if rev(i) == 0
                rxnNEW(i) = {sprintf('R%d_I',i)};
            else
                rxnNEW(i) = {sprintf('R%d_R',i)};
            end
            %%%%% SET ATP DEMAND AS AN UPTAKE REACTION TO REJECT ANY GENE 
            %%%%% REGULATION POSSIBILITY
            if ismember(model.rxns(i),uptake)
                rxnNEW(i) = {sprintf('Uptake%i',i)};
                upt(j) = i;
                j = j + 1;
            end
        end
    end
else
    for i=1:num_rxns
        if drains(i) 
            n = model.rxns(i);
            if ismember(n,uptake)
                rxnNEW(i) = {sprintf('Uptake%i',i)};
                upt(j) = i;
                j = j + 1;
            else
                if (rxnRev(i) > i)
                    rxnNEW(i) = {sprintf('E%d_f',ir2r(i))};
                elseif (rxnRev(i) == 0)
                    rxnNEW(i) = {sprintf('E%d',ir2r(i))};
                else 
                    rxnNEW(i) = {sprintf('E%d_b',ir2r(i))};
                end
            end
        else
            if rxnRev(i) == 0
                rxnNEW(i) = {sprintf('R%d',i)};
            elseif (rxnRev(i) > i )
                rxnNEW(i) = {sprintf('R%d_f',i)};
            else
                rxnNEW(i) = {sprintf('R%d_b',ir2r(i))};
            end
            
            %%%%% SET ATP DEMAND AS AN UPTAKE REACTION TO REJECT ANY GENE
            %%%%% REGULATION POSSIBILITY
            if ismember(model.rxns(i),uptake)
                rxnNEW(i) = {sprintf('Uptake%i',i)};
                upt(j) = i;
                j = j + 1;
            end
            
        end
    end
end
num_gene = size(model.genes,1);
for i=1:num_gene
    geneNEW(i,1) = {sprintf('G%d',i)};
end
tempg = model.genes;
model.genes = geneNEW;
model.genesName = tempg;
tempr = model.rxns;
model.rxns = rxnNEW';
model.reactions = tempr;

% upt = 'maytrix'

end