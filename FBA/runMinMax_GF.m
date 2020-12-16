% Function to get the maximum and minimum fluxes for each of the reactions
% of the model. 
% model = model we are using.
% end_rxn_index = No of the reaction we want the FVA algorithm to finish
% the iteration.
function [minmax] = runMinMax_GF(model,end_rxn_index)
[~, num_rxns,~] = size(model.S);
start_rxn_index=1;

if nargin < 2
    end_rxn_index=num_rxns;
end

count = 1;

max = zeros(end_rxn_index - start_rxn_index+1,1);
min = max;

for i=start_rxn_index:end_rxn_index 
    model.c=zeros(size(model.S,2),1);
    model.c(i)=1;
    sol = optimizeCbModel(model,'max'); 
    max(count,1) = sol.x(i);

    sol = optimizeCbModel(model,'min');
    min(count,1) = sol.x(i);

    count = count+1;
    
    if mod(count,1000) == 0
        sprintf("%d\n",count)
    end
    
end

minmax(:,1)=min;
minmax(:,2)=max;