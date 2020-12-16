%% Explanation
% Searching for all exhange reactions and returns 1 the reactions IDs that
% are EX_...

function isDrain=checkDrainRxns(model)

if isfield(model,'isDrain') % If theres is an isDrain in the model, remove it kalou kakou
    model = rmfield(model,'isDrain');
end

isDrain=false(length(model.rxns),1);

for i=1:length(model.rxns)
    %Full shows zeroes, otherwise just the ones and other coeffs.
    %full(model.S(:,i))
    %nnz(model.S(:,i))
    if nnz(model.S(:,i)) == 1 % Number of Non - Zero Elements, so if nnmz = 1
        isDrain(i,1)=true;       % then there is only one reaction which would drain the cell
    end                       % It creates a matrix with all the EX_.. reactions
end
