%% Checking if the model.rev is in fact correct.
% It compares if the model.lb is lower than 0 and
% checks whether model.ub > model.lb
%
% Returns a vector with 0 for one direction reactions 
% and 1 for bi direction reactions
%
%
% 10/04/2018 Athanasios Antonakoudis



%% Main Program
function rev = checkReversability(model)
lb = model.lb;
ub = model.ub;
imax = size(lb,1);
rev = zeros(imax,1);
for i = 1:imax
    if (lb(i) < 0 && ub(i) > 0)
        rev(i) = 1;
    elseif (lb(i) == 0 && ub(i) > 0)
        rev(i) = 0;
    elseif (lb(i) > 0 && ub(i) > 0)
        rev(i) = 0;
    elseif (lb(i) < 0 && ub(i) < 0)
        rev(i) = 0;
    elseif (lb(i) < 0 && ub(i) == 0)
        rev(i) = 0;
    end    
end
end


