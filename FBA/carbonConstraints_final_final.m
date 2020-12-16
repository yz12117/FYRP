function [carbonConst] = carbonConstraints_final_final(model,bnds,condition,carbonCount,cofPairs,minmax,runFinalMinMax,limEX,totalCarbon,skipRxns,direction)
[~,num_rxns] = size(model.S);

%% This script is designed to constrain the solution space based on the 
%   number of carbon atoms made available to the model by the exchange
%   reactions. This is done in order to avoid physiologically meaningless
%   futile cycles.
%
%   CREATED 2017-09-30 by Maximilian Lularevic
%   Edited by Athanasios Antonakoudis 10/12/2018
%
% INPUT
%   model       =   cobra model (must contain model.metFormulas)
%
%   bnds        =   bounds or constraints as cell or double. for double 
%                   first row corresponds to rxnIDs, second row to lower 
%                   and the third row to upper bnds
%
% OPTIONAL
%
%   condition   =   string containing information as to which condition is 
%                   run (e.g. exponential or stationary)
%
%   carbonCount =   carbonCount Structure as received form function 
%                   findMetCarbon.m (needs to contain carbonCount.carbonAtoms field)
%
%   cofParis    =   cofParis structure as received from function findCofactorPairs.m
%
%   minmax      =   FVA for model under experimental constraints. (needs 
%                   to be processed with fixMinMax.m)
%
%   runFinalMM  =   logical true or false
%
%   limEX       =   value to be used to constrain drain rxns. This is
%                   supposed to be a single carbon constraint (e.g. left
%                   over carbon can be set to CO2 rate which is a single
%                   carbon. If flux for CO2 (closing the carbon in excel)
%                   is 0.641 mmol/gDCW*h this would be the rate to fed into
%                   the limEX variable).
%
% OUTPUT
%   carbonConst =   structure containing the following fields
%                       - Title
%                       - Date and Time of job
%                       - constraintsName (given: e.g. 'LacProd')
%                       - model.description
%                       - total carbon consumed
%               *** optional ***
%                       - actual carbon used for constraining (only when
%                         passed into function)
%               *** optional ***
%                       - report
%                           - cell with the lenght of model.rxns which 
%                             states whther reaction was balanced or unbalanced  
%                           - number of balanced rxns
%                           - number of unbalanced rxns
%                           - number of set rxns (aka constraints)
%                           - number of rxn with no flux
%                       - bndAdjusted
%                           - rawData (logical 2 colum matrix)
%                           - number of lower bounds adjusted by cc
%                             algorithm
%                           - number of upper bnds adjusted by cc algorithm
%                       - carbonCount
%                           - dateAndTime
%                           - modelDescription
%                           - metName
%                           - metFormula
%                           - carbon atoms for each metabolite
%                           - table with all the above information compiled
%                       - cofactorPairs (cell with mets and met IDs and
%                         occurance of these mets with corresponding
%                         rxnIDs)
%                       - original minmax (either the one passed on to 
%                         function or herein calculated one)
%                       - newMinMax pased on carbon balance
%                       - solSpaceReduction size of new ssolution space
%                         compared to old one (calculatons based on minmax
%                         range comparison)
%                       - solutionNew (optimizeCbModel output with new bnds)
%                       - EXrxnsAdjusted -> title
%                       - limEX (single carbon limit passed into fucntion
%                         to constrain the Exchange reactions output)
%                       - number of EX bnds adjusted
%                       - solSpaceReductionEX in % compared to original
%                       - EXspaceReduction (how much has the space shrunk
%                         just comparing old EX ranges with new ones
%                       - post New MinMax -> title
%                       - postNewMinMax (new minmax run with runMinMax_GF
%                         and cc bnds)
%                       - solSpacereductionPostCC
%                       - solutionPostNew


%% check if minmax was passed into function and, run minmax if it hasn't
if exist('minmax','var')
    if isempty(minmax)
        max = 1000;
        model.rev=checkReversability(model);
        model = resetModel(model,-max,max,model.rev);
        [model,~] = setExperimentBounds_final(model,condition,max);
        minmax = runMinMax_GF(model);
    end
else
    max = 1000;
    model.rev=checkReversability(model);
    model = resetModel(model,-max,max,model.rev);
    [lb,ub,model] = setExperimentBounds_final(model,max,time);
    for i = 1:length(model.rxns)
        model.ub(i) = ub(i);
        model.lb(i) = lb(i);
    end
%     [lb,ub] = Susie_Bounds(1);
%     model = set_ub(model,ub);
%     model = set_lb(model,lb);
    minmax = runMinMax_GF(model); % runs the simple MinMax algorithm for the model;
%     minmax = fixMinMax(minmax);   % fixes the MinMax, if any lb>ub; usually there aren't any
end
%% check if condition was passed into function
if exist('condition','var')
    if isempty(condition)
        condition = 'n/a';
    end
else
    condition = 'n/a';
end

%% check if there are any rxns to be skipped
if exist('skipRxns','var')
    if isempty(skipRxns)
        skipRxns = [];
    else
        [s,z] = size(skipRxns);
        if z > s
            skipRxns = skipRxns';
        end
    end
else
    skipRxns = [];
end

%% extract lower and upper bnds from bnds variable
if iscell(bnds)
    k = 0;
    [num_vars,~] = size(bnds);
    for i=1:num_vars
        varName = bnds{i,1};
        var_index = find(ismember(model.rxns,varName));
        if isempty(var_index)
            fprintf('%s not found\n',bnds{i,1});
        else
            k = k+1;
            varID(k,1) = var_index;
            varMin(k,1) = bnds{i,2};
            varMax(k,1) = bnds{i,3};
        end
    end
elseif isnumeric(bnds)
    varID = bnds(:,1);
    varMin = bnds(:,2);
    varMax = bnds(:,3);
end
varID_tmp = varID;

%% find corresponding metIDs from varIDs
% first remove all bnd rxn that are not true drains (e.g. biomass, IgG formation)
isDrain = find(contains(model.rxns,'EX_'));
noExRxns = find(~ismember(varID,isDrain));
% adjust the previously retrieved variables
varID(noExRxns) = [];
varMin(noExRxns) = [];
varMax(noExRxns) = [];
% get all drains without the set bnds
isDrainNew = isDrain(~ismember(isDrain,varID));
% retrieve the actual metabolite IDs for carbon calculations
metID_tmp = find(model.S(:,varID));

[metID,~] = ind2sub(size(model.S),metID_tmp); % index to coordinates
metID = unique(metID,'stable');
for i = 1:length(metID)
    if ismember(model.mets(metID(i)),'fluxMeasure')
        metID(i) = [];
        break
    end
end
%% check if number of carbon atoms per metabolite was passed into function
if exist('carbonCount','var')
    if isempty(carbonCount)
        carbonCount = findMetCarbon(model);
    end
else
    carbonCount = findMetCarbon(model);
end
% extract the carbon vector (length of model.mets)
metCarb = carbonCount.carbonAtoms;
% metCarb = [40;40;40;40;45;45;45;45;50;50;50;50;5;5;20;5;5;6;6;20;6;9;20;9;20;9;26;9;21;9;21;31;21;33;21;37;37;39;39;26;39;26;21;39;21;47;21;21;26;35;20;10;3;20;3;20;20;7;3;7;20;3;3;7;3;3;3;20;3;21;21;3;5;3;9;20;20;20;4;9;20;20;26;10;26;4;10;25;26;29;3;3;20;12;4;29;4;20;20;9;21;6;9;6;6;21;21;6;21;3;4;2;9;8;8;8;8;8;8;8;8;9;9;9;9;9;9;10;10;10;10;9;5;5;4;4;4;4;4;4;10;11;11;19;19;17;17;6;57;58;57;18;31;31;7;25;25;25;25;31;31;31;33;33;37;37;23;23;25;26;4;21;21;20;28;27;28;27;27;27;7;7;7;5;5;5;5;6;5;5;27;9;9;6;3;28;28;28;27;27;27;6;6;13;28;28;28;27;27;27;3;2;4;4;9;9;57;57;58;58;5;5;4;4;5;8;7;7;24;24;3;35;35;9;4;9;9;7;9;4;9;4;9;5;4;5;4;9;26;27;4;6;9;4;20;6;20;4;5;6;10;19;10;19;8;19;5;19;11;5;0;25;7;25;28;25;0;10;8;19;19;9;8;10;8;6;5;5;5;9;24;8;24;39;39;8;24;8;78;78;8;39;39;5;6;4;5;9;6;5;20;10;10;5;10;20;6;20;20;24;24;8;24;6;6;17;24;24;6;7;17;10;6;10;7;24;6;24;10;6;10;31;6;31;11;6;11;14;12;14;23;11;12;10;21;29;6;9;28;6;16;6;6;19;6;6;11;19;19;6;9;11;9;20;9;9;9;20;20;20;18;20;12;20;20;12;20;5;8;20;8;20;5;9;6;6;39;20;20;20;20;20;4;20;20;39;39;12;39;20;11;19;6;10;19;6;44;19;30;18;44;30;44;18;30;18;44;13;22;19;19;19;5;18;19;5;19;18;23;18;6;23;18;23;19;6;20;6;5;20;5;6;20;36;20;20;20;9;20;20;20;10;9;20;44;10;20;10;44;18;44;20;18;20;44;18;27;20;18;20;20;20;27;49;49;20;49;20;16;20;63;49;20;49;16;20;15;49;20;49;20;20;20;49;4;21;4;21;4;11;27;16;20;16;20;15;20;15;8;15;33;10;33;10;39;10;10;43;3;43;10;35;10;41;37;45;10;20;7;7;45;11;24;7;35;7;24;20;24;7;17;7;24;27;29;7;7;27;7;29;27;27;7;27;27;20;4;27;8;9;19;27;35;1;27;1;37;33;1;27;1;1;39;27;37;22;35;27;20;7;5;27;18;20;27;20;18;27;18;27;18;20;27;20;20;27;20;20;27;20;20;27;20;20;27;20;20;27;11;20;27;20;7;20;20;5;20;20;4;20;20;20;20;11;9;5;16;20;5;20;16;6;14;20;6;20;11;9;11;20;45;20;11;45;9;20;45;20;10;20;9;9;9;10;20;18;37;33;18;9;18;11;19;15;19;17;20;19;20;19;9;12;9;11;10;27;28;16;15;27;10;27;27;27;54;54;60;60;27;16;27;12;10;8;8;27;27;20;27;27;46;20;57;27;20;19;57;21;20;21;21;26;21;15;22;15;43;43;28;41;15;20;15;45;20;43;15;27;15;41;27;39;21;6;21;39;16;39;21;11;21;45;26;21;26;14;20;14;14;14;28;17;25;20;25;23;7;20;20;23;20;6;23;7;23;18;20;74;20;20;20;44;0;44;20;0;20;35;35;20;0;20;69;0;69;20;38;20;38;20;29;20;29;20;30;60;30;60;30;28;20;28;20;26;20;28;41;28;20;41;29;28;41;28;29;18;41;18;28;45;28;19;43;19;28;28;45;18;6;21;43;18;45;18;21;43;9;41;18;45;18;18;18;41;41;45;40;20;18;45;18;20;43;20;18;43;18;18;20;20;45;20;18;45;18;20;43;41;18;45;18;41;18;45;18;41;41;41;18;41;18;41;41;41;41;8;41;8;41;43;43;8;43;8;43;43;43;9;43;20;43;43;43;18;45;43;45;43;18;43;43;45;18;45;18;45;43;45;18;45;18;45;45;45;41;43;41;41;41;45;39;39;45;39;41;45;41;39;39;27;39;39;27;48;27;20;39;20;39;20;20;20;41;20;41;9;41;20;41;9;0;28;41;28;0;28;41;0;20;28;0;28;20;8;20;28;8;28;20;11;41;41;28;28;28;41;8;41;41;28;28;28;28;8;28;41;8;39;9;39;39;8;41;19;41;11;19;27;41;41;41;20;18;20;20;41;18;20;20;19;20;20;19;20;20;18;20;20;20;20;18;20;20;7;6;20;10;20;20;20;20;20;41;20;41;20;20;30;41;20;39;30;20;30;39;20;39;30;20;9;18;11;18;18;9;12;41;10;41;0;41;20;41;20;41;30;20;30;41;20;30;41;20;28;20;41;28;41;20;18;20;20;20;20;20;30;20;28;8;23;30;30;23;23;1;30;20;23;20;23;23;23;41;23;41;10;10;23;10;41;42;41;41;10;42;18;41;42;41;22;42;22;39;42;39;22;26;22;39;26;18;22;42;22;18;42;18;56;41;22;56;22;41;41;41;22;22;41;22;41;50;41;22;50;22;42;44;52;44;22;52;22;44;12;44;22;12;21;40;42;22;40;22;42;32;11;22;32;22;28;26;42;26;42;42;22;22;22;40;20;19;20;22;22;20;22;20;11;20;22;18;22;20;30;29;30;22;29;22;30;29;20;18;29;9;20;29;20;20;29;20;30;29;30;30;28;29;28;30;29;30;28;22;28;30;22;30;28;30;20;22;20;21;20;21;20;21;20;20;12;20;20;9;20;20;9;20;20;20;20;20;20;20;30;9;25;29;9;29;29;20;11;20;29;7;20;20;9;20;20;20;20;0;20;16;0;40;0;0;20;20;20;0;20;0;0;20;9;20;0;0;0;28;0;20;0;0;12;0;20;0;20;0;0;0;20;13;20;24;0;24;20;20;0;0;0;0;0;20;4;20;43;20;20;0;43;0;0;43;0;45;0;45;20;45;20;20;45;20;0;0;45;0;20;45;20;0;47;0;20;47;20;0;0;0;20;0;20;0;20;0;0;0;0;0;20;20;0;0;12;0;0;0;0;0;3;16;0;0;0;3;0;6;6;7;0;0;0;6;7;10;6;0;6;5;0;5;6;0;9;6;11;10;12;9;0;3;0;16;3;0;16;0;6;4;3;0;6;4;0;4;6;0;11;6;11;0;0;10;20;20;20;6;10;15;6;12;15;12;9;24;9;9;24;9;12;24;12;13;13;26;6;6;26;6;6;26;6;12;2;26;16;6;0;12;1;0;3;9;0;12;3;0;22;48;0;0;22;0;26;22;13;26;22;22;22;13;26;2;22;22;26;1;0;0;1;2;0;2;22;0;22;2;0;2;8;0;8;4;20;4;1;20;1;4;20;4;22;20;22;4;20;2;8;20;8;2;20;2;1;1;20;1;2;20;2;23;1;20;22;23;20;23;22;22;20;8;23;20;23;8;20;0;23;3;20;0;20;3;20;20;3;0;9;0;3;7;7;3;53;11;5;30;53;53;9;5;9;8;5;8;8;14;8;4;8;59;79;79;0;4;0;25;79;93;10;93;25;10;25;93;10;73;3;10;3;59;8;10;8;5;72;5;8;72;8;72;8;21;8;10;10;53;53;2;10;4;53;10;67;10;53;67;4;10;4;59;10;39;6;7;6;19;4;27;9;7;27;7;16;7;1;16;6;15;8;15;8;8;15;15;15;8;15;73;15;15;15;73;22;73;14;22;10;10;73;73;73;10;22;10;67;67;43;67;10;43;10;73;43;73;10;10;73;29;11;6;29;19;11;9;11;19;9;19;11;10;56;19;17;19;70;17;70;70;19;5;19;25;56;19;62;7;25;25;25;9;25;9;7;25;7;9;14;11;14;11;14;20;14;10;10;4;14;4;0;9;0;4;9;6;9;6;4;9;4;16;9;16;8;10;10;6;10;10;10;10;11;20;10;10;6;10;8;29;21;10;21;29;10;29;6;5;10;29;5;29;5;29;10;5;62;25;5;25;62;5;5;3;5;5;3;5;3;20;15;3;20;3;20;5;3;5;20;3;1200000;3;41;3;9;3;4;41;3;41;4;41;3;27;4;3;45;27;3;20;45;6;45;20;39;6;39;39;33;33;33;33;16;16;16;16;20;20;25;25;10;10;20;20;20;41;41;41;41;27;27;27;27;6;6;6;6;10;6;6;6;4;4;4;4;10;10;21;4;4;21;4;5;28;28;7;37;37;37;5;37;7;7;5;31;37;5;37;31;37;24;39;31;24;31;39;24;39;31;24;31;39;24;41;17;45;17;43;45;43;17;32;45;35;43;45;43;33;48;33;43;48;29;33;48;33;29;48;31;33;48;33;33;48;33;31;5;11;11;34;37;5;35;35;35;11;11;11;11;5;26;35;35;12;5;35;27;35;12;27;27;27;35;27;37;13;27;13;37;37;27;13;27;28;27;29;27;29;15;24;15;14;24;15;6;15;21;6;30;15;15;6;15;9;26;15;6;15;9;36;6;15;36;15;6;15;0;15;0;30;0;36;0;4;0;0;4;13;4;0;22;0;19;22;19;19;4;43;10;19;43;19;10;43;10;19;29;19;9;29;9;19;11;7;7;7;9;40;40;9;10;7;9;7;4;10;9;5;21;9;21;5;9;62;21;9;21;62;20;21;21;20;21;62;20;1;20;1;21;14;22;1;23;7;23;22;1;45;31;23;1;9;37;1;45;9;1;9;45;1;37;9;14;1;14;31;31;1;23;14;1;14;45;1;45;11;1;5;37;37;21;45;5;21;1;45;21;45;111;21;27;31;45;21;37;5;21;5;31;21;31;23;37;21;45;17;45;23;22;45;45;16;45;45;24;37;45;16;45;37;16;31;45;16;45;45;31;23;16;31;45;27;9;37;45;45;27;45;9;9;45;10;45;9;10;1;37;50;31;1;50;1;31;95;31;23;3;95;3;31;95;37;45;3;2;39;2;45;26;45;31;7;26;7;26;9;4;10;4;9;9;10;9;45;10;45;45;45;9;10;10;10;10;45;45;10;10;10;10;10;10;4;41;4;10;10;5;30;10;5;5;10;5;30;10;16;8;10;10;16;16;10;32;10;10;10;32;10;32;10;3;10;18;18;3;18;10;3;10;5;3;31;3;5;3;11;31;9;21;10;21;10;27;10;21;27;21;10;27;22;27;22;31;48;31;22;48;22;14;9;48;9;22;48;43;9;19;24;43;43;6;43;24;24;45;95;6;45;6;45;101;101;19;95;45;19;9;9;19;9;95;19;101;9;19;22;101;19;95;22;19;43;95;8;19;43;19;43;8;19;8;29;8;29;8;0;9;9;6;9;9;0;8;8;9;9;9;48;9;21;5;21;9;10;9;21;10;5;33;10;33;5;31;36;23;31;12;27;31;12;27;31;10;33;31;33;10;31;10;10;10;10;10;32;10;32;10;10;10;10;10;32;9;10;4;10;4;10;6;20;10;10;20;10;20;41;20;41;27;27;18;10;9;9;18;20;25;9;20;9;25;41;55;9;9;41;9;55;41;53;9;27;9;55;27;41;9;2;9;0;9;10;9;48;0;4;18;10;18;5;5;8;4;18;1;24;24;24;17;17;1;19;32;18;1;18;32;1;32;24;6;24;22;24;22;18;8;18;10;18;45;39;15;33;24;24;15;24;18;18;33;45;6;15;6;18;6;0;6;6;6;65;0;51;2;79;65;73;2;65;65;2;59;73;65;65;2;51;2;65;51;65;2;51;22;45;6;57;57;22;6;6;6;57;71;6;71;71;10;6;6;57;65;97;27;65;27;97;65;97;27;62;27;55;27;55;62;55;27;62;27;51;62;31;31;31;27;27;1;62;1;25;25;62;25;1;57;1;25;57;25;1;6;57;0;25;71;6;71;0;71;0;91;0;6;6;91;6;0;91;8;7;77;10;42;77;17;45;77;45;17;85;56;45;71;45;0;85;0;2;0;2;57;26;13;57;0;57;0;26;51;26;42;67;51;67;19;51;19;71;19;1;67;67;85;67;67;79;67;6;65;67;51;61;63;53;6;53;63;10;6;63;6;10;57;10;6;57;10;57;10;1200026;4;5;16;4;16;16;5;7;5;16;6;16;1200026;1200026;64;5;5;10;48;64;6;5;6;70;5;70;5;6;12;70;5;5;70;5;18;76;24;8;76;8;6;171;6;6;25;6;25;6;6;6;25;82;25;6;82;6;6;26;177;2;88;2;88;2;89;183;89;2;3;2;89;8;89;89;2;2;6;2;89;6;89;3;23;3;5;5;10;10;10;23;10;5;23;3;10;37;78;78;3;37;3;78;37;78;3;37;3;3;78;78;37;72;3;37;64;37;64;23;20;3;23;3;20;27;20;10;3;18;2;10;18;10;2;18;2;10;10;26;10;4;26;4;10;10;47;10;66;47;66;5;47;5;66;33;42;33;18;12;3;3;6;8;6;8;6;8;11;6;6;6;2760000;11;6;7;2760024;6;2760026;7;5;5;2760026;5;1;5;6;56;1;40;0;56;27;56;0;27;0;0;50;27;50;0;27;0;42;17;10;0;4;0;10;4;10;0;9;0;10;9;10;0;0;3;0;10;5;22;0;17;0;30;17;100;10;10;38;0;100;38;100;0;24;0;100;24;100;0;3;0;100;3;17;57;0;0;3;0;57;20;0;0;20;28;79;55;28;22;14;73;55;79;14;1;55;1;87;55;79;1;57;9;79;49;79;9;49;4;47;43;16;43;16;16;47;5;43;16;79;16;5;16;45;16;5;37;5;16;37;16;16;31;16;31;2;23;2;79;0;118;110;81;104;0;73;25;104;67;96;67;241;6;67;90;69;90;61;82;79;76;76;79;68;6;79;6;62;79;62;10;79;10;54;79;48;233;6;6;6;40;31;34;222;10;37;10;222;45;210;10;51;10;210;59;194;188;0;65;0;129;47;129;135;143;0;6;6;143;55;149;6;61;5;157;69;157;69;163;75;171;5;83;5;79;83;171;5;177;10;10;185;10;185;10;10;68;10;191;0;199;199;10;0;11;205;0;213;6;0;6;213;0;87;6;219;54;227;6;54;10;10;10;227;233;241;17;241;10;0;0;247;5;0;5;93;101;0;5;0;26;101;107;0;115;5;0;10;115;0;121;10;0;10;81;3;81;81;26;3;26;70;3;70;12;3;12;58;3;58;50;0;3;0;3;0;9;44;3;44;9;3;247;36;3;30;247;3;247;81;33;30;39;188;30;180;47;12;53;174;61;12;174;12;166;61;67;12;75;75;160;160;152;6;146;89;89;6;89;146;241;6;241;78;78;8;66;138;8;20;132;66;89;20;36;132;124;20;118;20;20;20;8;20;8;20;63;20;0;63;20;4;0;0;20;0;4;20;0;14;20;0;20;29;14;20;4;18;20;18;4;20;4;39;4;30;24;39;30;25;24;30;24;25;0;25;0;24;25;24;29;25;18;24;24;23;12;12;18;18;18;23;39;12;42;23;42;39;23;39;36;28;36;39;28;25;25;30;12;30;25;12;25;18;24;18;18;24;18;18;24;24;18;24;24;6;18;39;45;39;6;45;6;25;45;25;6;31;6;18;31;18;6;13;6;18;18;13;39;13;10;61;39;10;39;61;14;47;1;39;25;25;8;117;9;1;3;1;9;3;0;1;3;10;27;25;10;0;3;27;27;0;11;6;0;6;6;6;3;3;0;26;6;6;50;6;5;42;5;5;5;58;5;123;53;58;45;20;58;20;26;58;26;6;58;129;50;59;50;51;135;6;0;39;46;0;141;6;58;6;141;0;52;6;6;6;52;52;6;147;58;6;6;66;6;58;74;58;6;9;6;58;0;153;6;0;6;64;64;64;6;6;0;6;64;159;0;70;6;6;0;6;70;165;6;63;6;6;6;21;6;21;7;21;6;0;21;0;21;0;21;0;6;21;20;9;21;10;20;1;21;25;10;21;0;8;21;8;25;21;25;8;21;8;14;21;8;21;24;21;24;45;45;21;45;5;21;2;45;21;31;31;2;21;2;21;2;5;21;6;5;21;5;0;103;0;5;0;5;0;0;0;5;8;0;0;6;0;0;8;8;6;8;0;6;0;12;0;12;0;0;0;11;0;11;0;4;0;4;11;8;11;4;11;0;11;55;0;55;11;0;11;55;66;0;66;11;0;11;66;0;19;11;11;11;19;11;29;11;11;11;11;29;29;11;29;11;17;11;18;17;18;11;0;11;18;0;18;11;42;11;18;54;18;11;11;18;11;18;11;11;11;18;11;29;11;11;18;11;5;11;25;11;25;11;0;27;8;8;0;8;11;8;11;11;10;10;29;39;39;39;8;39;10;8;10;39;29;39;10;29;10;10;39;9;25;10;9;10;9;10;25;13;0;0;0;13;6;17;10;17;7;34;4;10;10;34;10;4;34;10;10;6;20;10;9;10;20;9;5;27;31;8;5;31;5;41;5;7;7;20;7;10;20;5;7;20;7;41;8;41;9;5;0;7;7;30;2535;0;2559;2561;20;20;0;20;2561;0;3;20;0;20;3;20;0;20;0;0;0;20;0;20;0;37;0;20;37;20;10;20;37;20;10;37;10;20;37;20;10;23;20;23;20;23;20;23;24;9;24;20;9;5;20;9;5;3;5;20;20;3;20;8;3;8;20;3;20;8;10;8;20;24;20;8;24;3;41;24;26;9642006;0;9642032;26;0;26;9642032;0;9642032;5;0;5;15;0;15;17;0;17;36;36;6;22;11;22;6;4;4;5;6;4;6;4;5;34;90;11;34;11;11;90;34;90;8;0;8;84;0;84;84;8;0;8;5;8;84;0;84;8;21;8;76;76;21;8;21;8;0;21;0;8;21;8;0;21;0;10;10;21;3;0;0;21;0;3;19;3;0;19;0;3;59;19;59;0;19;0;59;40;59;0;40;0;59;5;59;0;5;0;7;5;5;0;5;0;5;5;5;5;0;5;0;0;5;5;5;7;5;5;21;25;3;19;11;3;3;21;20;20;20;3;20;3;3;3;20;3;20;29;3;29;20;0;20;6;0;6;20;0;20;14;0;14;0;15;35;0;0;35;15;0;10;25;23;26;10;23;31;26;18;26;31;18;0;35;24;18;0;18;3;24;45;18;45;7;18;4;24;18;24;4;18;4;45;45;24;45;5;24;5;31;24;31;5;24;0;24;18;24;0;18;3;45;45;3;18;3;45;7;3;7;3;31;31;24;3;3;7;3;24;10;45;45;25;10;25;45;10;31;25;30;31;6;9;10;9;10;9;0;39;9;0;39;12;48;39;27;39;27;12;25;15;27;25;27;15;66;27;18;37;37;15;18;15;3;18;19;15;18;11;19;18;19;11;39;10;19;39;19;27;39;19;19;25;12;19;19;25;25;12;12;8;12;25;8;25;12;4;12;19;4;19;12;4;12;0;4;0;12;4;12;0;25;20;14;6;6;25;4;4;14;14;12;4;21;4;21;4;31;5;21;31;20;5;6;5;20;10;6;10;20;6;8;10;6;10;8;6;8;15;2;15;9;2;9;9;15;2;9;9;17;47;17;20;47;20;17;26;9;20;26;20;9;26;9;20;41;1;41;9;1;9;41;1;27;9;1;9;27;1;12;17;26;17;12;0;17;26;17;0;26;15;35;15;27;15;15;15;9;9;9;27;9;15;27;15;9;18;27;18;15;15;15;39;39;5;5;5;28;28;28;27;27;27;20;20;20;12;12;16;16;19;19;19;14;14;5;5;10;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;28;28;28;28;28;27;27;27;27;27;27;27;27;27;27;27;27;27;27;27;5;5;10;10;5;5;5;5;5;5;5;5;5;5;5;28;28;28;15;15;14;14;14;9;9;9;9;9;9;9;40;4;4;5;5;5;5;6;1;1;1;1;9;9;9;9];
%% multiply the lower and upper bnds witht their respective number of carbon atoms
varMinCarb = varMin.*metCarb(metID);
if exist('totalCarbon','var')
    if isempty(totalCarbon)
        total_carbon_cons = sum(abs(varMinCarb(varMinCarb<0)));
    elseif ~isnumeric(totalCarbon)
        h = warndlg('Warning the totalCarbon input is not numeric. The totalCarbon will be calculated based on the lower bnds of the input variables. If you wish to continue press OK or else kill the code and input a numeric value for totalCarbon.','Problem with input');
        total_carbon_cons = sum(abs(varMinCarb(varMinCarb<0)));
    else
        total_carbon_cons = totalCarbon;
        actual_total_carbon_cons = sum(abs(varMinCarb(varMinCarb<0)));
    end
else
    total_carbon_cons = sum(abs(varMinCarb(varMinCarb<0)));
end
total_carbon_cons = total_carbon_cons;
% if ismember(condition,'exponential')
%     total_carbon_cons = total_carbon_cons * 0.05;
% elseif ismember(condition,'stationary')
%     total_carbon_cons = total_carbon_cons * 0.01;
% end
%% check if cofactorPairs have been passed on and set cofactors to 0 in model.S

% check for cofPari variable and/or calculate new
if exist('cofPairs','var')
    if isempty(cofPairs)
        smallMets = findSmallMets(model);
        metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
        cofPairs = findCofactorPairs_max(metPairs,110);
        cofPairs(5,:) = [];
        cofPairs(4,:) = [];
    end
else
    smallMets = findSmallMets(model);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% THIS NEEDS CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%
    metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
    %%%%%%%%%%%%%%%%%%%%%%%%% THIS NEEDS CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cofPairs = findCofactorPairs_max(metPairs,110);
    cofPairs(4,:)=[];
end

% initialize variable
cofMetIDs = [];
% retrieve all the metabolite IDs for cofactors
for i = 1:size(cofPairs,1)
    cof_tmp = cell2mat(cofPairs{i,6});
    cofMetIDs = [cofMetIDs;cof_tmp];
end

% set cofMets to 0 in model.S
model_tmp = model;
model_tmp.S(cofMetIDs,:) = 0;

% find all rxns containing coa mets on both sides of the rxn
coa = findCoAs(model_tmp,' ','C');

% find all rxns containing cytochromes, quinones, and hemes
lrgMlc = findLrgMlcs(model,metCarb);

%% calculating new bnds algorithm !!!!!
minmax_new = minmax;
mmAdj = zeros(num_rxns,2);
for i = 1:num_rxns
    % First, check if the rxn is an important one
    if (i == findRxnIDs(model,'BIOMASS_cho')) % Biomass non producing
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        bericht{i,1} ='Biomass nonproducing';
        constraint = abs(total_carbon_cons/subCarbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;
        end
        continue
    elseif (i == findRxnIDs(model,'BIOMASS_cho_producing')) %Biomass producing
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        bericht{i,1} ='Biomass producing';
        constraint = abs(total_carbon_cons/subCarbon);
            if abs(minmax(i,2)) > constraint
                minmax_new(i,2) = constraint;
                mmAdj(i,2) = 1;
            end
            continue
    elseif (i == findRxnIDs(model,'igg_hc')) %HC Synthesis
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        prodIDs = find(model.S(:,i) > 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        prodCarbon = sum(abs(full(model.S(prodIDs,i))).*metCarb(prodIDs,1));
        igg_hc_Carbon = subCarbon-prodCarbon;
        bericht{i,1} ='IgG Heavy Chain Synthesis';
        constraint = abs(total_carbon_cons/igg_hc_Carbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;
        end
        continue
    elseif (i == findRxnIDs(model,'DM_igg_g')) %igg demand reaction
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = 150e3/12;
        bericht{i,1} ='IgG Demand Reaction';
        constraint = abs(total_carbon_cons/subCarbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;
        end
        continue
    elseif (i == findRxnIDs(model,'igg_lc')) %LC Synthesis
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        prodIDs = find(model.S(:,i) > 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        prodCarbon = sum(abs(full(model.S(prodIDs,i))).*metCarb(prodIDs,1));
        igg_lc_Carbon = subCarbon-prodCarbon;
        bericht{i,1} ='IgG Heavy Chain Synthesis';
        constraint = abs(total_carbon_cons/igg_lc_Carbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;
        end
        continue
    end
    % evaluate if rxn is part of the measured constraints
    if ismember(i,varID_tmp)
        % If rxn is measured, then we don't set carbon constraints
        % as there are already lower and upper constraints to the reaction.
        bericht{i,1} = 'set_bnd'; 
        continue  
    % evaluate if rxn is part of the skip-rxns set       
    elseif ismember(i,skipRxns)
        bericht{i,1} = 'unconstr - oxPhos';
        continue
    % evaluate if rxn is balanced or unbalanced
    else
        % First, check if the reaction carries any flux under normal conditions
        if (minmax(i,1) == 0) && (minmax(i,2) == 0)
            bericht{i,1} = 'no_flux';
            continue
        % Check wheather the reaction is aprt fo the drain set/exchange
        % between the cell and the enviroment.            
        elseif ismember(i,isDrain)    
            bericht{i,1} = 'balanced - drain';
            subIDs = find(model_tmp.S(:,i) < 0);
            prodIDs = find(model_tmp.S(:,i) > 0);

            % account for stoichiometric coefficients and add up carbon
            % As we set all the lower bounds for the exchange reactions to
            % 0, the atom constrain will take effect only on the upper
            % bound
            subCarbon = sum(abs(full(model_tmp.S(subIDs,i))).*metCarb(subIDs,1));
            prodCarbon = sum(abs(full(model_tmp.S(prodIDs,i))).*metCarb(prodIDs,1));

            if subCarbon == 0 && prodCarbon == 0 
                continue
            elseif subCarbon == 0
                constraint = abs(total_carbon_cons/prodCarbon);
            else
                constraint = abs(total_carbon_cons/subCarbon);
            end
            if model_tmp.rev(i,1) == 1 % reversible
                if abs(minmax(i,1)) > constraint
                    minmax_new(i,1) = -constraint;
                    mmAdj(i,1) = 1;
                end
                if abs(minmax(i,2)) > constraint
                    minmax_new(i,2) = constraint;
                    mmAdj(i,2) = 1;
                end
            else % irreversible
                if abs(minmax(i,2)) > constraint
                    if constraint > abs(minmax(i,1))
                        minmax_new(i,2) = constraint;
                        mmAdj(i,2) = 1;
                    end
                end
            end            
        else
            % Meaning, if the reaction is not 0 or drain
            
            % find substrates and products 
            subIDs = find(model_tmp.S(:,i) < 0);
            prodIDs = find(model_tmp.S(:,i) > 0);            
            % account for stoichiometric coefficients and add up carbon
            subCarbon = sum(abs(full(model_tmp.S(subIDs,i))).*metCarb(subIDs,1));
            prodCarbon = sum(abs(full(model_tmp.S(prodIDs,i))).*metCarb(prodIDs,1));
            
            if ismember(i,coa.coaRxnIDs)
                if coa.missMatch
                    warning('Stoichiometry of CoA is not matching. Smaller number is subtracted from carbon balance. Check reaction reaction manually.')
                    if coa.subAtom(ismember(coa.coaRxnIDs,i)) > coa.subCarbon(ismember(coa.coaRxnIDs,i))
                        subCarbon = subCarbon - coa.subAtom(ismember(coa.coaRxnIDs,i));
                        prodCarbon = prodCarbon - coa.prdAtom(ismember(coa.coaRxnIDs,i));
                    else
                        subCarbon = subCarbon - coa.subAtom(ismember(coa.coaRxnIDs,i));
                        prodCarbon = prodCarbon - coa.prdAtom(ismember(coa.coaRxnIDs,i));
                    end
                else
                    subCarbon = subCarbon - coa.subAtom(ismember(coa.coaRxnIDs,i));
                    prodCarbon = prodCarbon - coa.prdAtom(ismember(coa.coaRxnIDs,i));
                end
            end     
            
            % Account for the large molecules cytochromes, hemes, and quinones
            % and remove the atoms taking part there (need to check why
            % though???????)                    
            if ismember(i,lrgMlc.lmRxnIDs)
                subCarbon = subCarbon - (lrgMlc.sub(ismember(lrgMlc.lmRxnIDs,i)));
                prodCarbon = prodCarbon - lrgMlc.prd(ismember(lrgMlc.lmRxnIDs,i));
            end         
            
            % Check whether the number of atoms in the substrates or
            % products are equal to 0; even if the minmax is different tha
            % 0 the molecules could not have any carbons or nitrogens etc.                 
            if subCarbon == 0 && prodCarbon == 0
                bericht{i,1} = 'no C in the reaction';
                continue
            end

            % If in fact, the reaction is balanced then do the following:
            if subCarbon == prodCarbon     
                bericht{i,1} = 'balanced'; 
                constraint = abs(total_carbon_cons/subCarbon);
                if model_tmp.rev(i,1) == 1
                    if abs(minmax(i,1)) > constraint
                        minmax_new(i,1) = -constraint;
                        mmAdj(i,1) = 1;
                    end
                    if abs(minmax(i,2)) > constraint
                       minmax_new(i,2) = constraint;
                       mmAdj(i,2) = 1;
                    end
                else
                    if abs(minmax(i,2)) > constraint
                        if constraint > abs(minmax(i,1))
                            minmax_new(i,2) = constraint;
                            mmAdj(i,2) = 1;
                        end
                    end
                end     
            else
                
                % If the reaction is not balanced                
                % First check for the cofactor pairs              
                if find(ismember(cofMetIDs,find(model.S(:,i))))
                    % if the sum of all carbon is even then it is balanced 
                    % (not bulletproof but if off by only one carbon it works)
                    % atoms for the substrates, products and the coenzymes
                    % - cofactors
                    if mod(sum([subCarbon;prodCarbon;metCarb(cofMetIDs(ismember(cofMetIDs,find(model.S(:,i)))))]),2) == 0 
                        bericht{i,1} = 'balanced - CofPair';
                        % if balanced then its go with the smaller one (constraining less).
                        % This is because there seems to be only 1 cofactor in this equation 
                        % (e.g. atp -> amp, this is not picked up as a cofactor pair but atp 
                        % is a cofactor) so the carbon of atp/amp should not be considered when constraining
                        if subCarbon > prodCarbon 
                            constraint = abs(total_carbon_cons/prodCarbon);
                        else
                            constraint = abs(total_carbon_cons/subCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                               minmax_new(i,2) = constraint;
                               mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    else
                        
                        bericht{i,1} = 'unbalanced - CofPair';
                        % if unbalanced then go with the smaller one (constraining less)
                        % usually contains AMP which isn't considered to be
                        % a cofactor with no CoA though                        
                        if subCarbon < prodCarbon 
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                               minmax_new(i,2) = constraint;
                               mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                        
                        % check whether there is a metabolite with X formula type.
                        flag = 0;    
                        f=1;
                        while (f<= size(subIDs,1) && flag == 0)
                            if (model.metFormulas{subIDs(f,1)} == 'X')
                                flag = 1;
                            end
                            f = f + 1;
                        end
                        f = 1;
                        while (f <= size(prodIDs,1) && flag == 0)
                            if (model.metFormulas{prodIDs(f,1)} == 'X')
                                flag = 1;
                            end
                            f = f + 1;
                        end
                        if (flag == 1)
                            bericht{i,1} = 'Contains X';
                            % Set as the constrain the larger one if there is an X
                            % metabolite; since normally the reaction would have
                            % been balanced - apprantly every time that's the case.
                            if subCarbon > prodCarbon
                                constraint = abs(total_carbon_cons/subCarbon);
                            else
                                constraint = abs(total_carbon_cons/prodCarbon);
                            end
                            if model_tmp.rev(i,1) == 1
                                if abs(minmax(i,1)) > constraint
                                    minmax_new(i,1) = -constraint;
                                    mmAdj(i,1) = 1;
                                end
                                if abs(minmax(i,2)) > constraint
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            else
                                if abs(minmax(i,2)) > constraint
                                    if constraint > abs(minmax(i,1))
                                        minmax_new(i,2) = constraint;
                                        mmAdj(i,2) = 1;
                                    end
                                end
                            end
                        end
                    end
                else
                    % When the number of atoms in the substrate and product
                    % are differe, usually in rxns with metabolite X; and
                    % there is no coa or anything

                    % check whether there is a metabolite with X formula type.
                    flag = 0;    
                    f=1;
                    while (f<= size(subIDs,1) && flag == 0)
                        if (model.metFormulas{subIDs(f,1)} == 'X')
                            flag = 1;
                        end
                        f = f + 1;
                    end
                    f = 1;
                    while (f <= size(prodIDs,1) && flag == 0)
                        if (model.metFormulas{prodIDs(f,1)} == 'X')
                            flag = 1;
                        end
                        f = f + 1;
                    end
                    if (flag == 1)
                        bericht{i,1} = 'Contains X';
                        % Set as the constrain the larger one if there is an X
                        % metabolite; since normally the reaction would have
                        % been balanced - apprantly every time that's the case.
                        if subCarbon > prodCarbon
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                                minmax_new(i,2) = constraint;
                                mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    else
                         % When there is no CofPair or X (there is probably
                         % EPO or HX but we didn't account for it yet :(
                        bericht{i,1} = 'unbalanced';
                        if subCarbon > prodCarbon
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                                minmax_new(i,2) = constraint;
                                mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% qunatify how much solution space was constraint compared to original minmax
mm_diff = minmax(:,2)-minmax(:,1);
mm_new_diff = minmax_new(:,2)-minmax_new(:,1);
diff = 100*(1-(sum(mm_new_diff)/sum(mm_diff)));
        
%% variable complier
unbRxns = sum(ismember(bericht,'unbalanced'))+sum(ismember(bericht,'unbalanced - CofPair'));
balRxns = sum(ismember(bericht,'balanced'))+sum(ismember(bericht,'balanced - X'))...
    +sum(ismember(bericht,'balanced - drain'))+sum(ismember(bericht,'unconstr - oxPhos'));
setRxns = sum(ismember(bericht,'set_bnd'));
nfRxns = sum(ismember(bericht,'no_flux'));
Unknowns = sum(ismember(bericht,'Contains X'));

carbonConst.title = '*** carbon closure constraining ***';
carbonConst.dateAndTime = datetime();
carbonConst.constraintsName = condition;
if isfield(model,'description')
    carbonConst.model = model.description;
end
carbonConst.totalCarbonConsumed = total_carbon_cons;
if exist('actual_total_carbon_cons','var')
    if ~isempty(actual_total_carbon_cons)
        carbonConst.actualTotalCarbonConsumed = actual_total_carbon_cons;
    end
end
carbonConst.report.summary = bericht;
carbonConst.report.unbalanced = unbRxns;
carbonConst.report.balanced = balRxns;
carbonConst.report.setRxns = setRxns;
carbonConst.report.noFlux = nfRxns;
carbonConst.report.Unknowns = Unknowns;
carbonConst.bndAjusted.rawData = mmAdj;
carbonConst.bndAjusted.numberOf_LB_adj = sum(mmAdj(:,1));
carbonConst.bndAjusted.numberOf_UB_adj = sum(mmAdj(:,2));
carbonConst.CoAevaluation = coa;
carbonConst.carbonCount = carbonCount;
carbonConst.cofactorPairs = cofPairs;
carbonConst.originalMinMax = minmax;
carbonConst.newMinMax = minmax_new;
carbonConst.solSpaceReduction = strcat(num2str(diff),'%');

%% Test if model solves under new conditions

model.lb = minmax_new(:,1);
model.ub = minmax_new(:,2);

% model.c(findRxnIDs(model,'BIOMASS_cho')) = 0;
% for i = 1 : 3561
%     model.lb(1:i) = minmax_new(1:i,1);
%     model.ub(1:i) = minmax_new(1:i,2);
%     sol = optimizeCbModel(model,'max')
% if isnumeric(sol.origStat)
%     if sol.origStat < 1
%         error('Model does not solve. Think about relaxing the bounds')
%     else
%         disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
%     end
% else
%     if ismember(cellstr(sol.origStat),{'INFEASIBLE'})
%         warning('Model does not solve. Think about relaxing the bounds')
%     else
%         disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
%     end
% end
% end

sol = optimizeCbModel(model,direction);

if isnumeric(sol.origStat)
    if sol.origStat < 1
        error('Model does not solve. Think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end
else
    if ismember(cellstr(sol.origStat),{'INFEASIBLE'})
        warning('Model does not solve. Think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end
end

carbonConst.solutionNew = sol;

%% constrain uptakes based on passed in limit (e.g. lost/closed carbon)

if exist('limEX','var')
    if ~isempty(limEX)
        carbonConst.EXrxnAdjust = '*** EX rxns adjustment ***';
        carbonConst.limEX = limEX;
        
        % get all the max values for the EX rxns not set
        varMaxN = minmax_new(isDrainNew,2);
        % get met indices
        metID_ind = find(model.S(:,isDrainNew));
        % get metIDs from Met indices
        [metID_EX,~] = ind2sub(size(model.S),metID_ind);
        % get rxnIDs where rxn actually has a flux and ignore any fluxes
        % that are set to 0
        tmpRxnID = isDrainNew(find(varMaxN));
        % how many rxns have been set to 0?
        numEXzero = size(isDrainNew,1)-size(tmpRxnID,1);
        % get metIDs from EX rxn IDs so carbon per met can be determined
        tmpMetID = metID_EX(find(varMaxN));
        % calculate new bnds
        limEXnew = limEX./metCarb(tmpMetID);
        % get rid of Inf
        limEXnew(~isfinite(limEXnew)) = 0;
        % find which reactions needf adjustment
        adjEX = tmpRxnID(limEXnew < varMaxN(ismember(isDrainNew,tmpRxnID)));
        %  check if there is any bnds to be adjusted
        if ~isempty(adjEX)
            % make calculations as to how much solution space shrinks
            % compared to original
            diffLim = 100*(1-(sum(limEXnew(ismember(tmpRxnID,adjEX)))/sum(minmax_new(adjEX,2))));
            limEXadj = limEXnew(ismember(tmpRxnID,adjEX));
            minmax_new(adjEX,2) = limEXnew(ismember(tmpRxnID,adjEX));
            
            mm_lim_diff = minmax_new(:,2)-minmax_new(:,1);
            diff_n = 100*(1-(sum(mm_lim_diff)/sum(mm_diff)));
            
            model.lb = minmax_new(:,1);
            model.ub = minmax_new(:,2);
            sol = optimizeCbModel(model,direction);

            if isempty(sol.f)
                warning('model does not solve. think about relaxing the bounds')
            else
                disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
            end
            
            carbonConst.numAdjusted = size(adjEX,1);
            carbonConst.solSpaceReductionEX  = strcat(num2str(diff_n),'%');
            carbonConst.EXspaceReduction = strcat(num2str(diffLim),'%');
        else
            warning('no adjustment of EXrxns was done as already constrained more than limEX')
            carbonConst.warning = 'no adjustment was done';
        end
    end
end

%% run minmax on new minmax bnds and see if solution space is reduced further
if runFinalMinMax
    model.lb = minmax_new(:,1);
    model.ub = minmax_new(:,2);
    minmax_post_new = runMinMax_GF(model);
    carbonConst.postNmm = '*** post new minmax ***';
    carbonConst.postNewMinMax = minmax_post_new;
    % qunatify how much solution space was constraint compared to original minmax
    mm_post_diff = minmax_post_new(:,2)-minmax_post_new(:,1);
    diff = 100*(1-(sum(mm_post_diff)/sum(mm_diff)));
    carbonConst.solSpaceReductionPostCC = strcat(num2str(diff),'%');

    % test if new minmax solves
    model.lb = minmax_post_new(:,1);
    model.ub = minmax_post_new(:,2);
    sol = optimizeCbModel(model,'max');

    if isempty(sol.f)
        warning('model does not solve. think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end

    carbonConst.solutionPostNew = sol;
end
end
