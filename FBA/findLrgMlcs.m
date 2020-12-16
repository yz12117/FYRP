%% This script is designed to find large molecule containing rxns
%
%   This script will go through a cobra model and identify all rxns
%   containing large molecules meaning cytochromes, hemes, and quinones
%   which do not contriobute to the carbon flux of a reaction but rather as
%   an electron carrier.
%   Only rxns where these molecules are present on both sides of the rxns, 
%   i.e. their oxidized and reduced will be identified. This will avoid 
%   flagging the reactions up where they get degraded or synthesiszed. 
%
%   atomNumber  =  the atom numbers from the findMetCarbon function,
%   we want to constrain
%
%
function [lrgMlcs] = findLrgMlcs(model,atomNumbers)

if exist('atomNumbers','var')
    if isempty(atomNumbers)
        atomNumbers = findMetCarbon(model);
        atomNumbers = atomNumbers.carbonAtoms
    elseif size(atomNumbers,1) < size(model.mets,1) || size(atomNumbers,1) > size(model.mets,1)
        atomNumbers = findMetCarbon(model);
        atomNumbers = atomNumbers.carbonAtoms
    end
else
    atomNumbers = findMetCarbon(model);
    atomNumbers = atomNumbers.carbonAtoms
end

% find all metIDs containing large molecules
quiCell = regexp(model.metNames,'quino');
cytCell = regexp(model.metNames,'cytoch');
hemCell = regexp(model.metNames,'heme');

% convert cell to metID dbl
quiIDs = find(~cellfun(@isempty,quiCell));
cytIDs = find(~cellfun(@isempty,cytCell));
hemIDs = find(~cellfun(@isempty,hemCell));

lrgIDs = [quiIDs;cytIDs;hemIDs];

rxnIDs = [];
for i = 1:size(lrgIDs,1)
    tmp = [];
    tmp = find(model.S(lrgIDs(i,1),:));
    rxnIDs = [rxnIDs;tmp'];
end

[C,~,~] = unique(rxnIDs);
valCount = hist(rxnIDs(:,1),C)';
trgtRxns = C(valCount>1);

% find logical vector for transportRxns
% isTrans  = checkTheTransport(model);
% isTrans  = findTrans(model);
%retrieve actual reactions containing electron transferring mlcls
% lrgMlcRxns= trgtRxns(~ismember(trgtRxns,find(isTrans)));
lrgMlcRxns= trgtRxns;
for i = 1:size(lrgMlcRxns,1)
    % find substrate and product IDs
    subIDs = find(model.S(:,lrgMlcRxns(i,1)) < 0);
    prdIDs = find(model.S(:,lrgMlcRxns(i,1)) > 0);
    
    % check if coa met is part of substarte and product list
    
    subLM = subIDs(ismember(subIDs,lrgIDs)); 
    prdLM = prdIDs(ismember(prdIDs,lrgIDs)); 
    
    if size(subLM,1) > 1 && size(prdLM,1) > 1
        for j = 1:size(subLM,1)
            sub(j,1) = atomNumbers(subLM(j,1),1)*...
                abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));

            prd(j,1) = atomNumbers(prdLM(j,1),1)*...
                abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
        end

        sub = sum(sub);
        prd = sum(prd);

    elseif size(subLM,1) > size(prdLM,1)
        if size(prdLM,1) > 1
            for j = 1:size(subLM,1)
                sub(j,1) = atomNumbers(subLM(j,1),1)*...
                    abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));
            end
            sub = sum(sub);

            for j = 1:size(prdLM,1)
                prd = atomNumbers(prdLM,1)*...
                    abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
            end
            prd = sum(prd);
        else
            if size(subLM,1) > 1
                for j = 1:size(subLM,1)
                    sub(j,1) = atomNumbers(subLM(j,1),1)*...
                        abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));
                end
                sub = sum(sub);
                prd = 0;
            else
                sub = atomNumbers(subLM(1,1),1)*...
                    abs(full(model.S(subLM(1,1),lrgMlcRxns(i,1))));
                prd = 0;
            end 
        end
    elseif size(subLM,1) < size(prdLM,1)
        if size(subLM,1) > 1
            for j = 1:size(prdLM,1)
                prd(j,1) = atomNumbers(prdLM(j,1),1)*...
                    abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
            end
            prd = sum(prd);

            for j = 1:size(subLM,1)
                sub = atomNumbers(subLM(j,1),1)*...
                    abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));
            end
            sub = sum(sub);
        else
            if size(prdLM,1) > 1
                for j = 1:size(prdLM,1)
                    prd(j,1) = atomNumbers(prdLM(j,1),1)*...
                        abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
                end
                prd = sum(sub);
                sub = 0;
            else
                prd = atomNumbers(prdLM(1,1),1)*...
                    abs(full(model.S(prdLM(1,1),lrgMlcRxns(i,1))));
                sub = 0;
            end 
        end
    elseif size(prdLM,1) == size(subLM,1)
        prd = atomNumbers(prdLM(1,1),1)*...
            abs(full(model.S(prdLM(1,1),lrgMlcRxns(i,1))));
        sub = atomNumbers(subLM,1)*...
            abs(full(model.S(subLM(1,1),lrgMlcRxns(i,1))));
    end

    subA(i,1) = sub;
    prdA(i,1) = prd;    
end

missMatches = find(subA - prdA);
lmRxnIDs = lrgMlcRxns;
lmRxnIDs(missMatches) = [];
subA(missMatches) = [];
prdA(missMatches) = [];

lrgMlcs.tiles = '*** Cyt/Qui/Heme Rxns ***';
lrgMlcs.time = date();
lrgMlcs.lmMetIDs = lrgIDs;
lrgMlcs.lmRxnIDs = lmRxnIDs;
lrgMlcs.sub = subA;
lrgMlcs.prd = prdA;
if ~isempty(missMatches)
    lrgMlcs.missMatch = true;
    lrgMlcs.numberOfMissMatches = size(missMatches,1);
    lrgMlcs.mM_IDs = lrgMlcRxns(missMatches);
else
    lrgMlcs.missMatch = false;
end  
