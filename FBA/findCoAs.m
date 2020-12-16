%% This script is designed to find CoA containing rxns
%
%   This script will go through a cobra model and identify all rxns
%   containing conenzyme A if it is present on both sides of the rxns. This
%   will avoid flagging the reactions up where coenzyme A gets degraded or
%   synthesiszed. CoA is not contributing to the actual carbon flux and
%   hence should not be taken into account for the carbon constraining
%   algorithm.
%

function [coa] = findCoAs(model,delimiter,atom)

if exist('atom','var')
    if isempty(atom)
        atom = 'C';
        warning('default atom C was chosen!')
    end
    if atom == 'C'
        Natom = 21;
    elseif atom == 'N'
        Natom = 7;
    elseif atom == 'O'
        Natom = 16;
    end
else
    atom = 'C';
    warning('default atom C was chosen!')
    Natom = 21;
end

if exist('delimiter','var')
    if isempty(delimiter)
        delimiter = '_';
        warning('default delimiter "_" was chosen!')
    end
else
    delimiter = '_';
    warning('default delimiter "_" was chosen!')
end

% concatenate pattern for regexp
pattern = strcat('coa',delimiter);
pattern2 = '-CoA';

% find all metIDs containing coa
coaCell = regexp(model.mets,pattern);
coaCell2 = regexp(model.metNames,pattern2);

% convert cell to metID dbl
coaIDs = find(~cellfun(@isempty,coaCell));
coaIDs2 = find(~cellfun(@isempty,coaCell2));

coaIDs = unique([coaIDs;coaIDs2]);

% find all rxns that contain a coa metabolite on both sides.
count = 0;
for i = 1:size(model.rxns,1)
    % find substrate and product IDs
    subIDs = find(model.S(:,i) < 0);
    prdIDs = find(model.S(:,i) > 0);
    
    % check if coa met is part of substarte and product list
    subCoA = sum(ismember(subIDs,coaIDs));
    prdCoA = sum(ismember(prdIDs,coaIDs));
    
    % if subCoA and prdCoA are both larger than one it means both sides of
    % the rxn contain the coenzyme and the ID is stored for later
    % processing
    if subCoA > 0 && prdCoA > 0
        count = count + 1;
        coaRxnIDs(count,1) = i;
        
        arr{count,1} = subIDs(ismember(subIDs,coaIDs)); % Places the CoA that take part in the rxn
        arr{count,2} = prdIDs(ismember(prdIDs,coaIDs)); % into the matrix arr column 1 and 2
        
        sub = 0;
        prd = 0;
               
        for j = 1:size(arr{count,1},1)
            sub = (abs(full(model.S(arr{count,1}(j,1),i)))*Natom) + sub; 
        end
        for j = 1:size(arr{count,2},1)
            prd = (full(model.S(arr{count,2}(j,1),i))*Natom) + prd; 
        end
        
        subAtom(count,1) = sub;
        prdAtom(count,1) = prd;
        
    end
end

missMatches = find(subAtom - prdAtom);



coa.tiles = '*** CoA Rxns ***';
coa.time = date();
coa.demilinter = delimiter;
coa.coaMetIDs = coaIDs;
coa.coaRxnIDs = coaRxnIDs;
coa.subAtom = subAtom;
coa.prdAtom = prdAtom;
if ~isempty(missMatches)
    coa.missMatch = true;
    coa.numberOfMissMatches = size(missMatches);
    coa.mM_IDs = coaRxnIDs(missMatches);
else
    coa.missMatch = false;
end  












