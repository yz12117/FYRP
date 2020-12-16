function cofactorPairs = findCofactorPairs_max(metPairs,MinNumOccurence)
%% A function to check metabolite pairs (found by function
% countMetPairs_max) and identify which ones are true co-factor pairs to
% be neglected in the carbon closure analysis. The criteria for exclusion
% of a pair are:
% Exclude a pair IF
% (1) The number of elements between products and substrates is different,
%     and the difference is larger than 1.
% (2) The number of elements between products and substrates is equal, but
%     the difference between participating atoms is larger than 5.
%
% INPUT
% - all_met_pairs*: This is a n-by-6 cell, where n is the number of met
%   pairs without repetition, found in the model. The 6 columns correspond
%   to:
%   (1) names of the metabolite pair separated by a colon punctuation
%   (2) number of occurence of this metabolite pair
%   (3) Formula of substrate
%   (4) Formula of product
%   (5) cell array with rxnIDs metPair participates in
%
%   *This is the output of countMetPairs_max.m
% - MinNumOccurence: Minimum number of occurences of each metabolite pair,
%   acting as a threshold for considering a pair to be a valid cofactor
%   pair. --> default 80
%
% OUTPUT
% - met_pairs_to_remove:


if nargin < 2
    MinNumOccurence = 80;
end

% Break down the input in separate cell arrays
MetPairNames = metPairs(:,1);
MetPairOccurence = metPairs(:,2);
MetPairOccurence = cell2mat(MetPairOccurence);
MetPairSubFormula = metPairs(:,3);
MetPairProdFormula = metPairs(:,4);
MetPairRxnIDs = metPairs(:,5);
MetPairMetIDs = metPairs(:,6);

% Find and keep only the pairs that have at least MinNumOccurence "hits"
PairsHighOccurence = find(MetPairOccurence >= MinNumOccurence);

MetPairNames = MetPairNames(PairsHighOccurence,1);
MetPairOccurence = MetPairOccurence(PairsHighOccurence,1);
MetPairSubFormula = MetPairSubFormula(PairsHighOccurence,1);
MetPairProdFormula = MetPairProdFormula(PairsHighOccurence,1);
MetPairRxnIDs = MetPairRxnIDs(PairsHighOccurence,1);
MetPairMetIDs = MetPairMetIDs(PairsHighOccurence,1);

% Split Metabolite Formulas to Elements and Number of participating atoms
for i=1:size(MetPairSubFormula,1)
    % Substrates
    % Extract the numerical entries. Use letters (small/capital) to split
    % the formula
    tmp_SubElNum = regexp(MetPairSubFormula{i},'[a-z A-Z]','split');
    % Extract the letters. Use numbers to split the formula
    tmp_SubElem = regexp(MetPairSubFormula{i},'[\d]','split');
    % Remove empty cells
    emptyNumCells = cellfun(@isempty,tmp_SubElNum);
    tmp_SubElNum(emptyNumCells) = [];
    emptyElemCells = cellfun(@isempty,tmp_SubElem);
    tmp_SubElem(emptyElemCells) = [];
        
    SubstrateElementsNum{i,:} = tmp_SubElNum;
    SubstrateElements{i,:} = tmp_SubElem;

    clear tmp_SubElNum tmp_SubElem
    
    % Products
    % Extract the numerical entries. Use letters (small/capital) to split
    % the formula
    tmp_SubElNum = regexp(MetPairProdFormula{i},'[a-z A-Z]','split');
    % Extract the letters. Use numbers to split the formula
    tmp_SubElem = regexp(MetPairProdFormula{i},'[\d]','split');
    % Remove empty cells
    emptyNumCells = cellfun(@isempty,tmp_SubElNum);
    tmp_SubElNum(emptyNumCells) = [];
    emptyElemCells = cellfun(@isempty,tmp_SubElem);
    tmp_SubElem(emptyElemCells) = [];
        
    ProductElementsNum{i,:} = tmp_SubElNum;
    ProductElements{i,:} = tmp_SubElem;

    clear tmp_SubElNum tmp_SubElem
end

NotCofactorPairs = zeros(size(PairsHighOccurence));
% Try to find those metabolite pairs that are not cofactor pairs
for i=1:size(PairsHighOccurence,1)
    % Discard a pair IF
    % (1) The number of elements between products and substrates is
    % different and the difference is larger than 1
    if ~isequal(cell2mat(SubstrateElements{i}),cell2mat(ProductElements{i})) &&...
            abs(size(SubstrateElements{i},2) - size(ProductElements{i},2)) > 1
        NotCofactorPairs(i) = 1;
    end
    % (2) The number of elements between products and substrates is
    % equal but the difference between participating atoms is larger than 4
    if abs(sum(str2double(SubstrateElementsNum{i})) - sum(str2double(ProductElementsNum{i}))) > 4
        NotCofactorPairs(i) = 1;
    end
    
end

cofactorPairs = metPairs(PairsHighOccurence,:);
cofactorPairs(NotCofactorPairs==1,:) = [];






