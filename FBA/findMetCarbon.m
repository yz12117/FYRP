function cStruc = findMetCarbon(model)
%% This script is designed to determine the number of carbon atoms per metabolite
%   This function goes through all metabolites of a cobra model, and generates
%   a list with the number of carbon atoms of that metabolite
%
%               *****************************************
%               *   created by - Maximilian Lularevic   *
%               *              2017-09-30               *
%               *****************************************
%
% INPUTS
%   model   =   cobra nodel (needs to contain model.metFormulas field)
%
% OUTPUT
%   cStruc  =   structure containing the following fields
%               
%
%   Version History:
%       Date:           Comment:                        Initials:
%       2017-11-17      Updated function to be more     MLU  
%                       accurate in carbon count 
%                       (COFULLR2COFULLR1 etc.)
%       Date:
%       2018-09-26      Added O, H, Na, Zn, Ca, K,      Author:
%                       S and P count for all the mets  Athanasios Antonakoudis    
%                       and the calculation of the 
%                       molar mass of each compound
%

%% Initialization
% check if metFormula is a field in model structure
if ~isfield(model,'metFormulas')
    error('model structure needs metFormula field')
end

% count the number of mets and store model.mets in mets variable
[met_num,~] = size(model.S);
mets = model.mets;

% remove the compartment information from mets (i.e. nadh_c -> nadh) and
% store new mets in metList
for i = 1:met_num
%     met_tmp = regexp(mets(i),'_[a-z]{1,2}$','split');
    met_tmp=regexp(mets(i),'\[','split');
    metList{i,1} = met_tmp{1,1}{1,1};
end

% find all fields with no metFormula in it
% NoFormulaStr = {'NA','None','NULL',''};
% NoFormulaStrID = find(ismember(model.metFormulas,NoFormulaStr));


% check for unique metabolites. before processing mets into metList
% metabolites such as nadh_c and nadh_m existed. since processing there are
% multiples in that list (i.e. nadh and nadh etc.)
% uList(ic) = metList; hence ic can be used to reconstruct the carbon
% vector
[uList,ia,ic] = unique(metList); % uList the uniq mets, ia the indexes of the 
                                % unique mets and ic the potition of
                                % uniques at original.

% extract metFormul;as into variable metForm
metForm = model.metFormulas(ia);

%% Detect Elements
% loop through all unique metabolites (ia = number of unique mets)
for i = 1:size(ia,1)
    %% Detect Carbons
    carbon_tmp = regexp(metForm(i,1),'(?<=C)\d+','match');
    if isempty(carbon_tmp{1,1})
        carbon_tmp{1,1} = 0;
    end
    isCarbon = regexpi(metForm(i,1),'(?-i)C[A-Z]');
    if isempty(isCarbon{1,1})
        carbon2 = 0;
    else
        carbon2 = size(isCarbon{1,1},2);
    end
    carbon_tmp = str2double(carbon_tmp{1,1});
    if isnan(carbon_tmp)
        carbon_tmp = 0;
    end
    carbon(i,1) = sum(carbon_tmp) + carbon2;
    
    %% Detect Hydrogens
    hydrogen_tmp = regexp(metForm(i,1),'(?<=H)\d+','match');
    isHydrogen = regexpi(metForm(i,1),'(?-i)H[A-Z]');
    SingleH = regexpi(metForm(i,1),'^\H{1}$');
    endH = regexpi(metForm(i,1),'\H$','match');    
    if isempty(endH{1,1})
        endH{1,1} = 0;
        endH = str2double(endH{1,1});
        if isnan(endH)
            endH = 0;
        end
    else    
        endH = size(endH,1);
    end    
    if isempty(SingleH{1,1})
        hydrogen3 = 0;
    else
        if (endH ~= 1)
        hydrogen3 = size(SingleH{1,1},2);
        end
    end    
    if isempty(hydrogen_tmp{1,1})
        hydrogen_tmp{1,1} = 0;
    end

    if isempty(isHydrogen{1,1})
        hydrogen2 = 0;
    else
        hydrogen2 = size(isHydrogen{1,1},2);
    end    
    if isempty(hydrogen_tmp{1,1})
        hydrogen_tmp = 0;
    end
    hydrogen_tmp = str2double(hydrogen_tmp{1,1});
    if isnan(hydrogen_tmp)
        hydrogen_tmp = 0;
    end    
    if isnan(carbon_tmp)
        hydrogen_tmp = 0;
    end
    hydrogen(i,1) = sum(hydrogen_tmp) + hydrogen2 + hydrogen3 + endH;
        
    %% Detect Oxygens
    oxygen_tmp = regexp(metForm(i,1),'(?<=O)\d+','match');
    isOxygen = regexpi(metForm(i,1),'(?-i)O[A-Z]');
    SingleO = regexpi(metForm(i,1),'^\O{1}$');
    endO = regexpi(metForm(i,1),'\O$','match');    
    if isempty(endO{1,1})
        endO{1,1} = 0;
        endO = str2double(endO{1,1});
        if isnan(endO)
            endO = 0;
        end
    else    
        endO = size(endO,1);
    end
   if isempty(isOxygen{1,1})
        oxygen2 = 0;
    else
        oxygen2 = size(isOxygen{1,1},2);
   end  
    if isempty(SingleO{1,1})
        oxygen3 = 0;
    else
        oxygen3 = size(SingleO{1,1},2);
    end
    if ((size(oxygen_tmp{1,1},2)) >= 2)
        a = oxygen_tmp{1,1};
        oxygen_tmp{1,1} = a{1,1} + a{1,2};
    end  
    if isempty(oxygen_tmp{1,1})
        oxygen_tmp{1,1} = 0;
    end    
    oxygen_tmp = str2double(oxygen_tmp{1,1});       
    if isnan(oxygen_tmp)
        oxygen_tmp = 0;
    end
    oxygen(i,1) = sum(oxygen_tmp) + oxygen2 + oxygen3 + endO;
    
    %% Detect Nitrogens
    nitrogen_tmp = regexp(metForm(i,1),'(?<=N)\d+','match');
    isNitrogen = regexpi(metForm(i,1),'(?-i)N[A-Z]');
    SingleN = regexpi(metForm(i,1),'^N{1}$');
    endN = regexpi(metForm(i,1),'(?-i)\N$','match');    
    if isempty(endN{1,1})
        endN{1,1} = 0;
        endN = str2double(endN{1,1});
        if isnan(endN)
            endN = 0;
        end
    else    
        endN = size(endN,1);
    end
    if isempty(isNitrogen{1,1})
        nitrogen2 = 0;
    else
        nitrogen2 = size(isNitrogen{1,1},2);
    end    
    if isempty(SingleN{1,1})
        nitrogen3 = 0;
    else
        nitrogen3 = size(SingleN{1,1},2);
    end
    if isempty(nitrogen_tmp{1,1})
        nitrogen_tmp{1,1} = 0;
    end
    nitrogen_tmp = str2double(nitrogen_tmp{1,1});    
    if isnan(nitrogen_tmp)
        nitrogen_tmp = 0;
    end
    nitrogen(i,1) =  sum(nitrogen_tmp) + nitrogen2 + nitrogen3 + endN;
    
    %% Detect Phosphorus
    phosphorus_tmp = regexp(metForm(i,1),'(?<=P)\d+','match');
    isPhosporus = regexpi(metForm(i,1),'(?-i)P[A-Z]');
    SingleP = regexpi(metForm(i,1),'^\P{1}$');
    endP = regexpi(metForm(i,1),'\P$','match');
     if isempty(endP{1,1})
        endP{1,1} = 0;
        endP = str2double(endP{1,1});
        if isnan(endP)
            endP = 0;
        end
    else    
        endP = size(endP,1);
     end         
    if isempty(isPhosporus{1,1})
        phosphorus2 = 0;
    else
        phosphorus2 = size(isPhosporus{1,1},2);
    end        
    if isempty(SingleP{1,1})
        phosphorus3 = 0;
    else
        phosphorus3 = size(SingleP{1,1},2);
    end
    if isempty(phosphorus_tmp{1,1})
        phosphorus_tmp{1,1} = 0;
    end    
    phosphorus_tmp = str2double(phosphorus_tmp{1,1});    
    if isnan(phosphorus_tmp)
        phosphorus_tmp = 0;
    end
    phosphorus(i,1) = sum(phosphorus_tmp) + phosphorus2 + phosphorus3 + endP;
    
    %% Detect Sulfur
    sulfur_tmp = regexp(metForm(i,1),'(?<=S)\d+','match');
    isSulfur = regexpi(metForm(i,1),'(?-i)S[A-Z]');
    SingleS = regexpi(metForm(i,1),'^S{1}$');
    endS = regexpi(metForm(i,1),'S$','match');    
    if isempty(endS{1,1})
        endS{1,1} = 0;
        endS = str2double(endS{1,1});
        if isnan(endS)
            endS = 0;
        end
    else    
        endS = size(endS,1);
    end       
    if isempty(isSulfur{1,1})
        sulfur2 = 0;
    else
        sulfur2 = size(isSulfur{1,1},2);
    end
        
    if isempty(SingleS{1,1})
        sulfur3 = 0;
    else
        sulfur3 = size(SingleS{1,1},2);
    end
    if isempty(sulfur_tmp{1,1})
        sulfur_tmp{1,1} = 0;
    end    
    sulfur_tmp = str2double(sulfur_tmp{1,1});    
    if isnan(sulfur_tmp)
        sulfur_tmp = 0;
    end    
    sulfur(i,1) = sum(sulfur_tmp) + sulfur2 + sulfur3 + endS;
    
    %% Detect Calcium
    SingleCa = regexpi(metForm(i,1),'^\Ca{1}$');     
    if isempty(SingleCa{1,1})
        calcium3 = 0;
    else
        calcium3 = size(SingleCa{1,1},2);
    end
    calcium(i,1) = calcium3;
    
    %% Detect Potassium
    SingleK = regexpi(metForm(i,1),'^\K{1}$');
    if isempty(SingleK{1,1})
        potassium3 = 0;
    else
        potassium3 = size(SingleK{1,1},2);
    end
    potassium(i,1) = potassium3;
    
    %% Detect Sodium
    SingleNa = regexpi(metForm(i,1),'^\Na{1}$');
    if isempty(SingleNa{1,1})
        sodium3 = 0;
    else
        sodium3 = size(SingleNa{1,1},2);
    end
    sodium(i,1) =sodium3;
    
    %% Detect Chlorine
    isChlorine = regexpi(metForm(i,1),'(?-i)Cl[A-Z]');
    SingleCl = regexpi(metForm(i,1),'^\Cl{1}$');
    if isempty(isChlorine{1,1})
        chlorine2 = 0;
    else
        chlorine2 = size(isChlorine{1,1},2);
    end    
    if isempty(SingleCl{1,1})
        chlorine3 = 0;
    else
        chlorine3 = size(SingleCl{1,1},2);
    end   
    chlorine(i,1) = chlorine2 + chlorine3;
    
    %% Detect Zinc
    SingleZn = regexpi(metForm(i,1),'^\Zn{1}$');         
    if isempty(SingleZn{1,1})
        zinc3 = 0;
    else
        zinc3 = size(SingleZn{1,1},2);
    end
    zinc(i,1) =zinc3;
    
    %% Iodine
    iodine_tmp = regexp(metForm(i,1),'(?<=I)\d+','match');
    isIodine = regexpi(metForm(i,1),'(?-i)I[A-Z]');
    SingleI = regexpi(metForm(i,1),'[A-Z]I{1}$');
    endI = regexpi(metForm(i,1),'I$','match');    
    if isempty(endI{1,1})
        endI{1,1} = 0;
        endI = str2double(endI{1,1});
        if isnan(endI)
            endI = 0;
        end
    else    
        endI = size(endI,1);
    end       
    if isempty(isIodine{1,1})
        iddine2 = 0;
    else
        iddine2 = size(isIodine{1,1},2);
    end
        
    if isempty(SingleI{1,1})
        iodine3 = 0;
    else
        iodine3 = size(SingleI{1,1},2);
    end
    if isempty(iodine_tmp{1,1})
        iodine_tmp{1,1} = 0;
    end    
    iodine_tmp = str2double(iodine_tmp{1,1});    
    if isnan(iodine_tmp)
        iodine_tmp = 0;
    end    
    iodine(i,1) = sum(iodine_tmp) + iddine2 + iodine3 + endI;
    
    %% Molar Weight Calculation
    
    C = 12.01;
    H = 1.01;
    O = 16;
    N = 14.01;
    P = 30.97;
    S = 32.07;
    Ca = 40.08;
    K = 39.10;
    Na = 23;
    Cl = 35.45;
    Zn = 65.38;
    I = 126.90;
    
    C_M = carbon(i,1) * C;
    H_M = hydrogen(i,1) * H;
    O_M = oxygen(i,1) * O;
    N_M = nitrogen(i,1) * N;
    P_M = phosphorus(i,1) * P;
    S_M = sulfur(i,1) * S;
    Ca_M = calcium(i,1) * Ca;
    K_M = potassium(i,1) * K;
    Na_M = sodium(i,1) * Na;
    Cl_M = chlorine(i,1) * Cl;
    Zn_M = zinc(i,1) * Zn;
    I_M = iodine(i,1) * I;

    mass(i,1) = C_M + H_M + O_M + N_M + P_M + S_M + Ca_M + K_M +Na_M + Cl_M + Zn_M + I_M;     

end

% metabolite name met formulas and number of carbon atoms in table
carbonList = table([1:met_num]',mets,metForm(ic),mass(ic), carbon(ic), hydrogen(ic), ...
    oxygen(ic), nitrogen(ic), phosphorus(ic), sulfur(ic), calcium (ic), ...
    potassium(ic), sodium(ic), chlorine(ic), zinc(ic), iodine(ic), 'VariableNames', ...
    {'metID','metName','chemFormula','mass', 'carbons', ' hydrogens' ...
    ' oxygens', 'nitrogens', 'phosphorus', 'sulfur', 'calcium', ' potassium', ...
    'sodium', 'chlorium', 'zinc', ' iodine' ...
    });

% constructing output structure cStruc
cStruc.title = '*** Met Carbon Count ***';
cStruc.dateAndTime = datetime();


if isfield(model,'description')
    cStruc.model = model.description;
end


cStruc.metName = mets;
cStruc.metForm = metForm(ic);
cStruc.Mass = mass(ic);
cStruc.carbonAtoms = carbon(ic);
cStruc.hydrogenAtoms = hydrogen(ic);
cStruc.oxygenAtoms = oxygen(ic);
cStruc.nitrogenAtoms = nitrogen(ic);
cStruc.phosporusAtoms = phosphorus(ic);
cStruc.sulfurAtoms = sulfur(ic);
cStruc.calciumAtoms = calcium(ic);
cStruc.potassiumAtoms = potassium(ic);
cStruc.sodiummAtoms = sodium(ic);
cStruc.chlorineAtoms = chlorine(ic);
cStruc.zincAtoms = zinc(ic);
cStruc.iodineAtoms = iodine(ic);

cStruc.table = carbonList;