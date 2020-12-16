function [model,bnds] = setExperimentBounds_final(model,condition,max,F)
if strcmp(F,'average')
    F = 6;
end
% condition: exponential or stationary
model.lb(find(model.ExchRxnBool)) = 0;
model.lb(find(model.DMRxnBool)) = 0;
model.lb(find(model.SinkRxnBool)) = -0.1;
model.ub(find(model.ExchRxnBool)) = max;
model.ub(find(model.DMRxnBool)) = max;
model.ub(find(model.SinkRxnBool)) = max;
lb_original = model.lb;
ub_original = model.ub;
data = readtable('Pavlos_FBA_Bounds_new.xlsx','sheet',condition,'HeaderLines',0); % read data from excel
metabolites_raw = data.Properties.VariableNames; % get the metabolites names from excel
metabolites = metabolites_raw(3:28); % remove extra elements
for i = 1:length(metabolites)
    metabolites{1,i} = lower(metabolites{1,i}); % change names to lower case
end
metabolites{1,8} = 'nh4'; % we assumed that the concentration of NH3 given in excel sheet equals to concentration of ammonium ion
ExcRxns = model.rxns(find(model.ExchRxnBool)); % find all exchange reactions
% printExcRxns(model) % print it out with reaction ids, equations, lower bounds and upper bounds
amino_acids = {'asp','ser','asn','glu','gly','gln','his','arg','thr','ala','leu','val','met','ile','pro','cys','tyr','lys','phe','trp'}; % a list of all amino acids in nature
amino_metabolites = {};
non_amino_metabolites = metabolites;
for i = 1:length(metabolites)
    for j = 1:length(amino_acids)
        if isequal(amino_acids{:,j},metabolites{:,i})
            amino_metabolites{end+1} = append(metabolites{:,i},'__L_e'); % change all names of amino acids involved to 'amino acid__L_e', as we only consider extracellular L-amino acid as the reactant of exchange reactions
        end
    end
end
for k = 1 : length(metabolites)
    for kk = 1 : length(amino_metabolites)
        if ismember(metabolites{k},amino_metabolites{kk})
            non_amino_metabolites{k} = []; % remove all metabolites that are amino acids
        end
    end
end
empties = find(cellfun(@isempty,non_amino_metabolites)); % identify the empty cells
non_amino_metabolites(empties) = []; % remove empty cells
for i = 1:length(non_amino_metabolites)
    non_amino_metabolites{:,i} = append(non_amino_metabolites{:,i},'_e'); % change all names of metabolites that are not amino acids involved to 'metabolites__L_e'
end
new_metabolites = [amino_metabolites,non_amino_metabolites]; % combine into 1 list
new_metabolites{:,5} = 'gly_e';
new_metabolites{:,24} = 'glc__D_e';
new_metabolites{:,25} = 'lac__L_e';
metabolites = sort(metabolites); % sort into alphabetic order
metabolites{1,17} = 'nh3'; % change it back to column names of excel
new_metabolites = sort(new_metabolites); % sort into alphabetic order
% https://doi.org/10.1016/j.cels.2016.10.020
model.lb(findRxnIDs(model,'EX_pydxn_e')) = -1.98e-5;
model.lb(findRxnIDs(model,'EX_gthox_e')) = 0.00205697;
model.lb(findRxnIDs(model,'EX_fol_e')) = -0.000337303;
% model.lb(findRxnIDs(model,'EX_btn_e')) = -3.38E-05;
model.lb(findRxnIDs(model,'EX_chol_e')) = -0.021895966;
model.lb(findRxnIDs(model,'EX_inost_e')) = -0.001988394;
for a = 1:length(metabolites)
    for b = 1:length(model.mets)
        if strcmp(model.mets{b,:},new_metabolites{:,a}) == 1
            if length(model.mets{b,:}) <= 8
                Rxns = findRxnsFromMets(model, new_metabolites{1,a}); % find the reactions involved extreacellular L-Asp amd print reaction ID and equation
                for i = 1:length(Rxns)
                    for j = 1:length(ExcRxns)
                        if isequal(ExcRxns{j,:},Rxns{i,:})
                           mets_ExcRxns = Rxns{i,:}; % find the exchange reactions involved this metabolite
                        else
                            continue
                        end
                    end
                end
                model.ub(findRxnIDs(model, mets_ExcRxns)) = table2array(data(2*F,upper(metabolites{1,a})));
                model.lb(findRxnIDs(model, mets_ExcRxns)) = table2array(data(2*F-1,upper(metabolites{1,a})));
            end
        end
    end
end
essential = {'EX_ac_e','EX_4hpro_LT_e','EX_thm_e','EX_ascb__L_e','EX_ribflv_e','EX_ncam_e','EX_sprm_e','EX_h_e','EX_h2o_e','EX_o2_e','EX_pi_e','EX_so4_e'};
essentialMets = findMetsFromRxns(model,essential);
essentialID = findMetIDs(model, essentialMets);
for i = 1:length(essentialID)
    metFormulas = model.metFormulas(essentialID(i));
    carbon_tmp = regexp(metFormulas,'(?<=C)\d+','match');
    if isempty(carbon_tmp{1,1})
        carbon_tmp{1,1} = 0;
    end
    isCarbon = regexpi(metFormulas(1,1),'(?-i)C[A-Z]');
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
    essentialRxns = findRxnsFromMets(model,model.mets(essentialID(i)));
    [Rxn,~,~] = intersect(essential,essentialRxns);
    if carbon(i,1) == 0
        model.lb(findRxnIDs(model,Rxn)) = -max;
    else
        model.lb(findRxnIDs(model,Rxn)) = model.lb(findRxnIDs(model,'EX_glc__D_e'))*6./carbon(i,1);
    end
end
lb_input = [];
ub_input = [];
rxnID_input = [];
lb = model.lb;
ub = model.ub;
rxnID = model.rxns;
for i = 1:length(model.rxns)
    if (lb(i) ~= lb_original(i)) || (ub(i) ~= ub_original(i))
        lb_input = vertcat(lb_input,lb(i));
        ub_input = vertcat(ub_input,ub(i));
        rxnID_input = vertcat(rxnID_input,rxnID(i));
    end
end
model.lb(findRxnIDs(model,'SK_Asn_X_Ser_Thr_r')) = -0.1;
model.lb(findRxnIDs(model,'SK_Ser_Thr_g')) = -0.1;
model.lb(findRxnIDs(model,'SK_Tyr_ggn_c')) = -0.1;
model.ub(findRxnIDs(model,'SK_Asn_X_Ser_Thr_r')) = 100;
model.ub(findRxnIDs(model,'SK_Ser_Thr_g')) = 100;
model.ub(findRxnIDs(model,'SK_Tyr_ggn_c')) = 100;
model.lb(findRxnIDs(model,'ATPM')) = 0;
model.ub(findRxnIDs(model,'ATPM')) = 100;
bnds = [rxnID_input,num2cell(lb_input),num2cell(ub_input)];
if F == 1
    if strcmp(condition,'exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.9319e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'stationary')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 0; %1.50656e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0.0077;%0.0065;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'early_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.9319E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.52465E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'late_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 2.58018E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.77862E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    end
elseif F == 2
    if strcmp(condition,'exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 0;%1.59272e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'stationary')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 0;%1.5146e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0.0077;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'early_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.59272E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.26762E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'late_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 2.45262E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.65874E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    end
elseif F == 3
    if strcmp(condition,'exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.76581e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'stationary')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.40666E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.10002e-5;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0.008;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'early_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.76581E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.31601E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'late_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 2.39289E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.62624E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    end
elseif F == 4
    if strcmp(condition,'exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.76098e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'stationary')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 0;%1.45915e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0.0079;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'early_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.76098E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.35723E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'late_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 2.56833E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.66714E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    end
elseif F==5
    if strcmp(condition,'exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.49578e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'stationary')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 0;%1.52683e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0.0077;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'early_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.49578E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.11858E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'late_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) =2.21044E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.71329E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    end
elseif F == 6
    if strcmp(condition,'exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.70944e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'stationary')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 0;%1.48276e-5;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 1000;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0.0078;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'early_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 1.70944E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.31682E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    elseif strcmp(condition,'late_exponential')
        model.lb(findRxnIDs(model,'DM_igg_g')) = 2.44089E-05;
        model.ub(findRxnIDs(model,'DM_igg_g')) = 2.6888E-05;
        model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
        model.ub(findRxnIDs(model,'BIOMASS_cho_producing')) = 1000;
    end
end
% igg_lb = model.lb(findRxnIDs(model,'DM_igg_g'));
% igg_ub = model.ub(findRxnIDs(model,'DM_igg_g'));
% biomass_lb = model.lb(findRxnIDs(model,'BIOMASS_cho_producing'));
% biomass_ub = model.ub(findRxnIDs(model,'BIOMASS_cho_producing'));
bnds = [bnds;{'SK_Asn_X_Ser_Thr_r',-0.1,100};{'SK_Ser_Thr_g',-0.1,100};{'SK_Tyr_ggn_c',-0.1,100};{'ATPM',8.6,100}];%{'DM_igg_g',igg_lb,igg_ub};{'BIOMASS_cho_producing',biomass_lb,biomass_ub}];
end