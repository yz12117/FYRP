function [modelIrrev] = PrintToGams(model,upt,FileName,feed,condition,solution)
[dateStr, timeStr] = getDateTimeStrings(date,clock);  

% doc = 'SUE'; % Susie FeedC
% doc = 'SAR_FeedC'; % Sarantos FeedC
% doc = 'SAR_FeedAll'; % Sarantos FeedAll
% doc = 'SAR_FeedAll40'; % Sarantos FeedAll plus 40%
doc = 'Pavlos'; % Combination of the previous set of data#

% doc = 'Min_Uptake'; % Minimal uptake scenario test

% MODEL = 'GENOME_SCALE';
MODEL = 'iCHO1766';
% MODEL = 'K1';
% MODEL = 'S';

% MODEL = 'GENOME_SCALE_NoEPO';
% MODEL = 'K1_NoEPO';

% doc = 'ECOLI_Core';
% MODEL = 'Core';

% doc = 'ECOLI_iJO1366';
% MODEL = 'iJO1366';

if exist('upt','var')
    mets_uptake = upt;
end

% FolderName = 'C:\Users\aa1018\Box Sync\PhD\Matlab Test Stuff\FVA Bounds\';
FolderName = 'C:\Users\dyycz\OneDrive - Imperial College London\Year 4\CENG97005 - Advanced Chemical Engineering Practice Research Project 2020-2021\Research project\MATLAB\Elemental Constrain';

dataset = FileName;
dataset = load(fullfile(FolderName, dataset));

%% IF E. Coli. then

% model.csense = model.csense';
% s1 = optimizeCbModel(model); s1.f 
% 
% % load FVA_EX_ECOLI_Core.mat;
% % load FVA_ST_ECOLI_Core.mat;
% % load FVA_EX_ECOLI_iJO1366.mat;
% % load FVA_ST_ECOLI_iJO1366.mat;

% % if (string(doc) == 'ECOLI_Core')
% %     phase = 2;
% %     atpm = 11;
% %     atpm_value = 5;
% % elseif string(doc) == 'ECOLI_iJO1366'
% %     phase = 1;
% %     atpm = 748;
% %     atpm_value = 10;
% % end
% % 
% % [lbg,ubg,rxng,model] = ecoli_Bounds(model,phase,max);
% % lbg = lbg';
% % ubg = ubg';
% % rxng = rxng';
% % mets_uptake = rxng;
% % 
% % % model.lb(atpm) = atpm_value;
% % % model = set_bounds(model,lbg,ubg,rxng);
% % model.lb(:,1) = FVA_EX(:,1);
% % model.ub(:,1) = FVA_EX(:,2);
% 
% % mets_uptake = {};
% a = checkDrainRxns(model);
% amax = size(a,1);
% k = 0;
% for i = 1:amax
%     model.lb(a(i));
%     if (a(i))==true
%         if model.lb(i)<0
%             k = k + 1;
%             mets_uptake(k) = model.rxns(i);
%         end
%     end
% end
% s2 = optimizeCbModel(model); s2.f 

%% IF CHO then

%%%% If CC bounds then -->
% load FVA_EX_SUE.mat
% load FVA_EX_SAR_FeedC.mat;
% load FVA_EX_SAR_FeedAll.mat;
% load FVA_EX_SAR_FeedAll40.mat;
% load FVA_EX_WID.mat

model.lb(:,1) = dataset.model_input.(feed).(condition).lb;
model.ub(:,1) = dataset.model_input.(feed).(condition).ub;

% [~,~,rxn,~] = FBA_Bounds(model,1,100,doc);
rxn = dataset.bnds.(feed).(condition)(:,1);
% rxn = { ...
% 'EX_asp_L_e_' ...;
% 'EX_asn_L_e_' ...;
% 'EX_gln_L_e_' ...;
% 'EX_glu_L_e_' ...;
% 'EX_nh4_e_' ...;
% 'EX_arg_L_e_' ...;
% 'EX_ser_L_e_' ...;
% 'EX_his_L_e_' ...;
% 'EX_pro_L_e_' ...;
% 'EX_gly_e_' ...;
% 'EX_ala_L_e_' ...;
% 'EX_thr_L_e_' ...;
% 'EX_tyr_L_e_' ...;
% 'EX_cys_L_e_' ...;
% 'EX_val_L_e_' ...;
% 'EX_met_L_e_' ...;
% 'EX_ile_L_e_' ...;
% 'EX_leu_L_e_' ...;
% 'EX_lys_L_e_' ...;
% 'EX_phe_L_e_' ...;
% 'EX_trp_L_e_' ...;
% 'EX_lac_L_e_' ...;
% 'EX_glc_e_' ...;
% 'EX_4hpro_e_' ...;
% 'EX_pydxn_e_' ...;
% 'EX_thm_e_' ...;
% 'EX_gthox_e_' ...;
% 'EX_so4_e_' ...;
% 'EX_ascb_L_e_' ...;
% 'EX_fol_e_' ...;
% 'EX_ribflv_e_' ...;
% 'EX_btn_e_' ...;
% 'EX_chol_e_' ...;
% 'EX_ncam_e_' ...;
% 'EX_inost_e_' ...;
% 'EX_pyr_e_' ...;
% 'EX_sprm_e_' ...;
% 'EX_hco3_e_' ...;
% 'EX_ac_e_' ...;
% 'EX_h_e_' ...;
% 'EX_h2o_e_' ...;
% 'EX_o2_e_' ...;
% 'EX_pi_e_' ...;
% 'EX_so4_e_' ...;
% 'sink_Asn_X_Ser/Thr[r]' ...;
% 'sink_Ser/Thr[g]' ...;
% 'sink_Tyr_ggn[c]' ...;
% 'DM_atp[c]' ...;
%     };

mets_uptake = rxn(1:end-2);


% [lbg,ubg,rxng,model] = FBA_Bounds(model,1,max,doc);
% lbg = lbg';
% ubg = ubg';
% rxng = rxng';
% model = set_bounds(model,lbg,ubg,rxng);
% mets_uptake = rxng';

% minmaxg = runMinMax_GF(model); 
% minmaxg = fixMinMax(minmaxg);  
% model.lb = minmaxg(:,1);
% model.ub = minmaxg(:,2);

%% Model Initialization and shit
% % Reversible model
% isDrain=checkDrainRxns(model);
% transports = checkTheTransport(model);
% rev = model.rev;
% [modelPseudo,upt] = PseudoRxns(model,isDrain,rev,1,mets_uptake);
% mets = model.mets;
% rxns = model.rxns;
% % genes_Pseudo = modelPseudo.genes;
% mets_Pseudo = modelPseudo.mets;
% rxns_Pseudo = modelPseudo.rxns;
% S_Pseudo = modelPseudo.S;
% % GeneS_Pseudo = modelPseudo.rxnGeneMat;
% ub = model.ub;
% lb = model.lb;
% uptakes = modelPseudo.rxns(upt);
% Drains = modelPseudo.rxns(isDrain);
% % exchanges = modelPseudo.rxns(~ isDrain);
% rxns_noeffect = NoRemoveRXN(modelPseudo)


% % Irreversible model
[modelIrrev, matchRev, ~, irrev2rev] = convertToIrreversible(model,'flipOrientation',true);
mets = modelIrrev.mets;
rxns = modelIrrev.rxns;
isDrain = contains(modelIrrev.rxns,'EX_');

[mets_uptake] = IrrevUptakes(model,mets_uptake,matchRev,modelIrrev);
% SubSrxn = ModelSubSystem(model);

% transports = checkTheTransport(modelIrrev);
[modelPseudo,upt] = PseudoRxns(modelIrrev,isDrain,matchRev,2,mets_uptake,irrev2rev);

% genes_Pseudo = modelPseudo.genes; 
mets_Pseudo = modelPseudo.mets;
rxns_Pseudo = modelPseudo.rxns;
S_Pseudo = modelPseudo.S;
% GeneS_Pseudo = modelPseudo.rxnGeneMat;
ub = modelPseudo.ub;
lb = modelPseudo.lb;
% Drains = modelPseudo.rxns(isDrain);
uptakes = modelPseudo.rxns(upt);
rxns_noeffect = findOrphanRxns(modelPseudo);%findOrphanRxns(model)
biomass = find(ismember(model.rxns,'BIOMASS_cho_producing'));
IgG = find(ismember(model.rxns,'DM_igg_g'));
mkdir(sprintf('GamsFiles/%s/%s/%s/%s','InputForGAMS',doc,feed,condition))
cd(sprintf('GamsFiles/%s/%s/%s/%s','InputForGAMS',doc,feed,condition))
% mkdir GamsFiles\InputForGAMS\Pavlos\      
% cd GamsFiles\InputForGAMS\Pavlos      
% path_f = sprintf('\\\\ce-samba.ce.ic.ac.uk\\aa1018\\gamsproject\\DoubleOptTest\\RKO_Dynamic\\RKO_Dynamic_%d\\InputForGAMS\\SarantosFeedC',day);
% cd(path_f)

%% Print METS to TXT, the cmp file
fileMET = fopen('ModelMETS.txt','w');
fprintf(fileMET,'/ \n');
imax = size(mets,1);
for i = 1:imax
    x = mets{i};
    fprintf(fileMET,'''%s'' \n',x);
end
fprintf(fileMET,'/');
fclose(fileMET);

disp("Stage 1")

%% Print RXNS to TXT
fileRXN = fopen('ModelRXNS.txt','w');
fprintf(fileRXN,'/ \n');
imax = size(rxns,1);
for i = 1:imax
    x = rxns{i};
    fprintf(fileRXN,'''%s'' \n',x);
end
fprintf(fileRXN,'/');
fclose(fileRXN);

disp("Stage 2")

%% Print RXNames to TXT, the rxnnames file
fileRXNames = fopen('ModelRXNames.txt','w');
fprintf(fileRXNames,'/ \n');
imax = size(rxns_Pseudo,1);
for i = 1:imax
    x = rxns_Pseudo{i};
    fprintf(fileRXNames,'''%s'' \n',x);
end
fprintf(fileRXNames,'/');
fclose(fileRXNames);

disp("Stage 3")

%% Print Uptakes to TXT, the rxnnames file
fileUPTAKEs = fopen('ModelUptakes.txt','w');
fprintf(fileUPTAKEs,'/ \n');
imax = size(uptakes,1);
for i = 1:imax
    x = uptakes{i};
    fprintf(fileUPTAKEs,'''%s'' \n',x);
end
fprintf(fileUPTAKEs,'/');
fclose(fileUPTAKEs);

disp("Stage 4")

%% Print RXNStype to TXT, the rxntype file (Irreversible model)
fileRXNtype = fopen('ModelRXNStype.txt','w');
fprintf(fileRXNtype,'/ \n');
imax = size(rxns_Pseudo,1);
for i = 1:imax
    x = rxns_Pseudo{i};
    if (isDrain(i))    
        % Exchange
        fprintf(fileRXNtype,'''%s'' 4 \n',x);
    elseif (matchRev(i) == 0)
        % IRREVERSIBLE
        fprintf(fileRXNtype,'''%s'' 0 \n',x); 
    elseif (matchRev(i) < i)
        % Reversible BACKWARD
        fprintf(fileRXNtype,'''%s'' 2 \n',x);
    elseif (matchRev(i) > i)
        % Reversible FORWARD
        fprintf(fileRXNtype,'''%s'' 1 \n',x); 
    end
end
fprintf(fileRXNtype,'/');
fclose(fileRXNtype);

disp("Stage 5")

%% Print RXNStype to TXT, the rxntype file (Reversible model)
% fileRXNtype = fopen('ModelRXNStype.txt','w');
% fprintf(fileRXNtype,'/ \n');
% 
% imax = size(rxns_Pseudo,1);
% for i = 1:imax
%     x = rxns_Pseudo{i};
%     if (ismember(i,upt))
%         % UPTAKE
%         fprintf(fileRXNtype,'''%s'' 5 \n',x); 
%     elseif (isDrain(i) == 1)
% %         Exchange
%         fprintf(fileRXNtype,'''%s'' 4 \n',x);  
% %     elseif (transports(i) == 1)
% %         %         Transport 
% %         fprintf(fileRXNtype,'''%s'' 3 \n',x);
%     elseif (rev(i) == 0)
% %         Irreversible
%         fprintf(fileRXNtype,'''%s'' 2 \n',x);     
%     else
% %         Reversible 
%         fprintf(fileRXNtype,'''%s'' 1 \n',x);
%     end 
% end
% fprintf(fileRXNtype,'/');
% fclose(fileRXNtype);
% 
% disp("Stage 5")

%% Print S to TXT, the sij file
fileS = fopen('ModelS.txt','w');
fprintf(fileS,'/ \n');
imax = size(S_Pseudo,1);
jmax = size(S_Pseudo,2);
for j = 1:jmax
    for i = 1:imax
        if (S_Pseudo(i,j) ~= 0)
            m = mets_Pseudo{i};
            r = rxns_Pseudo{j};
            s = full(S_Pseudo(i,j));       
            fprintf(fileS,'''%s''.''%s'' %12.8f \n',m,r,s);
        end        
    end
end
fprintf(fileS,'/');
fclose(fileS);

disp("Stage 6")

%% Print rxnUB to TXT
str = sprintf("ModelUpperBound_%s.txt",doc);
fileUBtype = fopen(str,'w');
% fileUBtype = fopen('ModelUpperBound.txt','w');
fprintf(fileUBtype,'/ \n');
imax = size(rxns_Pseudo,1);
for i = 1:imax
    x = rxns_Pseudo{i};
    y = ub(i);
    fprintf(fileUBtype,'''%s''   %12.8f \n',x,y);
end
fprintf(fileUBtype,'/');
fclose(fileUBtype);

disp("Stage 7")

%% Print rxnLB to TXT
str = sprintf("ModelLowerBound_%s.txt",doc);
fileLBtype = fopen(str,'w');
% fileLBtype = fopen('ModelLowerBound.txt','w');
fprintf(fileLBtype,'/ \n');
imax = size(rxns_Pseudo,1);
for i = 1:imax
    x = rxns_Pseudo{i};
    y = lb(i);
    fprintf(fileLBtype,'''%s''   %12.8f \n',x,y);
end
fprintf(fileLBtype,'/');
fclose(fileLBtype);

disp("Stage 8")

%% Print GeneS to TXT, the sij file
% fileGeneS = fopen('ModelGenesS.txt','w');
% fprintf(fileGeneS,'/ \n');
% rmax = size(GeneS_Pseudo,1);
% gmax = size(GeneS_Pseudo,2);
% for j = 1:gmax % loop genes
%     for i = 1:rmax % loop reactions
%         if (GeneS_Pseudo(i,j) ~= 0)
%             g = genes_Pseudo{j};
%             r = rxns_Pseudo{i};
%             s = full(GeneS_Pseudo(i,j));       
%             fprintf(fileGeneS,'''%s''.''%s'' %12.8f \n',g,r,s);
%         end
%         
%     end
% 
% end
% fprintf(fileGeneS,'/');
% fclose(fileGeneS);
% 
% disp("Stage 9")

%% Print GENENames to TXT, the rxnnames file
% fileGENEames = fopen('ModelGENEames.txt','w');
% fprintf(fileGENEames,'/ \n');
% imax = size(genes_Pseudo,1);
% for i = 1:imax
%     x = genes_Pseudo{i};
%     fprintf(fileGENEames,'''%s'' \n',x);
% end
% fprintf(fileGENEames,'/');
% fclose(fileGENEames);
% 
% disp("Stage 10")

%% Print ExchangeReation to TXT, the E_RXNS file
% fileE_RXN = fopen('ModelE_RXNS.txt','w');
% fprintf(fileE_RXN,'/ \n');
% imax = size(exchanges,1);
% for i = 1:imax
%     x = exchanges{i};
%     fprintf(fileE_RXN,'''%s'' \n',x);
% end
% fprintf(fileE_RXN,'/');
% fclose(fileE_RXN);
% 
% disp("Stage 11")

%% Print EssentialReactions to TXT, the EssentialRxns file
% str = sprintf("EssentialReactions_%s.txt",MODEL);
% FileName   = str;
% FolderName = 'C:\Users\aa1018\Box Sync\PhD\Matlab Test Stuff\EssentialReactions';
% % FolderName = 'C:\Users\prassa\Box\PhD\Matlab Test Stuff\EssentialReactions';
% File       = fullfile(FolderName, FileName);
% EssentialReactions = load(File);
% fileEssRxns = fopen('ModelEssentialRxns.txt','w');
% fprintf(fileEssRxns,'/ \n');
% imax = size(EssentialReactions,1);
% j = 0;
% if (exist('matchRev'))
%     for i = 1:imax
%         k = find(matchRev==EssentialReactions(i));
%         if (~isempty(k))
%             j = j + 1;
%             EssentialReactions(end+1) = k;
%         end
%     end
% else
%     for i = 1:imax    
%         x = rxns_Pseudo{EssentialReactions(i)};
%         fprintf(fileEssRxns,'''%s'' \n',x);
%     end
% end
% if (j~=0)
%     for i = 1:imax+j
%         x = rxns_Pseudo{EssentialReactions(i)};
%         fprintf(fileEssRxns,'''%s'' \n',x);
%     end
% end
% fprintf(fileEssRxns,'/');
% fclose(fileEssRxns);
% 
% disp("Stage 12")

%% Print Drains to TXT, the drainRXN file
% filedrainRXN = fopen('ModelDrainRXN.txt','w');
% 
% fprintf(filedrainRXN,'/ \n');
% 
% imax = size(Drains,1);
% % imax = size(SubSrxn,1);
% 
% for i = 1:imax
%     x = Drains{i};
%     fprintf(filedrainRXN,'''%s''\n',x);
% end
% fprintf(filedrainRXN,'/');
% fclose(filedrainRXN);
% disp("Stage 13")

%% Print Compartments to TXT, the CompRXN file
% fileCompRXN = fopen('ModelCompRXN.txt','w');
% 
% fprintf(fileCompRXN,'/ \n');
% 
% imax = size(rxns_Pseudo,1);
% % imax = size(SubSrxn,1);
% 
% for i = 1:imax
%     Lia = ismember(SubSrxn(1,:),modelPseudo.subSystems(i));
%     num = find(Lia);    
%     x = rxns_Pseudo{i};
%     fprintf(fileCompRXN,'''%s'' %d \n',x,num);
% end
% fprintf(fileCompRXN,'/');
% fclose(fileCompRXN);
% disp("Stage 13")

%% Print ReversabilityMatrices to TXT, the ForReac and BackReac file
fileBackReac = fopen('ModelBackwardRxns.txt','w');
fprintf(fileBackReac,'/ \n');
imax = size(matchRev,1);
for i = 1:imax
    if (matchRev(i) > i)
        x = rxns_Pseudo{matchRev(i)};
        fprintf(fileBackReac,'''%s'' %i \n',x,matchRev(i));
    else
        x = rxns_Pseudo{i};
        fprintf(fileBackReac,'''%s'' 0 \n', x);
    end
end
fprintf(fileBackReac,'/');
fclose(fileBackReac);
disp("Stage 14")

%% Print NoGeneticAlteration to TXT
str = sprintf("ModelNoGenes.txt");
filengenea = fopen(str,'w');
fprintf(filengenea,'/ \n');
imax = size(rxns_noeffect,1);
for i = 1:imax
    %fprintf(filengenea,'''%s'' \n',modelPseudo.rxns{rxns_noeffect(i)}); 
    fprintf(filengenea,'''%s'' \n',string(rxns_noeffect(i))); 
end
fprintf(filengenea,'/');
fclose(filengenea);

disp("Stage 15")

%% Print Model Description to TXT
str = sprintf("DescriptionOfTheSimulation.txt");
filendis = fopen(str,'w');
fprintf(filendis,'/ \n');
fprintf(filendis,'Model used: %s \n',MODEL);
fprintf(filendis,'Experimental data set used: %s \n',doc);
fprintf(filendis,'CHECK THE BOUNDS TO SEE IF IT IS CARBONS OR NORMAL ONES!!!! \n');
fprintf(filendis,'Biomass is: %s \n', rxns_Pseudo{biomass});
fprintf(filendis,'IgG is: %s \n',rxns_Pseudo{IgG});
fprintf(filendis,'Max growth rate is: %f \n',solution);
fprintf(filendis,'Date the files were created: on the %s at time %s \n',string(dateStr), string(timeStr));
fprintf(filendis,'/');
fclose(filendis);

disp("Stage 16")

%%

cd('C:\Users\dyycz\OneDrive - Imperial College London\Year 4\CENG97005 - Advanced Chemical Engineering Practice Research Project 2020-2021\Research project\MATLAB\Elemental Constrain')

% cd 'C:\Users\aa1018\Box Sync\PhD\Matlab Test Stuff'
