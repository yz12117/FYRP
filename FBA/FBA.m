clc;clear;close all;
solverName = 'gurobi';
solverType = 'LP'; 
changeCobraSolver(solverName, solverType);
model = load('iCHOv1.mat');
model = model.iCHOv1;
model = findSExRxnInd(model);
max = 1000;
model_original = model;
for F = 1:6
    if F == 6
        feed = 'average';
    else
        feed = append('F',string(F));
    end
    for c = 1:3
        if c == 1 
            condition = 'early_exponential';
            direction = 'max';
            model.c(findRxnIDs(model,'BIOMASS_cho_producing')) = 1;
            model.c(findRxnIDs(model,'BIOMASS_cho')) = 0;
            model.c(findRxnIDs(model,'DM_igg_g')) = 0;
        elseif c == 2
            condition = 'late_exponential';
            direction = 'max';
            model.c(findRxnIDs(model,'BIOMASS_cho_producing')) = 1;
            model.c(findRxnIDs(model,'BIOMASS_cho')) = 0;
            model.c(findRxnIDs(model,'DM_igg_g')) = 0;
        elseif c == 3
            condition = 'stationary';
            direction = 'max';
            obj = 'EX_met__L_e';
            model.c(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
            model.c(findRxnIDs(model,'BIOMASS_cho')) = 0;
            model.c(findRxnIDs(model,'DM_igg_g')) = 0;
            model.c(findRxnIDs(model,obj)) = 1;
        end
        disp(append(feed,' ',condition,' starts'))
        try
            [model,bnds.(feed).(condition)] = setExperimentBounds_final(model,condition,max,F);
            if c == 3
                model.lb(findRxnIDs(model,obj)) = -1000;
                model.ub(findRxnIDs(model,obj)) = 1000;
                for upts = 1:length(bnds.(feed).(condition))
                    if strcmp(bnds.(feed).(condition){upts},obj)
                        bnds.(feed).(condition){upts,1} = [];
                        bnds.(feed).(condition){upts,2} = [];
                        bnds.(feed).(condition){upts,3} = [];
                    end
                end
            end
            for column = 1:3
                bnds_new.(feed).(condition)(:,column) = bnds.(feed).(condition)(~cellfun(@isempty,bnds.(feed).(condition)(:,column)),column);
            end
            if c == 3
                bnds.(feed).(condition) = bnds_new.(feed).(condition);
                FBASolution_original.(feed).(condition) = optimizeCbModel(model,direction);
                disp(FBASolution.(feed).(condition).f)
            else
                bnds.(feed).(condition) = bnds_new.(feed).(condition);
                FBASolution_original.(feed).(condition) = optimizeCbModel(model,direction);
                disp(FBASolution_original.(feed).(condition).f)
                minmax.(feed).(condition) = runMinMax_GF(model);
                for i = 1:length(bnds.(feed).(condition))
                    if (bnds.(feed).(condition){i,2} == -1000) && (bnds.(feed).(condition){i,3} == 1000)
                        bnds.(feed).(condition){i,2} = minmax.(feed).(condition)(findRxnIDs(model,bnds.(feed).(condition)(i,1)),1);
                        bnds.(feed).(condition){i,3} = minmax.(feed).(condition)(findRxnIDs(model,bnds.(feed).(condition)(i,1)),2);
                    end
                end
                model_FVA.(feed).(condition) = model;
                model.rev = checkReversability(model);
                skipsubsystem = {'OXIDATIVE PHOSPHORYLATION'};
                skipRxns = find(ismember(model.subSystems,skipsubsystem));
                model_carbon.(feed).(condition) = model;
                carbonConst.(feed).(condition) = carbonConstraints_final_final(model,bnds.(feed).(condition),condition,[],[],minmax.(feed).(condition),false,[],[],skipRxns,direction);
                carbon_solution.(feed).(condition) = carbonConst.(feed).(condition).solutionNew.f;
                model_nitrogen.(feed).(condition) = model;
                nitrogenConst.(feed).(condition) = nitrogenConstraints_final(model,bnds.(feed).(condition),condition,[],[],minmax.(feed).(condition),false,[],[],skipRxns,direction);
                nitrogen_solution.(feed).(condition) = nitrogenConst.(feed).(condition).solutionNew.f;
                sulphurConst.(feed).(condition) = sulphurConstraints_final_final(model,bnds.(feed).(condition),condition,[],[],minmax.(feed).(condition),false,[],[],skipRxns,direction);
                sulphur_solution.(feed).(condition) = sulphurConst.(feed).(condition).solutionNew.f;
                NEWminmax.(feed).(condition) = selectMinMax(carbonConst.(feed).(condition),nitrogenConst.(feed).(condition),sulphurConst.(feed).(condition));
                model.lb = cell2mat(NEWminmax.(feed).(condition)(:,1)); 
                model.ub = cell2mat(NEWminmax.(feed).(condition)(:,2));
                FBASolution.(feed).(condition) = optimizeCbModel(model,direction);
                disp(FBASolution.(feed).(condition).f)
    %             pFBA_model.(feed).(condition) = model;
    %             [GeneClasses.(feed).(condition),RxnClasses.(feed).(condition),modelIrrevFM.(feed).(condition)] = pFBA(model,'geneoption',0);
    %             modelIrrevFM.(feed).(condition).c = zeros(length(modelIrrevFM.(feed).(condition).rxns),1);
    %             if c == 3
    %                 modelIrrevFM.(feed).(condition).c(findRxnIDs(modelIrrevFM.(feed).(condition),'DM_igg_g')) = 1;
    %             else
    %                 modelIrrevFM.(feed).(condition).c(findRxnIDs(modelIrrevFM.(feed).(condition),'BIOMASS_cho_producing')) = 1;
    %             end
    %             pFBASolution.(feed).(condition) = optimizeCbModel(modelIrrevFM.(feed).(condition));

    %             model_input.(feed).(condition) = convertToReversible(modelIrrevFM.(feed).(condition));
    %             model_input.(feed).(condition) = removeRxns(model_input.(feed).(condition),'netFlux');
    %             model_input.(feed).(condition) = model;
                model = model_original;
                disp(append(feed,' ',condition,' is done'))
            end
        catch
            disp(append(feed,' ',condition,' is infeasible'))
        end
    end
end

% date = datetime('now', 'Format','yyyy_MM_dd_HH_mm_ss');
% FileName = append('Output',string(date),'.mat');
% save(FileName)
% 
% for F = 3
%     feed = append('F',string(F));
%     if F == 6
%         feed = 'average';
%     end
%     for c = 3
%         if c == 1 
%             condition = 'early_exponential';
%         elseif c == 2
%             condition = 'late_exponential';
%         elseif c == 3
%             condition = 'stationary';
%         end
% %         try
%             modelGAMS = PrintToGams(model_input.(feed).(condition),[],FileName,feed,condition,FBASolution_original.(feed).(condition).f);
% %             continue
% %         catch
% %             continue
% %         end
%     end
% end
% %% gene intervention
% feed = 'F3';
% condition = 'stationary';
% c = 3;
% F = 3;
% obj = 'EX_met__L_e';
% [model,bnds.(feed).(condition)] = setExperimentBounds_final(model,condition,max,F);
% % igg_lb = 0:1e-5:3.3e-4;
% % igg_wild = zeros(length(igg_lb),1);
% % biomass_wild = zeros(length(igg_lb),1);
% % for i = 1:length(igg_lb)
% %     model.c(findRxnIDs(model,'BIOMASS_cho')) = 0;
% %     model.c(findRxnIDs(model,'BIOMASS_cho_producing')) = 1;
% %     model.lb(findRxnIDs(model,'DM_igg_g')) = igg_lb(i);
% %     model.lb(findRxnIDs(model,'BIOMASS_cho_producing')) = 0;
% %     originalSolution = optimizeCbModel(model,'max');
% %     try
% %     igg_wild(i) = originalSolution.v(findRxnIDs(model,'DM_igg_g'));
% %     biomass_wild(i) = originalSolution.v(findRxnIDs(model,'BIOMASS_cho_producing'));
% %     continue
% %     catch
% %         continue
% %     end
% % end
% %plot(biomass_wild,igg_wild)
% if c == 3
%     model.lb(findRxnIDs(model,obj)) = -1000;
%     model.ub(findRxnIDs(model,obj)) = 1000;
%     model.c = zeros(length(model.c),1);
%     model.c(findRxnIDs(model,obj)) = 1;
%     for upts = 1:length(bnds.(feed).(condition))
%         if strcmp(bnds.(feed).(condition){upts},obj)
%             met_ub = bnds.(feed).(condition){upts,2};
%             bnds.(feed).(condition){upts,1} = [];
%             bnds.(feed).(condition){upts,2} = [];
%             bnds.(feed).(condition){upts,3} = [];
%         end
%     end
% end
% orginalRevSolution = optimizeCbModel(model,'max');
% for column = 1:3
%     bnds_new.(feed).(condition)(:,column) = bnds.(feed).(condition)(~cellfun(@isempty,bnds.(feed).(condition)(:,column)),column);
% end
% bnds.(feed).(condition) = bnds_new.(feed).(condition);
% modelIrrev = convertToIrreversible(model,'flipOrientation',true);
% modelIrrev.c = zeros(length(modelIrrev.c),1);
% modelIrrev.c(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 1;
% modelIrrev.lb(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 0;%-1*orginalRevSolution.f;
% modelIrrev.ub(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 0.004110655;%1000;%-1*met_ub;
% % modelIrrev.ub(findRxnIDs(modelIrrev,'EX_met__L_e_f')) = 0.002226915;
% modelIrrev.ub(findRxnIDs(modelIrrev,'DM_igg_g')) = 1000;
% modelIrrev.ub(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')) = 1000;
% igg_lb = 0:1e-5:3.3e-4;
% biomass_lb = 0:0.01:0.1;
% igg_wild = [];
% biomass_wild = [];
% for i = 1:length(igg_lb)
%     for j =1:length(biomass_lb)
%     modelIrrev.lb(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')) = biomass_lb(j);
%     modelIrrev.lb(findRxnIDs(modelIrrev,'DM_igg_g')) = igg_lb(i);
%     originalSolution = optimizeCbModel(modelIrrev,'min');
%         try
% %             if originalSolution.f < -0.004110655*-1
%                 igg_wild = cat(1,igg_wild,originalSolution.v(findRxnIDs(modelIrrev,'DM_igg_g')));
%                 biomass_wild = cat(1,biomass_wild,originalSolution.v(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')));
% %             end
%             continue
%         catch
%             continue
%         end
%     end
% end
% plot(biomass_wild,igg_wild,'-')
% modelIrrev.lb(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 0;
% cd(sprintf('GamsFiles/InputForGAMS/Pavlos/%s/%s',feed,condition))
% fileID = fopen('ModelRXNS.txt','r');
% reactions = fscanf(fileID,'%c');
% reactions = split(reactions,"'");
% reaction = [];
% for i = 1:length(reactions)/2
%     reaction = cat(1,reaction,reactions(i*2));
% end
% cd('C:\Users\dyycz\OneDrive - Imperial College London\Year 4\CENG97005 - Advanced Chemical Engineering Practice Research Project 2020-2021\Research project\MATLAB\Elemental Constrain')
% e_up = 0.7;
% e_down = 0.3;
% %% txt file
% % %%%%% Increase the:
% % b = [];
% % for i=1:0
% % eval(sprintf('b = [b; r%d];',i));
% % end
% % r_up = b;
% % %%%%% Knock-out the:
% % r1=331;
% % r2=2513;
% % r3=2700;
% % r4=3635;
% % r5=4995;
% % r6=5190;
% % r7=6921;
% % r8=7561;
% % r9=7740;
% % r10=8456;
% % b = [];
% % for i=1:10
% % eval(sprintf('b = [b; r%d];',i));
% % end
% % r_null = b;
% % %%%%% Decrease the:
% % b = [];
% % for i=11:10
% % eval(sprintf('b = [b; r%d];',i));
% % end
% % r_down = b;
% %% txt file
% % %%%%% Increase the:
% % r1=8446;
% % b = [];
% % for i=1:1
% % eval(sprintf('b = [b; r%d];',i));
% % end
% % r_up = b;
% % %%%%% Knock-out the:
% % r2=934;
% % r3=2693;
% % r4=3228;
% % r5=3535;
% % r6=3537;
% % r7=3558;
% % r8=4816;
% % r9=5190;
% % r10=7569;
% % b = [];
% % for i=2:10
% % eval(sprintf('b = [b; r%d];',i));
% % end
% % r_null = b;
% % %%%%% Decrease the:
% % b = [];
% % for i=11:10
% % eval(sprintf('b = [b; r%d];',i));
% % end
% % r_down = b;
% %% txt
% %%%%% Increase the:
% b = [];
% for i=1:0
% eval(sprintf('b = [b; r%d];',i));
% end
% r_up = b;
% %%%%% Knock-out the:
% r1=1;
% r2=9;
% r3=215;
% r4=775;
% r5=2560;
% r6=3119;
% r7=3228;
% r8=3519;
% r9=4446;
% r10=6686;
% b = [];
% for i=1:10
% eval(sprintf('b = [b; r%d];',i));
% end
% r_null = b;
% %%%%% Decrease the:
% b = [];
% for i=11:10
% eval(sprintf('b = [b; r%d];',i));
% end
% r_down = b;
% %% reactions for intervention
% rxns_up =reaction(r_up);
% rxns_KO = reaction(r_null);
% rxns_down = reaction(r_down);
% %% FBA validation
% lb_up = modelIrrev.lb(findRxnIDs(modelIrrev,rxns_up));
% ub_up = modelIrrev.ub(findRxnIDs(modelIrrev,rxns_up));
% lb_down = modelIrrev.lb(findRxnIDs(modelIrrev,rxns_down));
% ub_down = modelIrrev.ub(findRxnIDs(modelIrrev,rxns_down));
% modelIrrev.lb(findRxnIDs(modelIrrev,rxns_up)) = lb_up + (ub_up-lb_up).*e_up;
% modelIrrev.ub(findRxnIDs(modelIrrev,rxns_down)) = ub_down - (ub_down-lb_down).*e_down;
% modelIrrev.lb(findRxnIDs(modelIrrev,rxns_KO)) = 0;
% modelIrrev.ub(findRxnIDs(modelIrrev,rxns_KO)) = 0;
% modelIrrev.ub(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')) = 1000;
% modelIrrev.lb(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')) = 0;
% modelIrrev.ub(findRxnIDs(modelIrrev,'DM_igg_g')) = 1000;
% modelIrrev.lb(findRxnIDs(modelIrrev,'DM_igg_g')) = 0;
% modelIrrev.c = zeros(length(modelIrrev.c),1);
% modelIrrev.c(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 1;
% geneIntervention = optimizeCbModel(modelIrrev,'min');
% modelIrrev.lb(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 0;
% modelIrrev.ub(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 0.004110655;%-1*met_ub;%0.0186
% % modelIrrev.ub(findRxnIDs(modelIrrev,'EX_met__L_e_f')) = 0.002226915;
% modelIrrev.c(findRxnIDs(modelIrrev,'EX_met__L_e_b')) = 1;
% igg_lb = 0:1e-5:3.3e-4;
% biomass = 0:0.01:0.1;
% igg_mutant = [];
% biomass_mutant = [];
% for i = 1:length(igg_lb)
%     for j = 1:length(biomass_lb)
%         modelIrrev.lb(findRxnIDs(modelIrrev,'DM_igg_g')) = igg_lb(i);
%         modelIrrev.lb(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')) = biomass_lb(j);
%         geneIntervention = optimizeCbModel(modelIrrev,'min');
% %     disp(num2str(geneIntervention.v(findRxnIDs(modelIrrev,'BIOMASS_cho_producing'))))
% %     disp(num2str(geneIntervention.v(findRxnIDs(modelIrrev,'DM_igg_g'))))
%         try
%             igg_mutant = cat(1,igg_mutant,geneIntervention.v(findRxnIDs(modelIrrev,'DM_igg_g')));
%             biomass_mutant = cat(1,biomass_mutant,geneIntervention.v(findRxnIDs(modelIrrev,'BIOMASS_cho_producing')));
%         catch
%             continue
%         end
%     end
% end
% hold on
% plot(biomass_mutant,igg_mutant,'--')
% legend('wild type CHO','mutant CHO')
% xlabel('Biomass growth rate')
% ylabel('IgG production rate')