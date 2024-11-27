%% MS results versus databases figures & MS versus Photometer figures
% Input for this script are the results from photometer and MS measurements
% of enzymatic activity when adding 79 potential effector molecules.
% First part of script converges the data of individual MS results from
% enzyme assays into 1 table, second part performs ROC curve analysis
% on MS data alone, third part imports photometer data and allows
% comparisons between the two datasets.
close all
clc
clear

%% INPUTS
analysis_path = "\\pasteur\SysBC-Home\chgruber\Desktop\==MANUSCRIPT WIP\==submission\=2=CODE\"
output_table = "existing" % "new" or "existing"
% "new" creates new table based on individual results of MS-based assays
% from output of =1_output=; advised if parts of analyses were changed and
% result tables differ from original ones.
% "existing" uses table created on individual results of MS-based assays
% with described parameters; advised if no changes were made to first part
% of analysis.


% create tables
template_table(1:20,1:96) = nan;
SCREENscores = struct;
% general info
SCREENscores.enzymelabel = string(template_table(:,1));
SCREENscores.ions_excl = cell(96,1);
SCREENscores.fit_function = string(template_table(:,1));
% parameter
SCREENscores.median = template_table;
SCREENscores.mean = template_table;
SCREENscores.std = template_table;
SCREENscores.eom = template_table;
SCREENscores.quant75 = template_table;
SCREENscores.quant75_raw_lb = template_table;
SCREENscores.quant75_raw_ub = template_table;
SCREENscores.quant75_DELTA_addcontrols = template_table;
SCREENscores.quant95 = template_table;
SCREENscores.quant95_raw_lb = template_table;
SCREENscores.quant95_raw_ub = template_table;
SCREENscores.quant95_raw = template_table;
SCREENscores.std_nearest = template_table;
SCREENscores.anyflag = template_table;
SCREENscores.flaglabellist = string(template_table);
% DELTA
SCREENscores.std_DELTA = template_table;
SCREENscores.quant75_DELTA = template_table;
SCREENscores.quant95_DELTA = template_table;
% confidence scores
SCREENscores.welchP_translACT_lumped_n4 = template_table;
SCREENscores.welchP_translACT_lumped_n200 = template_table;
SCREENscores.welchP_translACT_lumped_n4x8 = template_table;
SCREENscores.welchP_translACT_each_n4 = template_table;
SCREENscores.welchP_translACT_each_n200 = template_table;
SCREENscores.welchP_translACT_each_n4x8 = template_table;
SCREENscores.welchP_AUC_each_n4 = template_table;
SCREENscores.welchP_AUC_each_n200 = template_table;
SCREENscores.welchP_AUC_each_n4x8 = template_table;

% set parameters
bootstrap_err = "std";
transl_mode = "lumped";
lnACT_on = 2;
%
bypass_AUC = 0;
AUC2ac=0;
plot_r2_thresh = 0.7;
r2_threshold = plot_r2_thresh;
R2_cutoff = 0.7; % cutoff based on goodness of fit
split_on=1
excl_if_r2_lower_than_08_of_max = 0;
%
plot_only_relevant = 1; % show only relevant barplot
save_all_results = 0; % saving barplots in =1=, not here
% get info from screen overview
file = {analysis_path+"CODETABLES.xlsx"}; %path to overview
[overview_info_A,overview_info_B] = xlsread(file{1},"template","A1 : CS21"); % import screen overview table with datafile, substrate and product names
E_list = string(overview_info_B(2:end,1));
E_num = [1 : length(E_list)];

[ionselect_info_A,ionselect_info_B] = xlsread(file{1},"ions_excl_table","A1 : E20"); %import screen overview table with datafile, substrate and product names
% penalize 1 REACTANT ASSAYS
penalize1reactassay = 0;

if output_table == "new"
    % go through each enzyme and run barplot analysis script, save results
    % in one big table
    for i = 1 : 19
        clear num & enzymelabel & fiadata_effs & fiadata_scores & ions_excl & fit_function
        num = E_num(i);
        enzymelabel = E_list(i);
        ions_excl = ionselect_info_A(find(strcmp(string(ionselect_info_B(2:end,1)),enzymelabel)),:);
        if enzymelabel=="Gnd"
            enzymelabel="Gnd_reeval"
        elseif enzymelabel=="Ppc"
            enzymelabel="Ppc-repeat"
        elseif enzymelabel=="Acs"
            enzymelabel="Acs-repeat"
        elseif enzymelabel=="AckA"
            enzymelabel="AckA-repeat"
        elseif enzymelabel=="Icd"
            enzymelabel="Icd-reeval2"
        elseif enzymelabel=="MaeB"
            enzymelabel="MaeB-repipetrepeat"
        end

        load(analysis_path + "=1_output=STOREDonZENODO\" + enzymelabel + "_TABLEeffectors.mat");

        % ANALYSIS
        fit_function = fiadata_effs.fit_function;
        if sum(ions_excl == fiadata_effs.ions_excl) ~= 4
            enzymelabel
        end

        fiadata_scores = FUNCTION_ProgressCurveScreen_barplot_translated_bootstrap(analysis_path,fiadata_effs, enzymelabel, ions_excl,r2_threshold,AUC2ac,fit_function, split_on,excl_if_r2_lower_than_08_of_max,transl_mode,lnACT_on,num,penalize1reactassay,plot_only_relevant,save_all_results)
        % save into table
        SCREENscores.enzymelabel(num,1) = enzymelabel;
        SCREENscores.ions_excl(num,1) = {ions_excl};
        SCREENscores.fit_function(num,1) = fit_function;
        SCREENscores.median(num,:) = cell2mat(fiadata_scores.median(:,1)');
        SCREENscores.mean(num,:) = cell2mat(fiadata_scores.mean(:,1)');
        SCREENscores.std(num,:) = cell2mat(fiadata_scores.std(:,1)');
        SCREENscores.eom(num,:) = cell2mat(fiadata_scores.eom(:,1)');
        SCREENscores.quant75(num,:) = cell2mat(fiadata_scores.quant75(:,1)');

        SCREENscores.quant75_raw_lb(num,:) = cell2mat(fiadata_scores.quant75_raw(:,1)');
        SCREENscores.quant75_raw_ub(num,:) = cell2mat(fiadata_scores.quant75_raw(:,2)');

        SCREENscores.quant75_DELTA_addcontrols(num,:) = cell2mat(fiadata_scores.quant75_DELTA_addcontrols(:,1)');
        SCREENscores.quant95(num,:) = cell2mat(fiadata_scores.quant95(:,1)');

        SCREENscores.quant95_raw_lb(num,:) = cell2mat(fiadata_scores.quant95_raw(:,1)');
        SCREENscores.quant95_raw_ub(num,:) = cell2mat(fiadata_scores.quant95_raw(:,2)');

        SCREENscores.quant95_DELTA_addcontrols(num,:) = cell2mat(fiadata_scores.quant95_DELTA_addcontrols(:,1)');
        SCREENscores.std_nearest(num,:) = cell2mat(fiadata_scores.std_nearest(:,1)');
        SCREENscores.anyflag(num,:) = fiadata_scores.anyflag(1:96)';
        SCREENscores.flaglabellist(num,:) = fiadata_scores.flaglabellist(1:96);
        SCREENscores.std_DELTA(num,:) = cell2mat(fiadata_scores.std_DELTA(:,1)');
        SCREENscores.quant75_DELTA(num,:) = cell2mat(fiadata_scores.quant75_DELTA(:,1)');
        SCREENscores.quant95_DELTA(num,:) = cell2mat(fiadata_scores.quant95_DELTA(:,1)');
        SCREENscores.welchP_translACT_lumped_n4(num,:) = cell2mat(fiadata_scores.welch_translACT_lumped_n4(1:96,1)');
        SCREENscores.welchP_translACT_lumped_n200(num,:) = cell2mat(fiadata_scores.welch_translACT_lumped_n200(1:96,1)');
        SCREENscores.welchP_translACT_lumped_n4x8(num,:) = cell2mat(fiadata_scores.welch_translACT_lumped_n8x4(1:96,1)');
    end
elseif output_table == "existing"
    "loading existing datatable"
    load MSscreen_all_scores.mat
end

%% MS screen versus literature
SCREEN_excl = SCREENscores.anyflag(1:19,:); % labels in SCREENscores.flaglabellist
file = {analysis_path+"CODETABLES.xlsx"}; %path to screen overview
[TP_ecocyc_a,TP_ecocyc_b,TP_ecocyc_c] = xlsread(file{1},"TRUEPOSITIVES_ecocyc","B2 : CS20"); %import screen overview table with datafile, substrate and product names
[TP_brenda_a,TP_brenda_b,TP_brenda_c] = xlsread(file{1},"TRUEPOSITIVES_BRENDA","B2 : CS20");
[TP_altR_a,TP_altR_b,TP_altR_c] = xlsread(file{1},"TRUEPOSITIVES_altsubs","B2 : CS20");
[TP_smrn_a,TP_smrn_b,TP_smrn_c] = xlsread(file{1},"TRUEPOSITIVES_smrn","B2 : CS20");
% organize alternative reactants
TP_altR_a(isnan(TP_altR_a)) = 0;
TP_altR = abs(TP_altR_a);
% combine
TP_combined = abs(TP_ecocyc_a);
TP_combined(abs(TP_brenda_a)==1) = 1;
TP_combined(isnan(TP_combined)) = 0;
sum(sum(TP_combined))
% only overlap
TP_overlap = abs(TP_ecocyc_a + TP_brenda_a);
TP_overlap(isnan(TP_overlap)) = 0;
TP_overlap(TP_overlap==2) = 1;
sum(sum(TP_overlap))
SCORE_used = abs(SCREENscores.quant75); % USING FIRST QUARTILE AS SCORE
SCORE_used = SCORE_used(1:19,:);
SCORE_used(SCREEN_excl==1) = nan; % exclude all flagged ones
SCORE_used(TP_altR==1) = nan; % exclude known alternative reactants
SCORE_used(:,[94 95]) = nan;

% REremove glyoxylate from screen (is known to interact with other metabolites a lot, so we cant really derive any meaningful claims from our data)
SCORE_used(:,27) = nan;

%% ROC curve
figure(10)
hold on
xlabel('False positive rate')
ylabel('True positive rate')
title('ROCs')
plot([0 1],[0 1],'k','LineWidth',2)
hold on
%% ROC curve for Ecocyc database
TP_db = abs(TP_ecocyc_a); % SET DATABASE
TP_db(isnan(SCORE_used)) = nan;
TP_db(isnan(TP_db)) = 0;
TP_db(isnan(SCORE_used)) = nan;
nansum(nansum(TP_db)); % get number of possible literature finds
%
[X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(TP_db(:),SCORE_used(:),1);
figure(10)
plot(X,Y,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
"Ecocyc AUC: " + string(AUC)

%% ROC curve: add Brenda database
TP_db = abs(TP_brenda_a); % SET DATABASE
TP_db(isnan(SCORE_used)) = nan;
TP_db(isnan(TP_db)) = 0;
TP_db(isnan(SCORE_used)) = nan;
nansum(nansum(TP_db)); % get number of possible literature finds
%
[X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(TP_db(:),SCORE_used(:),1);
figure(10)
plot(X,Y,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
"BRENDA AUC: " + string(AUC)

%% ROC curve: add Overlap of Ecocyc and Brenda
TP_db = abs(TP_overlap); % SET DATABASE
TP_db(isnan(SCORE_used)) = nan;
TP_db(isnan(TP_db)) = 0;
TP_db(isnan(SCORE_used)) = nan;
nansum(nansum(TP_db)) % get number of possible literature finds
%
[X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(TP_db(:),SCORE_used(:),1);
figure(10)
plot(X,Y,'Color',[0.4940 0.1840 0.5560],'LineWidth',2)
"OVERLAP AUC: " + string(AUC)

label = append("ROC_curve_MS_versus_databases.png")
saveas(gca,label)





%% Validation experiments on Dehydrogenases: Photometer versus Mass Spectrometer
% import flagged TECAN assays (e.g. S or P, precipitation and more)
[~,TECAN_TP_flags,~] = xlsread(file{1},"photometer_results","B2 : CS5"); % results from =2=Photometer_analysis_script collected in excel sheet

tf = cellfun('isempty',TECAN_TP_flags); % true for empty cells
TECAN_TP_flags(tf) = {0} % all empty entries to zero
TECAN_TP_flags = string(TECAN_TP_flags);

TECAN_TP_flags(TECAN_TP_flags=="spac") = '';
TECAN_TP_flags(TECAN_TP_flags=="equ") = '';
TECAN_TP_flags(TECAN_TP_flags=="equ-strong") = '';
TECAN_TP_flags(TECAN_TP_flags=="equ-weak") = '';
TECAN_TP_flags(TECAN_TP_flags=="r2") = '';
TECAN_TP_flags(TECAN_TP_flags=="precip") = '';
TECAN_TP_flags(TECAN_TP_flags=="salt") = '';
TECAN_TP_flags(TECAN_TP_flags=="oor") = '';
TECAN_TP_flags(TECAN_TP_flags=="S") = '';
TECAN_TP_flags(TECAN_TP_flags=="P") = '';
TECAN_TP_flags(TECAN_TP_flags=="ctrl") = '';
TECAN_TP_flags(TECAN_TP_flags=="???") = '';
TECAN_TP_flags(TECAN_TP_flags=="altR") = '';
TECAN_TP_flags(TECAN_TP_flags=="null") = 0;

% import TECAN scores
clear TECAN_log2beta_mean & TECAN_log2beta_std4quant
TECAN_log2beta_mean(1:4,1:96) = nan;
TECAN_log2beta_std4quant(1:4,1:96) = nan;
[TECAN_log2beta_mean(1,1:96)] = xlsread(file{1},"photometer_results","B16:CS16");
[TECAN_log2beta_mean(2,1:96)] = xlsread(file{1},"photometer_results","B22:CS22");
[TECAN_log2beta_mean(3,1:96)] = xlsread(file{1},"photometer_results","B28:CS28");
[TECAN_log2beta_mean(4,1:96)] = xlsread(file{1},"photometer_results","B34:CS34");
[TECAN_log2beta_std4quant(1,1:96)] = xlsread(file{1},"photometer_results","B17:CS17");
[TECAN_log2beta_std4quant(2,1:96)] = xlsread(file{1},"photometer_results","B23:CS23");
[TECAN_log2beta_std4quant(3,1:96)] = xlsread(file{1},"photometer_results","B29:CS29");
[TECAN_log2beta_std4quant(4,1:96)] = xlsread(file{1},"photometer_results","B35:CS35");

tecandelta=0
flagcontrols=0
%% MS scores
screenscores_all = SCREENscores.quant75(1:19,:);
if tecandelta == 1
    MSscore_mean = SCREENscores.quant75(1:19,:);
else
    MSscore_mean = SCREENscores.mean(1:19,:);
end
MSscore_mean(SCREEN_excl==1) = nan;
MSscore_mean(TP_altR==1) = nan;
MSscore_mean(:,[94 95]) = nan;
MSscore_mean(:,27) = nan; % remove glyoxylate
% % % % % % % % %
MSscore_mean(:,[78:81]) = nan; %% exclude salt controls
%
MSscore_std = abs(SCREENscores.std(1:19,:));
MSscore_std(isnan(MSscore_mean)) = nan;


MSscore_quant50lb = abs(SCREENscores.quant75_raw_lb(1:19,:) - MSscore_mean);
MSscore_quant50lb(isnan(MSscore_mean)) = nan;
MSscore_quant50ub = abs(SCREENscores.quant75_raw_ub(1:19,:) - MSscore_mean);
MSscore_quant50ub(isnan(MSscore_mean)) = nan;
screenscores_all(isnan(MSscore_mean)) = nan;

% only DHs
clear TECANvsMS_mean & TECANvsMS_std & TECANvsMS_quant50raw_lb & TECANvsMS_quant50raw_ub & TECANvsMS_quant50score
TECANvsMS_mean(1,1:size(SCORE_used,2)) = MSscore_mean(1,:);
TECANvsMS_mean(2,1:size(SCORE_used,2)) = MSscore_mean(18,:);
TECANvsMS_mean(3,1:size(SCORE_used,2)) = MSscore_mean(2,:);
TECANvsMS_mean(4,1:size(SCORE_used,2)) = MSscore_mean(19,:);
TECANvsMS_std(1,1:size(SCORE_used,2)) = MSscore_std(1,:);
TECANvsMS_std(2,1:size(SCORE_used,2)) = MSscore_std(18,:);
TECANvsMS_std(3,1:size(SCORE_used,2)) = MSscore_std(2,:);
TECANvsMS_std(4,1:size(SCORE_used,2)) = MSscore_std(19,:);
% remove those that are invalid on TECAN
TECANvsMS_mean(TECAN_TP_flags=='') = nan;
TECANvsMS_std(TECAN_TP_flags=='') = nan;
%
TECANvsMS_quant50raw_lb(1,1:size(SCORE_used,2)) = MSscore_quant50lb(1,:);
TECANvsMS_quant50raw_lb(2,1:size(SCORE_used,2)) = MSscore_quant50lb(18,:);
TECANvsMS_quant50raw_lb(3,1:size(SCORE_used,2)) = MSscore_quant50lb(2,:);
TECANvsMS_quant50raw_lb(4,1:size(SCORE_used,2)) = MSscore_quant50lb(19,:);
TECANvsMS_quant50raw_ub(1,1:size(SCORE_used,2)) = MSscore_quant50ub(1,:);
TECANvsMS_quant50raw_ub(2,1:size(SCORE_used,2)) = MSscore_quant50ub(18,:);
TECANvsMS_quant50raw_ub(3,1:size(SCORE_used,2)) = MSscore_quant50ub(2,:);
TECANvsMS_quant50raw_ub(4,1:size(SCORE_used,2)) = MSscore_quant50ub(19,:);
%
TECANvsMS_quant50score(1,1:size(SCORE_used,2)) = screenscores_all(1,:);
TECANvsMS_quant50score(2,1:size(SCORE_used,2)) = screenscores_all(18,:);
TECANvsMS_quant50score(3,1:size(SCORE_used,2)) = screenscores_all(2,:);
TECANvsMS_quant50score(4,1:size(SCORE_used,2)) = screenscores_all(19,:);


% use nans in MS data on TECAN data
TECAN_log2beta_mean(isnan(TECANvsMS_mean)) = nan;
TECAN_log2beta_std4quant(isnan(TECANvsMS_mean)) = nan;
TECAN_log2beta_mean = reshape(TECAN_log2beta_mean',1,[]);
TECAN_log2beta_std4quant = reshape(TECAN_log2beta_std4quant',1,[]);

% calculate quantile from std for tecandata
TECAN_log2beta_quant = 0.675 * TECAN_log2beta_std4quant;

%
TECANvsMS_mean = reshape(TECANvsMS_mean',1,[]);
TECANvsMS_std = reshape(TECANvsMS_std',1,[]);
TECANvsMS_quant50raw_lb = reshape(TECANvsMS_quant50raw_lb',1,[]);
TECANvsMS_quant50raw_ub = reshape(TECANvsMS_quant50raw_ub',1,[]);
TECANvsMS_quant50score = reshape(TECANvsMS_quant50score',1,[]);

TP_db_new = [];
TP_db_new(1,1:size(SCORE_used,2)) = TP_db(1,:);
TP_db_new(2,1:size(SCORE_used,2)) = TP_db(18,:);
TP_db_new(3,1:size(SCORE_used,2)) = TP_db(2,:);
TP_db_new(4,1:size(SCORE_used,2)) = TP_db(19,:);
TP_db_new = reshape(TP_db_new',1,[])

TECAN_log2beta_mean(isnan(TECANvsMS_mean)) = nan;
TECANvsMS_mean(isnan(TECAN_log2beta_mean)) = nan;
TP_db_new(isnan(TECAN_log2beta_mean)) = nan;
TECANvsMS_quant50score(isnan(TECAN_log2beta_mean)) = nan;

log2transf=1;
if log2transf == 0
    TECAN_log2beta_mean=2.^TECAN_log2beta_mean
    TECANvsMS_mean = 2.^TECANvsMS_mean
end

figure
hold on
if tecandelta == 1
    plot(TECAN_log2beta_mean(:),TECANvsMS_mean(:),...
        'ok','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','white')
    hold on
    plot(TECAN_log2beta_mean(abs(TECANvsMS_mean)>1.1313),TECANvsMS_mean(abs(TECANvsMS_mean)>1.1313),...
        'ok','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor',	[0.3010, 0.7450, 0.9330])
    if flagcontrols == 1
        plot(TECAN_log2beta_mean(TP_db_new==1),TECANvsMS_mean(TP_db_new==1),...
            'ok','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',	'red')
    end
else
    errorbar(TECAN_log2beta_mean(:),TECANvsMS_mean(:),TECANvsMS_quant50raw_lb(:),TECANvsMS_quant50raw_ub(:),TECAN_log2beta_quant(:),TECAN_log2beta_quant(:),...
        'ok','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',	[1 1 1])

    hold on
    plot(TECAN_log2beta_mean(abs(TECANvsMS_quant50score)>1.1313),TECANvsMS_mean(abs(TECANvsMS_quant50score)>1.1313),...
        'ok','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor',	[0.3010, 0.7450, 0.9330])



    % Controls
    flagcontrols = 1
    if flagcontrols == 1

        plot(TECAN_log2beta_mean(TP_db_new==1),TECANvsMS_mean(TP_db_new==1),...
            'ok','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',	'red')


    end
end



hold on
plot([-9 9],[-9 9],'k--')
plot([-1.1313 1.1313],[-1.1313 -1.1313],'k-')
plot([-1.1313 1.1313],[1.1313 1.1313],'k-')
plot([-1.1313 -1.1313],[-1.1313 1.1313],'k-')
plot([1.1313 1.1313],[-1.1313 1.1313],'k-')
plot([-.5 .5],[-.5 -.5],'k--')
plot([-.5 .5],[.5 .5],'k--')
plot([-.5 -.5],[-.5 .5],'k--')
plot([.5 .5],[-.5 .5],'k--')


if tecandelta == 1
    xlabel("PHOTOMETER std distance")
    ylabel("MASS-SPECTROMETER quant50")
else
    xlabel("PHOTOMETER log2(activity)")
    ylabel("MASS-SPECTROMETER log2(activity)")
end

ylim([-5.2 2.5])
xlim([-5.2 2.5])
label = append("Scores_MS_versus_Photometer.png")
saveas(gca,label)





%% PLOT ALL SCORES / MEAN+QUANT OF MS screen (and mark literature (and more))

%% REMOVE ALL THAT ARE FLAGGED
MSscore_quant50ub(isnan(MSscore_mean))=nan;
MSscore_quant50lb(isnan(MSscore_mean))=nan;

% rearrange
MSscore_mean_t = reshape(MSscore_mean',1,[]);
MSscore_quant50ub_t = reshape(MSscore_quant50ub',1,[]);
MSscore_quant50lb_t  = reshape(MSscore_quant50lb',1,[]);
screenscores_all_t = reshape(screenscores_all',1,[]);
TP_ecocyc_a_t = reshape(TP_ecocyc_a',1,[]);
TP_brenda_a_t = reshape(TP_brenda_a',1,[]);
f = figure
errorbar(1:length(MSscore_mean_t(:)),MSscore_mean_t(:),MSscore_quant50lb_t(:),MSscore_quant50ub_t(:),...
    'ok','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','white')
hold on
plot([0 length(MSscore_mean_t(:))+1],[1.1313 1.1313],'k');
plot([0 length(MSscore_mean_t(:))+1],[-1.1313 -1.1313],'k');
plot([0 length(MSscore_mean_t(:))+1],[.5 .5],'k--');
plot([0 length(MSscore_mean_t(:))+1],[-.5 -.5],'k--');

% mark hits blue
markhits_x = 1:length(MSscore_mean_t(:));
plot(markhits_x(abs(screenscores_all_t)>1.1313),MSscore_mean_t(abs(screenscores_all_t)>1.1313),...
    'ok','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor',	[0.3010, 0.7450, 0.9330])

% mark ecocyc lit red
TP_db = abs(TP_ecocyc_a_t); % SET DATABASE
TP_db(isnan(MSscore_mean_t)) = nan;
TP_db(isnan(TP_db)) = 0;
TP_db(isnan(MSscore_mean_t)) = nan;

plot(markhits_x(TP_db==1),MSscore_mean_t(TP_db==1),...
    'ok','MarkerSize',10,'MarkerEdgeColor','red')
hold on

% MARK BRENDA
markbrenda = 1;
if markbrenda == 1
    TP_db = abs(TP_brenda_a_t); % SET DATABASE
    TP_db(isnan(MSscore_mean_t)) = nan;
    TP_db(isnan(TP_db)) = 0;
    TP_db(isnan(MSscore_mean_t)) = nan;
    plot(markhits_x(TP_db==1),MSscore_mean_t(TP_db==1),...
        'ok','MarkerSize',7,'MarkerEdgeColor','green')
    hold on
end
spacing96 = [96:96:96*length(E_list)];
hold on
for i = 1 : length(spacing96)
    plot([spacing96(i) spacing96(i)],[-10 10],'k--')
    if E_list(i) == "Gnd_reeval"
        E_str_split = "Gnd"
    else
        E_str = E_list(i);
        E_str_split = strsplit(E_str,"-");
    end
    text(spacing96(i)-70,-4,E_str_split(1))
end
ylim([-5 4])

f.Position = [100 300 2000 400];

% save figures
label = append("all_scores_MEAN_plus_QUANTILE.png")
saveas(gca,label)


% SCORES
figure
scoreplot_logon = 1;
if scoreplot_logon == 1
    screenscores_all_t_plot = screenscores_all_t;
else
    screenscores_all_t_plot = 2.^screenscores_all_t;
end

[out,idx] = sort(screenscores_all_t_plot) % sort scores and get indices
hold on
plot(1:size(out,2),out,'ok', 'MarkerSize',0.5)
plot(find(out>1.131),out(out>1.131),'ok', 'MarkerFaceColor',[0 0 0],'MarkerSize',1.5)
plot(find(out<-1.131),out(out<-1.131),'ok', 'MarkerFaceColor',[0 0 0],'MarkerSize',1.5)
sum(abs(screenscores_all_t_plot)>1.131)

% ecocyc
TP_db = abs(TP_ecocyc_a_t); % SET DATABASE
TP_db(isnan(screenscores_all_t_plot)) = nan;
TP_db(isnan(screenscores_all_t_plot)) = nan;
TP_db(isnan(TP_db)) = 0;
plot(find(TP_db(idx)==1),out(find(TP_db(idx)==1)),'o', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

% brenda
TP_db = abs(TP_brenda_a_t); % SET DATABASE
TP_db(isnan(screenscores_all_t_plot)) = nan;
TP_db(isnan(screenscores_all_t_plot)) = nan;
TP_db(isnan(TP_db)) = 0;
plot(find(TP_db(idx)==1),out(find(TP_db(idx)==1)),'o', 'MarkerSize',5,'Color',[0.6350 0.0780 0.1840],'LineWidth',1)
% count
"total possible BRENDA: " + string(sum(TP_db))
"recovered from BRENDA: " + string()


plot([1 size(out,2)],[1.1313 1.1313],'k-')
plot([1 size(out,2)],[-1.1313 -1.1313],'k-')
plot([1 size(out,2)],[1.1313/2 1.1313/2],'k--')
plot([1 size(out,2)],[-1.1313/2 -1.1313/2],'k--')

title("ranking of all scores (marked literature red/blue)")
ylabel("score: log2-fold change in kcat/kM")
xlim([0 1500])
label = append("ranking_and_marking_literature.png")
saveas(gca,label)











%% PIE charts, screen versus literature

% FOR SCORES
% add ecocyc TPs
TP_db = abs(TP_ecocyc_a_t);
% add brenda TPs
TP_db(abs(TP_brenda_a_t)==1) = 1;

TP_db(isnan(MSscore_mean_t)) = nan;
TP_db(isnan(TP_db)) = 0;
TP_db(isnan(MSscore_mean_t)) = nan;
% quantile-scores
lit_scores = screenscores_all_t(TP_db==1);
lit_total_n = length(lit_scores);
lit_recovered_n = sum(abs(lit_scores)>1.1313); % 19
lit_recovered_below05 = sum(abs(lit_scores)<1.1313/2); % 21
lit_between05andhit_n = lit_total_n-lit_recovered_n-lit_recovered_below05; % 7

figure
labels = ["recovered: "+string(lit_recovered_n)+" - "+string(round(lit_recovered_n/lit_total_n*100,0))+"%", ...
    "score below 0.5656: "+string(lit_recovered_below05)+" - "+string(round(lit_recovered_below05/lit_total_n*100,0))+"%", ...
    "score above 0.5656: "+string(lit_between05andhit_n)+" - "+string(round(lit_between05andhit_n/lit_total_n*100,0))+"%"];
pie([lit_recovered_n lit_recovered_below05 lit_between05andhit_n],labels)
title("pie chart SCORES")
subtitle(string(round(lit_recovered_below05/(lit_recovered_below05+lit_between05andhit_n)*100,0))+"% of missed hits below 0.5656")
label = append("how_many_missed_PIE_scores.png")
saveas(gca,label)

% pie for mean
litM_scores = MSscore_mean_t(TP_db==1);
litM_total_n = length(litM_scores);
litM_recovered_n = sum(abs(litM_scores)>1.1313); %
litM_recovered_below05 = sum(abs(litM_scores)<1.1313/2); %
litM_between05andhit_n = litM_total_n-litM_recovered_n-litM_recovered_below05; %

figure
labels = ["recovered: "+string(litM_recovered_n)+" - "+string(round(litM_recovered_n/litM_total_n*100,0))+"%", ...
    "score below 0.5656: "+string(litM_recovered_below05)+" - "+string(round(litM_recovered_below05/litM_total_n*100,0))+"%", ...
    "score above 0.5656: "+string(litM_between05andhit_n)+" - "+string(round(litM_between05andhit_n/litM_total_n*100,0))+"%"];
pie([litM_recovered_n litM_recovered_below05 litM_between05andhit_n],labels)
title("pie chart MEANS")
subtitle(string(round(litM_recovered_below05/(litM_recovered_below05+litM_between05andhit_n)*100,0))+"% of missed hits below 0.5656")
label = append("how_many_missed_PIE_means.png")
saveas(gca,label)