%% PROGRESS CURVE ANALYSIS from mass spectrometry data
% Input for this script is processed mass spectrometry data, i.e. full data
% data table after peak picking, annotation and quantification.
% It does processing (e.g. normalization and exclusion of faulty samples)
% and fits and plots all assays to extract information on change in enzyme
% activity when metabolite is added.
close all
clear
clc

%% Inputs
analysis_path = "\\pasteur\SysBC-Home\chgruber\Desktop\==MANUSCRIPT WIP\==submission\=2=CODE\"; % path with saved folders
enzyme_name = ["PykF"]; % enzyme name corresponding to CODETABLES.xlsx - overview sheet
inspect_controls = 1; % plot progress curves of controls
plot_on = 1; % plot progress curves of effectors
save_all_results = 1; % save all generated results (might overwrite previous ones)

%% Parameters
TIC_normalization = 1; % normalize ion counts to total ion count
ioncount_threshold = 3000; % exclude ions with count below threshold
errorbars_def = "std"; % standard deviation ("std") or error of mean ("eom")
errorbars_ioncount = "std";
bootstrap_err = "std";
TICfilter = "CUT50"; % "CUT50": exclude if TIC of sample 50% off from median TIC of timeseries
robustfit_on = 0;
plot_r2_thresh = 0.7; % exclude fit if Rsquared below 0.7
r2_threshold = plot_r2_thresh;
R2_cutoff = 0.7; % cutoff for CONTROLS based on goodness of fit
excl_if_r2_lower_than_08_of_max = 0;
split_bio_replicates = 0; % 0 = lumping replicates
ctrl_split_on = 1; % 
remove_ctrls_ratio_median = 0.7;
ctrl_collection_REscaled = 1;
plot_sep_scaled = "REscaled"; % plot progress curve scaled between 0 and 1
plot_OVERLAY_scaled = "REscaled"; % plot progress curve scaled between 0 and 1
get_transl_table = "overlay1to4and5to8"; % for each enzyme, data was generated in 2 MS measurements (with cleaning and recalibration in between), therefore analyzing first half and second half of tested effectors and controls separately
plot_only_95conf = 1; % plot 95% intervals of bootstrap instead of ALL
plot_only_relevant = 1;

%% Get information from overview
file = {analysis_path+"CODETABLES.xlsx"}; %path to screen overview
[overview_a,overview_b] = xlsread(file{1},"overview","A1 : AV20");
overview_labels = overview_b(1,:);
overview_b(1,:) = [];
Eidx = find(strcmp(overview_b(:,find(strcmp(overview_labels,"enzyme"))),enzyme_name));
enzyme = overview_b(Eidx,find(strcmp(overview_labels,"enzyme")));
enzymelabel = overview_b(Eidx,find(strcmp(overview_labels,"label")));
screen_idx = overview_a(Eidx,find(strcmp(overview_labels,"idx")));
setup = overview_a(Eidx,find(strcmp(overview_labels,"setup")));
fit_function = overview_b(Eidx,find(strcmp(overview_labels,"fitfunction")));

%% Load preprocessed mass spectrometry data
load(analysis_path+"=1_input=MSdata_preprocessed\"+screen_idx+"-"+enzymelabel+"_processed.mat")

%% Manual changes
% mainly excluding faulty injections based on visual inspection
if enzymelabel == "Edd"
    fiadata.S1(:,14) = [];
    fiadata.P1(:,14) = [];
    fiadata.S1_TICs(:,14) = [];
    fiadata.P1_TICs(:,14) = [];
elseif enzymelabel == "MaeA"
    fiadata.S1(318:319,7:end) = {nan};
    fiadata.S2(318:319,7:end) = {nan};
    fiadata.P1(318:319,7:end) = {nan};
    fiadata.P2(318:319,7:end) = {nan};
 elseif enzymelabel == ["MaeB-repipetrepeat"]
    fiadata.S1(192:193,7:end) = {nan};
    fiadata.S1(168:169,7:end) = {nan};
    fiadata.S2(192:193,7:end) = {nan};
    fiadata.S2(168:169,7:end) = {nan};    
    fiadata.P1(192:193,7:end) = {nan};
    fiadata.P1(168:169,7:end) = {nan};    
    fiadata.P2(192:193,7:end) = {nan};
    fiadata.P2(168:169,7:end) = {nan};    
    ioncount_threshold = 1000;
elseif enzymelabel == "Fbp"
    TICfilter = "CUT20"; % more stringent cutoff based on variance in TIC in cases where TIC was strongly influenced by reactant(s) of the reaction
elseif enzymelabel == "AceA"
    TICfilter = "CUT20"; % more stringent cutoff based on variance in TIC in cases where TIC was strongly influenced by reactant(s) of the reaction
elseif enzymelabel == "Icd-reeval2"
    fiadata.S1(48,7:end) = {nan}; 
    TICfilter = "CUT20"; % S determines TIC, therefore huge deviations especially when no effector is there
elseif enzymelabel == "Ppc-repeat" % product(s) start at very low counts, therefore no cutoff
    ioncount_threshold = 0
end

%% Fit (and plot) controls
fiadata_ctrls = FUNCTION_ProgressCurveScreen_fit_MSdata_controls(fiadata, ioncount_threshold, inspect_controls, TICfilter, TIC_normalization, R2_cutoff,fit_function,save_all_results);

%% Fit (and plot) effectors
[fiadata_effs,fiadata_SIMS] = FUNCTION_ProgressCurveScreen_fit_MSdata_effectors_bootstrap(enzymelabel, fiadata, fiadata_ctrls, ioncount_threshold, fit_function, robustfit_on, plot_on, plot_r2_thresh,TICfilter,TIC_normalization,plot_sep_scaled,plot_OVERLAY_scaled,errorbars_ioncount,ctrl_collection_REscaled,split_bio_replicates,ctrl_split_on,bootstrap_err,remove_ctrls_ratio_median,r2_threshold,get_transl_table,plot_only_95conf,analysis_path,plot_only_relevant,save_all_results)

%% Create overview to find most informative ion traces
[ions_excl_parameters,ions_excl] = FUNCTION_ProgressCurveScreen_ions_excl_parameters_overview(overview_a,overview_labels,Eidx,fiadata_ctrls,fiadata_effs);
% more information on chosen ion traces in CODETABLES.xlsx in the ions_excl_overview sheet

%% Step 7: Analysis
fiadata_effs.ions_excl_parameters = ions_excl_parameters;
fiadata_effs.ions_excl = ions_excl;
split_on=1;
bypass_AUC = 0
AUC2ac=0
enzymelabel = fiadata_effs.enzymelabel;
plot_r2_thresh = 0.7;
r2_threshold = plot_r2_thresh;
R2_cutoff = 0.7; % cutoff for CONTROLS
split_on=1
fit_function = fiadata_effs.fit_function;
lnACT_on = 2;
transl_mode = "lumped"; %"lumped" controls or "sep"-arate controls
num=0;
penalize1reactassay=0;
fiadata_scores = FUNCTION_ProgressCurveScreen_barplot_translated_bootstrap(analysis_path,fiadata_effs, enzymelabel, ions_excl,r2_threshold,AUC2ac,fit_function, split_on,excl_if_r2_lower_than_08_of_max,transl_mode,lnACT_on,num,penalize1reactassay,plot_only_relevant,save_all_results)

%% Save tables with results+
if save_all_results == 1
 save(analysis_path + enzymelabel + "_TABLEeffectors.mat", "fiadata_effs");
 save(analysis_path + enzymelabel + "_TABLEsimulations.mat", "fiadata_SIMS");
 save(analysis_path + enzymelabel + "_TABLEcontrols.mat", "fiadata_ctrls");
 save(analysis_path + enzymelabel + "_TABLEionexclusion.mat", "ions_excl_parameters");
end