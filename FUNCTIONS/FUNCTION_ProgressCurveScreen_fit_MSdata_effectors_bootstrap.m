function [fiadata_effs,fiadata_SIMS] = ProgressCurveScreen2_fit_fiadata_effectors_v14_bootstrap_resample(enzymelabel, fiadata, fiadata_ctrls, ioncount_threshold, fit_function, robustfit_on, plot_on, plot_r2_thresh,TICfilter,TIC_normalization,plot_sep_scaled,plot_OVERLAY_scaled,errorbars_ioncount,ctrl_collection_REscaled,split_bio_replicates,ctrl_split_on,bootstrap_err,remove_ctrls_ratio_median,r2_threshold,get_transl_table,plot_only_95conf,analysis_path,plot_only_relevant,save_all_results)
log_x_plot = 0;
extra_figures = 0;
fit_impossible = 0;
t0_cv_cutoff = 0.3; % maxCV of technical replicates of t0, otherwise excluded
exclude_all_if_t0_compromised = 0; % exclude whole bio replicate if t0 is excluded?
nboots = 200; % number of bootstraps

AUC2activity_activity = [0.01 : 0.001 : 10];

% get information on effectors from overview file and save in struct
file = {analysis_path+"CODETABLES.xlsx"}; %path to screen overview
[effsinfo_a,effsinfo_b,effsinfo_c] = xlsread(file{1},"effectors","A1 : D100"); %import screen overview table with datafile, substrate and product names
effsinfo_labels = effsinfo_b(1,:);
effsinfo_b(1,:) = [];
%
[effMASS_mass,effMASS_name] = xlsread(file{1},"effectors_masses","A1 : B88"); %import screen overview table with datafile, substrate and product names
% add controls for 8x12 plot labels
S = effMASS_name; % original sequence
ir = 11 ; % insert every x times
io = {'control'};
S = [S(:);nan(mod(-numel(S),ir),1)];
Sn = reshape(S,ir,[]);
Sn = [Sn;repmat(io,1,size(Sn,2))];
effMASS_name_plotting = Sn(:);
%
[overview_info_A,overview_info_B] = xlsread(file{1},"overview","A1 : BO34"); %import screen overview table with datafile, substrate and product names
Eindex = find(strcmp(string(overview_info_B(:,find(strcmp(overview_info_B(1,:),"label")))),enzymelabel));
reactant_masses = [overview_info_A(Eindex-1,find(strcmp(overview_info_B(1,:),"S1-m/z"))) overview_info_A(Eindex-1,find(strcmp(overview_info_B(1,:),"S2-m/z"))) overview_info_A(Eindex-1,find(strcmp(overview_info_B(1,:),"P1-m/z"))) overview_info_A(Eindex-1,find(strcmp(overview_info_B(1,:),"P2-m/z")))];
reactant_abbrev = string([overview_info_B(Eindex,find(strcmp(overview_info_B(1,:),"reactants-1"))) overview_info_B(Eindex,find(strcmp(overview_info_B(1,:),"reactants-2"))) overview_info_B(Eindex,find(strcmp(overview_info_B(1,:),"reactants-3"))) overview_info_B(Eindex,find(strcmp(overview_info_B(1,:),"reactants-4")))]);
precip_strong = string(overview_info_B(Eindex,find(strcmp(string(overview_info_B(1,:)),"precipitation_strong"))));
precip_weak = string(overview_info_B(Eindex,find(strcmp(string(overview_info_B(1,:)),"precipitation_weak"))));
%


fiadata_effs.enzymelabel = enzymelabel;
struct_template(1:size(effsinfo_b),5) = effsinfo_c(2:end,4);
struct_template(1:size(effsinfo_b),6) = {nan};
fiadata_effs.fit_function = fit_function;
fiadata_effs.analysis_legend = struct_template;
fiadata_effs.analysis_legend(:,1:4) = {1};

% exclude substrates and products from analysis
for r = 1 : length(reactant_abbrev)
    SandP_index = find(strcmp(string(fiadata_effs.analysis_legend(:,5)),reactant_abbrev(r)));
    if SandP_index ~= 0
        fiadata_effs.analysis_legend(SandP_index,1:4) = {0};
        if r < 3
            fiadata_effs.analysis_legend(SandP_index,6) = {"SUBSTRATE"};
        elseif r > 2
            fiadata_effs.analysis_legend(SandP_index,6) = {"PRODUCT"};
        end
    end
end

% mark precipitates
precip_strong_split = split(precip_strong,',');
for k = 1 : length(precip_strong_split)
    precip_str_index = find(strcmp(string(fiadata_effs.analysis_legend(:,5)),precip_strong_split(k)));
    if precip_str_index ~= 0
       fiadata_effs.analysis_legend(precip_str_index,1:4) = {0};
       fiadata_effs.analysis_legend(precip_str_index,6) = {"PRECIP"};
    end
end
precip_weak_split = split(precip_weak,',');
for m = 1 : length(precip_weak_split)
    precip_weak_index = find(strcmp(string(fiadata_effs.analysis_legend(:,5)),precip_weak_split(m)));
    if precip_weak_index ~= 0
       fiadata_effs.analysis_legend(precip_weak_index,1:4) = {0};
       fiadata_effs.analysis_legend(precip_str_index,6) = {"precip"};
    end
end

% fit options
clear opts
opts = statset('nlinfit');
if robustfit_on == 1
    opts.Robust = 'on';
    opts.RobustWgtFun = 'bisquare';
end

% preparing plots
dist_colors = distinguishable_colors(12)
dist_colors(4,:) = dist_colors(end,:)
dist_colors(end,:) = [0 0 0] % make effector 12 black (=> control)
for minim = 1
    % create figures
    if plot_on == 1
        figure(11)
        sgtitle(fiadata.enzymename+" set1 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1 : 4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(12)
        sgtitle(fiadata.enzymename+" set2 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(13)
        sgtitle(fiadata.enzymename+" set3 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(14)
        sgtitle(fiadata.enzymename+" set4 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(15)
        sgtitle(fiadata.enzymename+" set5 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(16)
        sgtitle(fiadata.enzymename+" set6 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(17)
        sgtitle(fiadata.enzymename+" set7 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(18)
        sgtitle(fiadata.enzymename+" set8 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(19)
        sgtitle(fiadata.enzymename+" set8 ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end
        
        figure(20)
        sgtitle(fiadata.enzymename+" CONTROLS ion counts")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("ion count")
        end



% ACTIVITY -> AUC plotting
        figure(40)
        sgtitle(fiadata.enzymename+" ACTIVITY -> AUC plotting")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
        end
        
        figure(41)
        sgtitle(fiadata.enzymename+" ACTIVITY -> AUC calibration")
        subplot(2,2,1)
        title(fiadata.S1_met)
        hold on
        subplot(2,2,2)
        title(fiadata.S2_met)
        hold on
        subplot(2,2,3)
        title(fiadata.P1_met)
        hold on
        subplot(2,2,4)
        title(fiadata.P2_met)
        hold on
        for i = 1:4
            subplot(2,2,i)
            xlabel("AUC")
            ylabel("activity (relative to ctrl)")
        end
        
        if extra_figures == 1
            figure(51)
            title(fiadata.enzymename+" set1 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(52)
            title(fiadata.enzymename+" set2 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(53)
            title(fiadata.enzymename+" set3 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(54)
            title(fiadata.enzymename+" set4 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(55)
            title(fiadata.enzymename+" set5 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(56)
            title(fiadata.enzymename+" set6 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(57)
            title(fiadata.enzymename+" set7 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(58)
            title(fiadata.enzymename+" set8 OVERLAY")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            movegui(figure(51),"northwest")
            movegui(figure(52),"north")
            movegui(figure(53),"northeast")
            movegui(figure(54),"west")
            movegui(figure(55),"center")
            movegui(figure(56),"east")
            movegui(figure(57),"southwest")
            movegui(figure(58),"south")
            
            figure(61)
            title(fiadata.enzymename+" set1 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(62)
            title(fiadata.enzymename+" set2 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(63)
            title(fiadata.enzymename+" set3 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(64)
            title(fiadata.enzymename+" set4 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(65)
            title(fiadata.enzymename+" set5 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(66)
            title(fiadata.enzymename+" set6 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(67)
            title(fiadata.enzymename+" set7 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            figure(68)
            title(fiadata.enzymename+" set8 singles")
            hold on
            xlabel("time [min]")
            ylabel("scaled ion count")
            ylim([0 1])
            
            movegui(figure(61),"northwest")
            movegui(figure(62),"north")
            movegui(figure(63),"northeast")
            movegui(figure(64),"west")
            movegui(figure(65),"center")
            movegui(figure(66),"east")
            movegui(figure(67),"southwest")
            movegui(figure(68),"south")
            
        end
    end
end




for minim = 2
    fig1 = figure(111)
    sgtitle(fiadata.enzymename+" BOOTSTRAP S1")
    fig1.WindowState = 'maximized';
    hold on
    subplot(8,12,1)
    for i = 1 : 8*12
        subplot(8,12,i)
        subtitle(string(effsinfo_b(i,4)))
        xlabel("time [min]")
        ylabel("scaled ion count")
    end
    hold on
    
    fig2 = figure(112)
    sgtitle(fiadata.enzymename+"BOOTSTRAP S2")
    fig2.WindowState = 'maximized';
    hold on
    subplot(8,12,1)
    for i = 1 : 8*12
        subplot(8,12,i)
        subtitle(string(effsinfo_b(i,4)))
        xlabel("time [min]")
        ylabel("scaled ion count")
    end
    hold on
    
    fig3 = figure(113)
    sgtitle(fiadata.enzymename+" BOOTSTRAP P1")
    fig3.WindowState = 'maximized';
    hold on
    subplot(8,12,1)
    for i = 1 : 8*12
        subplot(8,12,i)
        subtitle(string(effsinfo_b(i,4)))
        xlabel("time [min]")
        ylabel("scaled ion count")
    end
    hold on
    
    fig4 = figure(114)
    sgtitle(fiadata.enzymename+" BOOTSTRAP P2")
    fig4.WindowState = 'maximized';
    hold on
    subplot(8,12,1)
    for i = 1 : 8*12
        subplot(8,12,i)
        subtitle(string(effsinfo_b(i,4)))
        xlabel("time [min]")
        ylabel("scaled ion count")
    end
    hold on
end

% get indices of all effectors controls and breaking point
idx_effs = find(strcmp(fiadata.S1(:,6),"assay:"));
idx_break = min(find(strcmp(fiadata.S1(:,2),"set5")));
i_break = min(find(idx_effs>idx_break));

% check R2 of controls
if fit_function == "2SlambertW"
    ctrls_excl = fiadata_ctrls.TWOSUBSTRATElambertW_r2;
elseif fit_function == "lambertW"
    ctrls_excl = fiadata_ctrls.lambertW_r2;
elseif fit_function == "SIGMOID"
    ctrls_excl = fiadata_ctrls.SIGMOID_r2;
end

% create output table
fiadata_effs.Yfits = struct_template; % XVALUES x 200 fits
fiadata_effs.Yvalues = struct_template; % XVALUES x 200 fits
fiadata_effs.impossible_fits = struct_template; % 1 counter per eff and reactant
fiadata_effs.AUCs = struct_template; % 200 AUCs per eff and reactant
fiadata_effs.AUCsALL =  struct_template; % including the AUCs excluded because r2<0.7
fiadata_effs.stdAUCs = struct_template;
fiadata_effs.cvAUCs = struct_template;
fiadata_effs.meanAUCs = struct_template;
fiadata_effs.meanAUCsCURVE = struct_template;
fiadata_effs.act = struct_template; % 200 AUCs per eff and reactant
fiadata_effs.R2s = struct_template; % 200 R2s per eff and reactant
fiadata_effs.AUC2activity_all(1:size(AUC2activity_activity,2),1:5) = nan; % translation of AUC to activity change based on beta parameter simulations
fiadata_effs.AUC2activity_all(1:size(AUC2activity_activity,2),5) = AUC2activity_activity; % translation of AUC to activity change based on beta parameter simulations
fiadata_effs.AUC2activity_1st = fiadata_effs.AUC2activity_all; % translation of AUC to activity change based on beta parameter simulations
fiadata_effs.AUC2activity_2nd = fiadata_effs.AUC2activity_all; % translation of AUC to activity change based on beta parameter simulations
% v14
fiadata_effs.SIG_AUC2activity_all = fiadata_effs.AUC2activity_all; % SIGMOID-based translation of AUC to activity change based on beta parameter simulations
fiadata_effs.SIG_AUC2activity_1st = fiadata_effs.AUC2activity_all;  % SIGMOID-based translation of AUC to activity change based on beta parameter simulations
fiadata_effs.SIG_AUC2activity_2nd = fiadata_effs.AUC2activity_all;  % SIGMOID-based translation of AUC to activity change based on beta parameter simulations
fiadata_effs.MEANsepctrlAUC2activity_1st = fiadata_effs.AUC2activity_all; %
fiadata_effs.MEANsepctrlAUC2activity_2nd = fiadata_effs.AUC2activity_all; % 
fiadata_SIMS = struct;
fiadata_SIMS.S1 = nan(8,size(fiadata_effs.AUC2activity_all,1));
fiadata_SIMS.S2 = nan(8,size(fiadata_effs.AUC2activity_all,1));
fiadata_SIMS.P1 = nan(8,size(fiadata_effs.AUC2activity_all,1));
fiadata_SIMS.P2 = nan(8,size(fiadata_effs.AUC2activity_all,1));

% define fit functions
fun_lambertW_subs = @(x,t) x(3) * lambertw(0, x(1) * exp(x(1) - x(2)*t)); % a = x(1); b = x(2); k_m = x(3);
fun_lambertW_prod = @(x,t) x(3) * (x(1) - lambertw(0, x(1) * exp(x(1) - x(2)*t))); % a = x(1); b = x(2); k_m = x(3);
fun_TWOSUBSTRATElambertW_subs = @(x,t) x(3)/2 * (x(1)/x(3) - x(3)/x(1) - x(2) * t + sqrt((x(1)/x(3) - x(3)/x(1) - x(2) * t).^2 + 4)); % s0 = x(1); b = x(2); k_m = x(3);
fun_TWOSUBSTRATElambertW_prod = @(x,t) 1 - x(3)/2 * (x(1)/x(3) - x(3)/x(1) - x(2) * t + sqrt((x(1)/x(3) - x(3)/x(1) - x(2) * t).^2 + 4)); % s0 = x(1); b = x(2); k_m = x(3);
fun_SIGMOID_subs = @(x,t) 1 - x(1) * tanh(abs(x(2)/x(1))*(t-x(3))) ; % SIGMOID: x(1) height of curve; x(2) heat-parametersc; x(3) center point
fun_SIGMOID_prod = @(x,t) x(1) * tanh(abs(x(2)/x(1))*(t-x(3))) ; % SIGMOID: x(1) height of curve; x(2) heat-parametersc; x(3) center point

tic
%% LOOP THROUGH EFFECTORS
for i = 1 : 99
    if i ~= 97
    clear eff_abbr & eff_name & eff_set & eff_num & eff_mz
    
    % get effector info
    eff_abbr = effsinfo_b(i,find(strcmp(effsinfo_labels,"abbreviation")));
    eff_name = effsinfo_b(i,find(strcmp(effsinfo_labels,"metabolite")));
    eff_set = effsinfo_a(i,find(strcmp(effsinfo_labels,"set")));
    eff_num = effsinfo_a(i,find(strcmp(effsinfo_labels,"num")));
    eff_mz = effMASS_mass(find(strcmp(effMASS_name,eff_abbr)));
    
    % flag if  effector mass overlaps with known substrate/product
    mz_overlap_index = find((abs(eff_mz-reactant_masses) > 0.005) == 0);
    if ~isempty(mz_overlap_index)
        fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),mz_overlap_index) = {0};
        if string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "SUBSTRATE" & string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "PRODUCT"
            fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6) = {"OVERLAP with "+string(mz_overlap_index)};
        end
    end
    
    %% LOOP THROUGH REACTANTS
    for reactant_nr = 1 : 4
        clear bt_y_fit & bt_y_rescaled & impossible_fit_counter & bt_y_fit_AUC & bt_y_fit_r2 % clear parameters that become final output
        clear datatable & TICtable
        if reactant_nr == 1
            datatable = fiadata.S1;
            TICtable = fiadata.S1_TICs;
        elseif reactant_nr == 2
            datatable = fiadata.S2;
            TICtable = fiadata.S2_TICs;
        elseif reactant_nr == 3
            datatable = fiadata.P1;
            TICtable = fiadata.P1_TICs;
        elseif reactant_nr == 4
            datatable = fiadata.P2;
            TICtable = fiadata.P2_TICs;
        end
        % get TICmedian for TICnormalization
        TICall = cell2mat(TICtable(2:end,7:end));
        TICmedian = nanmedian(TICall(:));
        
        % continue if not all datatable entries are nan
        if ~isempty(datatable)
        Y_missing = sum(sum(ismissing(cell2mat(datatable(2:end,7:end)))));
        if Y_missing < size(cell2mat(datatable(2:end,7:end)),1) * size(cell2mat(datatable(2:end,7:end)),2)
            clear xopt & y_fit & fit_y_delta & y_fit_scaled & y_REscaled & y_raw & y_scaled & sampleidx & x & x_range & tic_raw
            %% GET SAMPLE INDICES
            if ~startsWith(eff_abbr,"ctrl_") % standard
                sampleidx = find(strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(eff_set))... % get indices of set
                    & strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"EFFECTOR")))),"eff"+string(eff_num))); % get indices of effector in set
            else % for lumping of controls
                if strcmp(eff_abbr,"ctrl_all") % all controls in sets 1 - 8
                    % get index of controls
                    sampleidx = find(strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"EFFECTOR")))),"eff"+string(eff_num)));
                    if remove_ctrls_ratio_median == 1
                        % remove series with bad R2
                        react_ctrls_excl = ctrls_excl(reactant_nr,:); % r2 of controls
                        react_ctrls_out_idx = find(abs(react_ctrls_excl / nanmedian(react_ctrls_excl)) < plot_r2_thresh); % find all that are off by 30% or more from median
                        sampleidx(react_ctrls_out_idx) = [];
                    elseif remove_ctrls_ratio_median == 0.7
                        react_ctrls_excl = ctrls_excl(reactant_nr,:); 
                        react_ctrls_out_idx = find(react_ctrls_excl < 0.7);
                        sampleidx(react_ctrls_out_idx) = [];
                    end
                elseif strcmp(eff_abbr,"ctrl_1st") % all controls in sets 1 - 4
                    % get index of controls
                    sampleidx = find((strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(1))...
                        | strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(2))...
                        | strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(3))...
                        | strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(4)))...
                        & strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"EFFECTOR")))),"eff"+string(12)));
                    if remove_ctrls_ratio_median == 1
                        % remove series with bad R2
                        react_ctrls_excl = ctrls_excl(reactant_nr,1:length(sampleidx)); % r2 of controls of first series
                        react_ctrls_out_idx = find(abs(react_ctrls_excl / nanmedian(react_ctrls_excl)) < plot_r2_thresh); % find all that are off by 30% or more from median
                        sampleidx(react_ctrls_out_idx) = [];
                    elseif remove_ctrls_ratio_median == 0.7
                        react_ctrls_excl = ctrls_excl(reactant_nr,1:length(sampleidx));
                        react_ctrls_out_idx = find(react_ctrls_excl < 0.7);
                        sampleidx(react_ctrls_out_idx) = [];
                        
                        %%% QC
                        fiadata_effs.QC_ctrl_1to4_Nexcl07(reactant_nr) = length(find(react_ctrls_excl < 0.7));
                        
                    end
                    if length(sampleidx) < 5
                        "ERROR: BAD CONTROLS FOR REACTANT " + string(reactant_nr)
                    end
                elseif strcmp(eff_abbr,"ctrl_2nd") % all controls in sets 5 - 8
                    sampleidx = find((strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(5))...
                        | strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(6))...
                        | strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(7))...
                        | strcmp(string(datatable(:,find(strcmp(string(datatable(1,:)),"SET")))),"set"+string(8)))...
                        & strcmp(string(datatable(:,strcmp(string(datatable(1,:)),"EFFECTOR"))),"eff"+string(12)));
                    if remove_ctrls_ratio_median == 1
                        % remove series with bad R2
                        react_ctrls_excl = ctrls_excl(reactant_nr,1+length(ctrls_excl)-length(sampleidx):end); % r2 of controls of second series
                        react_ctrls_out_idx = find(abs(react_ctrls_excl / nanmedian(react_ctrls_excl)) < plot_r2_thresh); % find all that are off by 30% or more from median
                        sampleidx(react_ctrls_out_idx) = [];
                    elseif remove_ctrls_ratio_median == 0.7
                        react_ctrls_excl = ctrls_excl(reactant_nr,1+length(ctrls_excl)-length(sampleidx):end);
                        react_ctrls_out_idx = find(react_ctrls_excl < 0.7);
                        sampleidx(react_ctrls_out_idx) = [];
                        
                        
                        %%% QC
                        fiadata_effs.QC_ctrl_5to8_Nexcl07(reactant_nr) = length(find(react_ctrls_excl < 0.7));
                        
                        
                    end
                    if length(sampleidx) < 5
                        "ERROR: BAD CONTROLS FOR REACTANT " + string(reactant_nr)
                    end
                end
            end
            
            if sampleidx

            
            
            %% GET X AND Y, PRESCALE
            x = cellfun(@str2num,fiadata.time(sampleidx,:)) / 60; % [s] to [min]
            x_range = min(min(x)) : 0.1 : max(max(x)); % for fits & plotting
            
            y_raw = cell2mat(datatable(sampleidx,7:end)); % measured ion count
            
            %% FILTER 1: low ion counts         
            if reactant_nr <= 2
                y_raw(y_raw<ioncount_threshold) = nan; % remove ions below threshold
            elseif reactant_nr > 2
                excl_idx = y_raw<ioncount_threshold;
                excl_idx(1) = 0; %  leave t0 if lower in case of products
                y_raw(excl_idx) = nan; % remove ions below threshold
            end
            
            %% FILTER 2: TIC outliers
            % get TICs and exclude in case of irregularities (e.g. bad injections)
            tic_raw = cell2mat(TICtable(sampleidx,7:end)); % getting TICs of each measurements
            if TICfilter == "CUT50"
                y_raw(abs(tic_raw ./ nanmedian(tic_raw,2)-1) > .5) = nan; % if TIC is less than 50% of median TIC = bad injection (-> remove from Y)
            elseif TICfilter == "CUT20"
                y_raw(abs(tic_raw ./ nanmedian(tic_raw,2)-1) > .8) = nan;
            end      
            
            %% Normalize to TIC
            if TIC_normalization == 1
                y_raw = (y_raw ./ tic_raw) * TICmedian; % normalized ion count = ion count / TIC * median(all TICs)
            else
                y_raw = y_raw;
            end

            %% FILTER 3: timepoint 1 OFF?
            % if techn replicate of first sample too different, exclude assay
            if ~startsWith(eff_abbr,"ctrl_")
                if exclude_all_if_t0_compromised == 1 % exclude whole bio replicate, if t0 is off
                    if  (abs(y_raw(1,1)-y_raw(2,1))/nanmean([max(y_raw(1,:)) max(y_raw(2,:))])) > t0_cv_cutoff %% normal CV doesnt work, instead DELTA of first injection vs highest counts
                        y_raw(1,:) = nan;
                        y_raw(2,:) = nan;
                    end
                    if (abs(y_raw(3,1)-y_raw(4,1))/nanmean([max(y_raw(3,:)) max(y_raw(4,:))])) > t0_cv_cutoff
                        y_raw(3,:) = nan;
                        y_raw(4,:) = nan;
                    end
                    % same if technical TICs very different
                    if (abs(tic_raw(1,1)-tic_raw(2,1))/nanmean(tic_raw(1:2,1))) > t0_cv_cutoff
                        y_raw(1,:) = nan;
                        y_raw(2,:) = nan;
                    end
                    if (abs(tic_raw(3,1)-tic_raw(4,1))/nanmean(tic_raw(3:4,1))) > t0_cv_cutoff
                        y_raw(3,:) = nan;
                        y_raw(4,:) = nan;
                    end
                else % only exclude t0 but keep rest if t0 is off
                    if (abs(y_raw(1,1)-y_raw(2,1))/nanmean([max(y_raw(1,:)) max(y_raw(2,:))])) > t0_cv_cutoff %% normal CV doesnt work, instead DELTA of first injection vs highest counts
                        y_raw(1,1) = nan;
                        y_raw(2,1) = nan;
                    end
                    if (abs(y_raw(3,1)-y_raw(4,1))/nanmean([max(y_raw(3,:)) max(y_raw(4,:))])) > t0_cv_cutoff
                        y_raw(3,1) = nan;
                        y_raw(4,1) = nan;
                    end
                    % same if technical TICs very different
                    if (abs(tic_raw(1,1)-tic_raw(2,1))/nanmean(tic_raw(1:2,1))) > t0_cv_cutoff
                        y_raw(1,1) = nan;
                        y_raw(2,1) = nan;
                    end
                    if (abs(tic_raw(3,1)-tic_raw(4,1))/nanmean(tic_raw(3:4,1))) > t0_cv_cutoff
                        y_raw(3,1) = nan;
                        y_raw(4,1) = nan;
                    end
                end
            end

            %% PLOT 1: normalized ion counts
            if ~startsWith(eff_abbr,"ctrl_")
                figure(10+eff_set)
            else
                figure(20) % control plots separately
            end
            subplot(2,2,reactant_nr)
            hold on
            if errorbars_ioncount == "std"
                y_raw_err = nanstd(y_raw);
            elseif errorbars_ioncount == "eom"
                y_raw_err = nanstd(y_raw)/sqrt(size(y_raw,1));
            end
            if log_x_plot == 1
                a_raw = plot(log2(x(:)+1),y_raw(:),".",'Color',dist_colors(eff_num,:)); % raw data points
                a_err = errorbar(nanmean(log2(x+1)),nanmean(y_raw),y_raw_err,'Color',dist_colors(eff_num,:)); % error bars
            else
                a_raw = plot(x(:),y_raw(:),".",'Color',dist_colors(eff_num,:)); % raw data points
                a_err = errorbar(nanmean(x),nanmean(y_raw),y_raw_err,'Color',dist_colors(eff_num,:)); % error bars
            end
            
            x = x(find(~all(isnan(y_raw),2)),:);
            y_raw = y_raw(find(~all(isnan(y_raw),2)),:);

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% BOOTSTRAP %%% BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%BOOTSTRAP %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % prepare tables and counters
            clear bt_y_rescaled & bt_y_fit & bt_y_fit_95lb & bt_y_fit_95ub & bt_y_fit_r2 & bt_y_fit_AUC & impossible_fit_counter
            bt_y_rescaled(1:nboots,1:size(x,2)) = nan;
            bt_y_fit(1:nboots,1:length(x_range)) = nan;
            bt_y_fit_95lb(1:nboots,1:length(x_range)) = nan;
            bt_y_fit_95ub(1:nboots,1:length(x_range)) = nan;
            bt_y_fit_r2(1:nboots) = nan;
            bt_y_fit_AUC(1:nboots) = nan;
            bt_y_fit_AUC_all(1:nboots) = nan;
            meanAUCidx = nan;
            impossible_fit_counter = 0;
            

            if ~isempty(y_raw)
            if split_bio_replicates == 0 % do we treat bio replicates separately or combine all 4 replicates?
                if ~all(isnan(y_raw(:,1))) % if not all t0 are nans, else FLAG
                    if length(y_raw(sum(isnan(y_raw)) == size(y_raw,1))) < 3 % 2 or less datapoints missing, else FLAG                    
                        % bootstrap
                        % --> OVERLAY NORMALIZED ION COUNTS (NO SCALING YET)
                        % get mean and error of raw ion counts
                        y_raw_mean = nanmean(y_raw,1); % mean
                        if bootstrap_err == "std" % either std or eom
                            y_raw_err = nanstd(y_raw);
                        elseif bootstrap_err == "eom"
                            for t_reps = 1 : size(y_raw,2)
                                y_raw_err(t_reps) = nanstd(y_raw(:,t_reps)) / sqrt(length(~isnan(y_raw(:,t_reps))));
                            end
                        end
                        
                        for boot = 1 : nboots
                            clear bt_y & tbl & mdl & xopt & ci & y_fit & y_fit_95lb & y_fit_95ub
                            

                            
                            % SAMPLE FROM MEAN + STD
                            bt_y = y_raw_mean + randn * y_raw_err; % get bootstrap sample based on standarddeviation

                            % scale each boot-sample drawn from raw ion counts between 0 and 1 before fitting (allows to choose same set of initial parameters for the fitting function across all experiments)
                            bt_y_scaled = rescale(bt_y);
                            tbl = table(nanmean(x,1)',bt_y_scaled'); % rearrange

                            % ONESUBSTRATE lambertW
                            if fit_function == "lambertW"
                                fit_impossible = 0;
                                try
                                    if reactant_nr <= 2 % =substrates
                                        mdl = fitnlm(tbl, fun_lambertW_subs, [1 0.01 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                        y_fit = xopt(3) * lambertw(0, xopt(1) * exp(xopt(1) - xopt(2)*x_range)); % = fun_lambertW_subs
                                        y_fit_95lb = xopt(3) * lambertw(0, xopt(1) * exp(xopt(1) - ci(2,1) * x_range));
                                        y_fit_95ub = xopt(3) * lambertw(0, xopt(1) * exp(xopt(1) - ci(2,2) * x_range));
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_lambertW_prod, [1 0.01 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                        y_fit = xopt(3) * (xopt(1) - lambertw(0, xopt(1) * exp(xopt(1) - xopt(2) * x_range))); % = fun_lambertW_prod
                                        y_fit_95lb = xopt(3) * (xopt(1) - lambertw(0, xopt(1) * exp(xopt(1) - ci(2,1) * x_range)));
                                        y_fit_95ub = xopt(3) * (xopt(1) - lambertw(0, xopt(1) * exp(xopt(1) - ci(2,2) * x_range)));
                                    end
                                catch
                                    fit_impossible = 1
                                end
                                % TWOSUSBTRATE lambertW
                            elseif fit_function == "2SlambertW"
                                fit_impossible = 0;
                                try
                                    if reactant_nr <= 2 % =substrates
                                        mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_subs, [0.01 0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                        y_fit = xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_subs
                                        y_fit_95lb =  xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,1) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,1) * x_range).^2 + 4));
                                        y_fit_95ub =  xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,2) * x_range).^2 + 4));
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_prod, [1 0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                        y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
                                        y_fit_95lb = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,1) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,1) * x_range).^2 + 4));
                                        y_fit_95ub = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,2) * x_range).^2 + 4));
                                    end
                                catch
                                    fit_impossible = 1
                                end
                                % SIGMOID
                            elseif fit_function == "SIGMOID"
                                fit_impossible = 0;
                                try
                                    if reactant_nr <= 2 % =substrates
                                        mdl = fitnlm(tbl, fun_SIGMOID_subs, [1 -0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                        y_fit = 1 - xopt(1) * tanh(abs(xopt(2)/xopt(1))*(x_range-xopt(3)));
                                        y_fit_95lb = 1 - xopt(1) * tanh(abs(ci(2,1)/xopt(1))*(x_range-xopt(3)))';
                                        y_fit_95ub = 1 - xopt(1) * tanh(abs(ci(2,2)/xopt(1))*(x_range-xopt(3)))';
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_SIGMOID_prod, [1 -0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                        y_fit = xopt(1) * tanh(abs(xopt(2)/xopt(1))*(x_range-xopt(3)));
                                        y_fit_95lb = xopt(1) * tanh(abs(ci(2,1)/xopt(1))*(x_range-xopt(3)))';
                                        y_fit_95ub = xopt(1) * tanh(abs(ci(2,2)/xopt(1))*(x_range-xopt(3)))';
                                    end
                                catch
                                    fit_impossible = 1
                                end
                            end
                            
                            % save parameters
                            if fit_impossible == 0
                                bt_y_rescaled(boot,:) = (bt_y_scaled - min(y_fit(:))) / (max(y_fit(:)) - min(y_fit(:))); % datapoints rescaled based on fit
                                bt_y_fit(boot,:) = rescale(y_fit); % rescale again between 0 and 1
                                bt_y_fit_95lb(boot,:) = rescale(y_fit_95lb);
                                bt_y_fit_95ub(boot,:) = rescale(y_fit_95ub);
                                bt_y_fit_r2(boot) = mdl.Rsquared.Adjusted;
                                if reactant_nr <= 2 % = substrate
                                    bt_y_fit_AUC(boot) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range,rescale(y_fit));
                                elseif reactant_nr >= 3 % = product
                                    bt_y_fit_AUC(boot) = trapz(x_range,rescale(y_fit));
                                end
                            elseif fit_impossible == 1
                                impossible_fit_counter = impossible_fit_counter + 1;                                %count and    %% FLAG
                            end
                        end
                        
                        if remove_ctrls_ratio_median == 1
                            bt_y_fit_AUC_all = bt_y_fit_AUC;
                            [excl_r2_idx] = find(bt_y_fit_r2<r2_threshold);
                            bt_y_fit_AUC(excl_r2_idx) = nan;
                            impossible_fit_counter = impossible_fit_counter + length(excl_r2_idx);
                        elseif remove_ctrls_ratio_median == 0.7
                            bt_y_fit_AUC_all = bt_y_fit_AUC;
                            [excl_r2_idx] = find(bt_y_fit_r2<r2_threshold);
                            bt_y_fit_AUC(excl_r2_idx) = nan;
                            impossible_fit_counter = impossible_fit_counter + length(excl_r2_idx);
                        end


                        % get 95% of population that remains
                        % fits and bad fits (r2<.7) are excluded)
                        [~,max5] = maxk(bt_y_fit_AUC,ceil(0.025*(nboots-impossible_fit_counter))); % max 2.5 percent
                        [~,min5] = mink(bt_y_fit_AUC,ceil(0.025*(nboots-impossible_fit_counter))); % min 2.5 percent
                        AUC_list = bt_y_fit_AUC; % save AUCs in separate table
                        AUC_list([max5;min5]) = nan; % remove top and low 2.5%
                        bt_y_fit_95 = bt_y_fit;
                        bt_y_fit_95(find(isnan(AUC_list)),:) = nan; % only keep 95%
                        [~,u95] = maxk(AUC_list,1); % find index of highest AUC within 95%
                        [~,l95] = mink(AUC_list,1); % find index of lowest AUC within 95%
                        bt_y_rescaled_95 = bt_y_rescaled;
                        bt_y_rescaled_95(find(isnan(AUC_list)),:) = nan; % only keep 95%
                        
                        % plot 95%
                        if ~startsWith(eff_abbr,"ctrl_")
                            figure(110+reactant_nr)
                            subplot(8,12,i)
                            hold on
                            %
                            markedd=0;
                            if ~cellfun(@isnumeric, fiadata_effs.analysis_legend(i,6)) % if effector is flagged as S, P or overlap?
                                if contains(string(fiadata_effs.analysis_legend(i,6)),"SUBSTRATE")
                                    text(max(x(:))/4,0.5,"SUBSTRATE");
                                    markedd = 1;
                                elseif contains(string(fiadata_effs.analysis_legend(i,6)),"PRODUCT")
                                    text(max(x(:))/4,0.5,"PRODUCT");
                                    markedd = 1;
                                end
                            end
                            if markedd == 0 & impossible_fit_counter > 0% all others
                                text(max(x(:))/4,0.5,"ex:"+string(impossible_fit_counter))
                            end

                           if  plot_only_95conf == 1
                            shadedplot(x_range, bt_y_fit_95(u95,:), bt_y_fit_95(l95,:), [0.3010 0.7450 0.9330],  [0.3010 0.7450 0.9330]);
                           else
                             [~,lb_idx]  = min(AUC_list);
                               [~,ub_idx]  = max(AUC_list);
                            shadedplot(x_range,  bt_y_fit(lb_idx,:),  bt_y_fit(ub_idx,:), [0.3010 0.7450 0.9330],  [0.3010 0.7450 0.9330]);
                           end
                            
                            hold on
                            plot(x(1,:),bt_y_rescaled_95,'.','Color',[0, 0.4470, 0.7410],'MarkerSize',5)
                            
                            [~,meanAUCidx] = min(abs(bt_y_fit_AUC-nanmean(bt_y_fit_AUC))); % find value closest to MEAN
                            hold on
                            plot(x_range,bt_y_fit(meanAUCidx,:),'b--') % plot

                        % plot controls on top of all effectors
                        elseif ctrl_split_on == 0
                            if strcmp(eff_abbr,"ctrl_all")
                                figure(110+reactant_nr)
                                hold on
                                for plotn = 1 : (8*12)
                                    subplot(8,12,plotn)
                                    hold on
                                    shadectrl = shadedplot(x_range, bt_y_fit_95(u95,:), bt_y_fit_95(l95,:), [0 0 0],  [0 0 0]);
                                end
                                %
                                figure(120+reactant_nr)
                                hold on
                                for plotn = 1 : (8*12)
                                    subplot(8,12,plotn)
                                    hold on
                                    shadectrl = shadedplot(x_range,SIG1fit_y_fit_95ub,SIG1fit_y_fit_95lb, [0 0 0],  [0 0 0]);
                                end
                                
                            end
                        elseif ctrl_split_on == 1
                            if strcmp(eff_abbr,"ctrl_1st")
                                figure(110+reactant_nr)
                                hold on
                                for plotn = 1 : (8*12)/2
                                    subplot(8,12,plotn)
                                    hold on
                                    shadectrl = shadedplot(x_range, bt_y_fit_95(u95,:), bt_y_fit_95(l95,:), [0 0 0],  [0 0 0]);
                                end
                                %
                                try % sigmoid plot, not important
                                    figure(120+reactant_nr)
                                    hold on
                                    for plotn = 1 : (8*12)/2
                                        subplot(8,12,plotn)
                                        hold on
                                        shadectrl = shadedplot(x_range,SIG1fit_y_fit_95ub,SIG1fit_y_fit_95lb, [0 0 0],  [0 0 0]);
                                    end
                                end
                            elseif strcmp(eff_abbr,"ctrl_2nd")
                                figure(110+reactant_nr)
                                hold on
                                for plotn = (8*12)/2+1 : (8*12)
                                    subplot(8,12,plotn)
                                    hold on
                                    shadectrl = shadedplot(x_range, bt_y_fit_95(u95,:), bt_y_fit_95(l95,:), [0 0 0],  [0 0 0]);
                                end
                                try
                                    figure(120+reactant_nr)
                                    hold on
                                    for plotn = (8*12)/2+1 : (8*12)
                                        subplot(8,12,plotn)
                                        hold on
                                        shadectrl = shadedplot(x_range,SIG1fit_y_fit_95ub,SIG1fit_y_fit_95lb, [0 0 0], [0 0 0]);
                                    end
                                end
                            end
                        end
                        
                        %% SIMULATE CHANGE IN BETA parameter (kcat/Km) and write down and plot translation
                        if startsWith(eff_abbr,"control")
                            clear vary_BETA & scaled_vary_BETA & AUC2activity_AUC & tbl & sepctrlsAUC2activity_AUC

                            vary_BETA(1:size(AUC2activity_activity,2),1:size(x_range,2)) = nan;
                            scaled_vary_BETA(1:size(AUC2activity_activity,2),1:size(x_range,2)) = nan;

                            auc_act_xvalues(1:size(bt_y_rescaled_95,1),1:length(nanmean(x))) = nan;
                            for rows_n = 1 : size(bt_y_rescaled_95,1)
                                auc_act_xvalues(rows_n,1:length(nanmean(x))) = nanmean(x); %
                            end

                            tbl = table(auc_act_xvalues(:),bt_y_rescaled_95(:)); % rearrange

                            % get optimal fit
                            xopt = nan;
                            try
                                if fit_function == "2SlambertW"
                                    if reactant_nr <= 2
                                        mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_subs, [0.01 0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_prod, [1 0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    end
                                elseif fit_function == "lambertW"
                                    if reactant_nr <= 2 % =substrates
                                        mdl = fitnlm(tbl, fun_lambertW_subs, [1 0.01 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_lambertW_prod, [1 0.01 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    end
                                end
                            catch
                                "CTRL ERROR"
                            end

                            if ~isnan(xopt)
                                % simulate change in activity based on beta
                                for f = 1 : size(AUC2activity_activity,2)
                                    if fit_function == "2SlambertW"
                                        if reactant_nr <= 2
                                            vary_BETA(f,:) = xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - (xopt(2)*AUC2activity_activity(f)) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) -     (xopt(2)*AUC2activity_activity(f))	* x_range).^2 + 4));
                                        elseif reactant_nr >= 3 % =products
                                            vary_BETA(f,:) = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2)*AUC2activity_activity(f))	* x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2)*AUC2activity_activity(f))	* x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
                                        end
                                    elseif fit_function == "lambertW"
                                        if reactant_nr <= 2
                                            vary_BETA(f,:) = xopt(3) * lambertw(0, xopt(1) * exp(xopt(1) -       (xopt(2)*AUC2activity_activity(f))      * x_range)); % xopt(2) == BETA
                                        elseif reactant_nr >= 3
                                            vary_BETA(f,:) = xopt(3) * (xopt(1) - lambertw(0, xopt(1) * exp(xopt(1) -       (xopt(2)*AUC2activity_activity(f))      * x_range))); % xopt(2) == BETA
                                        end
                                    elseif fit_function == "SIGMOID"
                                        "WRITE AUC->BETA FOR SIGMOID!"
                                    end
                                end
                            end
                            % rescale all simulated fits with varied BETA separately
                            scaled_vary_BETA(1:size(vary_BETA,1),1:size(vary_BETA,2)) = nan;
                            for betas = 1 : size(vary_BETA,1)
                                scaled_vary_BETA(betas,:) = rescale(vary_BETA(betas,:));
                            end
                            % extract AUC
                            sepctrlsAUC2activity_AUC(1:size(scaled_vary_BETA,1)) = nan;
                            if reactant_nr <= 2 % substrates
                                for sims = 1 : size(scaled_vary_BETA,1)
                                    sepctrlsAUC2activity_AUC(sims) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range,scaled_vary_BETA(sims,:));
                                end
                            elseif reactant_nr >= 3 % products
                                for sims = 1 : size(scaled_vary_BETA,1)
                                    sepctrlsAUC2activity_AUC(sims) = trapz(x_range,scaled_vary_BETA(sims,:));
                                end
                            end
                            clear allctrl_AUC2beta

                            if reactant_nr == 1
                                fiadata_SIMS.S1(i/12,:) = sepctrlsAUC2activity_AUC;
                            elseif reactant_nr == 2
                                fiadata_SIMS.S2(i/12,:) = sepctrlsAUC2activity_AUC;
                            elseif reactant_nr == 3
                                fiadata_SIMS.P1(i/12,:) = sepctrlsAUC2activity_AUC;
                            elseif reactant_nr == 4
                                fiadata_SIMS.P2(i/12,:) = sepctrlsAUC2activity_AUC;
                            end

                            figure(800+reactant_nr)
                            title("SEP CONTROLS: AUC to ß")
                            hold on
                            %fuse and save tables
                            allctrl_AUC2beta = [sepctrlsAUC2activity_AUC',AUC2activity_activity'];
                            fiadata_effs.AUC2activity_all(:,reactant_nr) = sepctrlsAUC2activity_AUC';
                            plot(sepctrlsAUC2activity_AUC,AUC2activity_activity,'.')
                            xlabel("AUC")
                            ylabel("change in ß")

                        end

                        %% LUMPING within ctrl1 and ctrl2
                        if startsWith(eff_abbr,"ctrl_")
                            clear vary_BETA & scaled_vary_BETA & AUC2activity_AUC & tbl
                            xopt = [nan nan nan];
                            vary_BETA(1:size(AUC2activity_activity,2),1:size(x_range,2)) = nan;
                            scaled_vary_BETA(1:size(AUC2activity_activity,2),1:size(x_range,2)) = nan;

                            if get_transl_table == "only_mean_AUC"
                                %% take only mean(AUC)
                                % X
                                auc_act_xvalues = nanmean(x)';
                                % Y
                                [~,mean_idx] = min(abs(bt_y_fit_AUC_all-nanmean(bt_y_fit_AUC_all)));
                                auc_act_yvalues = bt_y_rescaled(mean_idx,:);
                                tbl = table(auc_act_xvalues,auc_act_yvalues'); % rearrange

                            elseif get_transl_table == "rescale_together"
                                sim_y_input = rescale(y_raw);
                                sim_x_input = x;
                            elseif get_transl_table == "rescale_each"
                                auc_act_xvalues(1:size(sim_y_input,1),1:length(nanmean(x))) = nan;
                                for rows_n = 1 : size(sim_y_input,1)
                                    auc_act_xvalues(rows_n,1:length(nanmean(x))) = nanmean(x); %
                                end
                                tbl = table(sim_x_input(:),sim_y_input(:)); % rearrange
                            elseif get_transl_table ==  "overlay1to4and5to8"
                                clear sim_y_input & sim_y_r & sim_x_input & tbl
                                if strcmp(eff_abbr,"ctrl_all")
                                    sim_y_input = cell2mat(fiadata_effs.Yvalues(12:12:96,reactant_nr));
                                    sim_y_r2 = cell2mat(fiadata_effs.R2s(12:12:96,reactant_nr)');
                                elseif strcmp(eff_abbr,"ctrl_1st")
                                    sim_y_input = cell2mat(fiadata_effs.Yvalues(12:12:48,reactant_nr));
                                    sim_y_r2 = cell2mat(fiadata_effs.R2s(12:12:48,reactant_nr)');
                                elseif strcmp(eff_abbr,"ctrl_2nd")
                                    sim_y_input = cell2mat(fiadata_effs.Yvalues(60:12:96,reactant_nr));
                                    sim_y_r2= cell2mat(fiadata_effs.R2s(60:12:96,reactant_nr)');
                                end
                                sim_y_input(sim_y_r2<0.7,:)= nan;
                                sim_x_input(1:size(sim_y_input,1),1:length(nanmean(x,1))) = nan;
                                for rows_n = 1 : size(sim_y_input,1)
                                    sim_x_input(rows_n,1:length(nanmean(x,1))) = nanmean(x,1); %
                                end
                                tbl = table(sim_x_input(:),sim_y_input(:)); % rearrange
                            end
                            
                            
                            % get single optimal fit for visualisation
                            try
                            if fit_function == "2SlambertW"
                                    if reactant_nr <= 2
                                        if string(enzymelabel) == "Pta" % different startparameters for this one, otherwise fit fails
                                             mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_subs, [0.1 0.1 0.1],'Options',opts); % fitting procedure
                                        else
                                            mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_subs, [0.01 0.1 0.1],'Options',opts); % fitting procedure
                                        end
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_prod, [1 0.1 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    end
                                elseif fit_function == "lambertW"
                                    if reactant_nr <= 2 % =substrates
                                        mdl = fitnlm(tbl, fun_lambertW_subs, [1 0.01 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    elseif reactant_nr >= 3 % =products
                                        mdl = fitnlm(tbl, fun_lambertW_prod, [1 0.01 0.1],'Options',opts); % fitting procedure
                                        xopt = mdl.Coefficients{:, 'Estimate'};
                                        ci = coefCI(mdl);
                                    end
                                end
                            end

                            % simulate change in activity based on beta
                            for f = 1 : size(AUC2activity_activity,2)
                                if fit_function == "2SlambertW"
                                    if reactant_nr <= 2
                                        vary_BETA(f,:) = xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - (xopt(2)*AUC2activity_activity(f)) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) -     (xopt(2)*AUC2activity_activity(f))	* x_range).^2 + 4));
                                    elseif reactant_nr >= 3 % =products
                                        vary_BETA(f,:) = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2)*AUC2activity_activity(f))	* x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2)*AUC2activity_activity(f))	* x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
                                    end
                                elseif fit_function == "lambertW"
                                    if reactant_nr <= 2
                                        vary_BETA(f,:) = xopt(3) * lambertw(0, xopt(1) * exp(xopt(1) -       (xopt(2)*AUC2activity_activity(f))      * x_range)); % xopt(2) == BETA
                                    elseif reactant_nr >= 3
                                        vary_BETA(f,:) = xopt(3) * (xopt(1) - lambertw(0, xopt(1) * exp(xopt(1) -       (xopt(2)*AUC2activity_activity(f))      * x_range))); % xopt(2) == BETA
                                    end
                                elseif fit_function == "SIGMOID"
                                    "WRITE AUC->BETA FOR SIGMOID!"
                                end
                            end
                            
                            % rescale all simulated fits with varied BETA separately
                            scaled_vary_BETA(1:size(vary_BETA,1),1:size(vary_BETA,2)) = nan;
                            for betas = 1 : size(vary_BETA,1)
                                scaled_vary_BETA(betas,:) = rescale(vary_BETA(betas,:));
                            end
                            % extract AUC
                            AUC2activity_AUC(1:size(scaled_vary_BETA,1)) = nan;
                            if reactant_nr <= 2 % substrates
                                for sims = 1 : size(scaled_vary_BETA,1)
                                    AUC2activity_AUC(sims) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range,scaled_vary_BETA(sims,:));
                                end
                            elseif reactant_nr >= 3 % products
                                for sims = 1 : size(scaled_vary_BETA,1)
                                    AUC2activity_AUC(sims) = trapz(x_range,scaled_vary_BETA(sims,:));
                                end
                            end
                            clear allctrl_AUC2beta
                            if strcmp(eff_abbr,"ctrl_all")
                                figure(400)
                                title("ALL CONTROLS: AUC to ß")
                                hold on
                                %fuse and save tables
                                allctrl_AUC2beta = [AUC2activity_AUC',AUC2activity_activity'];
                                fiadata_effs.AUC2activity_all(:,reactant_nr) = AUC2activity_AUC';
                            elseif strcmp(eff_abbr,"ctrl_1st")
                                figure(401)
                                title("CONTROLS set1-set4: AUC to ß")
                                hold on
                                ctrl1st_AUC2beta  = [AUC2activity_AUC',AUC2activity_activity'];
                                fiadata_effs.AUC2activity_1st(:,reactant_nr) = AUC2activity_AUC';
                            elseif strcmp(eff_abbr,"ctrl_2nd")
                                figure(402)
                                title("CONTROLS set5-set8: AUC to ß")
                                hold on
                                ctrl2nd_AUC2beta = [AUC2activity_AUC',AUC2activity_activity'];
                                fiadata_effs.AUC2activity_2nd(:,reactant_nr) = AUC2activity_AUC';
                            end
                            plot(AUC2activity_AUC,AUC2activity_activity,'.')
                            xlabel("AUC")
                            ylabel("change in ß")
                            
                            % Same for sigmoid fit (obsolete)

                            clear SIG_vary_BETA & SIG_scaled_vary_BETA & SIG_AUC2activity_AUC & SIG_tbl &  SIG_AUC2activity_AUC

                            xopt = [nan nan nan];


                            SIGvary_BETA(1:size(AUC2activity_activity,2),1:size(x_range,2)) = nan;
                            SIGscaled_vary_BETA(1:size(AUC2activity_activity,2),1:size(x_range,2)) = nan;


                            auc_act_xvalues(1:size(bt_y_rescaled_95,1),1:length(nanmean(x))) = nan;
                            for rows_n = 1 : size(bt_y_rescaled_95,1)
                                auc_act_xvalues(rows_n,1:length(nanmean(x))) = nanmean(x); %
                            end
                            Sfit_tbl = table(auc_act_xvalues(:),bt_y_rescaled_95(:)); % rearrange

                            % get optimal SIG fit
                            try
                                if reactant_nr <= 2 % =substrates
                                    mdl = fitnlm(Sfit_tbl, fun_SIGMOID_subs, [1 -0.1 0.1],'Options',opts); % fitting procedure
                                    xopt = mdl.Coefficients{:, 'Estimate'};
                                    ci = coefCI(mdl);
                                elseif reactant_nr >= 3 % =products
                                    mdl = fitnlm(Sfit_tbl, fun_SIGMOID_prod, [1 -0.1 0.1],'Options',opts); % fitting procedure
                                    xopt = mdl.Coefficients{:, 'Estimate'};
                                    ci = coefCI(mdl);
                                end
                            catch
                                "SIG CTRL ERROR"
                            end
                            try
                                % simulate change in activity based on beta
                                for f = 1 : size(AUC2activity_activity,2)
                                    % SIGMOID FIT
                                    if reactant_nr <= 2
                                        SIGvary_BETA(f,:) = 1 - xopt(1) * tanh(abs((xopt(2)*AUC2activity_activity(f))/xopt(1))*(x_range-xopt(3))); % xopt2 = heat parameter
                                    elseif reactant_nr >= 3
                                        SIGvary_BETA(f,:) =  xopt(1) * tanh(abs((xopt(2)*AUC2activity_activity(f))/xopt(1))*(x_range-xopt(3))); % xopt2 = heat parameter
                                    end
                                end


                                % rescale all simulated fits with varied BETA separately
                                SIGscaled_vary_BETA(1:size(SIGvary_BETA,1),1:size(SIGvary_BETA,2)) = nan;
                                for betas = 1 : size(SIGvary_BETA,1)
                                    SIGscaled_vary_BETA(betas,:) = rescale(SIGvary_BETA(betas,:));
                                end
                                % extract AUC
                                SIG_AUC2activity_AUC(1:size(SIGscaled_vary_BETA,1)) = nan;
                                if reactant_nr <= 2 % substrates
                                    for sims = 1 : size(SIGscaled_vary_BETA,1)
                                        SIG_AUC2activity_AUC(sims) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range,SIGscaled_vary_BETA(sims,:));
                                    end
                                elseif reactant_nr >= 3 % products
                                    for sims = 1 : size(SIGscaled_vary_BETA,1)
                                        SIG_AUC2activity_AUC(sims) = trapz(x_range,SIGscaled_vary_BETA(sims,:));
                                    end
                                end


                                clear allctrl_AUC2beta & ctrl1st_SIG_AUC2beta & ctrl2nd_SIG_AUC2beta
                                if strcmp(eff_abbr,"ctrl_all")
                                    figure(900)
                                    title("SIGMOID ALL CONTROLS: AUC to ß")
                                    hold on
                                    %fuse and save tables
                                    allctrl_SIG_AUC2beta = [SIG_AUC2activity_AUC',AUC2activity_activity'];
                                    fiadata_effs.SIG_AUC2activity_all(:,reactant_nr) = SIG_AUC2activity_AUC';
                                elseif strcmp(eff_abbr,"ctrl_1st")
                                    figure(901)
                                    title("SIGMOID CONTROLS set1-set4: AUC to ß")
                                    hold on
                                    ctrl1st_SIG_AUC2beta  = [SIG_AUC2activity_AUC',AUC2activity_activity'];
                                    fiadata_effs.SIG_AUC2activity_1st(:,reactant_nr) = SIG_AUC2activity_AUC';
                                elseif strcmp(eff_abbr,"ctrl_2nd")
                                    figure(902)
                                    title("SIGMOID CONTROLS set5-set8: AUC to ß")
                                    hold on
                                    ctrl2nd_SIG_AUC2beta = [SIG_AUC2activity_AUC',AUC2activity_activity'];
                                    fiadata_effs.SIG_AUC2activity_2nd(:,reactant_nr) = SIG_AUC2activity_AUC';
                                end
                                plot(SIG_AUC2activity_AUC,AUC2activity_activity,'.')
                                xlabel("AUC")
                                ylabel("change in ß")
                            end
                        end
                        %%%% END OF SIGMOID FIT TEST



                        % FLAG IF MORE THAN 2 DATAPOINTS MISSING
                    else % -> more than 2 datapoints completely missing
                        if i < 97
                            figure(110+reactant_nr)
                            subplot(8,12,i)
                            hold on
                            text(max(x(:))/4,0.5,"3+ off")
                            xlim([0 max(x(:))])
                            ylim([0 1])
                        end
                        % find position in overview
                        fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),reactant_nr) = {0};
                        % if not already marked as reactant or overlap
                        if string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "SUBSTRATE" & string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "PRODUCT" & string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "OVERLAP"
                            if strcmp(num2str(fiadata_effs.analysis_legend{find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6}),'NaN')
                                fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6) = {append("reactant",string(reactant_nr),"-3+missing")};
                            else
                                fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6) = {append(string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)),"; reactant",string(reactant_nr),"-3+missing")};
                            end
                        end
                    end
                    % FLAG IF T0 IS COMPROMISED
                else % -> all t0s missing
                    % new
                    if i < 97
                        figure(110+reactant_nr)
                        subplot(8,12,i)
                        hold on
                        text(max(x(:))/4,0.5,"t0 off")
                        xlim([0 max(x(:))])
                        ylim([0 1])
                    end
                    % find position in overview
                    fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),reactant_nr) = {0};
                    % if not already marked as reactant or overlap
                    if string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "SUBSTRATE" & string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "PRODUCT" & string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)) ~= "OVERLAP"
                        if strcmp(num2str(fiadata_effs.analysis_legend{find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6}),'NaN') % is current field empty or already note inside?
                            fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6) = {append("reactant",string(reactant_nr),"-t0compromised")};
                        else
                            fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6) = {append(string(fiadata_effs.analysis_legend(find(strcmp(fiadata_effs.analysis_legend(:,5),eff_abbr)),6)),"; reactant",string(reactant_nr),"-t0compromised")};
                        end
                    end
                end
            end
            end % end if all y_raw nan
            end % end if sampleidx empty
        end
        if sampleidx
            % SAVE RESULTS IN STRUCT
            fiadata_effs.Yfits(i,reactant_nr) = {bt_y_fit}; % XVALUES x 200 fits
            fiadata_effs.Yvalues(i,reactant_nr)  = {bt_y_rescaled}; % XVALUES x 200 fits
            fiadata_effs.impossible_fits(i,reactant_nr)  = {impossible_fit_counter}; % 1 counter per eff and reactant
            fiadata_effs.AUCs(i,reactant_nr) = {bt_y_fit_AUC}; % 200 AUCs per eff and reactant after removing r2<0.7
            fiadata_effs.AUCsALL(i,reactant_nr) = {bt_y_fit_AUC_all}; % including AUCs excluded because r2<0.7
            fiadata_effs.meanAUCs(i,reactant_nr) = {nanmean(bt_y_fit_AUC)};
            if ~isnan(meanAUCidx)
                fiadata_effs.meanAUCsCURVE(i,reactant_nr) = {bt_y_fit(meanAUCidx,:)};
            else
                fiadata_effs.meanAUCsCURVE(i,reactant_nr) = {nan};
            end
            fiadata_effs.stdAUCs(i,reactant_nr) = {nanstd(bt_y_fit_AUC)};
            fiadata_effs.cvAUCs(i,reactant_nr) = {nanstd(bt_y_fit_AUC)/nanmean(bt_y_fit_AUC)*100};
            fiadata_effs.R2s(i,reactant_nr) = {bt_y_fit_r2}; % 200 R2s per eff and reactant

        end
        end
    end
  end
end
% save separately simulated series into overview
fiadata_effs.MEANsepctrlAUC2activity_1st(:,1) = nanmean(fiadata_SIMS.S1(1:4,:))';
fiadata_effs.MEANsepctrlAUC2activity_1st(:,2) = nanmean(fiadata_SIMS.S2(1:4,:))';
fiadata_effs.MEANsepctrlAUC2activity_1st(:,3) = nanmean(fiadata_SIMS.P1(1:4,:))';
fiadata_effs.MEANsepctrlAUC2activity_1st(:,4) = nanmean(fiadata_SIMS.P2(1:4,:))';
fiadata_effs.MEANsepctrlAUC2activity_2nd(:,1) = nanmean(fiadata_SIMS.S1(5:8,:))';
fiadata_effs.MEANsepctrlAUC2activity_2nd(:,2) = nanmean(fiadata_SIMS.S2(5:8,:))';
fiadata_effs.MEANsepctrlAUC2activity_2nd(:,3) = nanmean(fiadata_SIMS.P1(5:8,:))';
fiadata_effs.MEANsepctrlAUC2activity_2nd(:,4) = nanmean(fiadata_SIMS.P2(5:8,:))';

% legend plots
figure(400)
legend("P1","P2","S1","S2",'location','northwest')
if save_all_results == 1
    % savefig(append(enzymelabel,"-bterr-",bootstrap_err,"-AUC2activity-ctrls_all"))
end

figure(900)
legend("P1","P2","S1","S2",'location','northwest')

 for i = 1 : 8
     figure(800+i)
     legend("control1","control2","control3","control4","control5","control6","control7","control8")
     if save_all_results == 1
     % savefig(append(enzymelabel,"-bterr-",bootstrap_err,"sep-AUC2activity-reactant-",string(i)))
     end
 end
 for j = 11:18
     figure(j)
     if save_all_results == 1
     % savefig(append(enzymelabel,"-bterr-",bootstrap_err,"-ioncounts-set-",errorbars_ioncount,string(j-10)))
     
    % label = string(fiadata.enzymename) + "_EFFECTORS_ioncounts_set" + string(j-10) + ".png"
    % saveas(gca,label)
     end
     end
 for k = 111:114
     figure(k)
     if save_all_results == 1
     % savefig(append(enzymelabel,"-bterr-",bootstrap_err,"-",string(k-110)))
     label = string(fiadata.enzymename) + "_EFFECTORS_fitted_reactant" + string(k-110) + ".png"
     saveas(gca,label)
     end
 end
 
 if plot_only_relevant == 1
     figure(19)
     close
     figure(121)
     close
     figure(122)
     close
     figure(123)
     close
     figure(124)
     close
     figure(40)
     close
     figure(41)
     close
     figure(900)
     close
 end


 
 toc