%% PROGRESS CURVE ANALYSIS from photometer data
% Input for this script is photometry data in 96-well format. Tested
% enzymes are dehydrogenases, which produce NAD(P)H that is photometrically
% active.
% It imports the data, corrects blanks, plots raw data, does fits based on
% Haldane kinetics and based on simulations from control assays, translates
% area under the curve into approximate changes in enzyme activity when an
% effector is present.
% Parameters: 340 nm measurement wavelength, 9 nm bandwidth,
% 25 flashes, 9 ms settle time.
close all
clear
clc

%% Inputs
analysis_path = "\\pasteur\SysBC-Home\chgruber\Desktop\==MANUSCRIPT WIP\==submission\=2=CODE\"; % path with saved folders
enzyme = ["Zwf"]; % enzyme name corresponding to CODETABLES.xlsx - photometer_overview sheet
cutoffend = 1;
cutoffend_fitted = 1;
AUC2activity_activity = [0.01 : 0.001 : 10];
maxTECANsignal = 3.5; % excluded if absorption higher than this value
newequ_threshold = 1.2; % if raw absorption more than 1.2 times of control, then flag for visual inspection (could be new equilibrium that is reached)
lowequ_delta_cutoff = 1; % if abs change in last third of assay is larger than 1.5*of max control, then still increasing
%samesetttestCUT = 0.05/11;

%% Get information
file = {analysis_path+"CODETABLES.xlsx"};; %path to screen overview
% on TECAN experiment
[overview_info_A,overview_info_B] = xlsread(file{1},"photometer_overview","A1 : AF6"); %import screen overview table with datafile, substrate and product names
Eidx = find(strcmp(string(overview_info_B(:,find(strcmp(string(overview_info_B(1,:)),"enzyme")))),enzyme)); % idx of enzyme in overview table
TECAN_rawfile = analysis_path + "=2_input=PHOTOMETERdata\" + string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"TECAN-file")))); % data path
TECAN_setup = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"setup")))); % experimental setup
TECAN_date = string(overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"date")))); % experimental setup
s1_name = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"reactants-1"))));
s2_name = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"reactants-2"))));
p1_name = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"reactants-3"))));
p2_name = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"reactants-4"))));
n_replicates = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"n_replicates")));
precip_strong = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"precipitation strong"))));
precip_weak = string(overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"precipitation weak"))));
salt_calc = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"calc")));
salt_2licl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"2licl")));
salt_2kcl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"2kcl")));
salt_5nacl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"5nacl")));
salt_mgcl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"mgcl")));
rescale_ub = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"rescale upper bound")));
manual_lower_equ = overview_info_B(Eidx,find(strcmp(string(overview_info_B(1,:)),"new lower equ")));
% on effector names
[effsinfo_a,effsinfo_b,effsinfo_c] = xlsread(file{1},"effectors","A1 : D100"); %import screen overview table with datafile, substrate and product names
effsinfo_labels = effsinfo_b(1,:);
effsinfo_b(1,:) = [];
% replicates
replicates = [1:n_replicates];
% end of assay
assayend = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"assayend [min]"))); % normally 128 or 64
x_range = 0 : 0.1 : assayend;

%% Preparing plots
fig1 = figure(73642)
sgtitle(enzyme+"-TECAN raw")
fig1.WindowState = 'maximized';
hold on
subplot(8,12,1)
for i = 1 : 8*12
    subplot(8,12,i)
    subtitle(string(effsinfo_b(i,4)))
    xlabel("time [min]")
    ylabel("absorbtion")
end
hold on

fig2 = figure(73643)
sgtitle(enzyme+"-TECAN scaled & fitted")
fig2.WindowState = 'maximized';
hold on
subplot(8,12,1)
for i = 1 : 8*12
    subplot(8,12,i)
    subtitle(string(effsinfo_b(i,4)))
    xlabel("time [min]")
    ylabel("absorbtion")
end
hold on

fig3 = figure(73644)
sgtitle(enzyme+"-TECAN raw scaled")
fig3.WindowState = 'maximized';
hold on
subplot(8,12,1)
for i = 1 : 8*12
    subplot(8,12,i)
    subtitle(string(effsinfo_b(i,4)))
    xlabel("time [min]")
    ylabel("absorbtion")
end
hold on

%% Prepare fitting (all 2SlambertW)
% fitting function for substrate and product
fun_TWOSUBSTRATElambertW_subs = @(x,t) x(3)/2 * (x(1)/x(3) - x(3)/x(1) - x(2) * t + sqrt((x(1)/x(3) - x(3)/x(1) - x(2) * t).^2 + 4)); % s0 = x(1); b = x(2); k_m = x(3);
fun_TWOSUBSTRATElambertW_prod = @(x,t) 1 - x(3)/2 * (x(1)/x(3) - x(3)/x(1) - x(2) * t + sqrt((x(1)/x(3) - x(3)/x(1) - x(2) * t).^2 + 4)); % s0 = x(1); b = x(2); k_m = x(3);

%% Prepare table for metabolite-specific notes
met_notes = effsinfo_b;
met_notes(:,1) = {1};
met_notes(find(matches(met_notes(:,4),[s1_name,s2_name])),2) = {"S"};
met_notes(find(matches(met_notes(:,4),[p1_name,p2_name])),2) = {"P"};
met_notes(find(matches(met_notes(:,4),split(precip_strong,',')')),2) = {"PRECIP"}; % strong precipitation identified in well of 96-well plate
met_notes(find(matches(met_notes(:,4),split(precip_weak,',')')),2) = {"precip"}; % weak precipitation identified in well of 96-well plate
met_notes(find(matches(met_notes(:,4),["nadh","nadph","fad"])),2) = {"oor"}; % out of range

%% Prepare table to save control time series
CTRL_table_rawX = [];
CTRL_table_rawY = [];
CTRL_table_rawYscaled = [];
CTRL_table_scaledX = [];
CTRL_table_scaledY = [];

%% Prepare output tables
AUCs_raw(1:size(met_notes,1)-3,1:n_replicates+1) = {nan};
AUCs_raw(:,end) = met_notes(1:end-3,4);
AUCs_scaled = AUCs_raw;
AUCs_scaled_fitted = AUCs_raw;
fitr2_scaled = AUCs_raw; % r2 of fit
fit_imposs = AUCs_raw; % flag if fit was not possible
abs_equ = AUCs_raw; % absorption of last 5 RAW datapoints (~absorption at equilibrium)
equcheck_y = AUCs_raw; % save last third of y (and x) values to check if equilibrium was reached
equcheck_x = AUCs_raw; % save last third of y (and x) values to check if equilibrium was reached

%% Find common xvalue cutoff
for j = 1 : n_replicates
    blanktime = xlsread(TECAN_rawfile{1},strcat(enzyme,"-",string(replicates(j))),"B34"); % time [sec] passed between adding enzyme and TECAN starts measuring
    [~,~,data] = xlsread(TECAN_rawfile{1},strcat(enzyme,"-",string(replicates(j))),"A36 : ABS133");

    xvalues_sec_pre = [0 cell2mat(data(find(strcmp(string(data(:,1)),'Time [s]')),2:end))+blanktime]; % 0-point, assay factoring in delay in pipetting(=blanktime) and delay during measuring each well (scantime)
    xvalues_pre = xvalues_sec_pre/60; % convert from sec to min

    assayend_idx = find(xvalues_pre>assayend); % cut off end based on set end time
    xvalues_pre(min(assayend_idx):end) = [];

    assayend_cutoff_raw(j) = max(xvalues_pre); % for each replicate: highest value that is smaller than manually set assayend
end
assayend_cutoff = min(assayend_cutoff_raw); % earliest assayend across replicates

%% Loop and plot
% for each replicate
for j = 1 : n_replicates
    [~,~,blank] = xlsread(TECAN_rawfile{1},strcat(enzyme,"-",string(replicates(j)),"blank"),"A36 : D133");
    % get absorptions for blank (normally three consecutive measurements before adding enzyme)
    blanktime = xlsread(TECAN_rawfile{1},strcat(enzyme,"-",string(replicates(j))),"B34"); % time [sec] passed between adding enzyme and TECAN starts measuring
    if isempty(blanktime) == 1 % if blank time wasnt noted down, assume 30 seconds
        blanktime = 30;
    end
    [~,~,data] = xlsread(TECAN_rawfile{1},strcat(enzyme,"-",string(replicates(j))),"A36 : ABS133");
    scantime = cell2mat(data(find(strcmp(string(data(:,1)),'Time [s]')),3))/(size(data,1)-2); %% time between timepoint 1 and 2 (=delay in scanning all 96 wells)

    %% Manual removal
    % manually excluded last control (no E added?)
    if sum([enzyme == 'Gnd' enzyme == 'Eno' enzyme == 'Zwf' enzyme == 'Icd']) == 1
        data(98,2:end) = {nan};
        blank(98,2:end) = {nan};
    end
    % manually exclude dttp replicate from MaeB (was not pipetted)
    if enzyme == 'MaeB' & TECAN_date == '20210521'
        data(find(strcmp(data(:,1),'dttp')),2:end) = {nan};
    end

    % for each metabolite
    for i = 1 : size(met_notes,1)-3
        % find index of metabolites in data
        met_idx = find(strcmp(data(:,1),met_notes(i,4)));

        % get x-values
        xvalues_sec = [0 cell2mat(data(find(strcmp(string(data(:,1)),'Time [s]')),2:end))+blanktime+scantime*i]; %0-point, assay factoring in delay in pipetting(=blanktime) and delay during measuring each well (scantime)
        xvalues = xvalues_sec/60; % convert from sec to min
        assayend_idx = find(xvalues>assayend); % cut off end based on set end time
        xvalues(min(assayend_idx):end) = [];

        % change 'OVER' to nan
        ydata = data(met_idx,2:end);
        for y = 1 : length(ydata)
            if strcmp(ydata{y},'OVER')
                ydata{y} = NaN;
            end
        end
        blankdata = blank(met_idx,2:end);
        for y = 1 : length(blankdata)
            if strcmp(blankdata{y},'OVER')
                blankdata{y} = NaN;
            end
        end

        % create y-data array
        yvalues_raw = [nanmean(cell2mat(blankdata)) cell2mat(ydata)]; % first datapoint = mean of 3 blankreads; remaining ones from datatable
        yvalues_raw = yvalues_raw - nanmean(cell2mat(blankdata)); % subtract blank read
        % remove values higher than detection limit
        yvalues_raw(yvalues_raw>maxTECANsignal) = nan; % remove all with signal > max (out of range)
        yvalues_raw(assayend_idx) = [];
        yvalues = yvalues_raw;
        yvalues(min(assayend_idx):end) = [];
        % save absorbtion at end of assay (~ equilibrium)
        abs_equ(i,j) = {nanmean(yvalues(end-5:end))};
        equcheck_y(i,j) = {[yvalues(end-round(length(yvalues)/3):end)]}; % save last third of yvalues
        equcheck_x(i,j) = {[xvalues(end-round(length(yvalues)/3):end)]}; % save last third of xvalues

        % extract AUC
        if ~isempty(yvalues) %& yvalues~=0
            interp_xrange = x_range(x_range<=max(xvalues));
            y_interp = interp1(xvalues,yvalues,interp_xrange); % interpolate to allow to extract AUC from exactly same timespan in all assays
            y_interp_scaled = interp1(xvalues,rescale(yvalues,0,rescale_ub),interp_xrange); % same for rescaled y (between 0 and 1)
            if cutoffend == 1
                AUCs_raw(i,j) = {trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_interp(interp_xrange<=assayend_cutoff))}; % extract area under curve
                AUCs_scaled(i,j) = {trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_interp_scaled(interp_xrange<=assayend_cutoff))}; % extract area under curve
            elseif cutoffend == 0
                AUCs_raw(i,j) = {trapz(interp_xrange,y_interp)};
                AUCs_scaled(i,j) = {trapz(interp_xrange,y_interp_scaled)};
            end

            if max(interp_xrange) < assayend_cutoff % catch error
                "CHECK ASSAYEND CODING!!"
            end
        end

        % plot raw absorbtions
        figure(73642)
        subplot(8,12,i)
        hold on
        plot(xvalues,yvalues,'Color',[0.3010 0.7450 0.9330])
        plot(xvalues,yvalues,'.','MarkerSize',3,'Color',[0, 0.4470, 0.7410])
        if ~cellfun(@isempty,met_notes(i,2)) & j == 1
            text(assayend/4,0.25,string(met_notes(i,2)),'Color','black','FontSize',8)
        end

        figure(73644)
        subplot(8,12,i)
        hold on
        plot(xvalues,rescale(yvalues,0,rescale_ub),'Color',[0.3010 0.7450 0.9330])
        plot(xvalues,rescale(yvalues,0,rescale_ub),'.','MarkerSize',3,'Color',[0, 0.4470, 0.7410])
        if ~cellfun(@isempty,met_notes(i,2)) & j == 1
            text(assayend/4,0.25,string(met_notes(i,2)),'Color','black','FontSize',8)
        end






        %% FIT AND RESCALE based on Michaelis-Menten
        % remove nans
        yvalues(isnan(xvalues)) = [];
        xvalues(isnan(xvalues)) = [];
        xvalues(isnan(yvalues)) = [];
        yvalues(isnan(yvalues)) = [];

        % scale & rearrange
        y_scaled = rescale(yvalues,0,rescale_ub); % scale between 0 and 1 (same method that is used for MS data)
        tbl = table(xvalues',y_scaled'); % rearrange

        fit_impossible = 0;
        clear mdl & xopt & y_fit & y_fit_95lb & y_fit_95ub
        try
            mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_prod, [1 0.1 0.1]); % fitting procedure
            xopt = mdl.Coefficients{:, 'Estimate'};
            ci = coefCI(mdl);
            y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
            y_fit_95lb = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,1) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,1) * x_range).^2 + 4));
            y_fit_95ub = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - ci(2,2) * x_range).^2 + 4));
            fitr2_scaled(i,j) = {mdl.Rsquared.Ordinary};
            fit_imposs(i,j) = {fit_impossible};
        catch
            fit_impossible = 1;
            fit_imposs(i,j) = {fit_impossible};
        end

        if fit_impossible == 0


            y_backscaled = (y_scaled - min(y_fit(:))) / (max(y_fit(:)) - min(y_fit(:))); % datapoints rescaled based on fit
            y_fit_rescaled = rescale(y_fit,0,rescale_ub); % rescale again between 0 and 1
            y_fit_AUC = trapz(x_range,y_fit_rescaled);

            figure(73643)
            subplot(8,12,i)
            hold on
            plot(x_range,y_fit_rescaled,'Color',[0.3010 0.7450 0.9330])
            if ~cellfun(@isempty,met_notes(i,2)) & j == 1
                text(assayend/4,0.25,string(met_notes(i,2)),'Color','black','FontSize',8)
            end

            %% Extract area under curve
            % for raw absorption
            if ~isempty(yvalues)
                if cutoffend_fitted == 1
                    AUCs_scaled_fitted(i,j) = {trapz(x_range(x_range<=assayend_cutoff),y_fit_rescaled(x_range<=assayend_cutoff))}; % extract area under curve, but remove 2 min to allow exact same range across assays
                elseif cutoffend_fitted == 0
                    AUCs_scaled_fitted(i,j) = {trapz(x_range,y_fit_rescaled)};
                end
            end
        end

        if mod(i,12) == 0 % if i divisable by 12 (== CONTROL)
            CTRL_table_rawX = [CTRL_table_rawX xvalues];
            CTRL_table_rawY = [CTRL_table_rawY yvalues];
            CTRL_table_rawYscaled = [CTRL_table_rawYscaled rescale(yvalues,0,rescale_ub)];
            CTRL_collect_X(i/12,j) = {xvalues};
            CTRL_collect_Y(i/12,j) = {yvalues};
            % get raw AUC from control
            if cutoffend == 1
                beta_AUCs_raw(i/12,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_interp(interp_xrange<=assayend_cutoff)); % extract area under curve
            elseif cutoffend == 0
                beta_AUCs_raw(i/12,j) = trapz(interp_xrange,y_interp);
            end

            if fit_impossible == 0
                CTRL_table_scaledX = [CTRL_table_scaledX x_range];
                CTRL_table_scaledY = [CTRL_table_scaledY y_fit_rescaled];
            end
        end
    end
end

% add lines into plot for visual representation showing if threshold is surpassed
beta_up_weak = 1.4801; % 2^0.56565 = 1.4142
beta_down_weak = 1/1.4801;
beta_up_strong = 2.1906;  % 2^1.1313 = 2
beta_down_strong = 1/2.1906;


%% fit each control, vary BETA and extract AUC
figure(987)
hold on
figure(456)
hold on

%% calculate change in activity
% 1. SCALE FROM 0 TO 1
% 2. FIT AND EXTRACT AUC
% 3. CHANGE BETA AND EXTRACT AUC
% 4. RATIO OF 2 AND 3 * AUCRAW

for j = 1 : size(CTRL_collect_Y,2)
    for i = 1 : size(CTRL_collect_Y,1)
        x = cell2mat(CTRL_collect_X(i,j));
        y = cell2mat(CTRL_collect_Y(i,j));
        y_scaled = rescale(y,0,rescale_ub); % rescale_ub chosen manually and noted in excel table to optimize fit
        if ~isempty(y)
            interp_xrange = x_range(x_range<=max(x));
            y_interp = interp1(x,y,interp_xrange); % interpolate to allow to extract AUC from exactly same timespan in all assays
            y_interp_scaled = interp1(x,rescale(y,0,rescale_ub),interp_xrange); % same for rescaled y (between 0 and 1)
            if cutoffend == 1
                sim_AUCs_raw(i,j) = {trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_interp(interp_xrange<=assayend_cutoff))}; % extract area under curve
                sim_AUCs_scaled(i,j) = {trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_interp_scaled(interp_xrange<=assayend_cutoff))}; % extract area under curve
            elseif cutoffend == 0
                sim_AUCs_raw(i,j) = {trapz(interp_xrange,y_interp)};
                sim_AUCs_scaled(i,j) = {trapz(interp_xrange,y_interp_scaled)};
            end

            % raw plot
            hold on
            figure(987)
            plot(x,y_scaled,'k.')
            figure(456)
            plot(x,y,'k.')

            % DO BEST FIT AND EXTRACT AUC
            mdl = fitnlm(table(x',y_scaled'), fun_TWOSUBSTRATElambertW_prod, [1 0.1 0.1]); % fitting procedure
            xopt = mdl.Coefficients{:, 'Estimate'};
            ci = coefCI(mdl);
            y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * interp_xrange + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * interp_xrange).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
            y_fit = rescale(y_fit,0,rescale_ub)
            figure(987)
            plot(interp_xrange,y_fit,'k','LineWidth',1)
            figure(456)
            y_backscaled = rescale(y_fit,0,(y_fit(end)/y_scaled(end))*y(end));
            plot(interp_xrange,y_backscaled,'k','LineWidth',1)
            % extract
            if cutoffend == 1
                y_fit_interp = interp1(interp_xrange,y_fit,interp_xrange);
                AUCs_best_fit_scaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_fit_interp(interp_xrange<=assayend_cutoff));
                %
                y_rescaled_interp = interp1(interp_xrange,y_backscaled,interp_xrange);
                AUCs_best_fit_rescaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_rescaled_interp(interp_xrange<=assayend_cutoff));
            elseif cutoffend == 0
                AUCs_best_fit_scaled(i,j) = trapz(interp_xrange,y_fit); % extract area under curve, but remove 2 min to allow exact same range across assays
                AUCs_best_fit_rescaled(i,j) = trapz(interp_xrange,y_backscaled);
            end

            % up weak
            multiplier_crtrl = beta_up_weak;
            y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
            y_fit = rescale(y_fit,0,rescale_ub)
            figure(987)
            plot(interp_xrange,y_fit,'k','LineWidth',1,'Color',[0.6473 0.7456 0.4188])
            figure(456)
            y_backscaled = rescale(y_fit,0,(max(y_fit)/max(y_scaled))*max(y));
            plot(interp_xrange,y_backscaled,'LineWidth',1,'Color',[0.6473 0.7456 0.4188])
            % extract
            if cutoffend == 1
                y_fit_interp = interp1(interp_xrange,y_fit,interp_xrange);
                AUCs_up_weak_scaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_fit_interp(interp_xrange<=assayend_cutoff));
                %
                y_rescaled_interp = interp1(interp_xrange,y_backscaled,interp_xrange);
                AUCs_up_weak_rescaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_rescaled_interp(interp_xrange<=assayend_cutoff));
            elseif cutoffend == 0
                AUCs_up_weak_scaled(i,j) = trapz(interp_xrange,y_fit); % extract area under curve, but remove 2 min to allow exact same range across assays
                AUCs_up_weak_rescaled(i,j) = trapz(interp_xrange,y_backscaled);
            end
            sim_ratio_up_weak_0_1_scaled(i,j) = AUCs_up_weak_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_scaled(i,j));
            sim_ratio_up_weak_rawscaled(i,j) = AUCs_up_weak_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_raw(i,j));

            % down weak
            multiplier_crtrl = beta_down_weak;
            y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
            y_fit = rescale(y_fit,0,rescale_ub)
            hold on
            figure(987)
            plot(interp_xrange,y_fit,'k','LineWidth',1,'Color',[0.3010 0.7450 0.9330])
            figure(456)
            y_backscaled = rescale(y_fit,0,(max(y_fit)/max(y_scaled))*max(y));
            plot(interp_xrange,y_backscaled,'LineWidth',1,'Color',[0.3010 0.7450 0.9330])
            % extract
            if cutoffend == 1
                y_fit_interp = interp1(interp_xrange,y_fit,interp_xrange);
                AUCs_down_weak_scaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_fit_interp(interp_xrange<=assayend_cutoff));
                y_rescaled_interp = interp1(interp_xrange,y_backscaled,interp_xrange);
                AUCs_down_weak_rescaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_rescaled_interp(interp_xrange<=assayend_cutoff));
            elseif cutoffend == 0
                AUCs_down_weak_scaled(i,j) = trapz(interp_xrange,y_fit); % extract area under curve, but remove 2 min to allow exact same range across assays
                AUCs_down_weak_rescaled(i,j) = trapz(interp_xrange,y_backscaled);
            end
            sim_ratio_down_weak_0_1_scaled(i,j) = AUCs_down_weak_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_scaled(i,j));
            sim_ratio_down_weak_rawscaled(i,j) = AUCs_down_weak_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_raw(i,j));

            % up strong
            multiplier_crtrl = beta_up_strong;
            y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
            y_fit = rescale(y_fit,0,rescale_ub)
            figure(987)
            plot(interp_xrange,y_fit,'k','LineWidth',1,'Color',[0.3500 0.6740 0.1880])
            figure(456)
            y_backscaled = rescale(y_fit,0,(max(y_fit)/max(y_scaled))*max(y));
            plot(interp_xrange,y_backscaled,'LineWidth',1,'Color',[0.3500 0.6740 0.1880])
            % extract
            if cutoffend == 1
                y_fit_interp = interp1(interp_xrange,y_fit,interp_xrange);
                AUCs_up_strong_scaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_fit_interp(interp_xrange<=assayend_cutoff));
                y_rescaled_interp = interp1(interp_xrange,y_backscaled,interp_xrange);
                AUCs_up_strong_rescaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_rescaled_interp(interp_xrange<=assayend_cutoff));
            elseif cutoffend == 0
                AUCs_up_strong_scaled(i,j) = trapz(interp_xrange,y_fit); % extract area under curve, but remove 2 min to allow exact same range across assays
                AUCs_up_strong_rescaled(i,j) = trapz(interp_xrange,y_backscaled);
            end
            sim_ratio_up_strong_0_1_scaled(i,j) = AUCs_up_strong_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_scaled(i,j));
            sim_ratio_up_strong_rawscaled(i,j) = AUCs_up_strong_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_raw(i,j));

            % down strong
            multiplier_crtrl = beta_down_strong;
            y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - (multiplier_crtrl*xopt(2)) * interp_xrange).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
            y_fit = rescale(y_fit,0,rescale_ub)
            figure(987)
            plot(interp_xrange,y_fit,'k','LineWidth',1,'Color',[0 0.4470 0.7410])
            figure(456)
            y_backscaled = rescale(y_fit,0,(max(y_fit)/max(y_scaled))*max(y));
            plot(interp_xrange,y_backscaled,'LineWidth',1,'Color',[0 0.4470 0.7410])
            % extract
            if cutoffend == 1
                y_fit_interp = interp1(interp_xrange,y_fit,interp_xrange);
                AUCs_down_strong_scaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_fit_interp(interp_xrange<=assayend_cutoff));
                y_rescaled_interp = interp1(interp_xrange,y_backscaled,interp_xrange);
                AUCs_down_strong_rescaled(i,j) = trapz(interp_xrange(interp_xrange<=assayend_cutoff),y_rescaled_interp(interp_xrange<=assayend_cutoff));
            elseif cutoffend == 0
                AUCs_down_strong_scaled(i,j) = trapz(interp_xrange,y_fit); % extract area under curve, but remove 2 min to allow exact same range across assays
                AUCs_down_strong_rescaled(i,j) = trapz(interp_xrange,y_backscaled);
            end
            sim_ratio_down_strong_0_1_scaled(i,j) = AUCs_down_strong_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_scaled(i,j));
            sim_ratio_down_strong_rawscaled(i,j) = AUCs_down_strong_scaled(i,j)/AUCs_best_fit_scaled(i,j)*cell2mat(sim_AUCs_raw(i,j));
        end
    end
end


%% Plot controls (black) on top of effector data (blue)
for i = 1 : 96
    % raw controls
    figure(73642)
    subplot(8,12,i)
    hold on
    ctrl_plot1 = plot(CTRL_table_rawX,CTRL_table_rawY,'k.','MarkerSize',3);
    figure(73644)
    subplot(8,12,i)
    hold on
    ctrl_plot3 = plot(CTRL_table_rawX,CTRL_table_rawYscaled,'k.','MarkerSize',3);
    figure(73643)
    subplot(8,12,i)
    hold on
    ctrl_plot2 = plot(CTRL_table_scaledX,CTRL_table_scaledY,'k.','MarkerSize',3);
end

%% Flag (if new equilbrium is reached, r2 is bad, etc)
r2_flag = (sum(cell2mat(fitr2_scaled(:,1:n_replicates))' > 0.7) <= 1)'; % flagged if 2 or more r2 < 0.7

% flag if comment in exp notes
exp_flag = string(met_notes(1:end-3,2)) ~= '';

% flag if equilibrium absorbance (raw, not scaled) is far away from control
newequ_overview(1:size(abs_equ,1)+1,1:7) = {nan};
newequ_overview(1,:) = [{"mean(equ_abs)"} {"ctrl(equ_abs)"} {"ctrl_delta"} {"allctrl(equ_abs)"} {"allctrl_delta"} {"effector"} {"flag_newequ"}];
newequ_overview(2:end,6) = met_notes(1:end-3,4);
mean_equ_abs = nanmean(cell2mat(abs_equ(:,1:n_replicates)),2);
newequ_overview(2:end,1) = num2cell(mean_equ_abs);
ctrls_abs_idx = ([12 : 12 : 96])
if sum([enzyme == 'Gnd' enzyme == 'Eno' enzyme == 'Zwf' enzyme == 'Icd']) == 1 % for some enzymes no 8th control, therefore use closest (=7)
    ctrls_abs_idx(end) = [];
end
% ctrl absorption at end of assay
abs_equ_ctrl = mean_equ_abs(ctrls_abs_idx);
abs_equ_ctrl_all = nanmean(abs_equ_ctrl);
% ctrl slope at end of assay
abs_equ_ctrl = mean_equ_abs(ctrls_abs_idx);
ctrl_slope(1:length(ctrls_abs_idx),1:n_replicates)=nan;
for i = 1 : length(ctrls_abs_idx)
    for j = 1 : n_replicates
        % slope
        x = cell2mat(equcheck_x(ctrls_abs_idx(i),j));
        y = cell2mat(equcheck_y(ctrls_abs_idx(i),j));
        P = polyfit(x,y,1); % linear fit
        ctrl_slope(i,j) =  P(1);
        % abs change
        ctrl_abschange(i,j) = y(end)-y(1);
    end
end
mean_ctrl_slope_all = mean(ctrl_slope(:));
mean_ctrl_abschange_all = nanmean(ctrl_abschange(:));
median_ctrl_abschange_all = nanmedian(ctrl_abschange(:));
max_ctrl_abschange_all = max(ctrl_abschange(:));

for i = 1 : length(mean_equ_abs)
    absequ_ctrl_sameset = mean_equ_abs(min(ctrls_abs_idx(ctrls_abs_idx>=i)));
    if ~isempty(absequ_ctrl_sameset)
        newequ_overview(i+1,2) = {absequ_ctrl_sameset}; % get control of same set
    else
        [~,closestidx] = min(abs(i-ctrls_abs_idx));
        absequ_ctrl_sameset = mean_equ_abs(ctrls_abs_idx(closestidx));
        newequ_overview(i+1,2) = {absequ_ctrl_sameset}; % get control of same set
    end
    newequ_overview(i+1,4) = {abs_equ_ctrl_all};

    % calculate delta [%] from
    delta_set = mean_equ_abs(i) / absequ_ctrl_sameset;
    newequ_overview(i+1,3) = {delta_set};
    delta_all = mean_equ_abs(i) / abs_equ_ctrl_all;
    newequ_overview(i+1,5) = {delta_all};

    % check whether new 'stable' equilibrium was reached
    flag = 0;
    flag2 = 0;
    % higher equilibrium?
    if delta_all > 1.2 & delta_set > 1.2
        flag = 1;

        % smaller equ (check if endpoint away from control and not increasing anymore
    elseif delta_all < 0.8 & delta_set < 0.8 % endpoint 20% away from control?
        for j = 1:n_replicates % calculate absolute change in last third of assay
            y = cell2mat(equcheck_y(i,j));
            abschange(j) = y(end)-y(1);
        end

        if max(ctrl_abschange(:)) < 0
            lowequ_threshold = 1-(lowequ_delta_cutoff-1) * max_ctrl_abschange_all;
        elseif max(ctrl_abschange(:)) > 0
            lowequ_threshold = lowequ_delta_cutoff * max_ctrl_abschange_all;
        end

        if nanmean(abschange) > lowequ_threshold
            % still increasing
        else
            flag2 = 1;
        end
    end
    newequ_h_overview(i+1,7) = {flag}
    newequ_l_overview(i+1,7) = {flag2}
end
flag_equ_h = cell2mat(newequ_h_overview(2:end,7));

%% Transfer flags into overview
%% AUTOMATIC flags
met_notes_flagged(1:size(met_notes,1)-3+1,1:18) = {nan};
met_notes_flagged(1,:) = [{"metabolite"} {"met"} {"exp_notes"} {"exp_flag"} {"flag_r2"} {"flag_equ_h"} {"flag_notes"} {"flag_total"} {""} {"MANUAL_check_equ_l"} {"MANUAL_salt_calc"} {"MANUAL_salt_2licl"} {"MANUAL_salt_2kcl"} {"MANUAL_salt_5nacl"} {"MANUAL_salt_mgcl"} {"MANUAL_salt_total"} {""} {"any flag"}]; % MANUAL_salt = ["licl" "nacl" etc] from excel % MANUAL_else "comment"
met_notes_flagged(2:end,1:2) = met_notes(1:end-3,3:4);
met_notes_flagged(2:end,3) = met_notes(1:end-3,2);
met_notes_flagged(2:end,4) = num2cell(exp_flag);
met_notes_flagged(2:end,5) = num2cell(r2_flag);
met_notes_flagged(2:end,8) = num2cell((exp_flag+r2_flag)>=1);
for i = 2 : size(met_notes_flagged,1)
    if cell2mat(met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"flag_total")))) == 1
        if cell2mat(met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"exp_flag")))) == 1
            met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))) = met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"exp_notes")));
        elseif cell2mat(met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"flag_r2")))) == 1
            met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))) = {"r2"};
        else
            "ERROR"
        end
    else
        met_notes_flagged(i,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))) = {""};
    end
end
%% MANUAL flags (saved in excel overview file)
empty_flags = zeros(size(met_notes_flagged,1)-1,1); % empty table
% lower equilibrium (2do: maybe automatize?)
man_lowerequ_flags = empty_flags;
man_lowerequ_flags(find(matches(met_notes(:,4),split(manual_lower_equ,',')'))) = 1;
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_check_equ_l"))) = num2cell(man_lowerequ_flags);

% salt activators and inhibitors
salt_calc_flagged = empty_flags;
if or(salt_calc == -1, salt_calc == 1)
    salt_calc_flagged(find(matches(met_notes(:,4),["calc","panto"]))) = salt_calc_flagged(find(matches(met_notes(:,4),["calc","panto"]))) + salt_calc;
end
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_calc"))) = num2cell(salt_calc_flagged);
%
salt_2licl_flagged = empty_flags;
if or(salt_2licl == -1, salt_2licl == 1)
    salt_2licl_flagged(find(matches(met_notes(:,4),["kdpg","acp","2licl","carb-p","accoa","glyc3p"]))) = salt_2licl_flagged(find(matches(met_notes(:,4),["kdpg","acp","2licl","carb-p","accoa","glyc3p"]))) + salt_2licl;
end
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_2licl"))) = num2cell(salt_2licl_flagged);
%
salt_2kcl_flagged = empty_flags;
if or(salt_2kcl == -1, salt_2kcl == 1)
    salt_2kcl_flagged(find(matches(met_notes(:,4),["2kcl","acp","icit","gal1p","f1p"]))) = salt_2kcl_flagged(find(matches(met_notes(:,4),["2kcl","acp","icit","gal1p","f1p"]))) + salt_2kcl;
end
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_2kcl"))) = num2cell(salt_2kcl_flagged);
%
salt_5nacl_flagged = empty_flags;
if or(salt_5nacl == -1, salt_5nacl == 1)
    salt_5nacl_flagged(find(matches(met_notes(:,4),["5nacl","bpg","nadph"]))) = salt_5nacl_flagged(find(matches(met_notes(:,4),["5nacl","bpg","nadph"]))) + salt_5nacl;
end
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_5nacl"))) = num2cell(salt_5nacl_flagged);
%
salt_mgcl_flagged = empty_flags;
if or(salt_mgcl == -1, salt_mgcl == 1)
    salt_mgcl_flagged(find(matches(met_notes(:,4),["mgcl","dhap"]))) = salt_mgcl_flagged(find(matches(met_notes(:,4),["mgcl","dhap"]))) + salt_mgcl;
end
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_mgcl"))) = num2cell(salt_mgcl_flagged);
% ALL manuals
manual_salt_total_flagged = abs(salt_calc_flagged)+abs(salt_2licl_flagged)+abs(salt_2kcl_flagged)+abs(salt_5nacl_flagged)+abs(salt_mgcl_flagged)
met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_total"))) = num2cell(manual_salt_total_flagged);
% ALL manuals + automatics
any_flag = (abs(manual_salt_total_flagged)+abs(man_lowerequ_flags)+exp_flag+r2_flag+flag_equ_h)>=1;
met_notes_flagged(2:end,strcmp(string(met_notes_flagged(1,:)),"any flag")) = num2cell(any_flag);


%% Simulate change in beta and save
% X and Y
tbl = table(CTRL_table_scaledX',CTRL_table_scaledY');

% do one fit based on data from ALL controls
mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_prod, [1 0.1 0.1]); % fitting procedure
xopt = mdl.Coefficients{:, 'Estimate'};
ci = coefCI(mdl);
bestfit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2))	* x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2))	* x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
if cutoffend == 1
    bestfit_AUC = trapz(x_range(x_range<=assayend_cutoff),bestfit(x_range<=assayend_cutoff));
else
    bestfit_AUC = trapz(x_range,bestfit);
end
% rescale and extract AUC - also remove last two minutes to make comparable
AUC2activity_AUC(1:length(AUC2activity_activity)) = nan;
for sims = 1 : length(AUC2activity_activity)
    vary_BETA = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2)*AUC2activity_activity(sims))	* x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) -	(xopt(2)*AUC2activity_activity(sims))	* x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
    scaled_vary_BETA = rescale(vary_BETA,0,rescale_ub); % rescale between 0 and 1 to make comparable
    if cutoffend == 1
        AUC2activity_AUC(sims) = trapz(x_range(x_range<=assayend_cutoff),scaled_vary_BETA(x_range<=assayend_cutoff));
        AUC2activity_AUCratio(sims) = AUC2activity_AUC(sims)/bestfit_AUC;
    else
        AUC2activity_AUC(sims) = trapz(x_range,scaled_vary_BETA);
        AUC2activity_AUCratio(sims) = AUC2activity_AUC(sims)/bestfit_AUC;
    end
end
% plot
figure(737400)
title(enzyme + ": AUC to ß")
hold on
AUC2activity = [AUC2activity_AUC',AUC2activity_activity'];
plot(AUC2activity_AUC,AUC2activity_activity,'.')
xlabel("AUC")
ylabel("change in ß")

% plot
figure(737401)
title(enzyme + ": AUC/AUCctrl to ß")
hold on
AUCratio2activity = [AUC2activity_AUCratio',AUC2activity_activity'];
plot(AUC2activity_AUCratio,AUC2activity_activity,'.')
xlabel("AUC/AUCctrl")
ylabel("change in ß")

figure(45251) % control figure
plot(CTRL_table_scaledX',CTRL_table_scaledY','.')
hold on
% bestfit
plot(x_range,bestfit,'k')
plot(x_range,scaled_vary_BETA,'b')



%% Translate
AUCs_scaled_fitted_translated(1:size(AUCs_scaled_fitted,1),1:n_replicates) = nan; % empty table
for r = 1 : n_replicates % loop through replicates
    AUCs_collection = cell2mat(AUCs_scaled_fitted(:,r));
    for k = 1 : length(AUCs_collection) % translate each AUC
        [c index] = min(abs(AUC2activity(:,1)-AUCs_collection(k)));
        AUCs_scaled_fitted_translated(k,r) = AUC2activity(index,2);
    end
end
log2activities_scaled_fitted = log2(AUCs_scaled_fitted_translated);

%% Output table
TECANdata = struct;
TECANdata.enzyme = enzyme;
TECANdata.date = TECAN_date;
TECANdata.n_replicates = n_replicates;
TECANdata.AUCs_raw = AUCs_raw;
TECANdata.AUCs_scaled = AUCs_scaled;
TECANdata.AUCs_scaled_fitted = AUCs_scaled_fitted;
TECANdata.fitR2 = fitr2_scaled;
TECANdata.AUC2activity = AUC2activity;
TECANdata.activities = AUCs_scaled_fitted_translated;
TECANdata.log2act = log2(AUCs_scaled_fitted_translated);

%% Create output tables for bar plot analysis
OUTPUT_legend = string(effsinfo_b(1:end-3,4))'; % done
OUTPUT_anyflag = zeros(1,length(effsinfo_b(1:end-3,4))); % done
OUTPUT_flags = strings(1,length(effsinfo_b(1:end-3,4)));
% raw
OUTPUT_raw_effectsize(1:length(effsinfo_b(1:end-3,4))) = zeros(1,length(effsinfo_b(1:end-3,4)));
OUTPUT_raw_3xSTD(1:length(effsinfo_b(1:end-3,4))) = zeros(1,length(effsinfo_b(1:end-3,4)));
% scaled
OUTPUT_scaled_effectsize(1:length(effsinfo_b(1:end-3,4))) = zeros(1,length(effsinfo_b(1:end-3,4)));
OUTPUT_scaled_3xSTD(1:length(effsinfo_b(1:end-3,4))) = zeros(1,length(effsinfo_b(1:end-3,4)));

%% Bar plot 1 - raw TECAN data

% RAW DATA
datapoints = cell2mat(TECANdata.AUCs_raw(:,1:n_replicates))';

idx_flagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"flag_total"))))==1);
idx_notflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"flag_total"))))==0);
idx_saltflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_total"))))==1);
idx_lowequflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_check_equ_l"))))==1);

% assays
reps_auc_raw = datapoints; % datapoints(samples)
if n_replicates > 1
    mean_auc_raw = nanmean(reps_auc_raw); % mean(samples)
    std_auc_raw = nanstd(reps_auc_raw); % std(samples)
else
    mean_auc_raw = reps_auc_raw; % mean(samples)
    std_auc_raw = zeros(1,length(reps_auc_raw))
end


% controls
ctrls_abs_idx = ([12 : 12 : 96])
if sum([enzyme == 'Gnd' enzyme == 'Eno' enzyme == 'Zwf' enzyme == 'Icd']) == 1 % for some enzymes no 8th control, therefore use closest (=7)
    ctrls_abs_idx(end) = [];
end
ctrl_reps_auc_raw = reps_auc_raw(:,ctrls_abs_idx); % datapoints(ctrl)
ctrl_mean_auc_raw = nanmean(ctrl_reps_auc_raw(:)); % mean(ctrl)
ctrl_std_auc_raw = nanstd(ctrl_reps_auc_raw(:)); % std(ctrl)

figure(12345)
sgtitle(string(enzyme)+" raw AUCs")
hold on
% plot assays
salt_abbr = ["calc ","2licl ","2kcl ","5nacl ","mgcl "];



for b = 1 : length(mean_auc_raw)

    std_bounds = [mean_auc_raw(b)-std_auc_raw(b) mean_auc_raw(b)+std_auc_raw(b)];

    % new
    if or(isnan(mean_auc_raw(b)),isnan(std_auc_raw(b)))
        OUTPUT_raw_effectsize(b) = nan;
    else
        % inhibition or activation?
        if max(std_bounds) < nanmean(sim_ratio_down_strong_rawscaled(:))-nanstd(sim_ratio_down_strong_rawscaled(:)) % strong inhibition?
            OUTPUT_raw_effectsize(b) = -2;
        elseif max(std_bounds) < nanmean(sim_ratio_down_weak_rawscaled(:))-nanstd(sim_ratio_down_weak_rawscaled(:)) % weak inhibition?
            OUTPUT_raw_effectsize(b) = -1;
        end
        if min(std_bounds) > nanmean(sim_ratio_up_strong_rawscaled(:))+nanstd(sim_ratio_up_strong_rawscaled(:)) % strong activation?
            OUTPUT_raw_effectsize(b) = 2;
        elseif min(std_bounds) > nanmean(sim_ratio_up_weak_rawscaled(:))+nanstd(sim_ratio_up_weak_rawscaled(:)) % weak activation?
            OUTPUT_raw_effectsize(b) = 1;
        end
        % more than 3x std?
        if or(max(std_bounds) < ctrl_mean_auc_raw-(3*ctrl_std_auc_raw), min(std_bounds) > ctrl_mean_auc_raw+(3*ctrl_std_auc_raw))
            OUTPUT_raw_3xSTD(b) = 1;
        end
    end
    %

    if any(b == idx_flagged) % if flagged, exclude
        fnan = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
        uistack(fnan,'bottom')
        txt=text(b,15,met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))),'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        uistack(txt,'top')

        OUTPUT_anyflag(b) = 1;
        OUTPUT_flags(b) = string(met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))));

    elseif any(b == idx_lowequflagged)
        errorbar(b,mean_auc_raw(b),std_auc_raw(b),'.','Color',[0 0.4470 0.7410]);
        scatter(b,reps_auc_raw(:,b),50,'.k')
        plot(b,reps_auc_raw(:,b),'.k')


        fnan2 = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
        uistack(fnan2,'bottom')
        txt=text(b,15,"equ",'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        uistack(txt,'top')

        OUTPUT_anyflag(b) = 1;
        OUTPUT_flags(b) = "equ-l";

    else % if not flagged, plot mean and std
        errorbar(b,mean_auc_raw(b),std_auc_raw(b),'.','Color',[0 0.4470 0.7410]);
        scatter(b,reps_auc_raw(:,b),50,'.k')
        plot(b,reps_auc_raw(:,b),'.k')

        if any(b == idx_saltflagged) % if one of salt controls is flagged, add light grey bar
            fnan2 = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
            uistack(fnan2,'bottom')
            %
            salt_idx = cell2mat(met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_calc")):find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_mgcl"))));
            txt=text(b,15,salt_abbr(find(salt_idx)),'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
            uistack(txt,'top')

            OUTPUT_anyflag(b) = 1;
            OUTPUT_flags(b) = salt_abbr(find(salt_idx));
        end
    end
end

%upweak
errorbar(b+2.2,nanmean(sim_ratio_up_weak_rawscaled(:)),nanstd(sim_ratio_up_weak_rawscaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+2.2,sim_ratio_up_weak_rawscaled(:)','.','Color',[0.6473 0.7456 0.4188]);
plot([0 length(mean_auc_raw)+4],[nanmean(sim_ratio_up_weak_rawscaled(:))+nanstd(sim_ratio_up_weak_rawscaled(:)) nanmean(sim_ratio_up_weak_rawscaled(:))+nanstd(sim_ratio_up_weak_rawscaled(:))],'-','Color',[0.6473 0.7456 0.4188],'LineWidth',.5)
%downweak
errorbar(b+1.8,nanmean(sim_ratio_down_weak_rawscaled(:)),nanstd(sim_ratio_down_weak_rawscaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+1.8,sim_ratio_down_weak_rawscaled(:)','.','Color',[0.3010 0.7450 0.9330]);
plot([0 length(mean_auc_raw)+4],[nanmean(sim_ratio_down_weak_rawscaled(:))-nanstd(sim_ratio_down_weak_rawscaled(:)) nanmean(sim_ratio_down_weak_rawscaled(:))-nanstd(sim_ratio_down_weak_rawscaled(:))],'-','Color',[0.3010 0.7450 0.9330],'LineWidth',.5)
%upstrong
errorbar(b+3.2,nanmean(sim_ratio_up_strong_rawscaled(:)),nanstd(sim_ratio_up_strong_rawscaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+3.2,sim_ratio_up_strong_rawscaled(:)','.','Color',[0.3500 0.6740 0.1880]);
plot([0 length(mean_auc_raw)+4],[nanmean(sim_ratio_up_strong_rawscaled(:))+nanstd(sim_ratio_up_strong_rawscaled(:)) nanmean(sim_ratio_up_strong_rawscaled(:))+nanstd(sim_ratio_up_strong_rawscaled(:))],'-','Color',[0.3500 0.6740 0.1880],'LineWidth',.5)
%downstrong
errorbar(b+2.8,nanmean(sim_ratio_down_strong_rawscaled(:)),nanstd(sim_ratio_down_strong_rawscaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+2.8,sim_ratio_down_strong_rawscaled(:)','.','Color',[0 0.4470 0.7410]);
plot([0 length(mean_auc_raw)+4],[nanmean(sim_ratio_down_strong_rawscaled(:))-nanstd(sim_ratio_down_strong_rawscaled(:)) nanmean(sim_ratio_down_strong_rawscaled(:))-nanstd(sim_ratio_down_strong_rawscaled(:))],'-','Color',[0 0.4470 0.7410],'LineWidth',.5)

% plot controls
errorbar(length(mean_auc_raw)+1,ctrl_mean_auc_raw,ctrl_std_auc_raw,'.','Color',[0 0.4470 0.7410])
scatter(length(mean_auc_raw)+1,ctrl_reps_auc_raw(:),50,'.k')
% plot lines
plot([0 length(mean_auc_raw)+4],[ctrl_mean_auc_raw ctrl_mean_auc_raw],'k','LineWidth',3)
% std?
plot([0 length(mean_auc_raw)+4],[ctrl_mean_auc_raw-(ctrl_std_auc_raw) ctrl_mean_auc_raw-(ctrl_std_auc_raw)],'--','Color',[.6 .6 .6],'LineWidth',.5)
plot([0 length(mean_auc_raw)+4],[ctrl_mean_auc_raw+(ctrl_std_auc_raw) ctrl_mean_auc_raw+(ctrl_std_auc_raw)],'--','Color',[.6 .6 .6],'LineWidth',.5)
% 3x std?
plot([0 length(mean_auc_raw)+4],[ctrl_mean_auc_raw-(3*ctrl_std_auc_raw) ctrl_mean_auc_raw-(3*ctrl_std_auc_raw)],'--','Color',[.6 .6 .6],'LineWidth',.5)
plot([0 length(mean_auc_raw)+4],[ctrl_mean_auc_raw+(3*ctrl_std_auc_raw) ctrl_mean_auc_raw+(3*ctrl_std_auc_raw)],'--','Color',[.6 .6 .6],'LineWidth',.5)

% figure settings
set(gca, 'XTick', [1:size(met_notes_flagged,1)])
effectorlist = string(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"met"))));
set(gca,'xticklabel',[effectorlist' "ctrls"])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4]);
xtickangle(50)
ylim([0 max(max(reps_auc_raw(:,idx_notflagged)))+5])
ylabel(["AUC"])

%% Bar plot 2 - scaled TECAN data

% SCALED DATA
datapoints = cell2mat(TECANdata.AUCs_scaled(:,1:n_replicates))';

idx_flagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"flag_total"))))==1);
idx_notflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"flag_total"))))==0);
idx_saltflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_total"))))==1);
idx_lowequflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_check_equ_l"))))==1);

% assays
reps_auc_scaled = datapoints; % datapoints(samples)
if n_replicates > 1
    mean_auc_scaled = nanmean(reps_auc_scaled); % mean(samples)
    std_auc_scaled = nanstd(reps_auc_scaled); % std(samples)
else
    mean_auc_scaled = reps_auc_scaled; % mean(samples)
    std_auc_scaled = zeros(1,length(reps_auc_scaled))
end


% controls
ctrls_abs_idx = ([12 : 12 : 96])
if sum([enzyme == 'Gnd' enzyme == 'Eno' enzyme == 'Zwf' enzyme == 'Icd']) == 1 % for some enzymes no 8th control, therefore use closest (=7)
    ctrls_abs_idx(end) = [];
end
ctrl_reps_auc_scaled = reps_auc_scaled(:,ctrls_abs_idx); % datapoints(ctrl)
ctrl_mean_auc_scaled = nanmean(ctrl_reps_auc_scaled(:)); % mean(ctrl)
ctrl_std_auc_scaled = nanstd(ctrl_reps_auc_scaled(:)); % std(ctrl)

figure(12346)
sgtitle(string(enzyme)+" AUCs scaled")
hold on
% plot assays
salt_abbr = ["calc ","2licl ","2kcl ","5nacl ","mgcl "];
for b = 1 : length(mean_auc_scaled)

    std_bounds = [mean_auc_scaled(b)-std_auc_scaled(b) mean_auc_scaled(b)+std_auc_scaled(b)];

    %new
    if or(isnan(mean_auc_scaled(b)),isnan(std_auc_scaled(b)))
        OUTPUT_scaled_effectsize(b) = nan;
    else
        % inhibition or activation
        if max(std_bounds) < nanmean(sim_ratio_down_strong_0_1_scaled(:))-nanstd(sim_ratio_down_strong_0_1_scaled(:)) % strong inhibition?
            OUTPUT_scaled_effectsize(b) = -2;
        elseif max(std_bounds) < nanmean(sim_ratio_down_weak_0_1_scaled(:))-nanstd(sim_ratio_down_weak_0_1_scaled(:)) % weak inhibition?
            OUTPUT_scaled_effectsize(b) = -1;
        end
        if min(std_bounds) > nanmean(sim_ratio_up_strong_0_1_scaled(:))+nanstd(sim_ratio_up_strong_0_1_scaled(:)) % strong activation?
            OUTPUT_scaled_effectsize(b) = 2;
        elseif min(std_bounds) > nanmean(sim_ratio_up_weak_0_1_scaled(:))+nanstd(sim_ratio_up_weak_0_1_scaled(:)) % weak activation?
            OUTPUT_scaled_effectsize(b) = 1;
        end
        % more than 3x std?
        if or(max(std_bounds) < ctrl_mean_auc_scaled-(3*ctrl_std_auc_scaled), min(std_bounds) > ctrl_mean_auc_scaled+(3*ctrl_std_auc_scaled))
            OUTPUT_scaled_3xSTD(b) = 1;
        end
    end
    %

    if any(b == idx_flagged) % if flagged, exclude
        fnan = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
        uistack(fnan,'bottom')
        txt=text(b,15,met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))),'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        uistack(txt,'top')
    elseif any(b == idx_lowequflagged) % if low-equ, plot but make grey
        errorbar(b,mean_auc_scaled(b),std_auc_scaled(b),'.','Color',[0 0.4470 0.7410]);
        scatter(b,reps_auc_scaled(:,b),50,'.k')
        plot(b,reps_auc_scaled(:,b),'.k')


        fnan2 = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
        uistack(fnan2,'bottom')
        txt=text(b,15,"equ",'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        uistack(txt,'top')

    else % if not flagged, plot mean and std
        errorbar(b,mean_auc_scaled(b),std_auc_scaled(b),'.','Color',[0 0.4470 0.7410]);
        scatter(b,reps_auc_scaled(:,b),50,'.k')
        plot(b,reps_auc_scaled(:,b),'.k')

        if any(b == idx_saltflagged) % if one of salt controls is flagged, add light grey bar
            fnan2 = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
            uistack(fnan2,'bottom')
            %
            salt_idx = cell2mat(met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_calc")):find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_mgcl"))));
            txt=text(b,15,salt_abbr(find(salt_idx)),'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
            uistack(txt,'top')
        end
    end
end

%upweak
errorbar(b+2.2,nanmean(sim_ratio_up_weak_0_1_scaled(:)),nanstd(sim_ratio_up_weak_0_1_scaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+2.2,sim_ratio_up_weak_0_1_scaled(:)','.','Color',[0.6473 0.7456 0.4188]);
plot([0 length(mean_auc_scaled)+4],[nanmean(sim_ratio_up_weak_0_1_scaled(:))+nanstd(sim_ratio_up_weak_0_1_scaled(:)) nanmean(sim_ratio_up_weak_0_1_scaled(:))+nanstd(sim_ratio_up_weak_0_1_scaled(:))],'-','Color',[0.6473 0.7456 0.4188],'LineWidth',.5)

%downweak
errorbar(b+1.8,nanmean(sim_ratio_down_weak_0_1_scaled(:)),nanstd(sim_ratio_down_weak_0_1_scaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+1.8,sim_ratio_down_weak_0_1_scaled(:)','.','Color',[0.3010 0.7450 0.9330]);
plot([0 length(mean_auc_scaled)+4],[nanmean(sim_ratio_down_weak_0_1_scaled(:))-nanstd(sim_ratio_down_weak_0_1_scaled(:)) nanmean(sim_ratio_down_weak_0_1_scaled(:))-nanstd(sim_ratio_down_weak_0_1_scaled(:))],'-','Color',[0.3010 0.7450 0.9330],'LineWidth',.5)

%upstrong
errorbar(b+3.2,nanmean(sim_ratio_up_strong_0_1_scaled(:)),nanstd(sim_ratio_up_strong_0_1_scaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+3.2,sim_ratio_up_strong_0_1_scaled(:)','.','Color',[0.3500 0.6740 0.1880]);
plot([0 length(mean_auc_scaled)+4],[nanmean(sim_ratio_up_strong_0_1_scaled(:))+nanstd(sim_ratio_up_strong_0_1_scaled(:)) nanmean(sim_ratio_up_strong_0_1_scaled(:))+nanstd(sim_ratio_up_strong_0_1_scaled(:))],'-','Color',[0.3500 0.6740 0.1880],'LineWidth',.5)

%downstrong
errorbar(b+2.8,nanmean(sim_ratio_down_strong_0_1_scaled(:)),nanstd(sim_ratio_down_strong_0_1_scaled(:)),'.','Color',[0.3010 0.7450 0.9330]);
plot(b+2.8,sim_ratio_down_strong_0_1_scaled(:)','.','Color',[0 0.4470 0.7410]);
plot([0 length(mean_auc_scaled)+4],[nanmean(sim_ratio_down_strong_0_1_scaled(:))-nanstd(sim_ratio_down_strong_0_1_scaled(:)) nanmean(sim_ratio_down_strong_0_1_scaled(:))-nanstd(sim_ratio_down_strong_0_1_scaled(:))],'-','Color',[0 0.4470 0.7410],'LineWidth',.5)

% plot controls
errorbar(length(mean_auc_scaled)+1,ctrl_mean_auc_scaled,ctrl_std_auc_scaled,'.','Color',[0 0.4470 0.7410])
scatter(length(mean_auc_scaled)+1,ctrl_reps_auc_scaled(:),50,'.k')
% plot lines
plot([0 length(mean_auc_scaled)+4],[ctrl_mean_auc_scaled ctrl_mean_auc_scaled],'k','LineWidth',3)
% std?
plot([0 length(mean_auc_scaled)+4],[ctrl_mean_auc_scaled-(ctrl_std_auc_scaled) ctrl_mean_auc_scaled-(ctrl_std_auc_scaled)],'--','Color',[.6 .6 .6],'LineWidth',.5)
plot([0 length(mean_auc_scaled)+4],[ctrl_mean_auc_scaled+(ctrl_std_auc_scaled) ctrl_mean_auc_scaled+(ctrl_std_auc_scaled)],'--','Color',[.6 .6 .6],'LineWidth',.5)
% 3x std?
plot([0 length(mean_auc_scaled)+4],[ctrl_mean_auc_scaled-(3*ctrl_std_auc_scaled) ctrl_mean_auc_scaled-(3*ctrl_std_auc_scaled)],'--','Color',[.6 .6 .6],'LineWidth',.5)
plot([0 length(mean_auc_scaled)+4],[ctrl_mean_auc_scaled+(3*ctrl_std_auc_scaled) ctrl_mean_auc_scaled+(3*ctrl_std_auc_scaled)],'--','Color',[.6 .6 .6],'LineWidth',.5)
% figure settings
set(gca, 'XTick', [1:size(met_notes_flagged,1)])
effectorlist = string(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"met"))));
set(gca,'xticklabel',[effectorlist' "ctrls"])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4]);
xtickangle(50)
ylim([0 max(max(reps_auc_scaled(:,idx_notflagged)))+5])
ylabel(["AUC"])





%% 2nd way of translating
AUCs_scaled_translated(1:size(AUCs_scaled,1),1:n_replicates) = nan; % empty table
for r = 1 : n_replicates % loop through replicates
    AUCratio_collection = cell2mat(AUCs_scaled(:,r)) ./ ctrl_mean_auc_scaled;
    for k = 1 : length(AUCratio_collection) % translate each AUC
        [c index] = min(abs(AUCratio2activity(:,1)-AUCratio_collection(k)));
        AUCs_scaled_translated(k,r) = AUCratio2activity(index,2);
    end
end
log2activities_scaled = log2(AUCs_scaled_translated)





%
OUTPUToverview = cell(30,length(OUTPUT_legend)+1);
OUTPUToverview(1,1) = {"legend"};
OUTPUToverview(1,2:end) = cellstr(OUTPUT_legend);
OUTPUToverview(2,1) = {"anyflag"};
OUTPUToverview(2,2:end) = num2cell(OUTPUT_anyflag);
OUTPUToverview(3,1) = {"flag"};
OUTPUToverview(3,2:end) = cellstr(OUTPUT_flags);
%
OUTPUToverview(5,1) = {"raw_effectsize"};
OUTPUToverview(5,2:end) = num2cell(OUTPUT_raw_effectsize);
OUTPUToverview(6,1) = {"raw_3xSTD"};
OUTPUToverview(6,2:end) = num2cell(OUTPUT_raw_3xSTD);
%
OUTPUToverview(13,1) = {"scaled_effectsize"};
OUTPUToverview(13,2:end) = num2cell(OUTPUT_scaled_effectsize);
OUTPUToverview(14,1) = {"scaled_3xSTD"};
OUTPUToverview(14,2:end) = num2cell(OUTPUT_scaled_3xSTD);
%
OUTPUToverview(23,1) = {"mean(log2activities_scaled_fitted)"};
OUTPUToverview(23,2:end) = num2cell(log2(nanmean(AUCs_scaled_fitted_translated')));
OUTPUToverview(24,1) = {"std(log2activities_scaled_fitted)"};
OUTPUToverview(24,2:end) = num2cell(nanstd(AUCs_scaled_fitted_translated')./(nanmean(AUCs_scaled_fitted_translated')*log(2)));
%
OUTPUToverview(25,1) = {"OLD___mean(log2activities_scaled_fitted)"};
OUTPUToverview(25,2:end) = num2cell(nanmean(log2activities_scaled_fitted'));
OUTPUToverview(26,1) = {"OLD___std(log2activities_scaled_fitted)"};
OUTPUToverview(26,2:end) = num2cell(std(log2activities_scaled_fitted'));
OUTPUToverview(27,1) = {"mean(log2activities_scaled)"};
OUTPUToverview(27,2:end) = num2cell(nanmean(log2activities_scaled'));
OUTPUToverview(28,1) = {"std(log2activities_scaled)"};
OUTPUToverview(28,2:end) = num2cell(std(log2activities_scaled'));
% no log
OUTPUToverview(29,1) = {"mean(activities_scaled_fitted)"};
OUTPUToverview(29,2:end) = num2cell(nanmean(AUCs_scaled_fitted_translated'));
OUTPUToverview(30,1) = {"std(activities_scaled_fitted)"};
OUTPUToverview(30,2:end) = num2cell(std(AUCs_scaled_fitted_translated'));


%% Bar plot 3 - scaled and fitted TECAN data

% SCALED DATA
datapoints = log2(AUCs_scaled_fitted_translated');

idx_flagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"flag_total"))))==1);
idx_notflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"flag_total"))))==0);
idx_saltflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_total"))))==1);
idx_lowequflagged = find(cell2mat(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_check_equ_l"))))==1);



% assays
reps_auc_raw = datapoints; % datapoints(samples)
if n_replicates == 1
    mean_auc_raw = reps_auc_raw; % mean(samples)
    std_auc_raw = zeros(size(reps_auc_raw)); % std(samples)
else
    mean_auc_raw = nanmean(reps_auc_raw); % mean(samples)
    std_auc_raw = nanstd(reps_auc_raw); % std(samples)
end

% controls
ctrls_abs_idx = ([12 : 12 : 96])
if sum([enzyme == 'Gnd' enzyme == 'Eno' enzyme == 'Zwf' enzyme == 'Icd']) == 1 % for some enzymes no 8th control, therefore use closest (=7)
    ctrls_abs_idx(end) = [];
end
ctrl_reps_auc_raw = reps_auc_raw(:,ctrls_abs_idx); % datapoints(ctrl)
ctrl_mean_auc_raw = nanmean(ctrl_reps_auc_raw(:)); % mean(ctrl)
ctrl_std_auc_raw = nanstd(ctrl_reps_auc_raw(:)); % std(ctrl)

figure(1234567)
sgtitle(string(enzyme)+" log2(translated AUC-fit)")
hold on
% plot assays
for b = 1 : length(mean_auc_raw)
    if any(b == idx_flagged) % if flagged, exclude
        fnan = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[min(min(reps_auc_raw(:,idx_notflagged))) min(min(reps_auc_raw(:,idx_notflagged))) 500 500 min(min(reps_auc_raw(:,idx_notflagged)))],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
        uistack(fnan,'bottom')
        txt=text(b,15,met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"flag_notes"))),'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        uistack(txt,'top')
    elseif any(b == idx_lowequflagged)
        errorbar(b,mean_auc_raw(b),std_auc_raw(b),'.','Color',[0 0.4470 0.7410]);
        scatter(b,reps_auc_raw(:,b),50,'.k')
        plot(b,reps_auc_raw(:,b),'.k')

        fnan2 = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[1 1 500 500 1],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
        uistack(fnan2,'bottom')
        txt=text(b,15,"equ",'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        uistack(txt,'top')
    else % if not flagged, plot mean and std
        errorbar(b,mean_auc_raw(b),std_auc_raw(b),'.','Color',[0 0.4470 0.7410]);
        scatter(b,reps_auc_raw(:,b),50,'.k')
        plot(b,reps_auc_raw(:,b),'.k')
        if any(b == idx_saltflagged) % if one of salt controls is flagged, add light grey bar
            fnan2 = fill([b-0.4 b+0.4 b+0.4 b-0.4 b-0.4],[min(min(reps_auc_raw(:,idx_notflagged))) min(min(reps_auc_raw(:,idx_notflagged))) 500 500 min(min(reps_auc_raw(:,idx_notflagged)))],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
            uistack(fnan2,'bottom')
            %
            salt_idx = cell2mat(met_notes_flagged(b+1,find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_calc")):find(strcmp(string(met_notes_flagged(1,:)),"MANUAL_salt_mgcl"))));
            txt=text(b,15,salt_abbr(find(salt_idx)),'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
            uistack(txt,'top')
        end
    end
end
% plot controls
errorbar(length(mean_auc_raw)+1,ctrl_mean_auc_raw,ctrl_std_auc_raw,'.','Color',[0 0.4470 0.7410])
scatter(length(mean_auc_raw)+1,ctrl_reps_auc_raw(:),50,'.k')
% plot lines
plot([0 length(mean_auc_raw)+2],[ctrl_mean_auc_raw ctrl_mean_auc_raw],'k','LineWidth',3)
% 3x std?
plot([0 length(mean_auc_raw)+2],[ctrl_mean_auc_raw-(3*ctrl_std_auc_raw) ctrl_mean_auc_raw-(3*ctrl_std_auc_raw)],'--','Color',[.4 .4 .4],'LineWidth',0.5)
plot([0 length(mean_auc_raw)+2],[ctrl_mean_auc_raw+(3*ctrl_std_auc_raw) ctrl_mean_auc_raw+(3*ctrl_std_auc_raw)],'--','Color',[.4 .4 .4],'LineWidth',0.5)
plot([0 length(mean_auc_raw)+2],[0.56565 0.56565],'-','Color',[.6 .6 .6],'LineWidth',0.5)
plot([0 length(mean_auc_raw)+2],[-0.56565 -0.56565],'-','Color',[.6 .6 .6],'LineWidth',0.5)
plot([0 length(mean_auc_raw)+2],[1.1313 1.1313],'-','Color',[.6 .6 .6],'LineWidth',0.5)
plot([0 length(mean_auc_raw)+2],[-1.1313 -1.1313],'-','Color',[.6 .6 .6],'LineWidth',0.5)

set(gca, 'XTick', [1:size(met_notes_flagged,1)])
effectorlist = string(met_notes_flagged(2:end,find(strcmp(string(met_notes_flagged(1,:)),"met"))));
set(gca,'xticklabel',[effectorlist' "ctrls"])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4]);
xtickangle(50)
ylim([min(min(reps_auc_raw(:,idx_notflagged)))-0.5 max(max(reps_auc_raw(:,idx_notflagged)))+0.5])
xlabel(["effectors"])
ylabel(["log2(translated scaled AUC)"])

%% add manual grey for new lower equilibrium and/or SALTS... overview, done.
figure(73642)
label = append("TECAN-",enzyme,"-ProgressCurves-raw.png")
saveas(gca,label)

figure(73644)
label = append("TECAN-",enzyme,"-ProgressCurves-scaled.png")
saveas(gca,label)

figure(73643)
label = append("TECAN-",enzyme,"-ProgressCurves-scaled&fitted.png")
saveas(gca,label)

figure(12345)
label = append("TECAN-",enzyme,"-ProgressCurves-scaled&fitted.png")
saveas(gca,label)

figure(123456)
label = append("TECAN-",enzyme,"-BARPLOT-log2translatedAUCfit.png")
saveas(gca,label)

figure(737400)
label = append("TECAN-",enzyme,"-AUC2BETAplot.png")
saveas(gca,label)

figure(987)
label = append("TECAN-",enzyme,"-simulations_inhibition_activation-scaled&fitted.png")
saveas(gca,label)

save(append("TECAN-",enzyme,"-OUTPUToverview"),'OUTPUToverview')