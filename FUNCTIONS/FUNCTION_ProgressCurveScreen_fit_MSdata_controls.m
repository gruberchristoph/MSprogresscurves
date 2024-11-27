% fit and plot processed fiadata
function fiadata_ctrls = ProgressCurveScreen2_fit_fiadata_controls_ticnorm_v9(fiadata, ioncount_threshold, plot_on, TICfilter, TIC_normalization,R2_cutoff,fit_function,save_all_results)

opts = statset('nlinfit');
tic_excl_threshold = 0.5; % exclude if sample TIC more than 50% off from median TIC of timeseries

%prepare plots

for minim = 1
    if plot_on == 1
        % create figures
        figure(1)
        sgtitle(fiadata.enzymename+" ALL controls: LambertW")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(1),"northwest")
            ylim([-0.1 1.1])
            x1 = 0:128
            x2= zeros(1,129)
        end

        figure(2)
        sgtitle(fiadata.enzymename+" set 1 - 4 controls: LambertW")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(2),"north")
            ylim([-0.1 1.1])
        end

        figure(3)
        sgtitle(fiadata.enzymename+" set 5 - 8 controls: LambertW")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(3),"northeast")
            ylim([-0.1 1.1])
        end

        figure(4)
        sgtitle(fiadata.enzymename+" ALL controls: TWO-SUBSTRATE LambertW")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(4),"west")
            ylim([-0.1 1.1])
        end

        figure(5)
        sgtitle(fiadata.enzymename+" set 1 - 4 controls: TWO-SUBSTRATE LambertW")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(5),"center")
            ylim([-0.1 1.1])
        end

        figure(6)
        sgtitle(fiadata.enzymename+" set 5 - 8 controls: TWO-SUBSTRATE LambertW")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(6),"east")
            ylim([-0.1 1.1])
        end

        figure(7)
        sgtitle(fiadata.enzymename+" ALL controls: SIGMOID")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(7),"southwest")
            ylim([-0.1 1.1])
        end

        figure(8)
        sgtitle(fiadata.enzymename+" set 1 - 4 controls: SIGMOID")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(8),"south")
            ylim([-0.1 1.1])
        end

        figure(9)
        sgtitle(fiadata.enzymename+" set 5 - 8 controls: SIGMOID")
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
            plot(0:128,zeros(1,129),'--','Color',[.6 .6 .6])
            plot(0:128,zeros(1,129)+1,'--','Color',[.6 .6 .6])
            xlabel("time [min]")
            ylabel("scaled ion count")
            movegui(figure(9),"southeast")
            ylim([-0.1 1.1])
        end
    end
end

% get indices of all controls and breaking point
idx_controls = find(strcmp(fiadata.S1(:,4),"eff12"));
idx_break = min(find(strcmp(fiadata.S1(:,2),"set5")));
i_break = min(find(idx_controls>idx_break));

% output tables lambertW
fiadata_ctrls.lambertW_AUC(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_r2(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt1(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt2(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt3(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt1_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt1_u95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt2_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt2_u95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt3_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.lambertW_xopt3_u95(1:4,1:size(idx_controls)) = nan;
% output tables TWOSUBSTRATElambertW
fiadata_ctrls.TWOSUBSTRATElambertW_AUC(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_r2(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt1(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt2(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt3(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt1_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt1_u95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt2_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt2_u95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt3_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.TWOSUBSTRATElambertW_xopt3_u95(1:4,1:size(idx_controls)) = nan;
% output tables sigmoid
fiadata_ctrls.SIGMOID_AUC(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_r2(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt1(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt2(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt3(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt1_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt1_u95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt2_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt2_u95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt3_l95(1:4,1:size(idx_controls)) = nan;
fiadata_ctrls.SIGMOID_xopt3_u95(1:4,1:size(idx_controls)) = nan;

% OVERVIEW
ctrl_parameters_OVERVIEW(1:9,1:7) = {nan};
ctrl_parameters_OVERVIEW(2:9,1) = [{"1a"} {"1b"} {"2a"} {"2b"} {"3a"} {"3b"} {"4a"} {"4b"}];
ctrl_parameters_OVERVIEW(1,2:11) = [{"ctrls mean(r2)"} {"ctrls min(r2)"} {"maxstdY"} {"X"} {"ctrls 2S-lW CV (%)"} {"ctrls lW CV (%)"} {""} {"eff n_excluded"} {"eff CV mean"} {"eff CV median"}];

ctrl_parameters_OVERVIEWlump(1:5,1:6) = {nan};
ctrl_parameters_OVERVIEWlump(2:5,1) = [{"S1"} {"S2"} {"P1"} {"P2"}];
ctrl_parameters_OVERVIEWlump(1,1:6) = [{""} {"meanAUC"} {"medianAUC"} {"stdAUC"} {"meanR2"} {"minR2"}];

% loop through reactants
for reactant_nr = 1 : 4
    clear datatable & TICtable
    if reactant_nr == 1
        datatable = fiadata.S1
        TICtable = fiadata.S1_TICs
    elseif reactant_nr == 2
        datatable = fiadata.S2
        TICtable = fiadata.S2_TICs
    elseif reactant_nr == 3
        datatable = fiadata.P1
        TICtable = fiadata.P1_TICs
    elseif reactant_nr == 4
        datatable = fiadata.P2
        TICtable = fiadata.P2_TICs
    end

    % get TICmedian for TICnormalization
    TICall = cell2mat(TICtable(2:end,7:end))
    TICmedian = nanmedian(TICall(:))

    % check number of missing Y values
    Y_missing = sum(sum(ismissing(cell2mat(datatable(2:end,7:end)))));
    if Y_missing < size(cell2mat(datatable(2:end,7:end)),1) * size(cell2mat(datatable(2:end,7:end)),2)

        % output tables, save scaled y data
        lambertW_y_REscaled_collection(1:size(idx_controls),1:8) = nan;
        TWOSUBlambertW_y_REscaled_collection(1:size(idx_controls),1:8) = nan;
        SIGMOID_y_REscaled_collection(1:size(idx_controls),1:8) = nan;

        % fit and plot
        for i = 1:size(idx_controls,1)

            % get timepoints
            x = cellfun(@str2num,fiadata.time(idx_controls(i),:)) / 60;
            x_range = min(x) : 0.1 : max(x);

            % fitting MM-based lambert W (1S) function to data:
            clear xopt & y_fit & fit_y_delta & y_fit_scaled & y_REscaled & y_raw & y_scaled
            try
                % define function for substrate and product
                fun_lambertW_subs = @(x,t) x(3) * lambertw(0, x(1) * exp(x(1) - x(2)*t)); % a = x(1); b = x(2); k_m = x(3);
                fun_lambertW_prod = @(x,t) x(3) * (x(1) - lambertw(0, x(1) * exp(x(1) -x(2) * t ))); % a = x(1); b = x(2); k_m = x(3);

                % get y values (=measured ion count)
                y_raw = cell2mat(datatable(idx_controls(i),7:end)); % measured ion count

                if reactant_nr <= 2
                    y_raw(y_raw<ioncount_threshold) = nan; % remove ions below threshold
                elseif reactant_nr > 2
                    excl_idx = y_raw<ioncount_threshold;
                    excl_idx(1) = 0; %  leave t0 if lower in case of products
                    y_raw(excl_idx) = nan; % remove ions below threshold
                end

                % exclude outliers based on TIC
                tic_raw = cell2mat(TICtable(idx_controls(i),7:end)); % getting TICs corresponding to ion counts
                if TICfilter == "CUT50"
                	tic_raw(abs(tic_raw ./ nanmedian(tic_raw)-1) > tic_excl_threshold) = nan; % removing bad injections
                elseif TICfilter == "MAD"
                    [~,TF] = rmoutliers(tic_raw)
                    tic_raw(TF) = nan;
                end
                y_raw(isnan(tic_raw)) = nan;

                % normalize to TIC and scale from 0 to 1
                if TIC_normalization == 1
                    tic_norm = (y_raw ./ tic_raw);
                    y_scaled = rescale(tic_norm);  % scale between 0 and 1 before fitting (this allows us to choose 1 set of initial parameters for the fitting function across all experiments)
                else
                    y_scaled = rescale(y_raw);
                end

                % arrange data
                tbl = table(x', y_scaled');
                % fit
                if reactant_nr <= 2 % =substrates
                    mdl = fitnlm(tbl, fun_lambertW_subs, [1 0.01 0.1],'Options',opts); % fitting procedure
                    xopt = mdl.Coefficients{:, 'Estimate'};
                    y_fit = xopt(3) * lambertw(0, xopt(1) * exp(xopt(1) - xopt(2)*x_range)); % = fun_lambertW_subs
                elseif reactant_nr >= 3 % =products
                    mdl = fitnlm(tbl, fun_lambertW_prod, [1 0.01 0.1],'Options',opts); % fitting procedure
                    xopt = mdl.Coefficients{:, 'Estimate'};
                    y_fit = xopt(3) * (xopt(1) - lambertw(0, xopt(1) * exp(xopt(1) - xopt(2) * x_range))); % = fun_lambertW_prod
                end

                % rescaling to make S-infinity of fitted data = 0
                fit_y_delta =  max(y_fit(:)) - min(y_fit(:));
                y_fit_scaled = (y_fit - min(y_fit(:))) / fit_y_delta;
                y_REscaled =  (y_scaled - min(y_fit(:))) / fit_y_delta;
                if reactant_nr <= 2 % =substrates
                    fiadata_ctrls.lambertW_AUC(reactant_nr,i) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range, y_fit_scaled);
                elseif reactant_nr >= 3 % =products
                    fiadata_ctrls.lambertW_AUC(reactant_nr,i) = trapz(x_range, y_fit_scaled);
                end

                % extract parameters from 2nd fit
                fiadata_ctrls.lambertW_r2(reactant_nr,i) = mdl.Rsquared.Adjusted;
                fiadata_ctrls.lambertW_xopt1(reactant_nr,i) = xopt(1);
                fiadata_ctrls.lambertW_xopt2(reactant_nr,i) = xopt(2);
                fiadata_ctrls.lambertW_xopt3(reactant_nr,i) = xopt(3);
                ci = coefCI(mdl);
                fiadata_ctrls.lambertW_xopt1_l95(reactant_nr,i) = ci(1,1);
                fiadata_ctrls.lambertW_xopt1_u95(reactant_nr,i) = ci(1,2);
                fiadata_ctrls.lambertW_xopt2_l95(reactant_nr,i) = ci(2,1);
                fiadata_ctrls.lambertW_xopt2_u95(reactant_nr,i) = ci(2,2);
                fiadata_ctrls.lambertW_xopt3_l95(reactant_nr,i) = ci(3,1);
                fiadata_ctrls.lambertW_xopt3_u95(reactant_nr,i) = ci(3,2);

                % plot raw data point and fits
                if plot_on == 1
                    figure(1)
                    subplot(2,2,reactant_nr)
                    plot(x,y_REscaled,"k.")
                    hold on
                    plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                    hold on
                    if idx_controls(i) < idx_break
                        figure(2)
                        subplot(2,2,reactant_nr)
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x,y_REscaled,"k.")
                        else
                            plot(x,y_REscaled,"k.")
                        end
                        hold on
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                        else
                            plot(x_range, y_fit_scaled,'r');
                        end
                        hold on
                    elseif idx_controls(i)  >= idx_break
                        figure(3)
                        subplot(2,2,reactant_nr)
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x,y_REscaled,"k.")
                        else
                            plot(x,y_REscaled,"k.")
                        end
                        hold on
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                        else
                            plot(x_range, y_fit_scaled,'r');
                        end
                        hold on
                    end
                end

                lambertW_y_REscaled_collection(i,:) = y_REscaled;

            end

            % fitting TWO_SUBSTRATE lambert W function to data:
            clear xopt & y_fit & fit_y_delta & y_fit_scaled & y_REscaled & y_raw & y_scaled
            try
                % define functions
                fun_TWOSUBSTRATElambertW_subs = @(x,t) x(3)/2 * (x(1)/x(3) - x(3)/x(1) - x(2) * t + sqrt((x(1)/x(3) - x(3)/x(1) - x(2) * t).^2 + 4)); % s0 = x(1); b = x(2); k_m = x(3);
                fun_TWOSUBSTRATElambertW_prod = @(x,t) 1 - x(3)/2 * (x(1)/x(3) - x(3)/x(1) - x(2) * t + sqrt((x(1)/x(3) - x(3)/x(1) - x(2) * t).^2 + 4)); % s0 = x(1); b = x(2); k_m = x(3);

                % get y values and rescale
                y_raw = cell2mat(datatable(idx_controls(i),7:end)); % measured ion count

                if reactant_nr <= 2
                    y_raw(y_raw<ioncount_threshold) = nan; % remove ions below threshold
                elseif reactant_nr > 2
                    excl_idx = y_raw<ioncount_threshold;
                    excl_idx(1) = 0; %  leave t0 if lower in case of products
                    y_raw(excl_idx) = nan; % remove ions below threshold
                end
                % exclude outliers based on TIC
                tic_raw = cell2mat(TICtable(idx_controls(i),7:end)); % getting TICs corresponding to ion counts
                if TICfilter == "CUT50"
                	tic_raw(abs(tic_raw ./ nanmedian(tic_raw)-1) > tic_excl_threshold) = nan; % make bad injections = nan
                elseif TICfilter == "MAD"
                    [~,TF] = rmoutliers(tic_raw)
                    tic_raw(TF) = nan;
                end
                y_raw(isnan(tic_raw)) = nan;

                % normalize based on TIC
                if TIC_normalization == 1
                    tic_norm = (y_raw ./ tic_raw); % same as in fiaminer: ratio of ion count to TIC multiplied by median of all TICs
                    y_scaled = rescale(tic_norm);
                else
                    y_scaled = rescale(y_raw); % scale between 0 and 1 before fitting (this allows us to choose 1 set of initial parameters for the fitting function across all experiments)
                end

                % arrange data
                tbl = table(x', y_scaled');

                % fit
                if reactant_nr <= 2 % =substrates
                    mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_subs, [0.1 0.1 0.1],'Options',opts); % fitting procedure
                    xopt = mdl.Coefficients{:, 'Estimate'};
                    y_fit = xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_subs
                elseif reactant_nr >= 3 % =products
                    mdl = fitnlm(tbl, fun_TWOSUBSTRATElambertW_prod, [0.9 0.1 0.1],'Options',opts); % fitting procedure
                    xopt = mdl.Coefficients{:, 'Estimate'};
                    y_fit = 1 - xopt(3)/2 * (xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range + sqrt((xopt(1)/xopt(3) - xopt(3)/xopt(1) - xopt(2) * x_range).^2 + 4)); % = fun_TWOSUBSTRATElambertW_prod
                end

                % rescaling to make S-infinity of fitted data = 1
                fit_y_delta =  max(y_fit(:)) - min(y_fit(:));
                y_fit_scaled = (y_fit - min(y_fit(:))) / fit_y_delta;
                y_REscaled =  (y_scaled - min(y_fit(:))) / fit_y_delta;
                if reactant_nr <= 2 % =substrates
                    fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,i) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range, y_fit_scaled);
                elseif reactant_nr >= 3 % =products
                    fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,i) = trapz(x_range, y_fit_scaled);
                end

                % extract parameters from 2nd fit
                fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,i) = mdl.Rsquared.Adjusted;
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt1(reactant_nr,i) = xopt(1);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt2(reactant_nr,i) = xopt(2);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt3(reactant_nr,i) = xopt(3);
                ci = coefCI(mdl);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt1_l95(reactant_nr,i) = ci(1,1);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt1_u95(reactant_nr,i) = ci(1,2);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt2_l95(reactant_nr,i) = ci(2,1);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt2_u95(reactant_nr,i) = ci(2,2);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt3_l95(reactant_nr,i) = ci(3,1);
                fiadata_ctrls.TWOSUBSTRATElambertW_xopt3_u95(reactant_nr,i) = ci(3,2);

                % plot raw data point and fits
                if plot_on == 1
                    figure(4)
                    subplot(2,2,reactant_nr)
                    plot(x,y_REscaled,"k.")
                    hold on
                    plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                    hold on
                    if idx_controls(i) < idx_break
                        figure(5)
                        subplot(2,2,reactant_nr)
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x,y_REscaled,"k.")
                        else
                            plot(x,y_REscaled,"k.")
                        end
                        hold on
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                        else
                            plot(x_range, y_fit_scaled,'r');
                        end
                        hold on
                    elseif idx_controls(i)  >= idx_break
                        figure(6)
                        subplot(2,2,reactant_nr)
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x,y_REscaled,"k.")
                        else
                            plot(x,y_REscaled,"k.")
                        end
                        hold on
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                        else
                            plot(x_range, y_fit_scaled,'r');
                        end
                        hold on
                    end
                end
                TWOSUBlambertW_y_REscaled_collection(i,:) = y_REscaled;
            end

            % fitting SIGMOID function to data:
            clear xopt & y_fit & fit_y_delta & y_fit_scaled & y_REscaled & y_raw & y_scaled
            try
                % define functions
                fun_SIGMOID_subs = @(x,t) 1 - x(1) * tanh(abs(x(2)/x(1))*(t-x(3))) ; % SIGMOID: x(1) height of curve; x(2) heat-parametersc; x(3) center point
                fun_SIGMOID_prod = @(x,t) x(1) * tanh(abs(x(2)/x(1))*(t-x(3))) ; % SIGMOID: x(1) height of curve; x(2) heat-parametersc; x(3) center point

                % get y values and rescale
                y_raw = cell2mat(datatable(idx_controls(i),7:end)); % measured ion count

                if reactant_nr <= 2
                    y_raw(y_raw<ioncount_threshold) = nan; % remove ions below threshold
                elseif reactant_nr > 2
                    excl_idx = y_raw<ioncount_threshold;
                    excl_idx(1) = 0; %  leave t0 if lower in case of products
                    y_raw(excl_idx) = nan; % remove ions below threshold
                end

                % exclude outliers based on TIC
                tic_raw = cell2mat(TICtable(idx_controls(i),7:end)); % getting TICs corresponding to ion counts
                if TICfilter == "CUT50"
                	tic_raw(abs(tic_raw ./ nanmedian(tic_raw)-1) > tic_excl_threshold) = nan; % make bad injections = nan
                elseif TICfilter == "MAD"
                    [~,TF] = rmoutliers(tic_raw)
                    tic_raw(TF) = nan;
                end
                y_raw(isnan(tic_raw)) = nan;

                % normalize to TIC
                if TIC_normalization == 1
                    tic_norm = (y_raw ./ tic_raw);
                    y_scaled = rescale(tic_norm); % scale between 0 and 1 before fitting (this allows us to choose 1 set of initial parameters for the fitting function across all experiments)
                else
                    y_scaled = rescale(y_raw);
                end

                % arrange data
                tbl = table(x', y_scaled');

                % fit
                if reactant_nr <= 2 % =substrates
                    mdl = fitnlm(tbl, fun_SIGMOID_subs, [1 -0.1 0.1],'Options',opts); % fitting procedure
                    xopt = mdl.Coefficients{:, 'Estimate'};
                    y_fit = 1 - xopt(1) * tanh(abs(xopt(2)/xopt(1))*(x_range-xopt(3)));
                elseif reactant_nr >= 3 % =products
                    mdl = fitnlm(tbl, fun_SIGMOID_prod, [1 -0.1 0.1],'Options',opts); % fitting procedure
                    xopt = mdl.Coefficients{:, 'Estimate'};
                    y_fit = xopt(1) * tanh(abs(xopt(2)/xopt(1))*(x_range-xopt(3)));
                end

                % rescaling to make S-infinity of fitted data = 1
                fit_y_delta =  max(y_fit(:)) - min(y_fit(:));
                y_fit_scaled = (y_fit - min(y_fit(:))) / fit_y_delta;
                y_REscaled =  (y_scaled - min(y_fit(:))) / fit_y_delta;
                if reactant_nr <= 2 % =substrates
                    fiadata_ctrls.SIGMOID_AUC(reactant_nr,i) = trapz(x_range,ones(1,size(x_range,2))) - trapz(x_range, y_fit_scaled);
                elseif reactant_nr >= 3 % =products
                    fiadata_ctrls.SIGMOID_AUC(reactant_nr,i) = trapz(x_range, y_fit_scaled);
                end

                % extract parameters from 2nd fit
                fiadata_ctrls.SIGMOID_r2(reactant_nr,i) = mdl.Rsquared.Adjusted;
                fiadata_ctrls.SIGMOID_xopt1(reactant_nr,i) = xopt(1);
                fiadata_ctrls.SIGMOID_xopt2(reactant_nr,i) = xopt(2);
                fiadata_ctrls.SIGMOID_xopt3(reactant_nr,i) = xopt(3);
                ci = coefCI(mdl);
                fiadata_ctrls.SIGMOID_xopt1_l95(reactant_nr,i) = ci(1,1);
                fiadata_ctrls.SIGMOID_xopt1_u95(reactant_nr,i) = ci(1,2);
                fiadata_ctrls.SIGMOID_xopt2_l95(reactant_nr,i) = ci(2,1);
                fiadata_ctrls.SIGMOID_xopt2_u95(reactant_nr,i) = ci(2,2);
                fiadata_ctrls.SIGMOID_xopt3_l95(reactant_nr,i) = ci(3,1);
                fiadata_ctrls.SIGMOID_xopt3_u95(reactant_nr,i) = ci(3,2);

                % plot raw data point and fits
                if plot_on == 1
                    figure(7)
                    subplot(2,2,reactant_nr)
                    plot(x,y_REscaled,"k.")
                    hold on
                    plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                    hold on
                    if idx_controls(i) < idx_break
                        figure(8)
                        subplot(2,2,reactant_nr)
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x,y_REscaled,"k.")
                        else
                            plot(x,y_REscaled,"k.")
                        end
                        hold on
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                        else
                            plot(x_range, y_fit_scaled,'r');
                        end
                        hold on
                    elseif idx_controls(i)  >= idx_break
                        figure(9)
                        subplot(2,2,reactant_nr)
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x,y_REscaled,"k.")
                        else
                            plot(x,y_REscaled,"k.")
                        end
                        hold on
                        if mdl.Rsquared.Adjusted > R2_cutoff
                            plot(x_range, y_fit_scaled,'Color',[.4 .4 .4]);
                        else
                            plot(x_range, y_fit_scaled,'r');
                        end
                        hold on
                    end
                end
                SIGMOID_y_REscaled_collection(i,:) = y_REscaled;
            end


            % store parameters and add text to plots
            for minim = 1
                % when last control of reactant is reached
                if i == size(idx_controls,1)
                    fiadata_ctrls.ctrl_AUCmean_lambertW_all(reactant_nr,1) = nanmean(fiadata_ctrls.lambertW_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCmedian_lambertW_all(reactant_nr,1) =  nanmedian(fiadata_ctrls.lambertW_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCstd_lambertW_all(reactant_nr,1) = nanstd(fiadata_ctrls.lambertW_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2mean_lambertW_all(reactant_nr,1) = nanmean(fiadata_ctrls.lambertW_r2(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2median_lambertW_all(reactant_nr,1) = nanmedian(fiadata_ctrls.lambertW_r2(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2min_lambertW_all(reactant_nr,1) = min(fiadata_ctrls.lambertW_r2(reactant_nr,:))
                    if plot_on == 1
                        figure(1) % lambertW_all
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_lambertW_all(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_lambertW_all(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_lambertW_all(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_lambertW_all(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_all(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_lambertW_1to4(reactant_nr,1) = nanmean(fiadata_ctrls.lambertW_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCmedian_lambertW_1to4(reactant_nr,1) =  nanmedian(fiadata_ctrls.lambertW_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCstd_lambertW_1to4(reactant_nr,1) = nanstd(fiadata_ctrls.lambertW_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2mean_lambertW_1to4(reactant_nr,1) = nanmean(fiadata_ctrls.lambertW_r2(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2median_lambertW_1to4(reactant_nr,1) = nanmedian(fiadata_ctrls.lambertW_r2(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2min_lambertW_1to4(reactant_nr,1) = min(fiadata_ctrls.lambertW_r2(reactant_nr,1:i_break-1))
                    if plot_on == 1
                        figure(2) % lambertW_1to4
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_lambertW_1to4(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_lambertW_1to4(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_lambertW_1to4(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_lambertW_1to4(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_1to4(reactant_nr,1),2)))
                        text(70,0.3,"max(std(Y)): " + string(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_1to4(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_lambertW_5to8(reactant_nr,1) = nanmean(fiadata_ctrls.lambertW_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCmedian_lambertW_5to8(reactant_nr,1) =  nanmedian(fiadata_ctrls.lambertW_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCstd_lambertW_5to8(reactant_nr,1) = nanstd(fiadata_ctrls.lambertW_AUC(reactant_nr,i_break:end));
                    fiadata_ctrls.ctrl_AUCr2mean_lambertW_5to8(reactant_nr,1) = nanmean(fiadata_ctrls.lambertW_r2(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2median_lambertW_5to8(reactant_nr,1) = nanmedian(fiadata_ctrls.lambertW_r2(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2min_lambertW_5to8(reactant_nr,1) = min(fiadata_ctrls.lambertW_r2(reactant_nr,i_break:end))
                    if plot_on == 1
                        figure(3) % lambertW_5to8
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_lambertW_5to8(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_lambertW_5to8(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_lambertW_5to8(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_lambertW_5to8(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_5to8(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_all(reactant_nr,1) = nanmean(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_all(reactant_nr,1) =  nanmedian(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_all(reactant_nr,1) = nanstd(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_all(reactant_nr,1) = nanmean(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2median_TWOSUBSTRATElambertW_all(reactant_nr,1) = nanmedian(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_all(reactant_nr,1) = min(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,:))
                    if plot_on == 1
                        figure(4) % TWOSUBSTRATElambertW_all
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_all(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_all(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_all(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_all(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_all(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_1to4(reactant_nr,1) = nanmean(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_1to4(reactant_nr,1) =  nanmedian(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_1to4(reactant_nr,1) = nanstd(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_1to4(reactant_nr,1) = nanmean(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2median_TWOSUBSTRATElambertW_1to4(reactant_nr,1) = nanmedian(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_1to4(reactant_nr,1) = min(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,1:i_break-1))
                    if plot_on == 1
                        figure(5) % TWOSUBSTRATElambertW_1to4
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_5to8(reactant_nr,1) = nanmean(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_5to8(reactant_nr,1) =  nanmedian(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_5to8(reactant_nr,1) = nanstd(fiadata_ctrls.TWOSUBSTRATElambertW_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_5to8(reactant_nr,1) = nanmean(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2median_TWOSUBSTRATElambertW_5to8(reactant_nr,1) = nanmedian(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_5to8(reactant_nr,1) = min(fiadata_ctrls.TWOSUBSTRATElambertW_r2(reactant_nr,i_break:end))
                    if plot_on == 1
                        figure(6) % TWOSUBSTRATElambertW_5to8
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_SIGMOID_all(reactant_nr,1) = nanmean(fiadata_ctrls.SIGMOID_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCmedian_SIGMOID_all(reactant_nr,1) =  nanmedian(fiadata_ctrls.SIGMOID_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCstd_SIGMOID_all(reactant_nr,1) = nanstd(fiadata_ctrls.SIGMOID_AUC(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2mean_SIGMOID_all(reactant_nr,1) = nanmean(fiadata_ctrls.SIGMOID_r2(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2median_SIGMOID_all(reactant_nr,1) = nanmedian(fiadata_ctrls.SIGMOID_r2(reactant_nr,:))
                    fiadata_ctrls.ctrl_AUCr2min_SIGMOID_all(reactant_nr,1) = min(fiadata_ctrls.SIGMOID_r2(reactant_nr,:))
                    if plot_on == 1
                        figure(7) % SIGMOID_all
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_SIGMOID_all(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_SIGMOID_all(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_SIGMOID_all(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_SIGMOID_all(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_SIGMOID_all(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_SIGMOID_1to4(reactant_nr,1) = nanmean(fiadata_ctrls.SIGMOID_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCmedian_SIGMOID_1to4(reactant_nr,1) =  nanmedian(fiadata_ctrls.SIGMOID_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCstd_SIGMOID_1to4(reactant_nr,1) = nanstd(fiadata_ctrls.SIGMOID_AUC(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2mean_SIGMOID_1to4(reactant_nr,1) = nanmean(fiadata_ctrls.SIGMOID_r2(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2median_SIGMOID_1to4(reactant_nr,1) = nanmedian(fiadata_ctrls.SIGMOID_r2(reactant_nr,1:i_break-1))
                    fiadata_ctrls.ctrl_AUCr2min_SIGMOID_1to4(reactant_nr,1) = min(fiadata_ctrls.SIGMOID_r2(reactant_nr,1:i_break-1))
                    if plot_on == 1
                        figure(8) % SIGMOID_1to4
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_SIGMOID_1to4(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_SIGMOID_1to4(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_SIGMOID_1to4(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_SIGMOID_1to4(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_SIGMOID_1to4(reactant_nr,1),2)))
                    end
                    fiadata_ctrls.ctrl_AUCmean_SIGMOID_5to8(reactant_nr,1) = nanmean(fiadata_ctrls.SIGMOID_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCmedian_SIGMOID_5to8(reactant_nr,1) =  nanmedian(fiadata_ctrls.SIGMOID_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCstd_SIGMOID_5to8(reactant_nr,1) = nanstd(fiadata_ctrls.SIGMOID_AUC(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2mean_SIGMOID_5to8(reactant_nr,1) = nanmean(fiadata_ctrls.SIGMOID_r2(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2median_SIGMOID_5to8(reactant_nr,1) = nanmedian(fiadata_ctrls.SIGMOID_r2(reactant_nr,i_break:end))
                    fiadata_ctrls.ctrl_AUCr2min_SIGMOID_5to8(reactant_nr,1) = min(fiadata_ctrls.SIGMOID_r2(reactant_nr,i_break:end))
                    if plot_on == 1
                        figure(9) % SIGMOID_5to8
                        subplot(2,2,reactant_nr)
                        text(70,0.8,"mean: " + string(round(fiadata_ctrls.ctrl_AUCmean_SIGMOID_5to8(reactant_nr,1),2)))
                        text(70,0.7,"median: " + string(round(fiadata_ctrls.ctrl_AUCmedian_SIGMOID_5to8(reactant_nr,1),2)))
                        text(70,0.6,"std: " + string(round(fiadata_ctrls.ctrl_AUCstd_SIGMOID_5to8(reactant_nr,1),2)))
                        text(70,0.5,"mean(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2mean_SIGMOID_5to8(reactant_nr,1),2)))
                        text(70,0.4,"min(r2): " + string(round(fiadata_ctrls.ctrl_AUCr2min_SIGMOID_5to8(reactant_nr,1),2)))
                    end
                end
            end
        end

    end

    % max standarddeviation in Y
    lambertW_maxstdY(reactant_nr) = max(nanstd(lambertW_y_REscaled_collection))*100; % max std deviation
    TWOSUBlambertW_maxstdY(reactant_nr) = max(nanstd(TWOSUBlambertW_y_REscaled_collection))*100; % max std deviation
    SIGMOID_maxstdY(reactant_nr) = max(nanstd(SIGMOID_y_REscaled_collection))*100; % max std deviation
    lambertW_maxstdY_1to4(reactant_nr,:) = max(nanstd(lambertW_y_REscaled_collection(idx_controls < idx_break,:)))*100;
    lambertW_maxstdY_5to8(reactant_nr,:) = max(nanstd(lambertW_y_REscaled_collection(idx_controls >= idx_break,:)))*100;
    TWOSUBlambertW_maxstdY_1to4(reactant_nr,:) = max(nanstd(TWOSUBlambertW_y_REscaled_collection(idx_controls < idx_break,:)))*100;
    TWOSUBlambertW_maxstdY_5to8(reactant_nr,:) = max(nanstd(TWOSUBlambertW_y_REscaled_collection(idx_controls >= idx_break,:)))*100;
    SIGMOID_maxstdY_1to4(reactant_nr,:) = max(nanstd(SIGMOID_y_REscaled_collection(idx_controls < idx_break,:)))*100;
    SIGMOID_maxstdY_5to8(reactant_nr,:) = max(nanstd(SIGMOID_y_REscaled_collection(idx_controls >= idx_break,:)))*100;

    try % because problematic if less than 4 reactants
        if fit_function == "2SlambertW"
            % r2 parameters into overview table
            ctrl_parameters_OVERVIEW(1+reactant_nr*2-1,2) = num2cell(round(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2,2) = num2cell(round(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2-1,3) = num2cell(round(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_1to4(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2,3) = num2cell(round(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_5to8(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2-1,4) =  num2cell(TWOSUBlambertW_maxstdY_1to4(reactant_nr));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2,4) = num2cell(TWOSUBlambertW_maxstdY_5to8(reactant_nr));
            %
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_all,1),2) = num2cell(round(fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_all,1),3) = num2cell(round(fiadata_ctrls.ctrl_AUCmedian_TWOSUBSTRATElambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_all,1),4) = num2cell(round(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_all,1),5) = num2cell(round(fiadata_ctrls.ctrl_AUCr2mean_TWOSUBSTRATElambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_all,1),6) = num2cell(round(fiadata_ctrls.ctrl_AUCr2min_TWOSUBSTRATElambertW_all,2));
        elseif fit_function == "lambertW"
            ctrl_parameters_OVERVIEW(1+reactant_nr*2-1,2) = num2cell(round(fiadata_ctrls.ctrl_AUCr2mean_lambertW_1to4(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2,2) = num2cell(round(fiadata_ctrls.ctrl_AUCr2mean_lambertW_5to8(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2-1,3) = num2cell(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_1to4(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2,3) = num2cell(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_5to8(reactant_nr,1),2));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2-1,4) =  num2cell(lambertW_maxstdY_1to4(reactant_nr));
            ctrl_parameters_OVERVIEW(1+reactant_nr*2,4) = num2cell(lambertW_maxstdY_5to8(reactant_nr));
            %
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCmean_lambertW_all,1),2) = num2cell(round(fiadata_ctrls.ctrl_AUCmean_lambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCmedian_lambertW_all,1),3) = num2cell(round(fiadata_ctrls.ctrl_AUCmedian_lambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCstd_lambertW_all,1),4) = num2cell(round(fiadata_ctrls.ctrl_AUCstd_lambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCr2mean_lambertW_all,1),5) = num2cell(round(fiadata_ctrls.ctrl_AUCr2mean_lambertW_all,2));
            ctrl_parameters_OVERVIEWlump(2:1+size(fiadata_ctrls.ctrl_AUCr2min_lambertW_all,1),6) = num2cell(round(fiadata_ctrls.ctrl_AUCr2min_lambertW_all,2));
        end
    end
end

% add CVs to overview table
try %%% problematic if only 3 entries
    ctrl_parameters_OVERVIEW([2 4 6 8],6) = num2cell(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_1to4 ./ fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_1to4 * 100);
    ctrl_parameters_OVERVIEW([3 5 7 9],6) = num2cell(fiadata_ctrls.ctrl_AUCstd_TWOSUBSTRATElambertW_5to8 ./ fiadata_ctrls.ctrl_AUCmean_TWOSUBSTRATElambertW_5to8 * 100);
    ctrl_parameters_OVERVIEW([2 4 6 8],7) = num2cell(fiadata_ctrls.ctrl_AUCstd_lambertW_1to4 ./ fiadata_ctrls.ctrl_AUCmean_lambertW_1to4 * 100);
    ctrl_parameters_OVERVIEW([3 5 7 9],7) = num2cell(fiadata_ctrls.ctrl_AUCstd_lambertW_5to8 ./ fiadata_ctrls.ctrl_AUCmean_lambertW_5to8 * 100);
end

fiadata_ctrls.ctrl_parameters_OVERVIEW = ctrl_parameters_OVERVIEW;
fiadata_ctrls.ctrl_parameters_OVERVIEWlump = ctrl_parameters_OVERVIEWlump;

% maxstdY
fiadata_ctrls.lambertW_maxstdY=lambertW_maxstdY;
fiadata_ctrls.TWOSUBlambertW_maxstdY=TWOSUBlambertW_maxstdY;
fiadata_ctrls.SIGMOID_maxstdY=SIGMOID_maxstdY;
fiadata_ctrls.SIGMOID_maxstdY_1to4=SIGMOID_maxstdY_5to8;
fiadata_ctrls.SIGMOID_maxstdY_5to8=SIGMOID_maxstdY_5to8;
fiadata_ctrls.TWOSUBlambertW_maxstdY_1to4=TWOSUBlambertW_maxstdY_5to8;
fiadata_ctrls.TWOSUBlambertW_maxstdY_5to8=TWOSUBlambertW_maxstdY_5to8;
fiadata_ctrls.lambertW_maxstdY_1to4=lambertW_maxstdY_5to8;
fiadata_ctrls.lambertW_maxstdY_5to8=lambertW_maxstdY_5to8;

if save_all_results == 1
    if fit_function == "2SlambertW"
        figure(5)
        label = string(fiadata.enzymename) + "_CONTROLS_set1to4.png"
        saveas(gca,label)
        figure(6)
        label = string(fiadata.enzymename) + "_CONTROLS_set5to8.png"
        saveas(gca,label)
    elseif fit_function == "lambertW"
        figure(8)
        label = string(fiadata.enzymename) + "_CONTROLS_set1to4.png"
        saveas(gca,label)
        figure(9)
        label = string(fiadata.enzymename) + "_CONTROLS_set5to8.png"
        saveas(gca,label)
    end
end
end