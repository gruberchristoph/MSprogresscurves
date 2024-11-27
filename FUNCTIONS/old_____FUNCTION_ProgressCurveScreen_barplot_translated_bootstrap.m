function fiadata_scores = FUNCTION_ProgressCurveScreen_barplot_translated_bootstrap(analysis_path,fiadata_effs, enzymelabel, ions_excl,r2_threshold,AUC2ac,fit_function, bypass_AUCs, split_on,excl_if_r2_lower_than_08_of_max,transl_mode,lnACT_on,plot_only_relevant,save_all_results)
%% Options
set_each_ctrl_to_ONE = 1; % divide everything by mean of ctrl
lumped_CTRL_signif = 1; % significance based on fit across all controls
excl_r2 = 1; % exclude AUCs if R2 < 0.7

%% Create output table
fiadata_scores = struct;
fiadata_scores.ion_excl = ions_excl;
fiadata_scores.r2_threshold = r2_threshold;
fiadata_scores.AUC2ac = AUC2ac;
analysis_legend = fiadata_effs.analysis_legend;

%% Get information on salt effects (Mg2+,Mn2+,etc) from overview
file = {analysis_path+"CODETABLES.xlsx"}; %path to screen overview
[overview_info_A,overview_info_B] = xlsread(file{1},"overview","A1 : BO34"); % manually entered effects of salt controls
Eidx = find(strcmp(string(overview_info_B(:,find(strcmp(overview_info_B(1,:),"label")))),enzymelabel));
salt_calc = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"calc")));
salt_2licl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"2licl")));
salt_2kcl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"2kcl")));
salt_5nacl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"5nacl")));
salt_mgcl = overview_info_A(Eidx-1,find(strcmp(string(overview_info_B(1,:)),"mgcl")));
 
%% Flag if tested effector contains salt that was found to affect enzymatic activity
% flag salts
met_abbr = analysis_legend(:,5);
flags_SALT(1:length(met_abbr)+1,1:7) = {nan};
flags_SALT(1,:) = [{"met"} {"MANUAL_salt_calc"} {"MANUAL_salt_2licl"} {"MANUAL_salt_2kcl"} {"MANUAL_salt_5nacl"} {"MANUAL_salt_mgcl"} {"MANUAL_salt_total"}]; % MANUAL_salt = ["licl" "nacl" etc] from excel % MANUAL_else "comment"
flags_SALT(2:end,1) = met_abbr;
%
empty_flags = zeros(size(flags_SALT,1)-1,1); % empty table
salt_calc_flagged = empty_flags;
if or(salt_calc == -1, salt_calc == 1)
    salt_calc_flagged(find(matches(string(flags_SALT(2:end,1)),["calc","panto"]))) = salt_calc_flagged(find(matches(string(flags_SALT(:,1)),["calc","panto"]))) + salt_calc;
end
flags_SALT(2:end,find(strcmp(string(flags_SALT(1,:)),"MANUAL_salt_calc"))) = num2cell(salt_calc_flagged);
%
salt_2licl_flagged = empty_flags;
if or(salt_2licl == -1, salt_2licl == 1)
    salt_2licl_flagged(find(matches(string(flags_SALT(2:end,1)),["kdpg","acp","2licl","carb-p","accoa","glyc3p"]))) = salt_2licl_flagged(find(matches(string(flags_SALT(:,1)),["kdpg","acp","2licl","carb-p","accoa","glyc3p"]))) + salt_2licl;
end
flags_SALT(2:end,find(strcmp(string(flags_SALT(1,:)),"MANUAL_salt_2licl"))) = num2cell(salt_2licl_flagged);
%
salt_2kcl_flagged = empty_flags;
if or(salt_2kcl == -1, salt_2kcl == 1)
    salt_2kcl_flagged(find(matches(string(flags_SALT(2:end,1)),["2kcl","acp","icit","gal1p","f1p"]))) = salt_2kcl_flagged(find(matches(string(flags_SALT(:,1)),["2kcl","acp","icit","gal1p","f1p"]))) + salt_2kcl;
end
flags_SALT(2:end,find(strcmp(string(flags_SALT(1,:)),"MANUAL_salt_2kcl"))) = num2cell(salt_2kcl_flagged);
%
salt_5nacl_flagged = empty_flags;
if or(salt_5nacl == -1, salt_5nacl == 1)
    salt_5nacl_flagged(find(matches(string(flags_SALT(2:end,1)),["5nacl","bpg","nadph"]))) = salt_5nacl_flagged(find(matches(string(flags_SALT(:,1)),["5nacl","bpg","nadph"]))) + salt_5nacl;
end
flags_SALT(2:end,find(strcmp(string(flags_SALT(1,:)),"MANUAL_salt_5nacl"))) = num2cell(salt_5nacl_flagged);
%
salt_mgcl_flagged = empty_flags;
if or(salt_mgcl == -1, salt_mgcl == 1)
    salt_mgcl_flagged(find(matches(string(flags_SALT(2:end,1)),["mgcl","dhap"]))) = salt_mgcl_flagged(find(matches(string(flags_SALT(:,1)),["mgcl","dhap"]))) + salt_mgcl;
end
flags_SALT(2:end,find(strcmp(string(flags_SALT(1,:)),"MANUAL_salt_mgcl"))) = num2cell(salt_mgcl_flagged);
% ALL manuals
manual_salt_total_flagged = abs(salt_calc_flagged)+abs(salt_2licl_flagged)+abs(salt_2kcl_flagged)+abs(salt_5nacl_flagged)+abs(salt_mgcl_flagged)
flags_SALT(2:end,find(strcmp(string(flags_SALT(1,:)),"MANUAL_salt_total"))) = num2cell(manual_salt_total_flagged);
flags_SALT_any = cell2mat(flags_SALT(2:end,end))

% flag others
flags_SPprecip_any = or(string(analysis_legend(:,6)) == "PRECIP",or(string(analysis_legend(:,6)) == "SUBSTRATE",string(analysis_legend(:,6)) == "PRODUCT"))

% flag overview
anyflag = flags_SALT_any + flags_SPprecip_any;
anyflag(anyflag>=1) = 1;



%% INPUTS
% SCORING
data = fiadata_effs.AUCs; % 200 bootstrapped AUCs
r2 = fiadata_effs.R2s; % R2s corresponding to AUCs
impfit = fiadata_effs.impossible_fits; % number of failed fits

AUC2act_all = fiadata_effs.AUC2activity_all; % AUC translated to change in beta
if transl_mode == "lumped" % one fit on all ctrls lumped and then simulated
AUC2act_1st = fiadata_effs.AUC2activity_1st; % AUC translated to change in beta
AUC2act_2nd = fiadata_effs.AUC2activity_2nd; % AUC translated to change in beta
elseif transl_mode == "sep" % simulated for each control and then taken mean
AUC2act_1st = fiadata_effs.MEANsepctrlAUC2activity_1st; % AUC translated to change in beta
AUC2act_2nd = fiadata_effs.MEANsepctrlAUC2activity_2nd; % AUC translated to change in beta  
end

transl_act = fiadata_effs.act; % empty table for activities
transl_act_excl = fiadata_effs.act; % empty table for activities

% prepare tables
r2_exc(1:size(data,1),1:4) = {nan};
transOOR_exc(1:size(data,1),1:4) = {nan};
NUMEL_data = zeros(size(data,1),4);

% find split position
break_idx = find(strcmp(data(:,5),"control4"));
ctrl_break_idx = find(strcmp(data(:,5),"ctrl_all"));

% exclude bad ions from analysis (by removing from analysis legend completely and later using this to decide not to plot/fit)
for exc = 1 : size(ions_excl,2)
    if ions_excl(exc) == 0
        analysis_legend(:,exc) = {0};
    end
    if ions_excl(exc) == 11 % only first set
        analysis_legend(break_idx+1:ctrl_break_idx-1,exc) = {0};
    end
    if ions_excl(exc) == 12 % only second set
        analysis_legend(1:break_idx,exc) = {0};
    end
end
labels = string(analysis_legend(:,5))


fiadata_scores = struct;
scoringoverviewtemplate(1:ctrl_break_idx-1,1:4) = {nan};
scoringoverviewtemplate(1:ctrl_break_idx-1,5) = analysis_legend(1:ctrl_break_idx-1,5);

% SIGNIFICANCE TESTS: WELCH
% on each reactant separately: AUCs
 fiadata_scores.welch_AUCs_each_n4 = scoringoverviewtemplate;
 fiadata_scores.welch_AUCs_each_n200 = scoringoverviewtemplate;
 fiadata_scores.welch_AUCs_each_n8x4 = scoringoverviewtemplate;
% on each reactant separately: translated activities
 fiadata_scores.welch_translACT_each_n4 = scoringoverviewtemplate;
 fiadata_scores.welch_translACT_each_n200 = scoringoverviewtemplate;
 fiadata_scores.welch_translACT_each_n8x4 = scoringoverviewtemplate;
% LUMPED: AUCs
 fiadata_scores.welch_AUCs_lumped_n4 = scoringoverviewtemplate;
 fiadata_scores.welch_AUCs_lumped_n200 = scoringoverviewtemplate;
 fiadata_scores.welch_AUCs_lumped_n8x4 = scoringoverviewtemplate;
% LUMPED: translated activities
 fiadata_scores.welch_translACT_lumped_n4(1:99,1) = {nan};
 fiadata_scores.welch_translACT_lumped_n4(1:99,2) =  analysis_legend(:,5);
 fiadata_scores.welch_translACT_lumped_n200 = fiadata_scores.welch_translACT_lumped_n4;
 fiadata_scores.welch_translACT_lumped_n8x4 = fiadata_scores.welch_translACT_lumped_n4;
 
%% change all entries in data to NaN that (i) are excluded in ANALYSIS_LEGEND or (ii) correspond to low r2 (based on predefined cutoff)
%% also exclude all assays where less than 2 datapoints left
% then save exluded datapoints in new table
for j = 1 : 4 % for each rectant
    for i = 1 : size(data,1) % for each effector
        clear r2s & data_save & r2_exc_idx & r2_exc_save
        
        r2s_save = cell2mat(r2(i,j)); % get 200 R2s
        data_save = cell2mat(data(i,j)); % get 200 AUCs
        
        %% TRANSLATE AUC TO ACTIVITY AND SAVE
        if any(data_save)
            clear AUC2act_table & AUC2ac_save & CTRL_table & AUCs_signif
            AUC2ac_save(1:size(data_save)) = nan;
            % get correct translation table
            if split_on == 1
                if i <= break_idx % set 1-4
                    AUC2act_table = AUC2act_1st;
                elseif i > break_idx & i < ctrl_break_idx % set 5-8
                    AUC2act_table = AUC2act_2nd;
                elseif i >= ctrl_break_idx % CONTROL COLLECTIONS TRANSLATED WITH OWN TABLE
                    if strcmp(string(analysis_legend(i,5)),"ctrl_all")
                        AUC2act_table = AUC2act_all;
                    elseif strcmp(string(analysis_legend(i,5)),"ctrl_1st")
                        AUC2act_table = AUC2act_1st;
                    elseif strcmp(string(analysis_legend(i,5)),"ctrl_2nd")
                        AUC2act_table = AUC2act_2nd;
                    end
                end
                
            % get correct control AUCs for significance test
              if i <= break_idx % set 1-4
                  if lumped_CTRL_signif == 1
                    CTRL_table = cell2mat(data(find(strcmp(data(:,5),'ctrl_1st')),j)); % get set1-4 ctrl table
                    if excl_r2 == 1
                    CTRL_table(find(cell2mat(r2(find(strcmp(r2(:,5),'ctrl_1st')),j))<0.7)) = nan;
                    end
                  elseif lumped_CTRL_signif == 0
                      %%
                  end
                elseif i > break_idx & i < ctrl_break_idx % set 5-8
                  if lumped_CTRL_signif == 1
                    CTRL_table = cell2mat(data(find(strcmp(data(:,5),'ctrl_2nd')),j)); % get set1-4 ctrl table
                    if excl_r2 == 1
                    CTRL_table(find(cell2mat(r2(find(strcmp(r2(:,5),'ctrl_2nd')),j))<0.7)) = nan;
                    end
                  elseif lumped_CTRL_signif == 0
                      %%
                  end
               elseif i >= ctrl_break_idx % CONTROL COLLECTIONS TRANSLATED WITH OWN TABLE
                    if strcmp(string(analysis_legend(i,5)),"ctrl_all")
                        CTRL_table = cell2mat(data(find(strcmp(data(:,5),'ctrl_all')),j)); % get set1-4 ctrl table
                        if excl_r2 == 1
                        CTRL_table(find(cell2mat(r2(find(strcmp(r2(:,5),'ctrl_all')),j))<0.7)) = nan;
                        end
                    elseif strcmp(string(analysis_legend(i,5)),"ctrl_1st")
                        CTRL_table = cell2mat(data(find(strcmp(data(:,5),'ctrl_1st')),j)); % get set1-4 ctrl table
                        if excl_r2 == 1
                        CTRL_table(find(cell2mat(r2(find(strcmp(r2(:,5),'ctrl_1st')),j))<0.7)) = nan;
                        end
                    elseif strcmp(string(analysis_legend(i,5)),"ctrl_2nd")
                        CTRL_table = cell2mat(data(find(strcmp(data(:,5),'ctrl_2nd')),j)); % get set1-4 ctrl table
                        if excl_r2 == 1
                        CTRL_table(find(cell2mat(r2(find(strcmp(r2(:,5),'ctrl_2nd')),j))<0.7)) = nan;
                        end
                    end
              end 
              
              
 %% SIGNIFICANCE TESTS ON AUC, each reactant separately (RAW: before translation)
            % get raw AUCs for significance test
            AUCs_signif = data_save;
            AUCs_signif(r2s_save<0.7) = nan; %%%% remove bad fits before sign test?

            try
                clear n_free & WELCHdf_n4  & WELCHscore_n4 & WELCHp_n4 & WELCHdf_n200  & WELCHscore_n200 & WELCHp_n200 & WELCHdf_n8x4  & WELCHscore_n8x4 & WELCHp_n8x4
                
                n_free = 4; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
                WELCHscore_n4 = ( nanmean(CTRL_table) - nanmean(AUCs_signif) ) / ( sqrt( ( nanvar(CTRL_table) / n_free ) + ( nanvar(AUCs_signif) / n_free ) ) );
                WELCHdf_n4 = ( ( nanvar(CTRL_table) / n_free ) + ( nanvar(AUCs_signif) / n_free ) )^2 / ( ( ( nanvar(CTRL_table) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(AUCs_signif) / n_free )^2 / (n_free - 1) ) );
                WELCHp_n4 = 2*tcdf(-abs(WELCHscore_n4),WELCHdf_n4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
                fiadata_scores.welch_AUCs_each_n4(i,j) = {WELCHp_n4}
                
                n_free = 200; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
                WELCHscore_n200 = ( nanmean(CTRL_table) - nanmean(AUCs_signif) ) / ( sqrt( ( nanvar(CTRL_table) / n_free ) + ( nanvar(AUCs_signif) / n_free ) ) );
                WELCHdf_n200 = ( ( nanvar(CTRL_table) / n_free ) + ( nanvar(AUCs_signif) / n_free ) )^2 / ( ( ( nanvar(CTRL_table) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(AUCs_signif) / n_free )^2 / (n_free - 1) ) );
                WELCHp_n200 = 2*tcdf(-abs( WELCHscore_n200),WELCHdf_n200); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
                fiadata_scores.welch_AUCs_each_n200(i,j) = {WELCHp_n200};
                
                n_free = 8*4;
                WELCHscore_n8x4 = ( nanmean(CTRL_table) - nanmean(AUCs_signif) ) / ( sqrt( ( nanvar(CTRL_table) / n_free ) + ( nanvar(AUCs_signif) / n_free ) ) );
                WELCHdf_n8x4 = ( ( nanvar(CTRL_table) / n_free ) + ( nanvar(AUCs_signif) / n_free ) )^2 / ( ( ( nanvar(CTRL_table) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(AUCs_signif) / n_free )^2 / (n_free - 1) ) );
                WELCHp_n8x4 = 2*tcdf(-abs(WELCHscore_n8x4),WELCHdf_n8x4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
                fiadata_scores.welch_AUCs_each_n8x4(i,j) = {WELCHp_n8x4};
            end

            elseif split_on == 0
                AUC2act_table = AUC2act_all;
                "WRITE SPLIT CODE"
            end
            
            
            % translate each AUC
            % first for control AUCs % r2<0.7 already removed
            for ctrlauc_i = 1 : length(CTRL_table)
                if ~isnan(CTRL_table(ctrlauc_i))
                    % find closest AUC in AUC-translation table
                    [c ctrlindex] = min(abs(AUC2act_table(:,j)-CTRL_table(ctrlauc_i))); % 2do: is this working correctly in case of duplicates? SHORT ANSWER: there arent any duplicates(!?)
                    ctrlACT(ctrlauc_i) = AUC2act_table(ctrlindex,5);
                else
                    ctrlACT(ctrlauc_i) = nan;
                end
            end
            
            clear c & index & ACT_save
            for auc_i = 1 : length(data_save)
                if ~isnan(data_save(auc_i))
                    % find closest AUC in AUC-translation table
                    [c index] = min(abs(AUC2act_table(:,j)-data_save(auc_i))); % 2do: is this working correctly in case of duplicates? SHORT ANSWER: there arent any duplicates(!?)
                        ACT_save(auc_i) = AUC2act_table(index,5);
                else
                    ACT_save(auc_i) = nan;
                end
            end
        else
            ACT_save = nan;
            data_save = nan;
            ctrlACT = nan;
        end
        data_collection(i,j) = {data_save};
      
           if set_each_ctrl_to_ONE == 1
               ACT_save = ACT_save / nanmean(ctrlACT);
           end
           if lnACT_on == 0
               %
           elseif lnACT_on == 1
               ACT_save = log(ACT_save);
               ctrlACT = log(ctrlACT);
           elseif lnACT_on == 2
               ACT_save = log2(ACT_save);
               ctrlACT = log2(ctrlACT);
           end
       
        %% Remove bad fits (r2<0.7)
        if excl_r2 == 0
            transl_act(i,j) = {ACT_save};
        elseif sum(~isnan(ACT_save))>1 & length(ACT_save(r2s_save>=0.7)) > length(ACT_save(r2s_save<0.7)) % set zero if more than half are missing
            transl_act(i,j) = {ACT_save(r2s_save>=0.7)}; % all good fits
            transl_act_excl(i,j) = {ACT_save(r2s_save<0.7)}; % saves low r2 and nans into new table
        else
            transl_act(i,j) = {nan};
            transl_act_excl(i,j) = {nan};
        end
        
        
        %% Welch test of translated and log-transformed data here
        transl_act_signif = ACT_save;%(r2s_save>=0.7)
        transl_act_signif(r2s_save<0.7)=nan; % remove all with r2<0.7
        
        clear n_free & WELCHdf_ACTn4  & WELCHscore_ACTn4 & WELCHp_ACTn4 & WELCHdf_ACTn200  & WELCHscore_ACTn200 & WELCHp_ACTn200 & WELCHdf_ACTn8x4  & WELCHscore_ACTn8x4 & WELCHp_ACTn8x4 
        n_free = 4;
        WELCHscore_ACTn4 = ( nanmean(ctrlACT) - nanmean(transl_act_signif) ) / ( sqrt( ( nanvar(ctrlACT) / n_free ) + ( nanvar(transl_act_signif) / n_free ) ) );
        WELCHdf_ACTn4 = ( ( nanvar(ctrlACT) / n_free ) + ( nanvar(transl_act_signif) / n_free ) )^2 / ( ( ( nanvar(ctrlACT) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(transl_act_signif) / n_free )^2 / (n_free - 1) ) );
        WELCHp_ACTn4 = 2*tcdf(-abs(WELCHscore_ACTn4),WELCHdf_ACTn4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
        fiadata_scores.welch_translACT_each_n4(i,j) = {WELCHp_ACTn4}
        
        n_free = 200;
        WELCHscore_ACTn200 = ( nanmean(ctrlACT) - nanmean(transl_act_signif) ) / ( sqrt( ( nanvar(ctrlACT) / n_free ) + ( nanvar(transl_act_signif) / n_free ) ) );
        WELCHdf_ACTn200 = ( ( nanvar(ctrlACT) / n_free ) + ( nanvar(transl_act_signif) / n_free ) )^2 / ( ( ( nanvar(ctrlACT) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(transl_act_signif) / n_free )^2 / (n_free - 1) ) );
        WELCHp_ACTn200 = 2*tcdf(-abs(WELCHscore_ACTn200),WELCHdf_ACTn200); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
        fiadata_scores.welch_translACT_each_n200(i,j) = {WELCHp_ACTn200}    
        
        n_free = 8*4; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
        WELCHscore_ACTn8x4 = ( nanmean(ctrlACT) - nanmean(transl_act_signif) ) / ( sqrt( ( nanvar(ctrlACT) / n_free ) + ( nanvar(transl_act_signif) / n_free ) ) );
        WELCHdf_ACTn8x4 = ( ( nanvar(ctrlACT) / n_free ) + ( nanvar(transl_act_signif) / n_free ) )^2 / ( ( ( nanvar(ctrlACT) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(transl_act_signif) / n_free )^2 / (n_free - 1) ) );
        WELCHp_ACTn8x4 = 2*tcdf(-abs(WELCHscore_ACTn8x4),WELCHdf_ACTn8x4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
        fiadata_scores.welch_translACT_each_n8x4(i,j) = {WELCHp_ACTn8x4}    
    end
end

style = "std";
fiadata_scores.transl_act = transl_act;
fiadata_scores.transl_act_excl = transl_act_excl;

%% Translated
fig3 = figure(130)
hold on
max_ry = [];
min_ry = [];
for p = 1 : size(transl_act,1)    
%% randomize x for plotting datapoints
    r1 = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(transl_act(p,1))),1);
    r5 = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(transl_act(p,2))),1);
    r3 = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(transl_act(p,3))),1);
    r4 = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(transl_act(p,4))),1);
    r1e = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(transl_act_excl(p,1))),1);
    r2e = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(transl_act_excl(p,2))),1);
    r3e = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(transl_act_excl(p,3))),1);
    r4e = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(transl_act_excl(p,4))),1);
    
%% for barplot
    ry=[];
    rx=[];

    % get y of effector // and y of respective control for significance test
    if cell2mat(analysis_legend(p,1)) == 1
        ry = [ry cell2mat(transl_act(p,1))];
    end
    if cell2mat(analysis_legend(p,2)) == 1
        ry = [ry cell2mat(transl_act(p,2))];
    end
    if cell2mat(analysis_legend(p,3)) == 1
        ry = [ry cell2mat(transl_act(p,3))];
    end
    if cell2mat(analysis_legend(p,4)) == 1
        ry = [ry cell2mat(transl_act(p,4))];
    end
    ry(isnan(ry))= []; % remove nans

    % save max and min ry
    max_ry = [max_ry max(ry)];
    min_ry = [min_ry min(ry)];
    
    % create x
    rx = (p-0.4) + ((p+0.4)-(p-0.4)).*rand(length(ry),1);
    
    if isempty(ry)
        % do nothing
    elseif style == "box"
        boxchart(ones(1,length(ry))*p,ry,"BoxWidth",0.9,'BoxFaceColor',[.4 .4 .4],'MarkerColor',[.4 .4 .4])
    elseif style == "std"
        bar(p,nanmean(ry),'FaceColor',[.8 .8 .8])
        hold on
        errorbar(p,nanmean(ry),nanstd(ry),nanstd(ry),'.k')
        hold on
        
%% DATAPOINTS
    % plot y-values
    % plot excluded because of r2
    if excl_r2 == 1
        plot(r1e,cell2mat(transl_act_excl(p,1)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.4)
        plot(r2e,cell2mat(transl_act_excl(p,2)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.4)
        plot(r3e,cell2mat(transl_act_excl(p,3)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.4)
        plot(r4e,cell2mat(transl_act_excl(p,4)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.4)
    end

    if cell2mat(analysis_legend(p,1)) == 1
        plot(r1,cell2mat(transl_act(p,1)),'.','Color',[0.8500 0.3250 0.0980],'MarkerSize',.4) % orange
    else
%% 
    end
    if cell2mat(analysis_legend(p,2)) == 1
    plot(r5,cell2mat(transl_act(p,2)),'.','Color',[0.9290 0.6940 0.1250],'MarkerSize',.4) % yellow
    else
%%
    end
    if cell2mat(analysis_legend(p,3)) == 1
    plot(r3,cell2mat(transl_act(p,3)),'.','Color',[0.4660 0.6740 0.1880],'MarkerSize',.4) % green
    else
%%
    end
    if cell2mat(analysis_legend(p,4)) == 1
    plot(r4,cell2mat(transl_act(p,4)),'.','Color',[0 0.4470 0.7410],'MarkerSize',.4) % blue
    else
%%
    end
    hold on
    end
   end

% add flags
flaglabellist(1:size(transl_act,1)) = "";
for p = 1 : size(transl_act,1)
    if flags_SPprecip_any(p) == 1
        fnan = fill([p-0.45 p+0.45 p+0.45 p-0.45 p-0.45],[min(min_ry)-.5 min(min_ry)-.5 max(max_ry)+.5 max(max_ry)+.5 min(min_ry)-.5],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
        uistack(fnan,'bottom')
        if string(analysis_legend(p,6)) == "PRODUCT"
            flaglabel = "P";
            flaglabellist(p) = "P";
        elseif string(analysis_legend(p,6)) == "SUBSTRATE"
            flaglabel = "S";
            flaglabellist(p) = "S";
        elseif string(analysis_legend(p,6)) == "PRECIP"
            flaglabel = "precip";
            flaglabellist(p) = "precip";
        end
        txt=text(p,min(min_ry)+.2,flaglabel,'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black') 
    elseif flags_SALT_any(p) == 1
        fnan = fill([p-0.45 p+0.45 p+0.45 p-0.45 p-0.45],[min(min_ry)-.5 min(min_ry)-.5 max(max_ry)+.5 max(max_ry)+.5 min(min_ry)-.5],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
        uistack(fnan,'bottom')
        txt=text(p,min(min_ry)+.2,"salt",'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
        flaglabellist(p) = "salt";
    end
end


% Label
plot([0 size(transl_act,1)+1], [0 0],'-k','LineWidth',3)
plot([0 size(transl_act,1)+1], [log2(2) log2(2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(5) log2(5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(10) log2(10)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.5) log2(0.5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.2) log2(0.2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.1) log2(0.1)],'--','LineWidth',.2,'Color',[.6 .6 .6])
set(gca,'xtick',[1:size(transl_act,1)],'xticklabel',string(analysis_legend(:,5)))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4])
ylabel("log2(activity-change)")
ylim([min(min_ry)-.5 max(max_ry)+.5])
sgtitle(enzymelabel+" beta - lumped reactants")



% Plot figures
figure(260)
hold on
for p = 1 : size(transl_act,1)
    % create random numbers around p as x-values
    r1 = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(transl_act(p,1))),1);
    r5 = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(transl_act(p,2))),1);
    r3 = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(transl_act(p,3))),1);
    r4 = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(transl_act(p,4))),1);
    r1e = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(transl_act_excl(p,1))),1);
    r2e = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(transl_act_excl(p,2))),1);
    r3e = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(transl_act_excl(p,3))),1);
    r4e = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(transl_act_excl(p,4))),1);
    
    % Plot y-values
    % Plot excluded because of r2 in grey
    if excl_r2 == 1
        plot(r1e,cell2mat(transl_act_excl(p,1)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
        plot(r2e,cell2mat(transl_act_excl(p,2)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
        plot(r3e,cell2mat(transl_act_excl(p,3)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
        plot(r4e,cell2mat(transl_act_excl(p,4)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end

    if cell2mat(analysis_legend(p,1)) == 1
        plot(r1,cell2mat(transl_act(p,1)),'.','Color',[0.8500 0.3250 0.0980]) % orange
    else
        plot(r1,cell2mat(transl_act(p,1)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,2)) == 1
    plot(r5,cell2mat(transl_act(p,2)),'.','Color',[0.9290 0.6940 0.1250]) % yellow
    else
    plot(r5,cell2mat(transl_act(p,2)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,3)) == 1
    plot(r3,cell2mat(transl_act(p,3)),'.','Color',[0.4660 0.6740 0.1880]) % green
    else
        plot(r3,cell2mat(transl_act(p,3)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,4)) == 1
    plot(r4,cell2mat(transl_act(p,4)),'.','Color',[0 0.4470 0.7410]) % blue
    else
        plot(r4,cell2mat(transl_act(p,4)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    hold on
end
plot([0 size(transl_act,1)+1], [0 0],'-k','LineWidth',3)
plot([0 size(transl_act,1)+1], [log2(2) log2(2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(5) log2(5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(10) log2(10)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.5) log2(0.5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.2) log2(0.2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.1) log2(0.1)],'--','LineWidth',.2,'Color',[.6 .6 .6])
set(gca,'xtick',[1:size(transl_act,1)],'xticklabel',string(analysis_legend(:,5)))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4])
ylabel("log2(activity-change)")
sgtitle(enzymelabel+" beta - separate reactants")

%% Not Translated
figure(131)
hold on
for p = 1 : size(data_collection,1)
    ry=[];
    rx=[];
    if cell2mat(analysis_legend(p,1)) == 1
        ry = [ry cell2mat(data_collection(p,1))];
    end
    if cell2mat(analysis_legend(p,2)) == 1
        ry = [ry cell2mat(data_collection(p,2))];
    end
    if cell2mat(analysis_legend(p,3)) == 1
        ry = [ry cell2mat(data_collection(p,3))];
    end
    if cell2mat(analysis_legend(p,4)) == 1
        ry = [ry cell2mat(data_collection(p,4))];
    end
    rx = (p-0.4) + ((p+0.4)-(p-0.4)).*rand(length(ry),1);
    ry(isnan(ry))= []; % remove nans
    if isempty(ry)
        % do nothing
    elseif style == "box"
        boxchart(ones(1,length(ry))*p,ry,"BoxWidth",0.9,'BoxFaceColor',[.4 .4 .4],'MarkerColor',[.4 .4 .4])
    elseif style == "std"
        bar(p,nanmean(ry),'FaceColor',[.8 .8 .8])
        hold on
        errorbar(p,nanmean(ry),nanstd(ry),nanstd(ry),'.k')
        hold on
    end
end
set(gca,'xtick',[1:size(data_collection,1)],'xticklabel',string(analysis_legend(:,5)))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4])
ylabel("AUC")
sgtitle(enzymelabel+" AUC - lumped reactants")

figure(261)
hold on
for p = 1 : size(data_collection,1)
    % Create random numbers around p as x-values
    r1 = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(data_collection(p,1))),1);
    r5 = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(data_collection(p,2))),1);
    r3 = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(data_collection(p,3))),1);
    r4 = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(data_collection(p,4))),1);
    r1e = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(transl_act_excl(p,1))),1);
    r2e = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(transl_act_excl(p,2))),1);
    r3e = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(transl_act_excl(p,3))),1);
    r4e = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(transl_act_excl(p,4))),1);
    
    % Plot y-values
    % Plot excluded because of r2
    if excl_r2 == 1
        plot(r1e,cell2mat(transl_act_excl(p,1)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
        plot(r2e,cell2mat(transl_act_excl(p,2)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
        plot(r3e,cell2mat(transl_act_excl(p,3)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
        plot(r4e,cell2mat(transl_act_excl(p,4)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,1)) == 1
        plot(r1,cell2mat(data_collection(p,1)),'.','Color',[0.8500 0.3250 0.0980]) % orange
    else
        plot(r1,cell2mat(data_collection(p,1)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,2)) == 1
    plot(r5,cell2mat(data_collection(p,2)),'.','Color',[0.9290 0.6940 0.1250]) % yellow
    else
    plot(r5,cell2mat(data_collection(p,2)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,3)) == 1
    plot(r3,cell2mat(data_collection(p,3)),'.','Color',[0.4660 0.6740 0.1880]) % green
    else
        plot(r3,cell2mat(data_collection(p,3)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    if cell2mat(analysis_legend(p,4)) == 1
    plot(r4,cell2mat(data_collection(p,4)),'.','Color',[0 0.4470 0.7410]) % blue
    else
        plot(r4,cell2mat(data_collection(p,4)),'x','Color',[0.8 0.8 0.8],'MarkerSize',.5)
    end
    hold on
end
set(gca,'xtick',[1:size(data_collection,1)],'xticklabel',string(analysis_legend(:,5)))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4])
ylabel("AUC")
sgtitle(enzymelabel+" AUC - lumped reactants")


%% Save scores

% Create score output table
% fiadata_scores = struct;
 fiadata_scores.enzyme = enzymelabel;
 fiadata_scores.fit_function = fit_function;
 fiadata_scores.kolsmirn(1:ctrl_break_idx-1,1) = {nan};
 fiadata_scores.kolsmirn(1:ctrl_break_idx-1,2) = analysis_legend(1:ctrl_break_idx-1,5);
 fiadata_scores.median = fiadata_scores.kolsmirn;
 fiadata_scores.quant75 = fiadata_scores.kolsmirn;
 fiadata_scores.quant95 = fiadata_scores.kolsmirn;
 fiadata_scores.mean = fiadata_scores.kolsmirn;
 fiadata_scores.eom = fiadata_scores.kolsmirn;
 fiadata_scores.std = fiadata_scores.kolsmirn;
 fiadata_scores.std_nearest = fiadata_scores.kolsmirn;
 fiadata_scores.std_DELTA = fiadata_scores.kolsmirn;
 fiadata_scores.quant75_DELTA = fiadata_scores.kolsmirn;
 fiadata_scores.quant95_DELTA = fiadata_scores.kolsmirn;
 fiadata_scores.test_cohenD = fiadata_scores.kolsmirn;
 fiadata_scores.test_effsizeDIVstd = fiadata_scores.kolsmirn;
 fiadata_scores.welchLUMPED_n4 = fiadata_scores.kolsmirn;
 fiadata_scores.welchLUMPED_n200 = fiadata_scores.kolsmirn;
 include_salt_controls = 0

for p = 1 : (ctrl_break_idx-1) % all but combined controls
    clear fitgood & flagged & elig & effs & ctrls
    % Check whether more than 50% of fits good
    fitgood = [length(cell2mat(transl_act(p,1)))>=100 length(cell2mat(transl_act(p,2)))>=100 length(cell2mat(transl_act(p,3)))>=100 length(cell2mat(transl_act(p,4)))>=100];
    % Check if flagged in previous analysis (1 = okay)
    flagged = [cell2mat(analysis_legend(p,1)) cell2mat(analysis_legend(p,2)) cell2mat(analysis_legend(p,3)) cell2mat(analysis_legend(p,4))];
    % Combine, if 2 then eligible
    elig = fitgood + flagged;
    
    % Lump eff and ctrl data together IF eligible
    effs = []
    ctrls = []
    lumpedctrls = [] 
    effs_AUCs = []
    ctrls_AUCs = []
    lumpedctrls_AUCs = [] 
    if elig(1) == 2 % if reactant 1 eligible
        % lump effectors
        effs = [effs cell2mat(transl_act(p,1))]; % save tranlated activities in effs % r2<0.7 already removed
        ctrl_indices = find(contains(string(analysis_legend(:,5)),"control"));
        clear ctrl_act_set
        
        % put separate controls in same table (ctrl)
        if p <= break_idx % set 1-4
            ctrl_indices_sel = ctrl_indices(ctrl_indices<= break_idx);         
        elseif p > break_idx & p < ctrl_break_idx % set 5-8
            ctrl_indices_sel = ctrl_indices(ctrl_indices>break_idx & ctrl_indices < ctrl_break_idx);
        end
        for sel = 1 : length(ctrl_indices_sel)
            ctrls = [ctrls cell2mat(transl_act(ctrl_indices_sel(sel),1))];
        end 
        
        %lumped controls
        if p < break_idx % set1-4
            lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_1st")),1))];
         elseif p >= break_idx & p < 97 % set5-8
             lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_2nd")),1))];
        elseif p >= 97 % ctrls
            lumpedctrls = [lumpedctrls cell2mat(transl_act(p,1))];
        end
        
      clear AUCs_signif & AUCs_signif_lumpedctrls
       AUCs_signif = cell2mat(data(p,1));
       AUCs_signif(cell2mat(r2(p,1))<0.7) = nan; % removing low r2
        effs_AUCs = [effs_AUCs AUCs_signif];
        %lumped controls
        if p < break_idx % set1-4
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_1st")),1));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_1st")),1))<0.7) = nan;
        elseif p >= break_idx & p < 97 % set5-8
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_2nd")),1));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_2nd")),1))<0.7) = nan;
        elseif p >= 97
            AUCs_signif_lumpedctrls = cell2mat(data(p,1));
            AUCs_signif_lumpedctrls(cell2mat(r2(p,1))<0.7) = nan;
        end
        lumpedctrls_AUCs = [lumpedctrls_AUCs AUCs_signif_lumpedctrls];

    
    end
    
    
    
    
    if elig(2) == 2 % if reactant 1 eligible
        effs = [effs cell2mat(transl_act(p,2))]; % save tranlated activities in effs
        ctrl_indices = find(contains(string(analysis_legend(:,5)),"control"));
        clear ctrl_act_set
        if p <= break_idx % set 1-4
            ctrl_indices_sel = ctrl_indices(ctrl_indices<= break_idx);
            
        elseif p > break_idx & p < ctrl_break_idx % set 5-8
            ctrl_indices_sel = ctrl_indices(ctrl_indices>break_idx & ctrl_indices < ctrl_break_idx);
        end
        for sel = 1 : length(ctrl_indices_sel)
            ctrls = [ctrls cell2mat(transl_act(ctrl_indices_sel(sel),2))];
        end 
        
                %lumped controls
        if p < break_idx % set1-4
            lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_1st")),2))];
         elseif p >= break_idx & p < 97 % set5-8
             lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_2nd")),2))];
        elseif p >= 97 % ctrls
            lumpedctrls = [lumpedctrls cell2mat(transl_act(p,2))];
        end
        
      clear AUCs_signif & AUCs_signif_lumpedctrls
       AUCs_signif = cell2mat(data(p,2));
       AUCs_signif(cell2mat(r2(p,2))<0.7) = nan; % removing low r2
        effs_AUCs = [effs_AUCs AUCs_signif];
        %lumped controls
        if p < break_idx % set1-4
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_1st")),2));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_1st")),2))<0.7) = nan;
        elseif p >= break_idx & p < 97 % set5-8
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_2nd")),2));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_2nd")),2))<0.7) = nan;
        elseif p >= 97
            AUCs_signif_lumpedctrls = cell2mat(data(p,2));
            AUCs_signif_lumpedctrls(cell2mat(r2(p,2))<0.7) = nan;
        end
        lumpedctrls_AUCs = [lumpedctrls_AUCs AUCs_signif_lumpedctrls];

        
    end
    
    if elig(3) == 2 % if reactant 1 eligible
        effs = [effs cell2mat(transl_act(p,3))]; % save tranlated activities in effs
        ctrl_indices = find(contains(string(analysis_legend(:,5)),"control"));
        clear ctrl_act_set
        if p <= break_idx % set 1-4
            ctrl_indices_sel = ctrl_indices(ctrl_indices<= break_idx);         
        elseif p > break_idx & p < ctrl_break_idx % set 5-8
            ctrl_indices_sel = ctrl_indices(ctrl_indices>break_idx & ctrl_indices < ctrl_break_idx);
        end
        for sel = 1 : length(ctrl_indices_sel)
            ctrls = [ctrls cell2mat(transl_act(ctrl_indices_sel(sel),3))];
        end
        
                %lumped controls
        if p < break_idx % set1-4
            lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_1st")),3))];
         elseif p >= break_idx & p < 97 % set5-8
             lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_2nd")),3))];
        elseif p >= 97 % ctrls
            lumpedctrls = [lumpedctrls cell2mat(transl_act(p,3))];
        end
        
              clear AUCs_signif & AUCs_signif_lumpedctrls
       AUCs_signif = cell2mat(data(p,3));
       AUCs_signif(cell2mat(r2(p,3))<0.7) = nan; % removing low r2
        effs_AUCs = [effs_AUCs AUCs_signif];
        %lumped controls
        if p < break_idx % set1-4
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_1st")),3));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_1st")),3))<0.7) = nan;
        elseif p >= break_idx & p < 97 % set5-8
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_2nd")),3));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_2nd")),3))<0.7) = nan;
        elseif p >= 97
            AUCs_signif_lumpedctrls = cell2mat(data(p,3));
            AUCs_signif_lumpedctrls(cell2mat(r2(p,3))<0.7) = nan;
        end
        lumpedctrls_AUCs = [lumpedctrls_AUCs AUCs_signif_lumpedctrls];

        
        
    end  
    
    if elig(4) == 2 % if reactant 1 eligible
        effs = [effs cell2mat(transl_act(p,4))]; % save tranlated activities in effs
        ctrl_indices = find(contains(string(analysis_legend(:,5)),"control"));
        clear ctrl_act_set
        if p <= break_idx % set 1-4
            ctrl_indices_sel = ctrl_indices(ctrl_indices<= break_idx);        
        elseif p > break_idx & p < ctrl_break_idx % set 5-8
            ctrl_indices_sel = ctrl_indices(ctrl_indices>break_idx & ctrl_indices < ctrl_break_idx);
        end
        for sel = 1 : length(ctrl_indices_sel)
            ctrls = [ctrls cell2mat(transl_act(ctrl_indices_sel(sel),4))];
        end 
        
                %lumped controls
        if p < break_idx % set1-4
            lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_1st")),4))];
         elseif p >= break_idx & p < 97 % set5-8
             lumpedctrls = [lumpedctrls cell2mat(transl_act(find(strcmp(string(transl_act(:,5)),"ctrl_2nd")),4))];
        elseif p >= 97 % ctrls
            lumpedctrls = [lumpedctrls cell2mat(transl_act(p,4))];
        end
        
              clear AUCs_signif & AUCs_signif_lumpedctrls
       AUCs_signif = cell2mat(data(p,4));
       AUCs_signif(cell2mat(r2(p,4))<0.7) = nan; % removing low r2
        effs_AUCs = [effs_AUCs AUCs_signif];
        %lumped controls
        if p < break_idx % set1-4
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_1st")),4));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_1st")),4))<0.7) = nan;
        elseif p >= break_idx & p < 97 % set5-8
            AUCs_signif_lumpedctrls = cell2mat(data(find(strcmp(string(data(:,5)),"ctrl_2nd")),4));
            AUCs_signif_lumpedctrls(cell2mat(r2(find(strcmp(string(data(:,5)),"ctrl_2nd")),4))<0.7) = nan;
        elseif p >= 97
            AUCs_signif_lumpedctrls = cell2mat(data(p,4));
            AUCs_signif_lumpedctrls(cell2mat(r2(p,4))<0.7) = nan;
        end
        lumpedctrls_AUCs = [lumpedctrls_AUCs AUCs_signif_lumpedctrls];

        
    end  
        
    if isempty(effs)
        % do nothing
    else
        

% WELCH translACT LUMPED
try
    n_free = 4; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
    WELCHscore_n4 = ( nanmean(lumpedctrls) - nanmean(effs) ) / ( sqrt( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs) / n_free ) ) );
    WELCHdf_n4 = ( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs) / n_free ) )^2 / ( ( ( nanvar(lumpedctrls) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(effs) / n_free )^2 / (n_free - 1) ) );
    WELCHp_n4 = 2*tcdf(-abs(WELCHscore_n4),WELCHdf_n4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
    fiadata_scores.welch_translACT_lumped_n4(p,1) = {WELCHp_n4}
    
    n_free = 200; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
    WELCHscore_n200 = ( nanmean(lumpedctrls) - nanmean(effs) ) / ( sqrt( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs) / n_free ) ) );
    WELCHdf_n200 = ( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs) / n_free ) )^2 / ( ( ( nanvar(lumpedctrls) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(effs) / n_free )^2 / (n_free - 1) ) );
    WELCHp_n200 = 2*tcdf(-abs( WELCHscore_n200),WELCHdf_n200); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
    fiadata_scores.welch_translACT_lumped_n200(p,1) = {WELCHp_n200};
    
    n_free = 8*4;
    WELCHscore_n8x4 = ( nanmean(lumpedctrls) - nanmean(effs) ) / ( sqrt( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs) / n_free ) ) );
    WELCHdf_n8x4 = ( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs) / n_free ) )^2 / ( ( ( nanvar(lumpedctrls) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(effs) / n_free )^2 / (n_free - 1) ) );
    WELCHp_n8x4 = 2*tcdf(-abs(WELCHscore_n8x4),WELCHdf_n8x4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
    fiadata_scores.welch_translACT_lumped_n8x4(p,1) = {WELCHp_n8x4};
end


% WELCH LUMPED AUCs
try
    n_free = 4; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
    WELCHscore_n4 = ( nanmean(lumpedctrls) - nanmean(effs_AUCs) ) / ( sqrt( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs_AUCs) / n_free ) ) );
    WELCHdf_n4 = ( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs_AUCs) / n_free ) )^2 / ( ( ( nanvar(lumpedctrls) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(effs_AUCs) / n_free )^2 / (n_free - 1) ) );
    WELCHp_n4 = 2*tcdf(-abs(WELCHscore_n4),WELCHdf_n4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
    fiadata_scores.welch_AUCs_lumped_n4(p,1) = {WELCHp_n4}
    
    n_free = 200; %8*4; % n to calculate degrees of freedom; 8 datapoints x 4 replicates
    WELCHscore_n200 = ( nanmean(lumpedctrls) - nanmean(effs_AUCs) ) / ( sqrt( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs_AUCs) / n_free ) ) );
    WELCHdf_n200 = ( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs_AUCs) / n_free ) )^2 / ( ( ( nanvar(lumpedctrls) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(effs_AUCs) / n_free )^2 / (n_free - 1) ) );
    WELCHp_n200 = 2*tcdf(-abs( WELCHscore_n200),WELCHdf_n200); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
    fiadata_scores.welch_AUCs_lumped_n200(p,1) = {WELCHp_n200};
    
    n_free = 8*4;
    WELCHscore_n8x4 = ( nanmean(lumpedctrls) - nanmean(effs_AUCs) ) / ( sqrt( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs_AUCs) / n_free ) ) );
    WELCHdf_n8x4 = ( ( nanvar(lumpedctrls) / n_free ) + ( nanvar(effs_AUCs) / n_free ) )^2 / ( ( ( nanvar(lumpedctrls) / n_free )^2 / (n_free - 1) ) + ( ( nanvar(effs_AUCs) / n_free )^2 / (n_free - 1) ) );
    WELCHp_n8x4 = 2*tcdf(-abs(WELCHscore_n8x4),WELCHdf_n8x4); %tpdf(WELCHscore_n4,WELCHdf_n4) %tcdf(WELCHscore_n4,WELCHdf_n4)
    fiadata_scores.welch_AUCs_lumped_n8x4(p,1) = {WELCHp_n8x4};
end

        [~,kolsmirn_p] = kstest2(ctrls,effs);
        fiadata_scores.kolsmirn(p,1) = {kolsmirn_p};
        fiadata_scores.median(p,1) = {nanmedian(effs)};
        fiadata_scores.mean(p,1) = {nanmean(effs)};
        fiadata_scores.std(p,1) = {nanstd(effs)};
        fiadata_scores.eom(p,1) = {nanstd(effs)/sqrt(length(effs))};
        
        % std closest to 0
        mean_std = [nanmean(effs)+nanstd(effs) nanmean(effs)-nanstd(effs)];
        if sum(mean_std<0) == 1 % if one higher and other lower than 0, score = 0
            fiadata_scores.std_nearest(p,1) = {0};
             closeststdtozero = nan;
        else
            closeststdtozero = mean_std(find(abs(mean_std) == min(abs(mean_std)))); % save value
            fiadata_scores.std_nearest(p,1) = {mean_std(find(abs(mean_std) == min(abs(mean_std))))}; % value closer to zero saved
        end
        
        %% GET DELTA-STD SCORE
        mean_std_ctrls = [nanmean(lumpedctrls)+nanstd(lumpedctrls) nanmean(lumpedctrls)-nanstd(lumpedctrls)]; % get ctrl std boundaries
        % get delta of std closest to std of controls
        if closeststdtozero > min(lumpedctrls) & closeststdtozero < max(lumpedctrls) % within ctrl range?
            fiadata_scores.std_DELTA(p,1) = {0};
        else % outside of ctrl range
            [mindelta,closestIndex] = min(abs( mean_std_ctrls -(mean_std(find(abs(mean_std) == min(abs(mean_std)))))));
            fiadata_scores.std_DELTA(p,1) = {mindelta};
        end
        
        
        
        %% get stds of respective controls // set1-4 or 5-8
        % find the one closer to std_nearest
        
        % TESTING COHEN AND OTHERS
        fiadata_scores.test_cohenD(p,1) = {(nanmean(effs)-nanmean(ctrls))/sqrt((nanstd(effs)^2 + nanstd(ctrls)^2)/2)};
        fiadata_scores.test_effsizeDIVstd(p,1) = {(nanmean(effs)-nanmean(ctrls))/nanstd(ctrls)};
        
        % WELCH LUMPED (assuming normal distribution)
        
        fiadata_scores.welchLUMPED_n4(p,1) = {(nanmean(effs) - nanmean(ctrls)) / sqrt(((nanstd(effs)^2/  4 )+(nanstd(ctrls)^2/  4 )))};
        fiadata_scores.welchLUMPED_n200(p,1) = {(nanmean(effs) - nanmean(ctrls)) / sqrt(((nanstd(effs)^2/  200 )+(nanstd(ctrls)^2/  200 )))}; 
        
        
        clear mindelta & closestIndex & quant75_raw & quant75save & quant75_ctrl_raw
        % quantiles 25 75
        quant75_raw = quantile(effs,[0.25 0.75]);
        if sum(quant75_raw<0) == 1 % if one value below and one beyond 0, set score to 0
            fiadata_scores.quant75(p,1) = {0};
            quant75save = nan;
        else
            
            %if enzymelabel == "MaeB-repipetrepeat"
                quant75save = quant75_raw(find(abs(quant75_raw) == min(abs(quant75_raw))));
                fiadata_scores.quant75(p,1) = {quant75save(1)};
           % else
            %    fiadata_scores.quant75(p,1) = {quant75_raw(find(abs(quant75_raw) == min(abs(quant75_raw))))}; % value closer to zero saved
            %end
        end
        % GET DELTA-QUANT75
        quant75_ctrl_raw = quantile(lumpedctrls,[0.25 0.75]);
        % get delta of std closest to std of controls
        if quant75save > min(quant75_ctrl_raw) & quant75save < max(quant75_ctrl_raw) % within ctrl range?
            fiadata_scores.quant75_DELTA(p,1) = {0};
        else % outside of ctrl range
            [mindelta,closestIndex] = min(abs( quant75_ctrl_raw -(quant75_raw(find(abs(quant75_raw) == min(abs(quant75_raw)))))));
            fiadata_scores.quant75_DELTA(p,1) = {mindelta};
        end

        % quantiles 2.5 97.5
        quant95_raw = quantile(effs,[0.025 0.975]);
        if sum(quant95_raw<0) == 1 % if one value below and one beyond 0, set score to 0
            fiadata_scores.quant95(p,1) = {0};
            quant95save = nan;
        else
                quant95save = quant95_raw(find(abs(quant95_raw) == min(abs(quant95_raw))));
                fiadata_scores.quant95(p,1) = {quant95save(1)};
            %fiadata_scores.quant95(p,1) = {quant95_raw(find(abs(quant95_raw) == min(abs(quant95_raw))))}; % value closer to zero saved
        end
        
        % GET DELTA-QUANT95
        quant95_ctrl_raw = quantile(lumpedctrls,[0.025 0.975]);
        % get delta of std closest to std of controls
        if quant95save > min(quant95_ctrl_raw) & quant95save < max(quant95_ctrl_raw) % within ctrl range?
            fiadata_scores.quant95_DELTA(p,1) = {0};
        else % outside of ctrl range
            [mindelta,closestIndex] = min(abs( quant95_ctrl_raw -(quant95_raw(find(abs(quant95_raw) == min(abs(quant95_raw)))))));
            fiadata_scores.quant95_DELTA(p,1) = {mindelta};
        end
        
        
    end
end


figure(130)
for q = 1 : size(fiadata_scores.quant95,1)
    plot([q-0.45 q+0.45],[cell2mat(fiadata_scores.quant95(q,1)) cell2mat(fiadata_scores.quant95(q,1))],'-c','LineWidth',2)
    if enzymelabel == "MaeB-repipetrepeat"
        plot([q-0.45 q+0.45],[max(cell2mat(fiadata_scores.quant75(q,1))) max(cell2mat(fiadata_scores.quant75(q,1)))],'-b','LineWidth',2)
    else
        plot([q-0.45 q+0.45],[cell2mat(fiadata_scores.quant75(q,1)) cell2mat(fiadata_scores.quant75(q,1))],'-b','LineWidth',2)
    end
end


figure(136)
hold on
for q = 1 : size(fiadata_scores.quant95,1)
    plot([q-0.45 q+0.45],[cell2mat(fiadata_scores.quant95(q,1)) cell2mat(fiadata_scores.quant95(q,1))],'-c','LineWidth',2)
    if enzymelabel == "MaeB-repipetrepeat"
        plot([q-0.45 q+0.45],[max(cell2mat(fiadata_scores.quant75(q,1))) max(cell2mat(fiadata_scores.quant75(q,1)))],'-b','LineWidth',2)
    else
        plot([q-0.45 q+0.45],[cell2mat(fiadata_scores.quant75(q,1)) cell2mat(fiadata_scores.quant75(q,1))],'-b','LineWidth',2)
    end
end

% label
plot([0 size(transl_act,1)+1], [0 0],'-k','LineWidth',3)
plot([0 size(transl_act,1)+1], [log2(2) log2(2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(5) log2(5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(10) log2(10)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.5) log2(0.5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.2) log2(0.2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
plot([0 size(transl_act,1)+1], [log2(0.1) log2(0.1)],'--','LineWidth',.2,'Color',[.6 .6 .6])
set(gca,'xtick',[1:size(transl_act,1)],'xticklabel',string(analysis_legend(:,5)))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4])
ylabel("log2(activity-change)")
ylim([min(min_ry)-.5 max(max_ry)+.5])
sgtitle(enzymelabel+" beta - lumped reactants")


% add flags to overview
fiadata_scores.anyflag = anyflag;
fiadata_scores.flaglabellist = flaglabellist;












% 
% % SIGMOID TEST
% % data = fiadata_effs.SIGMOIDsingle_AUC; % 200 bootstrapped AUCs
% % data_lb =  fiadata_effs.SIGMOIDsingle_AUC95lb; % 200 bootstrapped AUCs
% % data_ub =  fiadata_effs.SIGMOIDsingle_AUC95ub; % 200 bootstrapped AUCs
% % r2 = fiadata_effs.SIGMOIDsingle_r2; % R2s corresponding to AUCs
% 
% 
% 
% 
% AUC2act_all = fiadata_effs.SIG_AUC2activity_all; % AUC translated to change in beta
% AUC2act_1st = fiadata_effs.SIG_AUC2activity_1st; % AUC translated to change in beta
% AUC2act_2nd = fiadata_effs.SIG_AUC2activity_2nd; % AUC translated to change in beta
% 
% transl_act = fiadata_effs.act; % empty table for activities
% transl_act_excl = fiadata_effs.act; % empty table for activities
% 
% 
% % find split position
% break_idx = find(strcmp(data(:,5),"control4"));
% ctrl_break_idx = find(strcmp(data(:,5),"ctrl_all"));
% 
% 
% 
% 
% %% change all entries in data to NaN that (i) are excluded in ANALYSIS_LEGEND or (ii) correspond to low r2 (based on predefined cutoff)
% %% also exclude all where <2 datapoints left
% % then save exluded datapoints in new table
% for j = 1 : 4 % for each rectant
%     for i = 1 : size(data,1) % for each effector
%         clear r2s & data_save & r2_exc_idx & r2_exc_save
%         
%         r2s_save = cell2mat(r2(i,j)); % get 200 R2s
%         data_save = cell2mat(data(i,j)); % get perfect SIGfit AUCs
%         lb_save = cell2mat(data_lb(i,j)); % get perfect SIGfit AUCs
%         ub_save = cell2mat(data_ub(i,j)); % get perfect SIGfit AUCs
%         
%         %% TRANSLATE AUC TO ACTIVITY AND SAVE
%         if any(data_save) & r2s_save>0.7
%             clear AUC2act_table & AUC2ac_save
%             AUC2ac_save(1:size(data_save)) = nan;
%             % get correct translation table
%             if split_on == 1
%                 if i <= break_idx % set 1-4
%                     AUC2act_table = AUC2act_1st;
%                 elseif i > break_idx & i < ctrl_break_idx % set 5-8
%                     AUC2act_table = AUC2act_2nd;
%                 elseif i >= ctrl_break_idx % CONTROL COLLECTIONS TRANSLATED WITH OWN TABLE
%                     if strcmp(string(analysis_legend(i,5)),"ctrl_all")
%                         AUC2act_table = AUC2act_all;
%                     elseif strcmp(string(analysis_legend(i,5)),"ctrl_1st")
%                         AUC2act_table = AUC2act_1st;
%                     elseif strcmp(string(analysis_legend(i,5)),"ctrl_2nd")
%                         AUC2act_table = AUC2act_2nd;
%                     end
%                 end
%             elseif split_on == 0
%                 AUC2act_table = AUC2act_all;
%             end
%             % translate each AUC
%             clear c & index & ACT_save & c_lb & index_lb & ACT_save_lb & c_ub & index_ub & ACT_save_ub
%             [c index] = min(abs(AUC2act_table(:,j)-data_save));
%             [c_lb index_lb] = min(abs(AUC2act_table(:,j)-lb_save));
%             [c_ub index_ub] = min(abs(AUC2act_table(:,j)-ub_save));
%             if ~isnan(data_save)
%                 if lnACT_on == 0
%                     ACT_save = AUC2act_table(index,5);
%                     ACT_lb = AUC2act_table(index_lb,5);
%                     ACT_ub = AUC2act_table(index_ub,5);
%                 elseif lnACT_on == 1
%                     ACT_save = log(AUC2act_table(index,5));
%                     ACT_lb = log(AUC2act_table(index_lb,5));
%                     ACT_ub = log(AUC2act_table(index_ub,5));
%                 elseif lnACT_on == 2
%                     ACT_save = log2(AUC2act_table(index,5));
%                     ACT_lb = log2(AUC2act_table(index_lb,5));
%                     ACT_ub = log2(AUC2act_table(index_ub,5));
%                 end
%             else
%                 ACT_save = nan;
%                 ACT_lb = nan;
%                 ACT_ub = nan;
%             end
%         else
%             ACT_save = nan;
%             ACT_lb = nan;
%             ACT_ub = nan;
%         end
%         SIG_transl_act(i,j) = {ACT_save}; % fit
%         SIG_transl_act_ub(i,j) = {ACT_ub}; % ub
%         SIG_transl_act_lb(i,j) = {ACT_lb}; % lb
%         
%     end
% end
% 
% style = "std"; % "std"
% fiadata_scores.SIG_transl_act = SIG_transl_act;
% fiadata_scores.SIG_transl_act_ub = SIG_transl_act_ub;
% fiadata_scores.SIG_transl_act_lb = SIG_transl_act_lb;
% 
% 
% 
% 

% %% Translated
% 
% fig4 = figure(930)
% hold on
% max_ry = [];
% min_ry = [];
% for p = 1 : size(SIG_transl_act,1)    % each effector
%     figure(930)
% %% x for datapoints
%     r1 = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(SIG_transl_act(p,1))),1);
%     r2 = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(SIG_transl_act(p,2))),1);
%     r3 = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(SIG_transl_act(p,3))),1);
%     r4 = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(SIG_transl_act(p,4))),1);
%     r1e = (p-0.4) + ((p-0.2)-(p-0.4)).*rand(length(cell2mat(transl_act_excl(p,1))),1);
%     r2e = (p-0.2) + ((p)-(p-0.2)).*rand(length(cell2mat(transl_act_excl(p,2))),1);
%     r3e = (p) + ((p+0.2)-(p)).*rand(length(cell2mat(transl_act_excl(p,3))),1);
%     r4e = (p+0.2) + ((p+0.4)-(p+0.2)).*rand(length(cell2mat(transl_act_excl(p,4))),1);
%     
% %% for barplot
%     ry=[];
%     ry_lb =[];
%     ry_ub = [];
%     rx=[];
%     
%     % get y
%     if cell2mat(analysis_legend(p,1)) == 1
%         ry = [ry cell2mat(SIG_transl_act(p,1))];
%         ry_lb =  [ry_lb cell2mat(SIG_transl_act_lb(p,1))];
%         ry_ub =  [ry_ub cell2mat(SIG_transl_act_ub(p,1))];  
%     end
%     if cell2mat(analysis_legend(p,2)) == 1
%         ry = [ry cell2mat(SIG_transl_act(p,2))];
%         ry_lb =  [ry_lb cell2mat(SIG_transl_act_lb(p,2))];
%         ry_ub =  [ry_ub cell2mat(SIG_transl_act_ub(p,2))];
%     end
%     if cell2mat(analysis_legend(p,3)) == 1
%         ry = [ry cell2mat(SIG_transl_act(p,3))];
%         ry_lb =  [ry_lb cell2mat(SIG_transl_act_lb(p,3))];
%         ry_ub =  [ry_ub cell2mat(SIG_transl_act_ub(p,3))];
%     end
%     if cell2mat(analysis_legend(p,4)) == 1
%         ry = [ry cell2mat(SIG_transl_act(p,4))];
%         ry_lb =  [ry_lb cell2mat(SIG_transl_act_lb(p,4))];
%         ry_ub =  [ry_ub cell2mat(SIG_transl_act_ub(p,4))];
%     end
%     ry_ub(isnan(ry))= []; % remove nans
%     ry_lb(isnan(ry))= []; % remove nans
%     ry(isnan(ry))= []; % remove nans
%     
% %     % save max and min ry
%      max_ry = [max_ry max(ry)];
%      min_ry = [min_ry min(ry)];
%     
%     % create x
%      rx = (p-0.4) + ((p+0.4)-(p-0.4)).*rand(length(ry),1);
%     
%     %%
%   %%  
% %% 2do: 
%      if isempty(ry)
%          % do nothing
%      elseif style == "box"
%          boxchart(ones(1,length(ry))*p,ry,"BoxWidth",0.9,'BoxFaceColor',[.4 .4 .4],'MarkerColor',[.4 .4 .4])
%      elseif style == "std"
%          bar(p,nanmean(ry),'FaceColor',[.8 .8 .8])
%          hold on
%          
%          %% DOUBLECHECK
%          errorbar(p,nanmean(ry),abs(nanmean(ry) - min([ry_ub ry_lb])),abs(nanmean(ry) - max([ry_ub ry_lb])),'.k') %%%%%%%% IS THIS CORRECT???????????????????
%          %%%%%%%%%%%%%%%%%%%%%
%          
%          
%          hold on
%      end
% end
%  

% 
% figure(930)
% % add flags
% flaglabellist(1:size(transl_act,1)) = "";
% for p = 1 : size(transl_act,1)
%     if flags_SPprecip_any(p) == 1
%         fnan = fill([p-0.45 p+0.45 p+0.45 p-0.45 p-0.45],[min(min_ry)-.5 min(min_ry)-.5 max(max_ry)+.5 max(max_ry)+.5 min(min_ry)-.5],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
%         uistack(fnan,'bottom')
%         if string(analysis_legend(p,6)) == "PRODUCT"
%             flaglabel = "P";
%             flaglabellist(p) = "P";
%         elseif string(analysis_legend(p,6)) == "SUBSTRATE"
%             flaglabel = "S";
%             flaglabellist(p) = "S";
%         elseif string(analysis_legend(p,6)) == "PRECIP"
%             flaglabel = "precip";
%             flaglabellist(p) = "precip";
%         end
%         txt=text(p,min(min_ry)+.2,flaglabel,'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black') 
%     elseif flags_SALT_any(p) == 1
%         fnan = fill([p-0.45 p+0.45 p+0.45 p-0.45 p-0.45],[min(min_ry)-.5 min(min_ry)-.5 max(max_ry)+.5 max(max_ry)+.5 min(min_ry)-.5],[.8 .8 .8],'EdgeColor',[.9 .9 .9]);
%         uistack(fnan,'bottom')
%         txt=text(p,min(min_ry)+.2,"salt",'VerticalAlignment', 'middle','rotation',90,'FontSize',8,'Color','black')
%         flaglabellist(p) = "salt";
%     end
% end
% 
% 
% % label
% plot([0 size(transl_act,1)+1], [0 0],'-k','LineWidth',3)
% plot([0 size(transl_act,1)+1], [log2(2) log2(2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
% plot([0 size(transl_act,1)+1], [log2(5) log2(5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
% plot([0 size(transl_act,1)+1], [log2(10) log2(10)],'--','LineWidth',.2,'Color',[.6 .6 .6])
% plot([0 size(transl_act,1)+1], [log2(0.5) log2(0.5)],'--','LineWidth',.2,'Color',[.6 .6 .6])
% plot([0 size(transl_act,1)+1], [log2(0.2) log2(0.2)],'--','LineWidth',.2,'Color',[.6 .6 .6])
% plot([0 size(transl_act,1)+1], [log2(0.1) log2(0.1)],'--','LineWidth',.2,'Color',[.6 .6 .6])
% set(gca,'xtick',[1:size(transl_act,1)],'xticklabel',string(analysis_legend(:,5)))
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1.0, 0.4])
% ylabel("log2(activity-change)")
% ylim([min(min_ry)-.5 max(max_ry)+.5])
% sgtitle(enzymelabel+" beta - lumped reactants")



% 
% % barplot deltascore
% figure
% bar(1:length(cell2mat(fiadata_scores.std_DELTA(:,1))),cell2mat(fiadata_scores.std_DELTA(:,1)))
% 
% %volcano plot
% figure
% volcano_x = cell2mat(fiadata_scores.mean(:,1));
% volcano_x(find(fiadata_scores.anyflag(1:96)))=nan;
% 
% volcano_y = -log10(cell2mat(fiadata_scores.welch_translACT_lumped_n4(1:96,1)))
% % plot mean of lumped vs -log10 of lumped pvalue (n=4?)
% plot(volcano_x,volcano_y,'.k','MarkerSize',15)
% %errorbar(volcano_x,volcano_y,cell2mat(fiadata_scores.std(:,1)),'horizontal','ok')
% title("volcano")
% labels = fiadata_scores.kolsmirn(:,2);
% text(volcano_x,volcano_y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
if save_all_results == 1
    figure(130)
    label = string(enzymelabel) + "_EFFECTORS_activitychange_barplot" + ".png"
    saveas(gca,label)
end

if plot_only_relevant == 1
    figure(260)
    close
    figure(131)
    close
    figure(261)
    close
    figure(136)
    close
end
end