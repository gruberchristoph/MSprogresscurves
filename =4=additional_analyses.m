% Histograms of Supplementary Figure 1 and 2
close all
clc
clear

%% INPUTS
analysis_path = "W:\~archive\Christoph\code_github_upload\"



%% importing absolute intracellular metabolite concentrations and plotting histogram
file = {analysis_path+"CODETABLES.xlsx"}; %path to overview
[max_conc_A,max_conc_B] = xlsread(file{1},"max_conc","C6 : C90");
figure(1)
histogram(max_conc_A,34,'FaceColor',[.5 .5 .5])
ylabel("Number of metabolites")
xlabel("Highest measured intracellular concentration [mM]")

label = append("maximum_intracellular_conc_mM.png")
saveas(gca,label)

%% importing reported Km values and plotting histogram
file = {analysis_path+"CODETABLES.xlsx"}; %path to overview
[Km_A,Km_B] = xlsread(file{1},"Km_values","I4 : I23");
figure(2)
histogram(Km_A,34,'FaceColor',[.5 .5 .5])
ylabel("Number of enzymes")
xlabel("Highest measured KM among substrates [uM]")

label = append("KMs_uM.png")
saveas(gca,label)





%% CHEMICAL SIMILARITY ANALYSIS FOR INHIBITORS
%% OPTIONS
only_new_hits = 0;

%% global or local similarity?
sim_matching = "local" % "local" or "global" (both SIMPCOMP2), "usrcat", "tanimoto" (both USRCAT tool)

% import data on metabolite similarity obtained from SIMCOMP2 online tool
if sim_matching == "global"
[simCOMP2_a,simCOMP2_b,~] = xlsread(analysis_path+"=4_inputs=\chemsimilaritymatrix_global_SIMCOMP2.csv","A3:C50000");
elseif sim_matching == "local"
[simCOMP2_a,simCOMP2_b,~] = xlsread(analysis_path+"=4_inputs=\chemsimilaritymatrix_local_SIMCOMP2.csv","A3:C50000");
elseif sim_matching == "usrcat"
[simCOMP2_a,simCOMP2_b,~] = xlsread(analysis_path+"=4_inputs=\chemsimilaritymatrix_usrcat.csv","A2:C50000");
elseif sim_matching == "tanimoto"
[simCOMP2_a,simCOMP2_b,~] = xlsread(analysis_path+"=4_inputs=\chemsimilaritymatrix_tanimoto.csv","A2:C50000");
end

% remove all nan rows
simCOMP2_b = simCOMP2_b(~isnan(simCOMP2_a),:);
simCOMP2_a = simCOMP2_a(~isnan(simCOMP2_a),:);
% combine to find duplicates
findDUPS_save = strcat(simCOMP2_b(:,1),simCOMP2_b(:,2));

% get info from excel
file = {analysis_path+"CODETABLES.xlsx"};
[screen_results_a,screen_results_b,screen_results_c] = xlsread(file{1},"overview_allhits","A1 : AB1861");

enzyme_labels = string(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"Enzyme"))));
eff_labels = string(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"Effectors"))));
in_ANY_database = cell2mat(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"ANY (database)"))));
SCORE = cell2mat(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"1st quantile [log2-fold activity change]"))));
mode = string(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"mode"))));
flagged = cell2mat(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"Any flag?"))));

% enzymes -> map reactants
file = {analysis_path+"CODETABLES.xlsx"};%path to screen overview
[map_reactant_a,map_reactant_b,map_reactant_c] = xlsread(file{1},"map_E_to_met","A1 : G21");

% effectors -> map keggs
file = {analysis_path+"CODETABLES.xlsx"};
[map_kegg_a,map_kegg_b,map_kegg_c] = xlsread(file{1},"map_eff_to_KEGG","A1 : E89");


%% CREATE TABLE WITH MAXIMUM SIMILARITIES
% create output tables
simmatrix_all = cell(length(enzyme_labels),1);
simmatrix_max = nan(length(enzyme_labels),1);
hit_matrix = zeros(length(enzyme_labels),1);
weakhit_matrix = zeros(length(enzyme_labels),1);
% for each enzymes
for i = 1 : length(enzyme_labels)
    clear enz & enz_idx & reactant_array & reactant_KEGGs & KEGG_id_save
    
    if flagged(i) ~= 0;
        simmatrix_all(i) = {nan};
        simmatrix_max(i) = nan;
        hit_matrix(i) = 0;
        weakhit_matrix(i) = 0; 
    else
        % enzyme name
        enz = enzyme_labels(i);
        if enz == "Zwf*"
            enz = "Zwf";
        end
        % enzyme idx in overview;
        enz_idx = find(strcmp(map_reactant_b(:,find(strcmp(string(map_reactant_b(1,:)),"enzyme"))),enz));
        % get reactants
        reactant_array = string(map_reactant_b(enz_idx(1),find(ismember(string(map_reactant_b(1,:)),{'reactants-1','reactants-2','reactants-3','reactants-4','reactants-5'}))));
    
        % get KEGG IDs
        reactant_KEGGs = ["" "" "" "" ""];
        KEGG_id_save = [map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(1))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG")))...
            map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(2))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG")))...
            map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(3))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG")))...
            map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(4))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG")))...
            map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(5))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG")))];
        reactant_KEGGs(1:length(KEGG_id_save)) = string(KEGG_id_save);
        
        %% check if any are not found
        if isempty(map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(1))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG"))))
            reactant_array(1);
        end
        if isempty(map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(2))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG"))))
            reactant_array(2);
        end
        if isempty(map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(3))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG"))))
            reactant_array(3);
        end
        if isempty(map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(4))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG"))))
            reactant_array(4);
        end
        if isempty(map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),reactant_array(5))),find(strcmp(map_kegg_b(1,:),"CAS/KEGG"))))
            reactant_array(5);
        end
    
    
       
        % for each effectors
        eff = eff_labels(i);
        eff_KEGG = string(map_kegg_b(find(strcmp(string(map_kegg_b(:,find(strcmp(map_kegg_b(1,:),"abbr.")))),eff)),find(strcmp(map_kegg_b(1,:),"CAS/KEGG"))));
        
        eff_react_KEGGcat = strcat(eff_KEGG,reactant_KEGGs);
        
        % get similarity for each reactant
        sim_SCORES_all = [simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(1))))...
            simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(2))))...
            simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(3))))...
            simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(4))))...
            simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(5))))];
        
           if strlength(eff_react_KEGGcat(1))> 6
             if isempty(simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(1)))))
                 strcat("ERROR 1: ",string(eff_react_KEGGcat(1)))
             end
           end
           if strlength(eff_react_KEGGcat(2))> 6
             if isempty(simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(2)))))
                 strcat("ERROR 2: ",string(eff_react_KEGGcat(2)))
             end
           end
           if strlength(eff_react_KEGGcat(3))> 6
             if isempty(simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(3)))))
                strcat("ERROR 3: ",string(eff_react_KEGGcat(3)))
             end
           end
           if strlength(eff_react_KEGGcat(4))> 6
             if isempty(simCOMP2_a(find(ismember(findDUPS_save,eff_react_KEGGcat(4)))))
                 strcat("ERROR 4: ",string(eff_react_KEGGcat(4)))
             end
           end
        
        % get max similarity
        sim_SCORES_max = nanmax(sim_SCORES_all);
        
        % save in overview
        simmatrix_all(i) = {sim_SCORES_all};
        simmatrix_max(i) = sim_SCORES_max;
          
        % mark if inhibitor (negative score) and if strong or not
        if SCORE(i) < -1.131
            hit_matrix(i) = 1;
            weakhit_matrix(i) = 1;

        elseif SCORE(i) < -0.566
                weakhit_matrix(i) = 1;
        end
    end
end

nbins = 20;

%% FIGURE
figure(13)
ONLYHITS_simmatrix_max = simmatrix_max(find(hit_matrix))
data = ONLYHITS_simmatrix_max(:); % Create sample data.
edges = linspace(0, 1, nbins+1); % Create 20 bins.
% Plot the histogram.
histogram(data, 'BinEdges',edges);
% Fancy up the graph.
grid on;
xlim([0, 1]);
xlabel('Data Value', 'FontSize', 14);
ylabel('Bin Count', 'FontSize', 14);
title('Histogram of Data', 'FontSize', 14);
%
"percent below 0.5: " + string(round(sum(ONLYHITS_simmatrix_max<0.5)/length(ONLYHITS_simmatrix_max)*100,1))
"percent below 0.25: " + string(round(sum(ONLYHITS_simmatrix_max<0.25)/length(ONLYHITS_simmatrix_max)*100,1))


label = append("chemsimilarity.png")
saveas(gca,label)









%% BARPLOTS HITS PER PATHWAY, HITS PER METABOLITE, HITS PER ENZYME

% get info from excel
file = {analysis_path+"CODETABLES.xlsx"};
[screen_results_a,screen_results_b,screen_results_c] = xlsread(file{1},"overview_allhits","A1 : AA1861");
file = {analysis_path+"CODETABLES.xlsx"};
[pathway_mapping_a,pathway_mapping_b,pathway_mapping_c] = xlsread(file{1},"pathway_mapping","A1 : C80");
[enzyme_mapping_a,enzyme_mapping_b,enzyme_mapping_c] = xlsread(file{1},"pathway_mapping","G2 : G21");

% extract results
enzyme_labels = string(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"Enzyme"))));
eff_labels = string(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"Effectors"))));
in_ANY_database = cell2mat(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"ANY (database)"))));
SCORE = cell2mat(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"1st quantile [log2-fold activity change]"))));
mode = string(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"mode"))));
flagged = cell2mat(screen_results_c(2:end,find(strcmp(string(screen_results_c(1,:)),"Any flag?"))));
% extract mapping
map_pw1 = string(pathway_mapping_c(2:end,find(strcmp(string(pathway_mapping_c(1,:)),"pathway1"))));
map_pw2 = string(pathway_mapping_c(2:end,find(strcmp(string(pathway_mapping_c(1,:)),"pathway2"))));
map_effs = string(pathway_mapping_c(2:end,find(strcmp(string(pathway_mapping_c(1,:)),"effectors"))));
map_enz = string(enzyme_mapping_b);

%% per metabolite
nhits = zeros(length(map_effs),1);
for i = 1 : length (map_effs)
   eff = map_effs(i)

   eff_idx = find(strcmp(eff_labels,eff)); % find index of all tests of this effector
   nhits(i) = sum(abs(SCORE(eff_idx))>1.131 & flagged(eff_idx)==0)
end

figure(123)
nbins = length(map_effs)
figure
h = bar(nhits)
h.FaceColor = [0.9290 0.6940 0.1250];
ylabel("Number of affected enzymes")
xticks([1:1:length(map_effs)])
xticklabels(map_effs)

%% per enzyme
Ehits = zeros(length(map_enz),1);
for i = 1 : length (map_enz)
   enz = map_enz(i);

   enz_idx = find(strcmp(enzyme_labels,enz)); % find index of all tests of this effector
   Ehits(i) = sum(abs(SCORE(enz_idx))>1.131 & flagged(enz_idx)==0)
end

figure(124)
nbins = length(map_enz)
figure
h = bar(Ehits)
h.FaceColor = [0.4940 0.1840 0.5560];
ylabel("Number of identified effectors")
xticks([1:1:length(map_enz)])
xticklabels(map_enz)

label = append("hits_per_enzyme.png")
saveas(gca,label)


%% per pathway1
pwlevel = map_pw1;

pw = unique(pwlevel,'stable');
PWhits = zeros(length(pw),1);
PWhitsratio = zeros(length(pw),1);
effs_per_pathway = zeros(length(pw),1);
for i = 1 : length (pw)
   pw_effs = map_effs(find(strcmp(pwlevel,pw(i))))


   for j = 1 : length (pw_effs)
       eff = pw_effs(j)
    
       eff_idx = find(strcmp(eff_labels,eff)); % find index of all tests of this effector
       
       PWhits(i) = PWhits(i) + sum(abs(SCORE(eff_idx))>1.131 & flagged(eff_idx)==0);
   end
   PWhitsratio(i) = PWhits(i) / length(pw_effs);
  effs_per_pathway(i) = length(pw_effs);
end

figure(125)
nbins = length(pw)
figure
h = bar(PWhitsratio)
h.FaceColor = [0.9290 0.6940 0.1250];
ylabel("Average number of effects per metabolite")
xticks([1:1:length(pw)])
xticklabels(pw)

label = append("hits_per_metabolite.png")
saveas(gca,label)

%% per pathway2
pwlevel = map_pw2;

pw = unique(pwlevel,'stable');
PWhits = zeros(length(pw),1);
PWhitsratio = zeros(length(pw),1);
effs_per_pathway = zeros(length(pw),1);
for i = 1 : length (pw)
   pw_effs = map_effs(find(strcmp(pwlevel,pw(i))))


   for j = 1 : length (pw_effs)
       eff = pw_effs(j)
    
       eff_idx = find(strcmp(eff_labels,eff)); % find index of all tests of this effector
       
       PWhits(i) = PWhits(i) + sum(abs(SCORE(eff_idx))>1.131 & flagged(eff_idx)==0);
   end
   PWhitsratio(i) = PWhits(i) / length(pw_effs);
   effs_per_pathway(i) = length(pw_effs);
end

figure(126)
nbins = length(pw)
figure
h = bar(PWhitsratio)
h.FaceColor = [0.6 0.6 0.6];
ylabel("Average number of effects per metabolite")
xticks([1:1:length(pw)])
xticklabels(pw)

label = append("avg_hits_per_pathway.png")
saveas(gca,label)




%% plot boxplot of number of interactions

% get info from excel
file = {analysis_path+"CODETABLES.xlsx"};
[int_num_a,int_num_b,int_num_c] = xlsread(file{1},"figure4d","A2 : C39");

legend = ["split","irrev","interm","other"];
split_array = int_num_a(find(strcmp(int_num_b(:,1),"split")));
irrev_array = int_num_a(find(strcmp(int_num_b(:,1),"irrev")));
interm_array = int_num_a(find(strcmp(int_num_b(:,1),"interm")));
other_array = int_num_a(find(strcmp(int_num_b(:,1),"other")));

x = [irrev_array; interm_array; split_array; other_array];
g = [zeros(length(irrev_array), 1); ones(length(interm_array), 1); 2*ones(length(split_array), 1); 3*ones(length(other_array), 1)];

figure(1234)
hold on
boxplot(x, g)
scatter(ones(length(irrev_array),1)-0.2+rand(length(irrev_array),1)*0.2,irrev_array,'k','MarkerFaceColor',[.6 .6 .6])
scatter(2*ones(length(interm_array),1)-0.2+rand(length(interm_array),1)*0.2,interm_array,'k','MarkerFaceColor',[.6 .6 .6])
scatter(3*ones(length(split_array),1)-0.2+rand(length(split_array),1)*0.2,split_array,'k','MarkerFaceColor',[.6 .6 .6])
scatter(4*ones(length(other_array),1)-0.2+rand(length(other_array),1)*0.2,other_array,'k','MarkerFaceColor',[.6 .6 .6])
ylabel("number of effectors")

label = append("boxplots_hits_per_enzyme_role.png")
saveas(gca,label)
