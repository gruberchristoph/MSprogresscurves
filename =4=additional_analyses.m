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


%% importing reported Km values and plotting histogram
file = {analysis_path+"CODETABLES.xlsx"}; %path to overview
[Km_A,Km_B] = xlsread(file{1},"Km_values","I4 : I23");
figure(2)
histogram(Km_A,34,'FaceColor',[.5 .5 .5])
ylabel("Number of enzymes")
xlabel("Highest measured KM among substrates [uM]")


