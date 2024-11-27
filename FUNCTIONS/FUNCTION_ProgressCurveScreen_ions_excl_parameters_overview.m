function [ions_excl_parameters,ions_excl] = FUNCTION_ProgressCurveScreen_ions_excl_parameters_overview(overview_a,overview_labels,Eidx,fiadata_ctrls,fiadata_effs)

%exclude if ctrl mean(r2) much worse than for others...
%exclude if too many assays are excluded because of >100
ions_excl_parameters = fiadata_ctrls.ctrl_parameters_OVERVIEW;
break_idx = find(strcmp(fiadata_effs.impossible_fits(:,5),'control4'));
fin_idx = find(strcmp(fiadata_effs.impossible_fits(:,5),'control8'));
ions_excl_parameters(2,9) = {sum(cell2mat(fiadata_effs.impossible_fits(1:break_idx,1)) > 100)};
ions_excl_parameters(3,9) = {sum(cell2mat(fiadata_effs.impossible_fits(break_idx+1:fin_idx,1)) > 100)};
ions_excl_parameters(4,9) = {sum(cell2mat(fiadata_effs.impossible_fits(1:break_idx,2)) > 100)};
ions_excl_parameters(5,9) = {sum(cell2mat(fiadata_effs.impossible_fits(break_idx+1:fin_idx,2)) > 100)};
ions_excl_parameters(6,9) = {sum(cell2mat(fiadata_effs.impossible_fits(1:break_idx,3)) > 100)};
ions_excl_parameters(7,9) = {sum(cell2mat(fiadata_effs.impossible_fits(break_idx+1:fin_idx,3)) > 100)};
ions_excl_parameters(8,9) = {sum(cell2mat(fiadata_effs.impossible_fits(1:break_idx,4)) > 100)};
ions_excl_parameters(9,9) = {sum(cell2mat(fiadata_effs.impossible_fits(break_idx+1:fin_idx,4)) > 100)};

ions_excl_parameters(2,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(1:break_idx,1)))};
ions_excl_parameters(3,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,1)))};
ions_excl_parameters(4,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(1:break_idx,2)))};
ions_excl_parameters(5,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,2)))};
ions_excl_parameters(6,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(1:break_idx,3)))};
ions_excl_parameters(7,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,3)))};
ions_excl_parameters(8,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(1:break_idx,4)))};
ions_excl_parameters(9,10) = {nanmean(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,4)))};

ions_excl_parameters(2,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(1:break_idx,1)))};
ions_excl_parameters(3,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,1)))};
ions_excl_parameters(4,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(1:break_idx,2)))};
ions_excl_parameters(5,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,2)))};
ions_excl_parameters(6,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(1:break_idx,3)))};
ions_excl_parameters(7,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,3)))};
ions_excl_parameters(8,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(1:break_idx,4)))};
ions_excl_parameters(9,11) = {nanmedian(cell2mat(fiadata_effs.cvAUCs(break_idx+1:fin_idx,4)))};

S1_eligible = overview_a(Eidx,find(strcmp(overview_labels,"S1-eligible")));
S2_eligible = overview_a(Eidx,find(strcmp(overview_labels,"S2-eligible"))); 
P1_eligible = overview_a(Eidx,find(strcmp(overview_labels,"P1-eligible"))); 
P2_eligible = overview_a(Eidx,find(strcmp(overview_labels,"P2-eligible")));

ions_excl = [S1_eligible S2_eligible P1_eligible P2_eligible];