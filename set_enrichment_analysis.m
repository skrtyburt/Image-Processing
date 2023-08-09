function [outputval1,outputval2,difference] = set_enrichment_analysis(group1,group2)
%% ---- Compute Set Enrichment Analysis on Regional SUVR Data Between Two Groups ---- %%
% INPUTS:
%           group1:         tier1 connectomic output from covnet_workflow.m, should contain control group data
%           group2:         tier1 connectomic output from covnet_workflow.m, should contain experimental group data
% 

%% load in data from connectomic outputs %%
g1 = load(group1);
g2 = load(group2);

%% Pull out the average value of every region for each cohort %%
[rw1,cl1] = size(g1.cellData());
[rw2,cl2] = size(g2.cellData());

data_in_g1 = g1.cellData(2:rw1,2:cl1);
data_in_g2 = g2.cellData(2:rw2,2:cl2);

[g1_factor1,g1_factor2] = size(data_in_g1);
[g2_factor1,g2_factor2] = size(data_in_g2);

%pull out specific regional SUVR value and either sum or average it across
%each cohort for group 1
Mean_Region_G1{10,10} = [] ;
for i=1:length(g1.cNames)
    Mean_Region_G1{1,i} = g1.cNames(1,i);
end

for i=1:length(g1.rNames)
    Mean_Region_G1{i,1} = g1.rNames(i,1);
end

for i=1:g1_factor1
    for j=1:g1_factor2
        cohort_avg = mean(data_in_g1{i,j},2);
        Mean_Region_G1{i+1,j+1} = cohort_avg;
    end
end
outputval1 = Mean_Region_G1;

% Do the same thing for group 2%
Mean_Region_G2{10,10} = [];
for i=1:length(g2.cNames)
    Mean_Region_G2{1,i} = g2.cNames(1,i);
end

for i=1:length(g2.rNames)
    Mean_Region_G2{i,1} = g2.rNames(i,1);
end

for i=1:g2_factor1
    for j=1:g2_factor2
        cohort_avg = mean(data_in_g2{i,j},2);
        Mean_Region_G2{i+1,j+1} = cohort_avg;
    end
end
outputval2 = Mean_Region_G2;

% %% Write out average region values of each cohort to Excel files (Optional) %%
% 
% fileout1 = 'Average Value of Regions for LEV';
% fileout2 = 'Average Value of Regions for 5XFAD';
% for i=1:g1_factor1
%     for j=1:g1_factor2
%         writematrix(Mean_Region_G1{i+1,j+1},[fileout1 '_' char(Mean_Region_G1{i+1,1}) '_' char(Mean_Region_G1{1,j+1})], 'FileType','spreadsheet')
%     end
% end
% for i=1:g2_factor1
%     for j=1:g2_factor2
%         writematrix(Mean_Region_G2{i+1,j+1},[fileout2 ' ' char(Mean_Region_G2{i+1,1}) ' ' char(Mean_Region_G2{1,j+1})], 'FileType', 'spreadsheet')
%     end
%% Create a ranked list of differential regional SUVR between two comparison groups %%
control_group_F = Mean_Region_G2{3,2};
control_group_M = Mean_Region_G2{3,3};
DiffSUVR{10,10} = [];
for i=2:rw1
    diff = Mean_Region_G1{i,2}-control_group_M;
    DiffSUVR{i,2} = diff;
end
for j=2:rw1
    diff = Mean_Region_G1{i,3}-control_group_F;
    DiffSUVR{j,3} = diff;
end

DiffSUVR{1,1} = 'Differential SUVR Output';
for i=2:length(g1.cNames)
    DiffSUVR{1,i} = g1.cNames(1,i);
end

for i=2:length(g1.rNames)
    DiffSUVR{i,1} = g1.rNames(i,1);
end

difference = DiffSUVR;
%% Take our differential scoring and organize in ascending order %%
% Need to Make sure that we still know what regions correspond with what
% values, or we can write out all nested vectors into Excel files%


end
