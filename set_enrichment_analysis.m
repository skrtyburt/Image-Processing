function set_enrichment_analysis(group1,group2,rois,annotations,filename)
%% ---- Compute Set Enrichment Analysis on Regional SUVR Data Between Two Groups ----
% INPUTS:
%           group1:         tier1 connectomic output from covnet_workflow.m, should contain experimental group data
%           group2:         tier1 connectomic output from covnet_workflow.m, should contain control group data
%           rois:           regions of interest, must have same row
%                           dimension as input data, stored in a .mat file
%                           with ROI variable named 'ROIs'
%           annotations:    functional annotations weighted for each region, stored in a .mat file with variable 'reg2func' 
%           filename:       name of the output file that you choose
%                           Ex: '<Study>_Functional_Determination.mat'                            

%% Version Control: Charles Burton, IU School of Medicine, 2024

%% Load in tier1 connectomic output data
exp = load(group1);
ctrl = load(group2);
load(annotations,annotations);
load(rois); %#ok<*LOAD>
beta = array2table(annotations);
location = pwd;
outlocation = strcat(location,'/',filename);
%% Pull out the average value of every region for each cohort
[rw1,cl1] = size(exp.cellData());
[rw2,cl2] = size(ctrl.cellData());

data_in_exp = exp.cellData(2:rw1,2:cl1);
data_in_ctrl = ctrl.cellData(2:rw2,2:cl2);

[exp_factor1,exp_factor2] = size(data_in_exp);
[ctrl_factor1,ctrl_factor2] = size(data_in_ctrl);

%pull out specific regional SUVR value and either sum or average it across
%each cohort for group 1
Mean_Region_Exp={};
for i=1:length(exp.cNames)
    Mean_Region_Exp{1,i} = exp.cNames(1,i); %#ok<*AGROW>
end

for i=1:length(exp.rNames)
    Mean_Region_Exp{i,1} = exp.rNames(i,1);
end

for i=1:exp_factor1
    for j=1:exp_factor2
        cohort_avg = mean(data_in_exp{i,j},2);
        Mean_Region_Exp{i+1,j+1} = cohort_avg;
    end
end

% Do the same thing for group 2
Mean_Region_Ctrl={};
for i=1:length(ctrl.cNames)
    Mean_Region_Ctrl{1,i} = ctrl.cNames(1,i);
end

for i=1:length(ctrl.rNames)
    Mean_Region_Ctrl{i,1} = ctrl.rNames(i,1);
end

for i=1:ctrl_factor1
    for j=1:ctrl_factor2
        cohort_avg = mean(data_in_ctrl{i,j},2);
        Mean_Region_Ctrl{i+1,j+1} = cohort_avg;
    end
end

%% Create a ranked list of differential regional SUVR between two comparison groups
% NOTE: We want to be able to specify this from the input in some capacity
% or else you will have to come in here and specify your control cell(s)
% for every run
control_group_F = Mean_Region_Ctrl{2,3};
control_group_M = Mean_Region_Ctrl{2,2};
difference={};
for i=2:rw1
    diff = Mean_Region_Exp{i,2}-control_group_M;
    difference{i,2} = diff;
end
for j=2:rw1
    diff = Mean_Region_Exp{j,3}-control_group_F;
    difference{j,3} = diff;
end
[rw,~] = size(difference);

%% Take our differential scoring, attach MRCC data, region label, z-score, beta coefficient, and organize in ascending order
%compute the z-score of each regional average
Z = {};
for i=2:rw
    zm = zscore(Mean_Region_Exp{i,2});
    zf = zscore(Mean_Region_Exp{i,3});
    Z{i,2} = zm;
    Z{i,3} = zf;
end

%add in regional values to each element in the vectors of DiffSUVR
sorted_diff = {};
diff_text = Mean_Region_Exp{1,1} + " REF " + Mean_Region_Ctrl{1,1};
sorted_diff{1,1} = diff_text;
for i=2:length(exp.cNames)
    sorted_diff{1,i} = exp.cNames(1,i);
end

for i=2:length(exp.rNames)
    sorted_diff{i,1} = exp.rNames(i,1);
end
% Add in MRCC Module Data to the tables for all Cohorts, sort SUVR in ascending order
header = ["ROIs" "Differential Value" "Module" "Z-Score" "Auditory" "Learning" "Motor" "Perception" "Sensory" "Integration" "Visual"];
for i=2:rw
    diff_table1 = table(ROIs,difference{i,2},exp.mrccPartition{i,2},Z{i,2});
    diff_table2 = table(ROIs,difference{i,3},exp.mrccPartition{i,3},Z{i,3});
    diff_table1 = [diff_table1 beta];
    diff_table2 = [diff_table2 beta];
    diff_table1.Properties.VariableNames = header;
    diff_table2.Properties.VariableNames = header;
    sortit1 = sortrows(diff_table1,2,'descend');
    sortit2 = sortrows(diff_table2,2,'descend');
    sorted_diff{i,2} = sortit1;
    sorted_diff{i,3} = sortit2;
end

%% Optional (For validation of automated quantification): Write out cohort tables into Excel spreadsheets
% [rw,cl] = size(sorted_diff);
% for i=2:rw
%     for j=2:cl
%         fileout = ['regional ontology data ' char(sorted_diff{1,1}) ' ' char(sorted_diff{i,1}) ' ' char(sorted_diff{1,j}) '.xlsx'];
%         writetable(sorted_diff{i,j},fileout, 'FileType','spreadsheet');
%     end
% end

%% Separate modules within each cohort
% generate nested cell array for every cohort containing sorted vectors by
% module

[rw,cl]=size(sorted_diff);
module_sorted = {};
complement_sorted = {};

% Title generation
module_text = "Sorted Modules of " + Mean_Region_Exp{1,1};   
complement_text = "Complement of " + Mean_Region_Exp{1,1};
module_sorted{1,1} = module_text;
complement_sorted{1,1} = complement_text;

% Naming convention for nested cell array
names = {};
for i=2:rw
    for j=2:cl
        row_name = sorted_diff{i,1};
        col_name = sorted_diff{1,j};
        label = string(row_name) + " " + string(col_name);
        names{end+1} = label;
    end
end
names = names.';
[rw_names,~] = size(names);
for i=1:rw_names
    module_sorted{i+1,1} = names{i,1};
    complement_sorted{i+1,1} = names{i,1};
end

rois_in_module = [];
rois_in_complement = [];
num = 0;
for i=2:rw
    for j=2:cl
        num = num + 1; 
        cohort=sorted_diff{i,j};
        cohort_data = table2array(cohort(:,2:11));
        [rois,~] = size(cohort);
        modules = table2array(range(cohort(:,3))+1);
        for l=1:modules
            for k=1:rois
                if cohort_data(k,2) == l
                    region = cohort(k,:);
                    rois_in_module = [rois_in_module; region];
                else
                    region = cohort(k,:);
                    rois_in_complement = [rois_in_complement; region];
                end
            end
            module_sorted{num+1,l+1} = rois_in_module;
            complement_sorted{num+1,l+1} = rois_in_complement;
            rois_in_module = [];
            rois_in_complement = [];
        end
    end
end
[~,cl] = size(module_sorted);

% column label generation
for i=2:cl
    name = "Module " + string(i-1);
    module_sorted{1,i} = name;
    complement_sorted{1,i} = name;
end

%% Compute the functional score of each module for all modules in a cohort
[rw,cl] = size(module_sorted);
[rw_complement,cl_complement] = size(complement_sorted);

fcn_scores = {};
fcn_scores_complement = {};
module_text = "Functional Scoring of " + Mean_Region_Exp{1,1} + " REF " + Mean_Region_Ctrl{1,1};
fcn_scores{1,1} = module_text;
fcn_scores_complement{1,1} = module_text;

% row label generation %
for i=1:rw_names
    fcn_scores{i+1,1} = names{i,1};
    fcn_scores_complement{i+1,1} = names{i,1};
end

% column label generation
for i=2:cl
    name = "Module " + string(i-1);
    fcn_scores{1,i} = name;
    fcn_scores_complement{1,i} = name;
end

%function score computation for the module
for i=2:rw
    for j=2:cl
        %num = num+1;
        auditory = 0;
        learning = 0;
        motor = 0;
        perception = 0;
        sensory = 0;
        integration = 0;
        visual = 0;
        tmod = module_sorted{i,j};
        if isempty(tmod) == true
            continue
        else
            size(tmod);
            module = table2cell(tmod);
            module = cell2mat(module(:,2:11));
            [mrw,~] = size(module);

            if mrw == 0
                continue
            else
                for k=1:mrw
                    auditory_score = ((module(k,1)*module(k,3))^2)*module(k,4);
                    learning_score = ((module(k,1)*module(k,3))^2)*module(k,5);
                    motor_score = ((module(k,1)*module(k,3))^2)*module(k,6);
                    perception_score = ((module(k,1)*module(k,3))^2)*module(k,7);
                    sensory_score = ((module(k,1)*module(k,3))^2)*module(k,8);
                    integration_score = ((module(k,1)*module(k,3))^2)*module(k,9);
                    visual_score = ((module(k,1)*module(k,3))^2)*module(k,10);
                    auditory = auditory+auditory_score;
                    learning = learning+learning_score;
                    motor = motor+motor_score;
                    perception = perception+perception_score;
                    sensory = sensory+sensory_score;
                    integration = integration+integration_score;
                    visual = visual+visual_score;
                end
            end
        set_of_scores = ["auditory" auditory; "learning" learning; "motor" motor; "perception" perception; "sensory" sensory; "integration" integration; "visual" visual];
        fcn_scores{i,j} = set_of_scores;
        auditory = 0; %#ok<*NASGU>
        learning = 0;
        motor = 0;
        perception = 0;
        sensory = 0;
        integration = 0;
        visual = 0;
        end
    end
end

%function score computation for the module's complement
for i=2:rw_complement
    for j=2:cl_complement
        %num = num+1;
        auditory = 0;
        learning = 0;
        motor = 0;
        perception = 0;
        sensory = 0;
        integration = 0;
        visual = 0;
        tmod = complement_sorted{i,j};
        if isempty(tmod) == true
            continue
        else
            size(tmod);
            complement = table2cell(tmod);
            complement = cell2mat(complement(:,2:11));
            [mrw,~] = size(complement);

            if mrw == 0
                continue
            else
                for k=1:mrw
                    auditory_score = ((complement(k,1)*complement(k,3))^2)*complement(k,4);
                    learning_score = ((complement(k,1)*complement(k,3))^2)*complement(k,5);
                    motor_score = ((complement(k,1)*complement(k,3))^2)*complement(k,6);
                    perception_score = ((complement(k,1)*complement(k,3))^2)*complement(k,7);
                    sensory_score = ((complement(k,1)*complement(k,3))^2)*complement(k,8);
                    integration_score = ((complement(k,1)*complement(k,3))^2)*complement(k,9);
                    visual_score = ((complement(k,1)*complement(k,3))^2)*complement(k,10);
                    auditory = auditory+auditory_score;
                    learning = learning+learning_score;
                    motor = motor+motor_score;
                    perception = perception+perception_score;
                    sensory = sensory+sensory_score;
                    integration = integration+integration_score;
                    visual = visual+visual_score;
                end
            end
        set_of_scores_complement = ["auditory" auditory; "learning" learning; "motor" motor; "perception" perception; "sensory" sensory; "integration" integration; "visual" visual];
        fcn_scores_complement{i,j} = set_of_scores_complement;
        auditory = 0;
        learning = 0;
        motor = 0;
        perception = 0;
        sensory = 0;
        integration = 0;
        visual = 0;
        end
    end
end

%% Significance Testing
% Take the differential list distribution as control distribution,
% differential of each region should be attached to all of it's
% corresponding functions, create sub-distributions for all functions

%Pass 1: Ungraded assignment of differential value to function according to
%region
[rw,cl]=size(sorted_diff);
functional = ["Auditory" "Learning" "Motor" "Perception" "Sensory" "Integration" "Visual"];
globalfunc = {};
globalval = [];
for i=2:rw 
    for j=2:cl
        module = sorted_diff{i,j}; %for every cell element in sorted_diff
        module2 = table2array(module(:,2:11)); %change the table to an array, omitting regional labels 
        %store all differential values according to function
        for k=4:10 %for all functional groups
            for m=1:27 %for every region
                if module2(m,k) == 0 %if regional value is zero (no functional value in a region)
                    diff_val = NaN; %diff_val variable is NaN
                else
                    diff_val = module2(m,1); %otherwise, diff_val is the nonzero regional value for a function
                end
                globalval(m,k-3) = diff_val; %globalval stores all of the diff_val values
            end
        end
        globalval = array2table(globalval);
        globalval.Properties.VariableNames = functional;
        globalfunc{i,j} = globalval;
        globalval = [];
    end
end

for i=2:length(exp.cNames)
    globalfunc{1,i} = exp.cNames(1,i);
end

for i=2:length(exp.rNames)
    globalfunc{i,1} = exp.rNames(i,1);
end

[rw,cl] = size(globalfunc);
vectorized = [];
house_of_vectors={};

n=0;
for i=2:rw %for every row in the globalfunc array
    for j=2:cl %for every column in the globalfunc array
        func_dist = globalfunc{i,j}; %opens the table corresponding to functional allocation of some cohort
        func_dist = table2array(func_dist); %converts table to array
        [funcrw,funccl] = size(func_dist); %grabs size of array
        n=n+1; %counter
        for k=1:funccl %for every function
            single_func = func_dist(:,k); %grabs the column containing functional allocation data
            for l=1:funcrw %for every row in the column
                if isnan(single_func(l,1)) == true %if the value is NaN
                    continue %skip
                else
                    single_func1 = single_func(l,1); % pulls nonzero functional value
                    vectorized = [vectorized single_func1]; % appends nonzero functional value to a vector called vectorized
                end
            end
            house_of_vectors{n+1,k+1} = vectorized;
            vectorized = [];
        end
    end
end

% concatenate the differential to the end of each cohort (@ the end of the
% row)
[rw,cl] = size(difference);
cohort_diffs = [];
for i=2:rw
    for j=2:cl
        cohort = difference{i,j};
        cohort_diffs = [cohort_diffs cohort];
    end
end

[rw,cl] = size(house_of_vectors);
for i=2:rw
    cohort_diff = cohort_diffs(:,i-1);
    house_of_vectors{i,9} = cohort_diff;
end


% naming the column, row headers of house of vectors
for i=2:cl
    house_of_vectors{1,i} = functional(i-1);
end
for i=2:rw
    rowlabel = module_sorted{i,1};
    house_of_vectors{i,1} = rowlabel;
end
house_of_vectors{1,9} = "Cohort Differentials";
house_of_vectors{1,1} = "Differential Mapping of Global Function";

%% Use module_sorted and complement_sorted to run 2-Sample KS-Test with module = set, complement = subset 
% Run for every functional distribution for every module
mod_fcnl_dist = {};
comp_fcnl_dist = {};

%% Compute the functional distribution of each module for all modules in a cohort
[rw,cl] = size(module_sorted);
[rw_complement,cl_complement] = size(complement_sorted);

module_text = "Functional Distribution of " + Mean_Region_Exp{1,1} + " REF " + Mean_Region_Ctrl{1,1};
mod_fcnl_dist{1,1} = module_text;
module_complement_text = "Complment Functional Distribution of " + Mean_Region_Exp{1,1} + " REF " + Mean_Region_Ctrl{1,1};
comp_fcnl_dist{1,1} = module_text;

% row label generation %
for i=1:rw_names
    mod_fcnl_dist{i+1,1} = names{i,1};
    comp_fcnl_dist{i+1,1} = names{i,1};
end

% column label generation
for i=2:cl
    name = "Module " + string(i-1);
    mod_fcnl_dist{1,i} = name;
    comp_fcnl_dist{1,i} = name;
end

%compute the function distribution for each module (including zeros)
%for every module
for i=2:rw
    for j=2:cl
        %create lists of each functional group
        auditory = [];
        learning = [];
        motor = [];
        perception = [];
        sensory = [];
        integration = [];
        visual = [];
        tmod = module_sorted{i,j};

        %if the module doesn't exist, pass on to the next cell
        if isempty(tmod) == true
            continue
        else
            % pull data from the module excluding region names
            module = table2cell(tmod);
            module = cell2mat(module(:,2:11));
            % determine the number of regions within the module
            [mrw,~] = size(module);
            % if there are no regions, move on to the next module
            if mrw == 0
                continue
            else
                % for every region
                for k=1:mrw
                    auditory_score = ((module(k,1)*module(k,3))^2)*module(k,4);
                    learning_score = ((module(k,1)*module(k,3))^2)*module(k,5);
                    motor_score = ((module(k,1)*module(k,3))^2)*module(k,6);
                    perception_score = ((module(k,1)*module(k,3))^2)*module(k,7);
                    sensory_score = ((module(k,1)*module(k,3))^2)*module(k,8);
                    integration_score = ((module(k,1)*module(k,3))^2)*module(k,9);
                    visual_score = ((module(k,1)*module(k,3))^2)*module(k,10);
                    auditory(k,1) = auditory_score;
                    learning(k,1) = learning_score;
                    motor(k,1) = motor_score;
                    perception(k,1) = perception_score;
                    sensory(k,1) = sensory_score;
                    integration(k,1) = integration_score;
                    visual(k,1) = visual_score;
                end
            end
        set_of_distributions = [auditory learning motor perception sensory integration visual];
        mod_fcnl_dist{i,j} = set_of_distributions;
        auditory = [];
        learning = [];
        motor = [];
        perception = [];
        sensory = [];
        integration = [];
        visual = [];
        end
    end
end

%function score computation for the module's complement
for i=2:rw_complement
    for j=2:cl_complement
        %create lists of each functional group
        auditory = [];
        learning = [];
        motor = [];
        perception = [];
        sensory = [];
        integration = [];
        visual = [];
        tmod = complement_sorted{i,j};

        %if the module doesn't exist, pass on to the next cell
        if isempty(tmod) == true
            continue
        else
            % pull data from the module excluding region names
            complement = table2cell(tmod);
            complement = cell2mat(complement(:,2:11));
            [mrw,~] = size(complement);

            if mrw == 0
                continue
            else
                for k=1:mrw
                    auditory_score = ((complement(k,1)*complement(k,3))^2)*complement(k,4);
                    learning_score = ((complement(k,1)*complement(k,3))^2)*complement(k,5);
                    motor_score = ((complement(k,1)*complement(k,3))^2)*complement(k,6);
                    perception_score = ((complement(k,1)*complement(k,3))^2)*complement(k,7);
                    sensory_score = ((complement(k,1)*complement(k,3))^2)*complement(k,8);
                    integration_score = ((complement(k,1)*complement(k,3))^2)*complement(k,9);
                    visual_score = ((complement(k,1)*complement(k,3))^2)*complement(k,10);
                    auditory(k,1) = auditory_score;
                    learning(k,1) = learning_score;
                    motor(k,1) = motor_score;
                    perception(k,1) = perception_score;
                    sensory(k,1) = sensory_score;
                    integration(k,1) = integration_score;
                    visual(k,1) = visual_score;
                end
                set_of_distributions_complement = [auditory learning motor perception sensory integration visual];
                comp_fcnl_dist{i,j} = set_of_distributions_complement;
                auditory = 0;
                learning = 0;
                motor = 0;
                perception = 0;
                sensory = 0;
                integration = 0;
                visual = 0;
            end 
        end
    end
end

%% Use functional arrays to run 2-Sample KS-Test with set and subset at input variables 
% for every element in the house_of_vectors with n>2 elements, compare against its
% respective differential distribution

% compare the distribution of the functional scores in the enriched set
% (module) compared to the complement set (everything outside of the module
[rw,cl] = size(mod_fcnl_dist);
ks_complement = {};
for i=2:rw
    for j = 2:cl
        module = mod_fcnl_dist{i,j};
        module_comp = comp_fcnl_dist{i,j};
        [~,functions] = size(module);
        if isempty(module) == true
            continue
        else
            for k=1:functions
                mod = module(:,k);
                comp = module_comp(:,k);
                [h,p,ks2stat] = kstest2(mod,comp);
                output = [h p ks2stat];
                ks_complement{i,j} = output;
            end
        end
    end
end
% Assign header labels to rows/columns

[~,cl] = size(module_sorted);

ks_text = "Two-Sample Complement KS Test";
ks_complement{1,1} = ks_text;

% row label generation %
for i=1:rw_names
   ks_complement{i+1,1} = names{i,1};
end

% column label generation
for i=2:cl
    name = "Module " + string(i-1);
    ks_complement{1,i} = name;
end

% run 2-sample KS-Test
[rw,cl] = size(module_sorted);
ks_global = {};
for i=2:rw
    for j=2:cl
        if length(house_of_vectors{i,j}) > 2 %this should be where the thresholding is applied, 2 is arbitrary and must be replaced by a more robust computation
            [h,p,ks2stat] = kstest2(house_of_vectors{i,j},house_of_vectors{i,9});
            output = [h p ks2stat];
            ks_global{i,j} = output;
        else
            ks_global{i,j} = "Not enough data";
        end
    end
end

% naming the column, row headers of ks_output
for i=2:cl
    ks_global{1,i} = functional(i-1);
end
for i=2:rw
    rowlabel = module_sorted{i,1};
    ks_global{i,1} = rowlabel;
end
ks_global{1,1} = "Two-sample Global KS Test";

%% Write out to .mat
save(outlocation, "ks_global","ks_complement","exp","ctrl")
end