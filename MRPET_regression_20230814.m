%% Regression with PET data

clc;clear

IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];


% data = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPsmoothed_and_Memory_bothSession_03-Aug-2023_imputed_outlierRemoved.xlsx');
data = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPsmoothed_and_Memory_bothSession_09-Aug-2023_imputed_outlierRemoved.xlsx');

load('/Users/alex/Dropbox/paperwriting/MRPET/data/meanTval_hipp_20230809.mat');
load('/Users/alex/Dropbox/paperwriting/MRPET/data/AppliedActivity.mat');

% New variable to append
% BPchange_Hippocampus_highDA = BPchange_hipp(:,1);
% BPchange_Hippocampus_lowDA = BPchange_hipp(:,2);
% BPchange_HPCcluster_highDA = BPchange_hipp_cluster(:,1);%BPchange_hipp(:,1);
% BPchange_HPCcluster_lowDA = BPchange_hipp_cluster(:,2);%BPchange_hipp(:,2);
% BPlpnt_HPCcluster_highDA = BP_lpnt_hipp_cluster(:,1);%BPchange_hipp(:,1);
% BPlpnt_HPCcluster_lowDA = BP_lpnt_hipp_cluster(:,2);%BPchange_hipp(:,2);
meanTval_HPCcluster_highDA = meanTval_hipp(:,1);
meanTval_HPCcluster_lowDA = meanTval_hipp(:,2);
AppliedActivity_highDA = AppliedActivity(:,2);
AppliedActivity_lowDA = AppliedActivity(:,3);

% data = addvars(data, AppliedActivity_highDA, 'NewVariableNames', 'AppliedActivity_highDA');
% data = addvars(data, AppliedActivity_lowDA, 'NewVariableNames', 'AppliedActivity_lowDA');


for foldcomments=1
% Set outlier detection threshold
threshold = 3; % Number of standard deviations from the mean

% Loop over each column
for i = 1:size(data, 2)
    col = data{:, i};
    
    % Calculate mean and standard deviation
    mu = nanmean(col);
    sigma = nanstd(col);
    
    % Find indices of outliers
    outliers = abs(col - mu) > threshold * sigma;
    
    % Replace outliers with NaN
    col(outliers) = NaN;
    
    % Update column in data table
    data{:, i} = col;
end
% Save cleaned data to new sheet in Excel file
%  writetable(ModelData, ['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_RegressionData_BPchangeLPNTsmoothed_and_Memory_bothSession_' date '_imputed_outlierRemoved.xlsx'], 'Sheet', 'Cleaned Data');
end

data = removevars(data, 'BPchange_SRTM_SubcorticalROIs_highDA');
data = removevars(data, 'BPchange_SRTM_SubcorticalROIs_lowDA');
data = removevars(data, 'BPchange_lpnt_SubcorticalROIs_highDA');
data = removevars(data, 'BPchange_lpnt_SubcorticalROIs_lowDA');

%% variables list

% Get variable names and number of samples
var_names = data.Properties.VariableNames';

AGEs=2;

% AppliedAct_highDA=787;
% AppliedAct_lowDA=788;

% OCCs_highDA=[453:463];
% OCCs_lowDA =[498:508];

BPchanges_SRTM_highDA=[3:12];
BPchanges_SRTM_lowDA=[227:236];

BPchanges_lpnt_highDA=[693:702];
BPchanges_lpnt_lowDA=[737:746];

BPchange_SRTM_allROIs_highDA = [3:46];
BPchange_SRTM_allROIs_lowDA = [227:270];


% BPchange_HPC_highDA=[5 697];
% BPchange_HPC_lowDA=[230 742];
% 
% BP_task_lpntPET_highDA=[138:148];
% BP_task_lpntPET_lowDA=[363:373];
% 
% BP_task_SRTM_highDA = [93:103];
% BP_task_SRTM_lowDA = [318:328];
% 
% BPbsl_SRTM_highDA = [48:58];
% BPbsl_SRTM_lowDA = [273:280 282:283];
% 
% Gamma_highDA=[183:193];
% Gamma_lowDA=[408:418];
% 
Dprimes_rew_highDA=[581:582];
Dprimes_rew_lowDA=[583:584];
Dprimes_rew_highconf_highDA=[589:590];
Dprimes_rew_highconf_lowDA=[591:592];

Dprimes_neu_highDA=[585:586];
Dprimes_neu_lowDA=[587:588];
Dprimes_neu_highconf_highDA=[593:594];
Dprimes_neu_highconf_lowDA=[595:596];

% Dprimes_highDA = [579:580];
% Dprimes_lowDA  = [581:582];
% 
% hit_highDA = [555:556];
% hit_lowDA  = [557:558];
% FAs_highDA = [559:560];
% FAs_lowDA  = [561:562];
% 
% hits_rew_highDA = [563:564];
% hits_rew_lowDA  = [565:566];
% FAs_rew_highDA  = [567:568];
% FAs_rew_lowDA   = [569:570];
% 

% ctrl_vars = [var_names(AGEs); var_names(AppliedAct_highDA)]; % for baseline measures
ctrl_vars = [var_names(AGEs)];
var_selected = [var_names(BPchanges_SRTM_highDA); var_names(BPchanges_SRTM_lowDA)];

n = height(data);

%% run PLS regression: BPchange

for FoldPLSscripts=1

ROInames_highDA = cellfun(@(x) x(15:(end-7)),var_selected(1:10),'Uniformoutput',false);
ROInames_lowDA = cellfun(@(x) x(15:(end-6)),var_selected(11:end),'Uniformoutput',false);

% Define ROIs and data
Data_BP_all_highDA= data(:, var_selected(1:10));
Data_BP_all_lowDA= data(:, var_selected(11:end));

Data_BP_all_highDA.Properties.VariableNames = cellfun(@(x) x(1:(end-7)),var_selected(1:10),'Uniformoutput',false);
Data_BP_all_lowDA.Properties.VariableNames = Data_BP_all_highDA.Properties.VariableNames;

Data_all = [Data_BP_all_highDA;Data_BP_all_lowDA];
SessionCondition = [ repmat({'HighDA'}, 33, 1);  repmat({'lowDA'}, 33, 1)];

roi = cellfun(@(x) x(1:(end-7)),var_selected(1:10),'Uniformoutput',false); % ... Add all ROIs here

% Create a sample data table (you will replace this with your real data)
ModelData = addvars(Data_all, SessionCondition, 'NewVariableNames', 'SessionCondition');
% ModelData = table({
%     'High' ; 'Low' ; 'High' ; 'Low'},... 
%     rand(4,1), rand(4,1), rand(4,1), 'VariableNames', {'SessionCondition', roi{:}});

alpha = 0.05;  % significance level

% Loop through each ROI and run the regression analysis
all_pvalues = [];  % to collect all the p-values across ROIs

for idx = 1:length(roi)
    disp('****************')
    % Extract the response variable and predictor from the data table
    response = ModelData.(roi{idx});
    SessionCondition = double(strcmp(ModelData.SessionCondition, 'HighDA'));
    
    % Create the linear regression model
    tbl = table(response, SessionCondition, 'VariableNames', {'response', 'SessionCondition'});
    mdl = fitlm(tbl, 'response ~ SessionCondition');
    
    
    % Collect p-values
    all_pvalues = [all_pvalues; mdl.Coefficients.pValue(2)];  % Assuming that SessionCondition is the second coefficient
    
    % Display the result
    disp(['Regression result for ' roi{idx} ':']);
    disp(mdl);
    disp('****************')
end

% FDR correction
adj_p = mafdr(all_pvalues, 'BHFDR', true);

% Decide which tests are significant after correction
significant_tests = adj_p < alpha;

% Display results after FDR correction
disp('Results after FDR correction:');
for idx = 1:length(roi)
    disp([roi{idx} ': ' num2str(significant_tests(idx))]);
end

%% paired-samples t-test

all_pvalues = []; % to collect all the p-values across ROIs

for idx = 1:length(roi)
    highDA_data = Data_BP_all_highDA.(roi{idx});
    lowDA_data = Data_BP_all_lowDA.(roi{idx});
    
    % Perform paired t-test
    [~, p] = ttest(highDA_data, lowDA_data);
    
    % Collect p-values
    all_pvalues = [all_pvalues; p];
end

% FDR correction
adj_p = mafdr(all_pvalues, 'BHFDR', true);

% Decide which tests are significant after correction
alpha = 0.05;  % significance level
significant_tests = adj_p < alpha;

% Display results after FDR correction
disp('Results after FDR correction:');
for idx = 1:length(roi)
    disp([roi{idx} ': p-value = ' num2str(adj_p(idx)) ', Significant: ' num2str(significant_tests(idx))]);
end

%% PLS attempts
%% partial least squares regression
%  doing this because our ROIs are super colinear

close all

% --- create dataset for the analyses

% convert SessionCondition to a numerical vector for PLSR
Y = double(strcmp(ModelData.SessionCondition, 'HighDA')); % HighDA as 1 and LowDA as 0

% extract ROI data as matrix for PLSR
X = table2array(ModelData(:, roi));

% standardize (mean-zero, unit variance) the data before PLSR
X_mean = nanmean(X);
X_std = nanstd(X);

% standardize the data
X_zscored = (X - X_mean) ./ X_std;

% remove rows with NaNs both X and Y
rowsToRemove = any(isnan(X_zscored), 2);
X_zscored(rowsToRemove, :) = [];
Y(rowsToRemove) = [];


% --- now run: determine the number of PLS components let's use 4
numComponents = 2;

[XL, YL, XS, YS, BETA, PCTVAR] = plsregress(X_zscored, Y, numComponents);

% PCTVAR contains the percent variance explained by each component for X and Y
explainedVarianceX = PCTVAR(1,:);  % Variance explained in X by each component
explainedVarianceY = PCTVAR(2,:);  % Variance explained in Y by each component

% predict Y using the PLS model
yfit = [ones(size(X_zscored,1),1) X_zscored] * BETA;


%% Create the boxplot
figure; grid on

% Create the boxplot
h = boxplot(yfit(:,1), Y);

% Set color properties
colors = {[0 0 1], [1 0 0]}; % Blue for Low DA, Red for High DA

for i = 1:2  % Assuming you have two groups (Low DA and High DA)
    patch(get(h(5,i), 'XData'), get(h(5,i), 'YData'), colors{i}, 'FaceAlpha',0.5);
end

% Adjust width of median line
set(h(6,:), 'LineWidth', 2.5);

% Label the median values
medians = findobj(h, 'tag', 'Median');
medianValues = get(medians, 'YData');
numMedians = length(medianValues);

for i = 1:numMedians
    text(i, medianValues{i}(1), num2str(medianValues{i}(1), '%.2f'), ...
         'VerticalAlignment','bottom', 'HorizontalAlignment','center');
end

% Adjust the rest of the plot aesthetics
set(gca, 'XTick', [1, 2], 'XTickLabel', {'LowDA', 'HighDA'});
ylabel('Predicted Values');
title('Spread of Predicted Values for Each Actual Condition');

% Examine which ROIs are driving the effect
disp('X-loadings (ROIs vs PLS components):');
disp(XL);


%% variance explained by each component

% Display the variance explained:
fprintf('Variance explained in X by each component: \n');
disp(explainedVarianceX);
fprintf('Variance explained in Y by each component: \n');
disp(explainedVarianceY);


% plotting component scores vs conditions
figure;

% Assign colors based on condition
colors = zeros(length(Y), 3); % RGB values
colors(Y == 0, :) = repmat([0 0 1], sum(Y == 0), 1); % Blue for Low DA
colors(Y == 1, :) = repmat([1 0 0], sum(Y == 1), 1); % Red for High DA

h(1) = scatter(Y(Y == 0), XS(Y == 0,1), 50, 'b', 'filled', 'DisplayName', 'Low DA');
hold on; % Ensure the next scatter plot is added to the same figure
h(2) = scatter(Y(Y == 1), XS(Y == 1,1), 50, 'r', 'filled', 'DisplayName', 'High DA');

xlabel('Actual Session Condition');
ylabel('Scores for Component 1');
title('Scores for Component 1 vs. Actual Session Condition');
xlim([-1 2]);
legend(h, 'Location', 'Best'); % Use the handles to define the legend
grid on


% plot against each compotnents
figure;
gscatter(YS(:,1), YS(:,2), Y, 'br', 'oo', 10);
xlabel('PLS Component 1 Scores');
ylabel('PLS Component 2 Scores');
title('Scores on PLS Component 1 vs Component 2');
legend({'Low DA', 'High DA'}, 'Location', 'Best');
xlim([-2.5 2.5])
grid on;


%% plot for visual inspection of component scores in each ROIs
figure;

% Component 1
subplot(1,2,1);
bar(XL(:,1));
title('Loadings for PLS Component 1');
xlabel('ROIs');
ylabel('Loading Value');
xticks(1:length(roi));
xticklabels(roi);
xtickangle(45);
ylim([-6 7])

% Add variance explained as text box for Component 1
varExplained1 = explainedVarianceY(1);  % Example value. Replace with your actual value.
eq_str1 = ['Variance explained by Y: ' num2str(varExplained1*100) '%'];
dim1 = [0.1, 0.8, 0.3, 0.1];
annotation('textbox', dim1, 'String', eq_str1, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'black', 'FaceAlpha', 0.5);

% Component 2
subplot(1,2,2);
bar(XL(:,2));
title('Loadings for PLS Component 2');
xlabel('ROIs');
ylabel('Loading Value');
xticks(1:length(roi));
xticklabels(roi);
xtickangle(45);
ylim([-6 7])

% Add variance explained as text box for Component 2
varExplained2 = explainedVarianceY(2);  % Example value. Replace with your actual value.
eq_str2 = ['Variance explained by Y: ' num2str(varExplained2*100) '%'];
dim2 = [0.6, 0.8, 0.3, 0.1];
annotation('textbox', dim2, 'String', eq_str2, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'black', 'FaceAlpha', 0.5);

       
% 
% % Component 3
% subplot(2,2,3);
% bar(XL(:,3));
% title('Loadings for PLS Component 3');
% xlabel('ROIs');
% ylabel('Loading Value');
% xticks(1:length(roi));
% xticklabels(roi);
% xtickangle(45);
% ylim([-6 7])
% 
% % Component 4
% subplot(2,2,4);
% bar(XL(:,4));
% title('Loadings for PLS Component 4');
% xlabel('ROIs');
% ylabel('Loading Value');
% xticks(1:length(roi));
% xticklabels(roi);
% xtickangle(45);
% ylim([-6 7])


%% Distribution of scores on PLS Component 1 for each condition
% Create jitter for better visualization

figure;

jitterAmount = 0.3;
jitterValues = jitterAmount * (2 * rand(size(YS,1), 1) - 1);  % Values between [-jitterAmount, jitterAmount]

% Split scores based on the condition
scores_highDA = YS(Y == 1,1);
scores_lowDA = YS(Y == 0,1);

% Set transparency level (0 is completely transparent, 1 is opaque)
alphaValue = 0.5; 

% Add jitter to the Y-values for visualization
scatter(repmat(1, numel(scores_highDA), 1) + jitterValues(Y == 1), scores_highDA, 50, 'r', 'filled', 'MarkerFaceAlpha', alphaValue, 'MarkerEdgeAlpha', alphaValue);
hold on;
scatter(repmat(2, numel(scores_lowDA), 1) + jitterValues(Y == 0), scores_lowDA, 50, 'b', 'filled', 'MarkerFaceAlpha', alphaValue, 'MarkerEdgeAlpha', alphaValue);

% Compute mean for each condition and plot them
mean_highDA = mean(scores_highDA);
mean_lowDA = mean(scores_lowDA);

scatter([1, 2], [mean_highDA, mean_lowDA], 100, 'k', 'filled');

% Set appropriate labels and title
set(gca, 'XTick', [1, 2], 'XTickLabel', {'High DA', 'Low DA'});
ylabel('PLS Component 1 Scores');
title('Distribution of scores on PLS Component 1 for each condition');
legend({'High DA Subjects', 'Low DA Subjects', 'Mean'}, 'Location', 'Best');
text(1, mean_highDA + 0.2, sprintf('Mean = %.2f', mean_highDA), 'HorizontalAlignment', 'center');
text(2, mean_lowDA + 0.2, sprintf('Mean = %.2f', mean_lowDA), 'HorizontalAlignment', 'center');
hold off;
xlim([0.5 2.5])
ylim([-2.5 2.5])


%% among ROIs

close all

% Preparing data
% Convert table data to matrix for easier manipulation
mat_highDA = table2array(Data_BP_all_highDA);
mat_lowDA = table2array(Data_BP_all_lowDA);

% Remove rows with NaN or Inf for each dataset
mat_highDA(any(isnan(mat_highDA), 2), :) = [];
mat_lowDA(any(isnan(mat_lowDA), 2), :) = [];

% Concatenate for combined analysis
mat_combined = [mat_highDA; mat_lowDA];


% PLS for highDA condition
[XL_highDA,~,XS_highDA,~,~,PCTVAR_highDA] = plsregress(mat_highDA, mat_highDA, 3);

% PLS for lowDA condition
[XL_lowDA,~,XS_lowDA,~,~,PCTVAR_lowDA] = plsregress(mat_lowDA, mat_lowDA, 3);

% PLS for combined data
[XL_combined,~,XS_combined,~,~,PCTVAR_combined] = plsregress(mat_combined, mat_combined, 3);

% Visualization
% Visualize Loadings for High DA
figure;
subplot(1,3,1);
bar(XL_highDA(:,1));
title('Loadings for PLS Component 1 (High DA)');
xlabel('ROIs');
ylabel('Loading Value');
xticks(1:length(ROInames_highDA));
xticklabels(ROInames_highDA);
xtickangle(45);

% Visualize Loadings for Low DA
subplot(1,3,2);
bar(XL_lowDA(:,1));
title('Loadings for PLS Component 1 (Low DA)');
xlabel('ROIs');
ylabel('Loading Value');
xticks(1:length(ROInames_lowDA));
xticklabels(ROInames_lowDA);
xtickangle(45);

% Visualize Loadings for Combined data
subplot(1,3,3);
bar(XL_combined(:,1));
title('Loadings for PLS Component 1 (Combined)');
xlabel('ROIs');
ylabel('Loading Value');
xticks(1:length(ROInames_highDA));  % Assuming same ROIs for both conditions
xticklabels(ROInames_highDA);
xtickangle(45);

end

%% hierarchical clustering

for FoldHierarchicalClusteringscripts=1
 
%% set variables
var_names = data.Properties.VariableNames';

AGEs = 2;

BPchange_SRTM_allROIs_highDA = [3:46];
BPchange_SRTM_allROIs_lowDA = [227:270];

BPbsl_SRTM_allROIs_highDA = [47:90];
BPbsl_SRTM_allROIs_lowDA = [270:313];

Dprimes_rew_highDA=[581:582];
Dprimes_rew_lowDA=[583:584];
Dprimes_rew_highconf_highDA=[589:590];
Dprimes_rew_highconf_lowDA=[591:592];

Dprimes_neu_highDA=[585:586];
Dprimes_neu_lowDA=[587:588];
Dprimes_neu_highconf_highDA=[593:594];
Dprimes_neu_highconf_lowDA=[595:596];

ctrl_vars = [var_names(AGEs)];
    
%% not controlling for age

var_selected = [var_names(BPbsl_SRTM_allROIs_highDA); var_names(BPbsl_SRTM_allROIs_lowDA)];
    
close all

% Get variable names and number of samples
var_names = data.Properties.VariableNames';

AGEs=2;

BPchange_SRTM_allROIs_highDA = [3:46];
BPchange_SRTM_allROIs_lowDA = [227:270];

ctrl_vars = [var_names(AGEs)];
var_selected = [var_names(BPchange_SRTM_allROIs_highDA); var_names(BPchange_SRTM_allROIs_lowDA)];

n = height(data);

ROInames_highDA = cellfun(@(x) x(15:(end-7)),var_selected(1:44),'Uniformoutput',false);
ROInames_lowDA = cellfun(@(x) x(15:(end-6)),var_selected(45:end),'Uniformoutput',false);

% Define ROIs and data
Data_BP_all_highDA= data(:, var_selected(1:44));
Data_BP_all_lowDA= data(:, var_selected(45:end));

Data_BP_all_highDA.Properties.VariableNames = cellfun(@(x) x(1:(end-7)),var_selected(1:44),'Uniformoutput',false);
Data_BP_all_lowDA.Properties.VariableNames = Data_BP_all_highDA.Properties.VariableNames;

Data_all = [Data_BP_all_highDA;Data_BP_all_lowDA];

% Extracting data from your table
dataMatrix = Data_all.Variables;

% Calculate the pairwise distance between the ROIs (columns)
distanceMatrix = pdist(dataMatrix', 'euclidean');

% Use hierarchical clustering on the distance matrix
linkageTree = linkage(distanceMatrix, 'average');

% Create labels for the 44 ROIs
allLabels = Data_all.Properties.VariableNames;

% Plot the dendrogram
figure;
[H, T, outperm] = dendrogram(linkageTree, 0); % 0 ensures all nodes are shown
set(H, 'LineWidth', 1);
xticks(1:length(allLabels));
xticklabels(allLabels);
xtickangle(45);
xlabel('ROIs');
ylabel('Distance');
title('Hierarchical clustering of ROIs');


%% controlling for age: Among ROIs

var_selected = [var_names(BPbsl_SRTM_allROIs_highDA); var_names(BPbsl_SRTM_allROIs_lowDA)];

close all

% Get variable names and number of samples
ROInames_highDA = cellfun(@(x) x(15:(end-7)),var_selected(1:44),'Uniformoutput',false);
ROInames_lowDA = cellfun(@(x) x(15:(end-6)),var_selected(45:end),'Uniformoutput',false);

% Define ROIs and data
Data_BP_all_highDA = data(:, var_selected(1:44));
Data_BP_all_lowDA = data(:, var_selected(45:end));

Data_BP_all_highDA.Properties.VariableNames = cellfun(@(x) x(1:(end-7)),var_selected(1:44),'Uniformoutput',false);
Data_BP_all_lowDA.Properties.VariableNames = Data_BP_all_highDA.Properties.VariableNames;

Data_all = [Data_BP_all_highDA;Data_BP_all_lowDA];

% Extracting data from your table
dataMatrix = Data_all.Variables;

% Partial out the effect of age from each column of dataMatrix
ages = repmat(data{:, AGEs}, 2, 1); % Repeating the age values for both conditions

for i = 1:size(dataMatrix, 2)
    [b,~,residuals] = regress(dataMatrix(:,i), [ones(length(ages), 1), ages]);
    dataMatrix(:,i) = residuals;
end

% Calculate the pairwise distance between the ROIs (columns)
distanceMatrix = pdist(dataMatrix', 'euclidean');

% Use hierarchical clustering on the distance matrix
linkageTree = linkage(distanceMatrix, 'average');

% Create labels for the 44 ROIs
allLabels = cellfun(@(x) x(12:(end-7)),var_selected(1:44),'Uniformoutput',false);%Data_all.Properties.VariableNames;

% Plot the dendrogram
figure;
[H, T, outperm] = dendrogram(linkageTree, 0);
set(H, 'LineWidth', 1);
xticks(1:length(allLabels));
xticklabels(allLabels);
set(gca, 'TickLabelInterpreter', 'none');
set(H, 'LineWidth', 2);  % Set the line width to 2
ax = gca;  % Get handle to current axis
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;% Set the font size to 12 (or your desired size)
xtickangle(45);
xlabel('ROIs');
ylabel('Distance');
titleHandle=title(['Hierarchical clustering of ROIs (controlled for Age):' var_selected{1}(1:10) ],'FontSize',25);
set(titleHandle, 'Interpreter', 'none');


% labels of the distances
for k = 1:length(linkageValues)
    idx = find(T == k);
    if isempty(idx)
        continue;
    end
    x = mean(get(H(idx), 'XData'));
    y = linkageValues(k) - 0.0005; % subtracting a small value to shift the label down
    text(x, y, sprintf('%.2f', linkageValues(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8);
end


%% controlling for age: ROIs and Memory

close all

ctrl_vars = [var_names(AGEs)];
var_selected = [var_names(Dprimes_rew_highconf_lowDA); var_names(BPbsl_SRTM_allROIs_lowDA)];

% Define ROIs and data
Data_BPbslAll_DprimeRew_lowDA = data(:, var_selected);
Data_all = Data_BPbslAll_DprimeRew_lowDA;

% Extracting data from your table
dataMatrix = Data_all.Variables;

% Partial out the effect of age from each column of dataMatrix
ages = data{:, AGEs}; % Repeating the age values for both conditions

for i = 1:size(dataMatrix, 2)
    [b,~,residuals] = regress(dataMatrix(:,i), [ones(length(ages), 1), ages]);
    dataMatrix(:,i) = residuals;
end

% Calculate the pairwise distance between the ROIs (columns)
distanceMatrix = pdist(dataMatrix', 'euclidean');

% Use hierarchical clustering on the distance matrix
linkageTree = linkage(distanceMatrix, 'average');

% Create labels for the ROIs
MemoryLabels={'Dprime HighConf Reward Immediate',' Dprime HighConf Reward Delayed'};
ROIlabels = cellfun(@(x) x(12:(end-6)), var_selected(3:end),'Uniformoutput',false);
allLabels = [MemoryLabels ROIlabels'];

% Plot the dendrogram
figure;
[H, T, outperm] = dendrogram(linkageTree, 0);
set(H, 'LineWidth', 1);
xticks(1:length(allLabels));
xticklabels(allLabels);
set(gca, 'TickLabelInterpreter', 'none');
set(H, 'LineWidth', 2);  % Set the line width to 2
ax = gca;  % Get handle to current axis
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;% Set the font size to 12 (or your desired size)
xtickangle(45);
xlabel('ROIs');
ylabel('Distance');
titleHandle=title(['Hierarchical clustering of ROIs (controlled for Age): BP Baseline SRTM and Dprime HighConfidence Reward - High DA' ],'FontSize',25);
set(titleHandle, 'Interpreter', 'none');

% Extract linkage values
linkageValues = linkageTree(:, 3);

% labels of the distances
for k = 1:length(linkageValues)
    idx = find(T == k);
    if isempty(idx)
        continue;
    end
    x = mean(get(H(idx), 'XData'));
    y = linkageValues(k) - 0.0005; % subtracting a small value to shift the label down
    text(x, y, sprintf('%.2f', linkageValues(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8);
end


%% controlling for age: ROIs (BPchange) and Memory

close all

ctrl_vars = [var_names(AGEs)];
var_selected = [var_names(Dprimes_rew_highDA); var_names(BPbsl_SRTM_allROIs_lowDA)];

% Define ROIs and data
Data_BPbslAll_DprimeRew_lowDA = data(:, var_selected);
Data_all = Data_BPbslAll_DprimeRew_lowDA;

% Extracting data from your table
dataMatrix = Data_all.Variables;

% Partial out the effect of age from each column of dataMatrix
ages = data{:, AGEs}; % Repeating the age values for both conditions

for i = 1:size(dataMatrix, 2)
    [b,~,residuals] = regress(dataMatrix(:,i), [ones(length(ages), 1), ages]);
    dataMatrix(:,i) = residuals;
end

% Calculate the pairwise distance between the ROIs (columns)
distanceMatrix = pdist(dataMatrix', 'euclidean');

% Use hierarchical clustering on the distance matrix
linkageTree = linkage(distanceMatrix, 'average');

% Create labels for the ROIs
MemoryLabels={'Dprime HighConf Reward Immediate',' Dprime HighConf Reward Delayed'};
ROIlabels = cellfun(@(x) x(12:(end-6)), var_selected(3:end),'Uniformoutput',false);
allLabels = [MemoryLabels ROIlabels'];

% Plot the dendrogram
figure;
[H, T, outperm] = dendrogram(linkageTree, 0);
set(H, 'LineWidth', 1);
xticks(1:length(allLabels));
xticklabels(allLabels);
set(gca, 'TickLabelInterpreter', 'none');
set(H, 'LineWidth', 2);  % Set the line width to 2
ax = gca;  % Get handle to current axis
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;% Set the font size to 12 (or your desired size)
xtickangle(45);
xlabel('ROIs');
ylabel('Distance');
titleHandle=title(['Hierarchical clustering of ROIs (controlled for Age): BP Baseline SRTM and Dprime HighConfidence Reward - High DA' ],'FontSize',25);
set(titleHandle, 'Interpreter', 'none');

% Extract linkage values
linkageValues = linkageTree(:, 3);

% labels of the distances
for k = 1:length(linkageValues)
    idx = find(T == k);
    if isempty(idx)
        continue;
    end
    x = mean(get(H(idx), 'XData'));
    y = linkageValues(k) - 0.0005; % subtracting a small value to shift the label down
    text(x, y, sprintf('%.2f', linkageValues(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8);
end

end