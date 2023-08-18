% Load Excel file with data

clc;clear

% IDs : try 40141 again!!!!!
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];


% data = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPsmoothed_and_Memory_bothSession_03-Aug-2023_imputed_outlierRemoved.xlsx');
data = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPsmoothed_and_Memory_bothSession_09-Aug-2023_imputed_outlierRemoved.xlsx');
data_memoryCollapsed = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_memoryTable_collapsed_loglinear_18-Aug-2023_imputed.xlsx');

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

% Append the new variable to the existing table
% data = addvars(data, BPchange_HPCcluster_highDA, 'NewVariableNames', 'BPchange_HPCcluster_highDA');
% data = addvars(data, BPchange_HPCcluster_lowDA, 'NewVariableNames', 'BPchange_HPCcluster_lowDA');

% data = addvars(data, meanTval_HPCcluster_highDA, 'NewVariableNames', 'meanTval_HPCcluster_highDA_2');
% data = addvars(data, meanTval_HPCcluster_lowDA, 'NewVariableNames', 'meanTval_HPCcluster_lowDA_2');

% data = addvars(data, AppliedActivity_highDA, 'NewVariableNames', 'AppliedActivity_highDA');
% data = addvars(data, AppliedActivity_lowDA, 'NewVariableNames', 'AppliedActivity_lowDA');

% data = addvars(data, BPlpnt_HPCcluster_highDA, 'NewVariableNames', 'BP_lpnt_HPCcluster_highDA');
% data = addvars(data, BPlpnt_HPCcluster_lowDA, 'NewVariableNames', 'BP_lpnt_HPCcluster_lowDA');

% data = addvars(data, BPchange_Hippocampus_highDA, 'NewVariableNames', 'BPchange_Hippocampus_highDA');
% data = addvars(data, BPchange_Hippocampus_lowDA, 'NewVariableNames', 'BPchange_Hippocampus_lowDA');


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
%  writetable(data, ['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPsmoothed_and_Memory_bothSession_' date '_imputed_outlierRemoved.xlsx'], 'Sheet', 'Cleaned Data');
end


%% fMRI

% load fMRI data
% hipp=spm_read_vols(spm_vol([paths_mask 'mni_icbm152_AAL3_hippocampus_bi.nii']));
hipp=spm_read_vols(spm_vol(['/Volumes/ALEX3/MRPET/analysis/eachSession/2ndLvL/STIM_rew_v_neu_p005_c370_HippClusterMask_conj.nii']));
fMRI=[];
for id=1:length(IDs)
    for d=1:2
        if days(id,d)==0
            fMRI.stim.rew_v_neu{id,d}=[];
        else
            fMRI.stim.rew_v_neu{id,d}=spm_read_vols(spm_vol(['/Volumes/ALEX3/MRPET/analysis/eachSession/' num2str(IDs(id)) '_' num2str(d) '/con_0002_mni.nii']));
        end
    end
end

% extract t-values within the hippocampus mask
idx_hipp=hipp>0;
Tval_hipp=[];
for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            Tval_hipp{id,d} = [];
            meanTval_hipp(id,d) = NaN;
        else
            Tval_hipp{id,d}     = fMRI.stim.rew_v_neu{id,d}(idx_hipp);
            meanTval_hipp(id,d) = nanmean(Tval_hipp{id,d}(:));
        end
    end
    
end
%% choose your variables

% Get variable names and number of samples
var_names = data.Properties.VariableNames';
var_names_mem = data_memoryCollapsed.Properties.VariableNames';

%% select variables

AGEs=2;

AppliedAct_highDA=787;
AppliedAct_lowDA=788;

OCCs_highDA=[453:463];
OCCs_lowDA =[498:508];

BPchanges_SRTM_highDA=[3:13];
BPchanges_SRTM_lowDA=[228:238];

BPchanges_lpnt_highDA=[695:705];
BPchanges_lpnt_lowDA=[740:750];

BPchange_HPC_highDA=[5 697];
BPchange_HPC_lowDA=[230 742];

BP_task_lpntPET_highDA=[138:148];
BP_task_lpntPET_lowDA=[363:373];

BP_task_SRTM_highDA = [93:103];
BP_task_SRTM_lowDA = [318:328];

BPbsl_SRTM_highDA = [48:58];
BPbsl_SRTM_lowDA = [273:280 282:283];

Gamma_highDA=[183:193];
Gamma_lowDA=[408:418];

Dprimes_rew_highDA=[583:584];
Dprimes_rew_lowDA=[585:586];
Dprimes_rew_highconf_highDA=[591:592];
Dprimes_rew_highconf_lowDA=[593:594];

Dprimes_highDA = [579:580];
Dprimes_lowDA  = [581:582];

DprimesAll_highDA = 2;
DprimesAll_lowDA  = 3;

DprimesAll_highconf_highDA = 4;
DprimesAll_highconf_lowDA = 5;

DprimesAll_rew_highDA = 6;
DprimesAll_neu_highDA = 7;

DprimesAll_rew_lowDA = 8;
DprimesAll_neu_lowDA = 9;

DprimesAll_highconf_rew_highDA = 10;
DprimesAll_highconf_neu_highDA = 11;

DprimesAll_highconf_rew_lowDA = 12;
DprimesAll_highconf_neu_lowDA = 13;

Dprimes_neu_highDA=[587 588];
Dprimes_neu_lowDA=[589 590];
Dprimes_neu_highconf_highDA=[595 596];
Dprimes_neu_highconf_lowDA=[597 598];

hit_highDA = [555:556];
hit_lowDA  = [557:558];
FAs_highDA = [559:560];
FAs_lowDA  = [561:562];

hits_rew_highDA = [563:564];
hits_rew_lowDA  = [565:566];
FAs_rew_highDA  = [567:568];
FAs_rew_lowDA   = [569:570];


Tstat_highDA_HPCclusters_Scene_Rew_v_Neu_p005=785;
Tstat_lowDA_HPCclusters_Scene_Rew_v_Neu_p005=786;

Tstat2_highDA_HPCclusters_Scene_Rew_v_Neu_p005=789;
Tstat2_lowDA_HPCclusters_Scene_Rew_v_Neu_p005=790;

% Specify variables to control for in partial correlation analysis
% for both memory tests collapsed, use var_names_mem

% ctrl_vars = [var_names(AGEs); var_names(AppliedAct_highDA)]; % for baseline measures
ctrl_vars = [var_names(AGEs)];
var_selected = [var_names(BPchanges_SRTM_lowDA); var_names(Dprimes_rew_lowDA); var_names(Dprimes_neu_lowDA)];%[var_names(OCCs);var_names(Dprimes)];

n = height(data);
    
%% now plot, single plots

close all
set(0, 'DefaultAxesFontSize', 14,'DefaultLineLineWidth',2);

sigmoidalFunc = @(b, x) 1 ./ (1 + exp(-b(1) * x + b(2)));

% Loop through all pairs of variables
for i = 1:length(var_selected)
    for j = (i+1):length(var_selected)
        
        % Create subset data table for partial correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}, ctrl_vars{:}});
        
        subset_data{:,2}=subset_data{:,2}.*-1;
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Perform partial correlation analysis
        [r, p] = partialcorr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), table2array(subset_data(:, 3:end)), 'Type', 'Spearman');
        
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Check if any correlations are significant after multiple comparisons correction
        if any(p_adj < 0.05)
            
            % Select only significant correlations
            sig_corr = r(p_adj < 0.05);
            
            % Plot scatter plot with fit line and equation
            figure;
            scatter(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 200,'MarkerFaceColor','b','MarkerFaceAlpha',0.5);
            hold on;
            pfit = polyfit(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 1);
            plot(table2array(subset_data(:, 1)), polyval(pfit, table2array(subset_data(:, 1))), 'r--');
            clear str1 str2
            str1 = var_selected{i};
            str1 = strrep(str1, '_', ' ');
            str2 = var_selected{j};
            str2 = strrep(str2, '_', ' ');
            xlabel(str1);
            ylabel(str2);
            
            eqn = sprintf('y = %.2fx + %.2f (p_F_D_R_c = %.5f)', pfit(1), pfit(2), min(p_adj(p_adj < 0.05)));
            text(min(table2array(subset_data(:, 1))), max(table2array(subset_data(:, 2))), eqn,'FontSize',20);
            
            title({'Partial Correlation Plots: Spearman''s Rho, FDR-corrected','Controlled: \it{AGE, Applied Activity}\rm'},'FontSize',24,'FontWeight','bold')
            
        end
        
    end
end

disp('done')

%% alternative plotting, Partial correlation controlling for age : i prefer this one

for FoldScripts=1

close all
set(0, 'DefaultAxesFontSize', 20,'DefaultLineLineWidth',2);

% Initialize a matrix to store partial correlation coefficients
num_vars = length(var_selected);
partial_corr_matrix = zeros(num_vars, num_vars);
corr_coeff_matrix = zeros(num_vars, num_vars);

% Loop through all pairs of variables
num_sig = 0;
for i = 1:length(var_selected)
    for j = (i+1):length(var_selected)
        
        % Create subset data table for partial correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}, ctrl_vars{:}});
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Perform partial correlation analysis
        [r, p] = partialcorr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), table2array(subset_data(:, 3:end)), 'Type', 'Spearman');
        
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Check if any correlations are significant after multiple comparisons correction
        if any(p_adj < 0.05)
            
            % Select only significant correlations
            sig_corr = r(p_adj < 0.05);
            num_sig = num_sig + 1;
            
            % Store the partial correlation coefficient in the matrix
            partial_corr_matrix(i, j) = sig_corr;
            corr_coeff_matrix(i, j) = r;
            
        end
        
    end
end

% Set number of rows and columns for subplot
num_rows = ceil(sqrt(num_sig));
num_cols = ceil(num_sig/num_rows);

% Loop through all pairs of variables again
subplot_idx = 1;
for i = 1:length(var_selected)
    for j = (i+1):length(var_selected)
        
        % Create subset data table for partial correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}, ctrl_vars{:}});
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Perform partial correlation analysis
        [r, p] = partialcorr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), table2array(subset_data(:, 3:end)), 'Type', 'Spearman');
        
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Check if any correlations are significant after multiple comparisons correction
        if any(p_adj < 0.05)
            
            % Select only significant correlations
            sig_corr = r(p_adj < 0.05);
            
            % Plot scatter plot with fit line and equation
            subplot(num_rows, num_cols, subplot_idx);
            scatter(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)),'MarkerFaceColor','b','MarkerFaceAlpha',0.5);
            hold on;
            pfit = polyfit(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 1);
            plot(table2array(subset_data(:, 1)), polyval(pfit, table2array(subset_data(:, 1))), 'r--');
            clear str1 str2
            str1 = var_selected{i};
            str1 = strrep(str1, '_', ' ');
            str2 = var_selected{j};
            str2 = strrep(str2, '_', ' ');
            xlabel(str1);
            ylabel(str2);
            set(gca,'FontSize',10,'LineWidth',1);
            
            eqn = sprintf('y = %.2fx + %.2f (p_F_D_R_c = %.5f , rho = %.3f )', pfit(1), pfit(2), min(p_adj(p_adj < 0.05)), r);
            text(min(table2array(subset_data(:, 1))), max(table2array(subset_data(:, 2))), eqn,'FontSize',14);
            
            subplot_idx = subplot_idx + 1;
        end
        
    end
end

% Set figure title
sgtitle({['Partial Correlation Plots: Spearman''s Rho, FDR-corrected','Controlled: \it{AGE}\rm']},'FontSize',24,'FontWeight','bold')


% --------------------------------------------------------- %

% Compute adjusted p-values for each correlation coefficient pair
for i = 1:num_vars
    for j = (i+1):num_vars
        % Create subset data table for partial correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}, ctrl_vars{:}});
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Perform partial correlation analysis
        [r, p] = partialcorr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), table2array(subset_data(:, 3:end)), 'Type', 'Spearman');
        
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Store the FDR-adjusted p-values in the partial_corr_matrix
        partial_corr_matrix(i, j) = p_adj;
        partial_corr_matrix(j, i) = p_adj; % Since it's symmetric
        
        corr_coeff_matrix(i, j) = r;
        corr_coeff_matrix(j, i) = r;

    end
end

% Create a figure for the correlation matrix
figure;

% Plot the correlation matrix as a heatmap without colors for cells below the diagonal
imagesc(triu(corr_coeff_matrix, 1),[-1, 1]);
caxis([-1, 1]); % Set the color axis limits based on the threshold

% Add red asterisks for significant correlations in the upper triangle
[num_rows, num_cols] = size(partial_corr_matrix);
for row = 1:num_rows
    for col = row:num_cols
        if partial_corr_matrix(row, col) ~= 1 && partial_corr_matrix(row, col) ~= 0
            % Check if the correlation coefficient is significant
            if partial_corr_matrix(row, col) < 0.05
                text(col, row, '*', 'Color', 'red', 'FontSize', 40, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
    end
end

% Add colorbar and labels
colorbar;
xticks(1:num_vars);
yticks(1:num_vars);
xticklabels(var_selected);
yticklabels(var_selected);
xtickangle(45);
xlabel('Variables in X');
ylabel('Variables in X');
title('Partial Correlation Matrix');

% Adjust figure size to accommodate x-axis labels
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), pos(3), pos(4) + 50]);
% --------------------------------------------------------- %

end

%% no partial, just correlation in spearman's rho

for FoldScripts=1
close all

clear num_rows num_cols

set(0, 'DefaultAxesFontSize', 20,'DefaultLineLineWidth',2);
 
% Initialize a matrix to store partial correlation coefficients
num_vars = length(var_selected);
corr_matrix = zeros(num_vars, num_vars);
corr_rho_coeff_matrix = zeros(num_vars, num_vars);
 
% Loop through all pairs of variables
num_sig = 0;
for i = 1:length(var_selected)
    for j = (i+1):length(var_selected)
        
        % Create subset data table for correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}});
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Calculate Spearman's correlation coefficient
        [r, p] = corr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 'Type', 'Spearman');
        
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Check if any correlations are significant after multiple comparisons correction
        if any(p_adj < 0.05)
            
          % Select only significant correlations
            sig_corr = r(p_adj < 0.05);
            num_sig = num_sig + 1;
            
            % Store the partial correlation coefficient in the matrix
            corr_matrix(i, j) = sig_corr;
            corr_rho_coeff_matrix(i, j) = r;

        end
        
    end
end
 
% Set number of rows and columns for subplot
num_rows = ceil(sqrt(num_sig));
num_cols = ceil(num_sig/num_rows);
 
% Loop through all pairs of variables again
subplot_idx = 1;
for i = 1:length(var_selected)
    for j = (i+1):length(var_selected)
        
       % Create subset data table for correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}});
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Calculate Spearman's correlation coefficient
        [r, p] = corr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 'Type', 'Spearman');
                
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Check if any correlations are significant after multiple comparisons correction
        if any(p_adj < 0.05)
            
            % Select only significant correlations
            sig_corr = r(p_adj < 0.05);
            
            % Plot scatter plot with fit line and equation
            subplot(num_rows, num_cols, subplot_idx);
            scatter(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)),'MarkerFaceColor','b','MarkerFaceAlpha',0.5);
            hold on;
            pfit = polyfit(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 1);
            plot(table2array(subset_data(:, 1)), polyval(pfit, table2array(subset_data(:, 1))), 'r--');
            clear str1 str2
            str1 = var_selected{i};
            str1 = strrep(str1, '_', ' ');
            str2 = var_selected{j};
            str2 = strrep(str2, '_', ' ');
            xlabel(str1);
            ylabel(str2);
            set(gca,'FontSize',10,'LineWidth',1);
            
            eqn = sprintf('y = %.2fx + %.2f (p_F_D_R_c = %.5f , rho = %.3f )', pfit(1), pfit(2), min(p_adj(p_adj < 0.05)), r);
            text(min(table2array(subset_data(:, 1))), max(table2array(subset_data(:, 2))), eqn,'FontSize',14);
            
            subplot_idx = subplot_idx + 1;
        end
        
    end
end
 
% Set figure title
sgtitle({'Non-parametric Correlation Plots: Spearman''s Rho, FDR-corrected'},'FontSize',24,'FontWeight','bold')
 
 
% --------------------------------------------------------- %
 
% Compute adjusted p-values for each correlation coefficient pair
for i = 1:num_vars
    for j = (i+1):num_vars
        
        % Create subset data table for correlation analysis
        subset_data = data(:, {var_selected{i}, var_selected{j}});
        
        % Remove any rows with missing values
        subset_data = rmmissing(subset_data);
        
        % Calculate Spearman's correlation coefficient
        [r, p] = corr(table2array(subset_data(:, 1)), table2array(subset_data(:, 2)), 'Type', 'Spearman');
        
        
        % Adjust p-values for multiple comparisons
        p_adj = mafdr(p, 'BHFDR', true);
        
        % Store the FDR-adjusted p-values in the partial_corr_matrix
        corr_matrix(i, j) = p_adj;
        corr_matrix(j, i) = p_adj; % Since it's symmetric
        
        corr_rho_coeff_matrix(i, j) = r;
        corr_rho_coeff_matrix(j, i) = r;

 
    end
end
 
% Create a figure for the correlation matrix
figure;
 
% Plot the correlation matrix as a heatmap without colors for cells below the diagonal
imagesc(triu(corr_rho_coeff_matrix, 1),[-1, 1]);
caxis([-1, 1]); % Set the color axis limits based on the threshold
 
% Add red asterisks for significant correlations in the upper triangle
[num_rows, num_cols] = size(corr_matrix);
for row = 1:num_rows
    for col = row:num_cols
        if corr_matrix(row, col) ~= 1 && corr_matrix(row, col) ~= 0
            % Check if the correlation coefficient is significant
            if corr_matrix(row, col) < 0.05
                text(col, row, '*', 'Color', 'red', 'FontSize', 40, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
    end
end
 
% Add colorbar and labels
colorbar;
xticks(1:num_vars);
yticks(1:num_vars);
xticklabels(var_selected);
yticklabels(var_selected);
xtickangle(45);
xlabel('Variables in X');
ylabel('Variables in X');
title('Non-parametric Correlation Matrix');
 
% Adjust figure size to accommodate x-axis labels
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), pos(3), pos(4) + 50]);
% --------------------------------------------------------- %

end

%% violin plots

close all

% Calculate the number of variables and the total number of subplots
num_vars = length(var_selected);

% Calculate the number of rows and columns for the subplot grid
num_rows = ceil(sqrt(num_vars));
num_cols = ceil(num_vars / num_rows);

% Create a new figure for violin plots
figure;

% Adjust subplot spacing
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

% Loop through each variable and plot violin plots
for i = 1:num_vars
    % Create the subplot
    subplot(num_rows, num_cols, i);
    x = data{:, var_selected{i}};
    
    % Compute the kernel density estimate for the violin plot
    [f, xi] = ksdensity(x);
    
    % Plot the left half of the violin plot
    fill([fliplr(xi) xi], [fliplr(-f) f], [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    hold on;
    
    % Plot the median line
    median_x = nanmedian(x);
    plot([median_x, median_x], [-max(f), max(f)], 'r-', 'LineWidth', 2);
    
    % Plot the quartile lines
    quartile1_x = prctile(x, 25);
    quartile3_x = prctile(x, 75);
    plot([quartile1_x, quartile3_x], [0, 0], 'k-', 'LineWidth', 2);
    
    hold off;
    
    % Set the x-axis label
%     xlabel(var_selected{i}, 'Interpreter', 'none');

    % Set the title
    clear str1
    str1 = var_selected{i};
    str1 = strrep(str1, '_', ' ');
    title(str1);
    
    % Set the y-axis label
    ylabel('kernel density');
    
    
    % Rotate x-axis labels if needed
    xtickangle(45);
    
%     xlim([-0.5 0.5])
end

% Adjust figure title
sgtitle('Violin Plots', 'FontSize', 24, 'FontWeight', 'bold');


%%
% Create a folder to save the PDF file
pdf_folder = '/Users/alex/Dropbox/paperwriting/MRPET/figures';

% Create a PDF file to save all violin plots
pdf_filename = fullfile(pdf_folder, 'violin_plots.ps');

% Loop through each variable and plot the violin plots
for i = 1:num_vars
    
    % Create a new figure for each violin plot
    figure;
    
    % Set the title
    clear str1
    str1 = var_selected{i};
    str1 = strrep(str1, '_', ' ');
    xlabel(str1);
    
    % Create a violin plot using the violin function
    violin(data{:, i}, 'xlabel', {str1}, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'k');
    
    title(str1)
    
    % Save the figure to the PDF file
%     if i == 1
%         % For the first plot, create a new PDF file
%         print('-dpsc2','-append','-bestfit', pdf_filename);
%     else
        % For subsequent plots, append to the existing PDF file
        print('-dpsc2','-append','-bestfit', pdf_filename);
%     end

    % Close the figure to clear memory
    close gcf;
end

%% remove outliers

% % Read in data from Excel file
% data = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BP_and_Memory_bothSession_05-May-2023_v2.xlsx');
% 
% % Loop through each variable
% for i = 1:554%size(data, 2)
%     
%     % Calculate mean and standard deviation of current variable
%     var_mean = mean(data{:, i}, 'omitnan');
%     var_std = std(data{:, i}, 'omitnan');
%     
%     % Find any values that are more than 3 standard deviations from the mean
%     outliers = abs(data{:, i} - var_mean) > 3 * var_std;
%     
%     % Replace outliers with NaN
%     data{outliers, i} = NaN;
%     
% end
% 
% % Save cleaned data to new sheet in Excel file
% writetable(data, '/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BP_and_Memory_bothSession_05-May-2023_outlierRemoved_Z3.xlsx');

