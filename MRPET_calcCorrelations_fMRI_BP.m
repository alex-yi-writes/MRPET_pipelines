%% calculate correlation with fMRI and BPND

%% prepare

clc;clear

% set up paths
paths_fMRI = '/Volumes/ALEX3/MRPET/analysis/eachSession/';
paths_mask = '/Users/alex/Dropbox/Masks/Masks/mni_icbm152/';

% IDs : try 40141 again!!!!!
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

% load PET data
% load('/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/MRPET_BPpackage_04-May-2023.mat')
load('/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/MRPET_BPpackage_All_smoothed3mm_07-Aug-2023.mat')
fMRIclusterBP=load('/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/MRPET_BPpackage_condensed_12-Jul-2023_fMRIclusters.mat');

% load fMRI data
% hipp=spm_read_vols(spm_vol([paths_mask 'mni_icbm152_AAL3_hippocampus_bi.nii']));
hipp=spm_read_vols(spm_vol([paths_fMRI '2ndLvL/STIM_rew_v_neu_p005_c370_HippClusterMask.nii']));
fMRI=[];
for id=1:length(IDs)
    for d=1:2
        if days(id,d)==0
            fMRI.stim.rew_v_neu{id,d}=[];
        else
            fMRI.stim.rew_v_neu{id,d}=spm_read_vols(spm_vol([paths_fMRI num2str(IDs(id)) '_' num2str(d) '/con_0002_mni.nii']));
        end
    end
end

%% extract t-values within the hippocampus mask

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

%% calc BP change in hippocampus
% BPchange_SRTM_save{id,d}.(FSlabels1{labels,1})=(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl - BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm) / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);

BPchange_hipp=[];
for id=1:length(IDs)
    for d=1:2
        if days(id,d)==0
            BPchange_hipp(id,d)=NaN;
            BP_lpnt_hipp(id,d) =NaN;
        else
            BPchange_hipp(id,d) = (BPdataSave{id,d}.Hipp.BP_srtm_bl - BPdataSave{id,d}.Hipp.BP_srtm) / (BPdataSave{id,d}.Hipp.BP_srtm_bl);
            BP_lpnt_hipp(id,d)  = BP_lp_save{id,d}.Hipp; 
        end
    end
end

%% calc BP change in the hippocampus activation cluster mask

BPchange_hipp_cluster=[];
for id=1:length(IDs)
    for d=1:2
        if days(id,d)==0
            BPchange_hipp_cluster(id,d)=NaN;
            BP_lpnt_hipp_cluster(id,d) =NaN;
        else
            BPchange_hipp_cluster(id,d) = (fMRIclusterBP.BPdataSave{id,d}.fMRI_HPC_STIM_rew_neu.BP_srtm_bl - fMRIclusterBP.BPdataSave{id,d}.fMRI_HPC_STIM_rew_neu.BP_srtm) / (fMRIclusterBP.BPdataSave{id,d}.fMRI_HPC_STIM_rew_neu.BP_srtm_bl);
            BP_lpnt_hipp_cluster(id,d)  = fMRIclusterBP.BP_lp_save{id,d}.fMRI_HPC_STIM_rew_neu; 
        end
    end
end

%% correlation : in HPC cluster

close all
figure

subplot(1,2,1)
hold on

% Sample data
y = BPchange_hipp_cluster(:,1);
x = meanTval_hipp(:,1);

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f**x + %.2f', p(1), p(2));
elseif P<=0.1
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('∆BP_{ND}', 'FontSize', 14);
xlabel('T_{HPC}', 'FontSize', 14);
title({'high DA session: T_{HPC} and ∆BP_{ND}','in fMRI actvation cluster', '(Stim:Rew>Neu)'}, 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

% correlation lowDA session

subplot(1,2,2)


% Sample data
y = BPchange_hipp_cluster(:,2);
x = meanTval_hipp(:,2);

% Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('∆BP_{ND}', 'FontSize', 14);
xlabel('T_{HPC}', 'FontSize', 14);
title({'low DA session: T_{HPC} and ∆BP_{ND}','in fMRI actvation cluster', '(Stim:Rew>Neu)'}, 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

%% correlation : highDA session

close all
figure

subplot(1,2,1)
hold on

% Sample data
y = BPchange_hipp(:,1);
x = meanTval_hipp(:,1);

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f**x + %.2f, p = %.3f, Rho = %.3f', p(1), p(2), P, Rho);
    
elseif P<0.1
    eq_str = sprintf('y = %.2f*x + %.2f, p = %.3f, Rho = %.3f', p(1), p(2), P, Rho);
else
    eq_str = sprintf('y = %.2fx + %.2f, p = %.3f, Rho = %.3f', p(1), p(2), P, Rho);
end

dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('∆BP_N_D', 'FontSize', 14);
xlabel('T_H_P_C', 'FontSize', 14);
title('high DA session: T_H_P_C and ∆BP_N_D', 'FontSize', 16);

% Set axis limits
% xlim([min(x)-0.5, max(x)+0.5]);
% ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');



% correlation lowDA session

subplot(1,2,2)


% Sample data
y = BPchange_hipp(:,2);
x = meanTval_hipp(:,2);

% Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f**x + %.2f, p = %.3f, Rho = %.3f', p(1), p(2), P, Rho);
    
elseif P<0.1
    eq_str = sprintf('y = %.2f*x + %.2f, p = %.3f, Rho = %.3f', p(1), p(2), P, Rho);
else
    eq_str = sprintf('y = %.2fx + %.2f, p = %.3f, Rho = %.3f', p(1), p(2), P, Rho);
end

dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('∆BP_N_D', 'FontSize', 14);
xlabel('T_H_P_C', 'FontSize', 14);
title('low DA session: T_H_P_C and ∆BP_N_D', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');


%% correlation both sessions

figure

% Sample data
y = [BPchange_hipp(:,1); BPchange_hipp(:,2)];
x = [meanTval_hipp(:,1); meanTval_hipp(:,2)];

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('∆BP_N_D', 'FontSize', 14);
xlabel('T_H_P_C', 'FontSize', 14);
title('both sessions: T_H_P_C and ∆BP_N_D', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');


%% highDA session correlation (lpntPET)

close all
figure

% Sample data
x = BP_lpnt_hipp(:,1);
y = meanTval_hipp(:,1);

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
eq_x = nanmean(x); % x coordinate for the equation label
eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
xlabel('BP_l_p_n_t', 'FontSize', 14);
ylabel('T_H_P_C', 'FontSize', 14);
title('high DA session: T_H_P_C and BP_l_p_n_t', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

%% lowDA session correlation (lpntPET)

figure

% Sample data
x = BP_lpnt_hipp(:,2);
y = meanTval_hipp(:,2);

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
eq_x = nanmean(x); % x coordinate for the equation label
eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
xlabel('BP_l_p_n_t', 'FontSize', 14);
ylabel('T_H_P_C', 'FontSize', 14);
title('low DA session: T_H_P_C and BP_l_p_n_t', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

%% both session correlation (lpntPET)

figure

% Sample data
x = [BP_lpnt_hipp(:,1); BP_lpnt_hipp(:,2)];
y = [meanTval_hipp(:,1); meanTval_hipp(:,2)];

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
eq_x = nanmean(x); % x coordinate for the equation label
eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
xlabel('BP_l_p_n_t', 'FontSize', 14);
ylabel('T_H_P_C', 'FontSize', 14);
title('both sessions: T_H_P_C and BP_l_p_n_t', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

%% both session correlation (lpntPET, fMRI Clusters)

figure

% Sample data
x = [BP_lpnt_hipp_cluster(:,1); BP_lpnt_hipp_cluster(:,2)];
y = [meanTval_hipp(:,1); meanTval_hipp(:,2)];

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
eq_x = nanmean(x); % x coordinate for the equation label
eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
xlabel('BP_l_p_n_t', 'FontSize', 14);
ylabel('T_H_P_C', 'FontSize', 14);
title('both sessions: T_H_P_C and BP_l_p_n_t', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

%% high DA
close all
figure

% Sample data
x = BP_lpnt_hipp_cluster(:,1);
y = meanTval_hipp(:,1);

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
eq_x = nanmean(x); % x coordinate for the equation label
eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
xlabel('BP_l_p_n_t', 'FontSize', 14);
ylabel('T_H_P_C', 'FontSize', 14);
title('high DA session: T_H_P_C and BP_l_p_n_t', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Data', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

%% for blitz talk

%% correlation : highDA session

close all
figure

subplot(1,2,1)
hold on

% Sample data
y = -BPchange_hipp(:,1);
x = meanTval_hipp(:,1);

% % Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('Estimate of endogenous DA release (-∆BP_{ND})', 'FontSize', 14);
xlabel('T_{HPC}', 'FontSize', 14);
title('high DA session: T_{HPC} and DA release', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
% ylim([min(y)-0.5, max(y)+0.5]);
ylim([-0.-0.04, 0.1]);

% Show legend
legend({'Subject', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');

% correlation lowDA session

subplot(1,2,2)


% Sample data
y = -BPchange_hipp(:,2);
x = meanTval_hipp(:,2);

% Remove NaN values
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% Create scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');

% Fit a linear regression line
p = polyfit(x, y, 1);

% Evaluate the line equation at x values
yfit = p(1)*x + p(2);

% Compute correlation coefficient and p-value
[Rho, P] = corr(x, y, 'type', 'Spearman');

% Plot the line
hold on;
plot(x, yfit, 'r-', 'LineWidth', 2);

% Add equation label
if P<=0.05
    eq_str = sprintf('y = %.2f*x + %.2f', p(1), p(2));
else
    eq_str = sprintf('y = %.2fx + %.2f', p(1), p(2));
end
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', eq_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);
% eq_x = nanmean(x); % x coordinate for the equation label
% eq_y = p(1)*eq_x + p(2); % y coordinate for the equation label
% text(eq_x, eq_y, eq_str, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

% Set axis labels and title
ylabel('Estimate of endogenous DA release (-∆BP_{ND})', 'FontSize', 14);
xlabel('T_{HPC}', 'FontSize', 14);
title('low DA session: T_{HPC} and DA release', 'FontSize', 16);

% Set axis limits
xlim([min(x)-0.5, max(x)+0.5]);
ylim([-0.-0.04, 0.1]);
% ylim([min(y)-0.5, max(y)+0.5]);

% Show legend
legend({'Subject', 'Fit Line'}, 'FontSize', 12, 'Location', 'northwest');