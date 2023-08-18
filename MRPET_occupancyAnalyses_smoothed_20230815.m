%% Occupancy analyses

clc;clear;

paths.parent  = '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/original/';
paths.TACs    = '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/Smoothed/';
paths.TACs_new= '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/Smoothed/';
paths.ROImask = '/Volumes/ALEX3/MRPET/coreg_roi/';
paths.figure  = '/Users/alex/Dropbox/paperwriting/MRPET/figures/modellingFigures/';
paths.rawImg  = '/Volumes/ALEX3/MRPET/img/';
paths.exports = '/Users/alex/Dropbox/paperwriting/MRPET/data/';

addpath(genpath('/Users/alex/Dropbox/paperwriting/MRPET/scripts/modelling'))

% set envs
% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

load([paths.exports 'Occupancy_DBP_20230815.mat'])

%% Robust Regression or LOESS
% there's an underlying trend but it's hidden behind noise... so i'm
% performing LOESS (locally estimated scatterplot smoothing) to get a 
% clearer picture of the trend over time. 

close all;

% Extract list of ROIs using fieldnames
rois_all = fieldnames(Occupancy{1, 1});
rois = rois_all([1:7 44:53]);

% Create a time vector
time = 1:111;

% Determine grid dimensions for subplots
numROIs = length(rois);
rows = floor(sqrt(numROIs));
cols = ceil(numROIs/rows);

% Create a single figure for all subplots
figure;

% Loop through each ROI
for r = 1:length(rois)
    roi = rois{r};

    % Initialize matrices to store values for each subject and condition
    data_highDA = zeros(111, length(Occupancy));
    data_lowDA = zeros(111, length(Occupancy));

    % Initialize count for clipped subjects
    clippedSubjects_highDA = 0;
    clippedSubjects_lowDA = 0;

    % Aggregate data for each subject and condition
    for subj = 1:length(Occupancy)
        wasClipped_highDA = false;
        wasClipped_lowDA = false;

        if ~isempty(Occupancy{subj, 1})
            data = Occupancy{subj, 1}.(roi);
            if subj == 31
                data = [zeros(14,1); data];
            end

            if any(data > 100) || any(data < 0)
                wasClipped_highDA = true;
            end

            data(data > 100) = 100;
            data(data < 0) = 0;
            data_highDA(:, subj) = data;
        end
        
        if ~isempty(Occupancy{subj, 2})
            data = Occupancy{subj, 2}.(roi);

            if any(data > 100) || any(data < 0)
                wasClipped_lowDA = true;
            end

            data(data > 100) = 100;
            data(data < 0) = 0;
            data_lowDA(:, subj) = data;
        end

        if wasClipped_highDA
            clippedSubjects_highDA = clippedSubjects_highDA + 1;
        end

        if wasClipped_lowDA
            clippedSubjects_lowDA = clippedSubjects_lowDA + 1;
        end
    end

    avg_highDA = mean(data_highDA, 2);
    avg_lowDA = mean(data_lowDA, 2);

    spanValue = 0.3;
    smoothed_avg_highDA = smooth(time, avg_highDA, 'loess', spanValue);
    smoothed_avg_lowDA = smooth(time, avg_lowDA, 'loess', spanValue);

    residuals_highDA = avg_highDA - smoothed_avg_highDA;
    residuals_lowDA = avg_lowDA - smoothed_avg_lowDA;

    Rsq_highDA = 1 - sum(residuals_highDA.^2)/sum((avg_highDA - mean(avg_highDA)).^2);
    Rsq_lowDA = 1 - sum(residuals_lowDA.^2)/sum((avg_lowDA - mean(avg_lowDA)).^2);

    subplot(rows, cols, r);
    plot(time, smoothed_avg_highDA, 'r-', 'LineWidth', 2);
    hold on;
    plot(time, smoothed_avg_lowDA, 'b-', 'LineWidth', 2);
    hold on;
    xline(100,'--k','LineWidth',2);
    xlim([80 111]);
    ylim([-10 50]);
    xlabel('Time');
    ylabel('Occupancy (%)');
    title([roi ' | R^2 HD: ' num2str(Rsq_highDA, '%.2f') ' LD: ' num2str(Rsq_lowDA, '%.2f')]);
    text(85, 45, ['Span: ' num2str(spanValue) ' | Clipped HD: ' num2str(clippedSubjects_highDA) ' LD: ' num2str(clippedSubjects_lowDA)], 'FontSize', 8);
end

sgtitle('LOESS smoothed Average Occupancy');
legend({'High DA', 'Low DA','Task Onset'},'Location','northwest');

%% without clipping

close all;

% Extract list of ROIs using fieldnames
rois_all = fieldnames(Occupancy{1, 1});
rois = rois_all([1:7 44:53]);

% Create a time vector
time = 1:111;

% Determine grid dimensions for subplots
numROIs = length(rois);
rows = floor(sqrt(numROIs));
cols = ceil(numROIs/rows);

% Create a single figure for all subplots
figure;

% Loop through each ROI
for r = 1:length(rois)
    roi = rois{r};

    % Initialize matrices to store values for each subject and condition
    data_highDA = zeros(111, length(Occupancy));
    data_lowDA = zeros(111, length(Occupancy));

    % Aggregate data for each subject and condition
    for subj = 1:length(Occupancy)
        if ~isempty(Occupancy{subj, 1})
            data = Occupancy{subj, 1}.(roi);
            if subj == 31
                data = [zeros(14,1); data];
            end
            data_highDA(:, subj) = data;
        end
        
        if ~isempty(Occupancy{subj, 2})
            data = Occupancy{subj, 2}.(roi);
            data_lowDA(:, subj) = data;
        end
    end

    avg_highDA = mean(data_highDA, 2);
    avg_lowDA = mean(data_lowDA, 2);

    spanValue = 0.3;
    smoothed_avg_highDA = smooth(time, avg_highDA, 'loess', spanValue);
    smoothed_avg_lowDA = smooth(time, avg_lowDA, 'loess', spanValue);

    residuals_highDA = avg_highDA - smoothed_avg_highDA;
    residuals_lowDA = avg_lowDA - smoothed_avg_lowDA;

    Rsq_highDA = 1 - sum(residuals_highDA.^2)/sum((avg_highDA - mean(avg_highDA)).^2);
    Rsq_lowDA = 1 - sum(residuals_lowDA.^2)/sum((avg_lowDA - mean(avg_lowDA)).^2);

    subplot(rows, cols, r);
    plot(time, smoothed_avg_highDA, 'r-', 'LineWidth', 2);
    hold on;
    plot(time, smoothed_avg_lowDA, 'b-', 'LineWidth', 2);
    hold on;
    xline(100,'--k','LineWidth',2);
    xlim([80 111]);
    ylim([-10 100]);
    xlabel('Time');
    ylabel('Occupancy (%)');
    title([roi ' | R^2 HD: ' num2str(Rsq_highDA, '%.2f') ' LD: ' num2str(Rsq_lowDA, '%.2f')]);
    text(85, 45, ['Span: ' num2str(spanValue)], 'FontSize', 8);
end

sgtitle('LOESS smoothed Average Occupancy');
legend({'High DA', 'Low DA','Task Onset'},'Location','northwest');

