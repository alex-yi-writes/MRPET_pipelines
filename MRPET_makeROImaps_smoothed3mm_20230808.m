%% make ROI maps

clc;clear

load('/Users/alex/Dropbox/paperwriting/MRPET/data/FS_labels_LCSNVTA.mat')
load('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_table.mat')

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      SRTM                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% t-test with FDRc

clear significant_results_table p_values stats t_stats fdr_corrected_p_values summary_table significant_results_table ...
    significant_ROIs allVariableNames indices_of_significant_variables

% Extract the subset of variables
selectedVariables = MRPET_highSession_Table(:, 2:47);

% Value you want to test against
testValue = 0; % practically 0 in cerebellar cortex

% Initialize p-values and T statistics arrays
p_values = zeros(1, width(selectedVariables));
t_stats = zeros(1, width(selectedVariables));

% Perform one-sample t-test for each variable and collect p-values and T statistics
for i = 1:width(selectedVariables)
    [~, p_values(i), ~, stats] = ttest(selectedVariables{:, i}, testValue);
    t_stats(i) = stats.tstat;
end

% Correct p-values using the False Discovery Rate (FDR) method
fdr_corrected_p_values = mafdr(p_values, 'BHFDR', true);

% Create a summary table with variable names, p-values, and T statistics
summary_table = array2table([p_values; fdr_corrected_p_values; t_stats]', ...
    'VariableNames', {'Original_P_Value', 'Corrected_P_Value', 'T_Stat'}, ...
    'RowNames', selectedVariables.Properties.VariableNames);

% Filter the table to only include rows where the corrected p-value is less than 0.05
significant_results_table = summary_table(fdr_corrected_p_values < 0.05, :);

% Convert the numeric values to strings with four decimal places for p-values
significant_results_table.Original_P_Value = cellstr(num2str(significant_results_table.Original_P_Value, '%.4f'));
significant_results_table.Corrected_P_Value = cellstr(num2str(significant_results_table.Corrected_P_Value, '%.4f'));

% You can also format T_Stat if needed
% significant_results_table.T_Stat = cellstr(num2str(significant_results_table.T_Stat, '%.4f'));

% Display the table containing only significant results
disp('Significant Results:');
disp(significant_results_table);



%% significance map!

% high-reward session
clear ZScoreImg baseimage FSlabel existingLabels hdr indices_of_significant_variables

FSlabels2=FSlabels1;
FSlabels2(1,:)=[];

allVariableNames=cellfun(@char,FSlabels2(:,1),'UniformOutput',false);

% Get the names of the significant variables from the table
significant_variable_names = significant_results_table.Properties.RowNames;
significant_ROIs           = cellfun(@(x) x(15:(end-7)),significant_variable_names,'Uniformoutput',false);

% Iterate over the significant_ROIs and find the matching index for each one
for i = 1:numel(significant_ROIs)
    indices_of_significant_variables(i) = find(strcmp(allVariableNames, significant_ROIs{i}), 1);
end

% Display the indices
disp('Indices of Significant Variables:');
disp(indices_of_significant_variables);

baseimage = spm_read_vols(spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']));
TstatImg = zeros(size(baseimage));
existingLabels = unique(baseimage);

cnt=0;
for labels=indices_of_significant_variables
    if FSlabels2{labels,2}==0
    else
        
        fprintf(['**' FSlabels2{labels,1}{1} '**\n'])
        
        cnt=cnt+1;
        
        % left
        disp('left')
        clear indL
        disp(['label:' num2str(FSlabels2{labels,2})])
        indL=find(baseimage==FSlabels2{labels,2});
        TstatImg(indL)=significant_results_table.T_Stat(cnt);
        significant_results_table.T_Stat(cnt)
        
        % right
        disp('right')
        clear indR
        disp(['label:' num2str(FSlabels2{labels,3})])
        indR=find(baseimage==FSlabels2{labels,3});
        TstatImg(indR)=significant_results_table.T_Stat(cnt);
        
    end
end


hdr = spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
hdr.fname = ['/Users/alex/Dropbox/paperwriting/MRPET/data/TstatMaps_highDAsession_BPchange_SRTM.nii'];
hdr.dim = size(TstatImg);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,TstatImg);

fprintf(['\n\n\n'])

%% lowDA session

clear significant_results_table p_values stats t_stats fdr_corrected_p_values summary_table significant_results_table ...
    significant_ROIs allVariableNames

% Extract the subset of variables
selectedVariables = MRPET_lowSession_Table(:, 2:47);

% Value you want to test against
testValue = 0; % practically 0 in cerebellar cortex

% Initialize p-values and T statistics arrays
p_values = zeros(1, width(selectedVariables));
t_stats = zeros(1, width(selectedVariables));

% Perform one-sample t-test for each variable and collect p-values and T statistics
for i = 1:width(selectedVariables)
    [~, p_values(i), ~, stats] = ttest(selectedVariables{:, i}, testValue);
    t_stats(i) = stats.tstat;
end

% Correct p-values using the False Discovery Rate (FDR) method
fdr_corrected_p_values = mafdr(p_values, 'BHFDR', true);

% Create a summary table with variable names, p-values, and T statistics
summary_table = array2table([p_values; fdr_corrected_p_values; t_stats]', ...
    'VariableNames', {'Original_P_Value', 'Corrected_P_Value', 'T_Stat'}, ...
    'RowNames', selectedVariables.Properties.VariableNames);

% Filter the table to only include rows where the corrected p-value is less than 0.05
significant_results_table = summary_table(fdr_corrected_p_values < 0.05, :);

% Convert the numeric values to strings with four decimal places for p-values
significant_results_table.Original_P_Value = cellstr(num2str(significant_results_table.Original_P_Value, '%.4f'));
significant_results_table.Corrected_P_Value = cellstr(num2str(significant_results_table.Corrected_P_Value, '%.4f'));

% You can also format T_Stat if needed
% significant_results_table.T_Stat = cellstr(num2str(significant_results_table.T_Stat, '%.4f'));

% Display the table containing only significant results
disp('Significant Results:');
disp(significant_results_table);



% significance map 2

% high-reward session
clear ZScoreImg baseimage FSlabel existingLabels hdr indices_of_significant_variables

FSlabels2=FSlabels1;
FSlabels2(1,:)=[];


allVariableNames=cellfun(@char,FSlabels2(:,1),'UniformOutput',false);

% Get the names of the significant variables from the table
significant_variable_names = significant_results_table.Properties.RowNames;
significant_ROIs           = cellfun(@(x) x(15:(end-6)),significant_variable_names,'Uniformoutput',false);

% Iterate over the significant_ROIs and find the matching index for each one
indices_of_significant_variables=[];
for i = 1:numel(significant_ROIs)
    indices_of_significant_variables(i) = find(strcmp(allVariableNames, significant_ROIs{i}), 1);
end

% Display the indices
disp('Indices of Significant Variables:');
disp(indices_of_significant_variables);

baseimage = spm_read_vols(spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']));
TstatImg = zeros(size(baseimage));
existingLabels = unique(baseimage);

cnt=0;
for labels=indices_of_significant_variables
%     if FSlabels2{labels,2}==0
%     else
        
        fprintf(['**' FSlabels2{labels,1}{1} '**\n'])
        
        cnt=cnt+1;
        
        % left
        disp('left')
        clear indL
        disp(['label:' num2str(FSlabels2{labels,2})])
        indL=find(baseimage==FSlabels2{labels,2});
        TstatImg(indL)=significant_results_table.T_Stat(cnt);
        significant_results_table.T_Stat(cnt)
        
        % right
        disp('right')
        clear indR
        disp(['label:' num2str(FSlabels2{labels,3})])
        indR=find(baseimage==FSlabels2{labels,3});
        TstatImg(indR)=significant_results_table.T_Stat(cnt);
        
%     end
end


hdr = spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
hdr.fname = ['/Users/alex/Dropbox/paperwriting/MRPET/data/TstatMaps_lowDAsession_BPchange_SRTM.nii'];
hdr.dim = size(TstatImg);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,TstatImg);

fprintf(['\n\n\n'])

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      LPNT                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% t-test with FDRc

clear significant_results_table p_values stats t_stats fdr_corrected_p_values summary_table significant_results_table ...
    significant_ROIs allVariableNames

% Extract the subset of variables
selectedVariables = MRPET_highSession_Table(:, 48:93);

% Value you want to test against
testValue = 0; % practically 0 in cerebellar cortex

% Initialize p-values and T statistics arrays
p_values = zeros(1, width(selectedVariables));
t_stats = zeros(1, width(selectedVariables));

% Perform one-sample t-test for each variable and collect p-values and T statistics
for i = 1:width(selectedVariables)
    [~, p_values(i), ~, stats] = ttest(selectedVariables{:, i}, testValue);
    t_stats(i) = stats.tstat;
end

% Correct p-values using the False Discovery Rate (FDR) method
fdr_corrected_p_values = mafdr(p_values, 'BHFDR', true);

% Create a summary table with variable names, p-values, and T statistics
summary_table = array2table([p_values; fdr_corrected_p_values; t_stats]', ...
    'VariableNames', {'Original_P_Value', 'Corrected_P_Value', 'T_Stat'}, ...
    'RowNames', selectedVariables.Properties.VariableNames);

% Filter the table to only include rows where the corrected p-value is less than 0.05
significant_results_table = summary_table(fdr_corrected_p_values < 0.05, :);

% Convert the numeric values to strings with four decimal places for p-values
significant_results_table.Original_P_Value = cellstr(num2str(significant_results_table.Original_P_Value, '%.4f'));
significant_results_table.Corrected_P_Value = cellstr(num2str(significant_results_table.Corrected_P_Value, '%.4f'));

% You can also format T_Stat if needed
% significant_results_table.T_Stat = cellstr(num2str(significant_results_table.T_Stat, '%.4f'));

% Display the table containing only significant results
disp('Significant Results:');
disp(significant_results_table);



% significance map!

% high-reward session
clear ZScoreImg baseimage FSlabel existingLabels hdr indices_of_significant_variables

FSlabels2=FSlabels1;
FSlabels2(1,:)=[];

allVariableNames=cellfun(@char,FSlabels2(:,1),'UniformOutput',false);

% Get the names of the significant variables from the table
significant_variable_names = significant_results_table.Properties.RowNames;
significant_ROIs           = cellfun(@(x) x(15:(end-7)),significant_variable_names,'Uniformoutput',false);

% Iterate over the significant_ROIs and find the matching index for each one
for i = 1:numel(significant_ROIs)
    indices_of_significant_variables(i) = find(strcmp(allVariableNames, significant_ROIs{i}), 1);
end

% Display the indices
disp('Indices of Significant Variables:');
disp(indices_of_significant_variables);

baseimage = spm_read_vols(spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']));
TstatImg = zeros(size(baseimage));
existingLabels = unique(baseimage);

cnt=0;
for labels=indices_of_significant_variables
    if FSlabels2{labels,2}==0
    else
        
        fprintf(['**' FSlabels2{labels,1}{1} '**\n'])
        
        cnt=cnt+1;
        
        % left
        disp('left')
        clear indL
        disp(['label:' num2str(FSlabels2{labels,2})])
        indL=find(baseimage==FSlabels2{labels,2});
        TstatImg(indL)=significant_results_table.T_Stat(cnt);
        significant_results_table.T_Stat(cnt)
        
        % right
        disp('right')
        clear indR
        disp(['label:' num2str(FSlabels2{labels,3})])
        indR=find(baseimage==FSlabels2{labels,3});
        TstatImg(indR)=significant_results_table.T_Stat(cnt);
        
    end
end


hdr = spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
hdr.fname = ['/Users/alex/Dropbox/paperwriting/MRPET/data/TstatMaps_highDAsession_BPchange_lpnt.nii'];
hdr.dim = size(TstatImg);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,TstatImg);

fprintf(['\n\n\n'])

%% lowDA session

clear significant_results_table p_values stats t_stats fdr_corrected_p_values summary_table significant_results_table ...
    significant_ROIs allVariableNames

% Extract the subset of variables
selectedVariables = MRPET_lowSession_Table(:, 48:93);

% Value you want to test against
testValue = 0; % practically 0 in cerebellar cortex

% Initialize p-values and T statistics arrays
p_values = zeros(1, width(selectedVariables));
t_stats = zeros(1, width(selectedVariables));

% Perform one-sample t-test for each variable and collect p-values and T statistics
for i = 1:width(selectedVariables)
    [~, p_values(i), ~, stats] = ttest(selectedVariables{:, i}, testValue);
    t_stats(i) = stats.tstat;
end

% Correct p-values using the False Discovery Rate (FDR) method
fdr_corrected_p_values = mafdr(p_values, 'BHFDR', true);

% Create a summary table with variable names, p-values, and T statistics
summary_table = array2table([p_values; fdr_corrected_p_values; t_stats]', ...
    'VariableNames', {'Original_P_Value', 'Corrected_P_Value', 'T_Stat'}, ...
    'RowNames', selectedVariables.Properties.VariableNames);

% Filter the table to only include rows where the corrected p-value is less than 0.05
significant_results_table = summary_table(fdr_corrected_p_values < 0.05, :);
%%% nothing significant after multiple-comparison correction

% 
% % Convert the numeric values to strings with four decimal places for p-values
% significant_results_table.Original_P_Value = cellstr(num2str(significant_results_table.Original_P_Value, '%.4f'));
% significant_results_table.Corrected_P_Value = cellstr(num2str(significant_results_table.Corrected_P_Value, '%.4f'));
% 
% % You can also format T_Stat if needed
% % significant_results_table.T_Stat = cellstr(num2str(significant_results_table.T_Stat, '%.4f'));
% 
% % Display the table containing only significant results
% disp('Significant Results:');
% disp(significant_results_table);
% 
% 
% 
% % significance map 2
% 
% % high-reward session
% clear ZScoreImg baseimage FSlabel existingLabels hdr
% 
% FSlabels2=FSlabels1;
% FSlabels2(1,:)=[];
% 
% 
% allVariableNames=cellfun(@char,FSlabels2(:,1),'UniformOutput',false);
% 
% % Get the names of the significant variables from the table
% significant_variable_names = significant_results_table.Properties.RowNames;
% significant_ROIs           = cellfun(@(x) x(15:(end-6)),significant_variable_names,'Uniformoutput',false);
% 
% % Iterate over the significant_ROIs and find the matching index for each one
% indices_of_significant_variables=[];
% for i = 1:numel(significant_ROIs)
%     indices_of_significant_variables(i) = find(strcmp(allVariableNames, significant_ROIs{i}), 1);
% end
% 
% % Display the indices
% disp('Indices of Significant Variables:');
% disp(indices_of_significant_variables);
% 
% baseimage = spm_read_vols(spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']));
% TstatImg = zeros(size(baseimage));
% existingLabels = unique(baseimage);
% 
% cnt=0;
% for labels=indices_of_significant_variables
% %     if FSlabels2{labels,2}==0
% %     else
%         
%         fprintf(['**' FSlabels2{labels,1}{1} '**\n'])
%         
%         cnt=cnt+1;
%         
%         % left
%         disp('left')
%         clear indL
%         disp(['label:' num2str(FSlabels2{labels,2})])
%         indL=find(baseimage==FSlabels2{labels,2});
%         TstatImg(indL)=significant_results_table.T_Stat(cnt);
%         significant_results_table.T_Stat(cnt)
%         
%         % right
%         disp('right')
%         clear indR
%         disp(['label:' num2str(FSlabels2{labels,3})])
%         indR=find(baseimage==FSlabels2{labels,3});
%         TstatImg(indR)=significant_results_table.T_Stat(cnt);
%         
% %     end
% end
% 
% 
% hdr = spm_vol(['/Users/alex/Dropbox/paperwriting/MRPET/data/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
% hdr.fname = ['/Users/alex/Dropbox/paperwriting/MRPET/data/TstatMaps_lowDAsession_BPchange_lpnt.nii'];
% hdr.dim = size(TstatImg);
% hdr = rmfield(hdr,'pinfo');
% hdr.nii = spm_write_vol(hdr,TstatImg);

fprintf(['\n\n\n'])

