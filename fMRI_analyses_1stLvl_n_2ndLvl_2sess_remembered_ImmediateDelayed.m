%% fMRI 1st-level analysis pipeline

%% work log

%   05-03-2019    created the script

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/alex/Documents/angela/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.physio  = [paths.parent 'E_data/EA_raw/EAE_physio/physio/3Tpilot/'];
paths.funx    = ['/Users/alex/Dropbox/literatures_IKND/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent];
paths.analyses= [paths.parent 'model_memory_immdel_epi2t1/'];
paths.behav   = '/Users/alex/Documents/angela/behav/';
load('/Users/alex/Dropbox/literatures_IKND/length_scan_allsubj.mat')
length_scan1=length_scan1(:,2)


% IDs
IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130 2131 2132 2233 2234 2235 2236]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

% IDs  = [2202 2204 2205 2206 2208 2109 2110 2112 2113 2114 2115 2116 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130]; %% 2203 2207 2217 delayed memory test lost
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
% d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

TR = 3.6;

% load experimental details
expdat = []; onsets_temp=[]; onsets=[];
for id = 1:length(IDs)
    
    % set up workspace
    mkdir([paths.analyses num2str(IDs(id))]);
    
    %% gather info
    
    for d = 1:2
        if days(id,d) == 0
            fname_beh{id,d} = {NaN};
            expdat{id,d} = {NaN};
            contingency{id,d} = [];
        else
            
            fname_beh{id,d}     = [ num2str(IDs(id)) '_' num2str(days(id,d)) '.mat' ];
            expdat{id,d} = load([paths.behav fname_beh{id,d}]);
            
            % define contingency
            eval(['contingency{id,d}(1) = str2num(expdat{id,d}.dat.day' num2str(d) '.RewardCategory(9))']);
            if contingency{id,d}(1) == 1
                contingency{id,d}(2) = 2;
            elseif contingency{id,d}(1) == 3
                contingency{id,d}(2) = 4;
            elseif contingency{id,d}(1) == 2
                contingency{id,d}(2) = 1;
            elseif contingency{id,d}(1) == 4
                contingency{id,d}(2) = 3;
            end
            
            % task info
            clear stim_task stim_del stim_imm stim_del_rem stim_imm_rem accuracies_imm accuracies_del stim_rewards stim_neutrals
            stim_task = eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.config.stim.fname']);
            stim_imm  = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl(:,[1 4])']);            
            accuracies_imm= eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.accu']);
            stim_imm_rem  = stim_imm(cell2mat(stim_imm(:,2))==1 & accuracies_imm==1,1); % old items and remembered
            stim_imm_for  = stim_imm(cell2mat(stim_imm(:,2))==1 & accuracies_imm==0,1); % old items and forgotten
            
            stim_task = eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.config.stim.fname']);
            if eval(['d' num2str(d) 'm(id,2)==0'])
            disp('no task data')
            else
            stim_del  = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.trl(:,[1 4])']);            
            accuracies_del= eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.accu']);
            stim_del_rem  = stim_del(cell2mat(stim_del(:,2))==1 & accuracies_del==1,1); % old items and remembered
            stim_del_for  = stim_del(cell2mat(stim_del(:,2))==1 & accuracies_del==0,1); % old items and forgotten
            end
            
            
            for q = 1:length(stim_task)
%                 idx = find(strcmp([C{:}], 'a'))
                ind_remembered_imm(q,1) = {find(strcmp(stim_task{q,1},stim_imm_rem))};
            end            
            tmp = cellfun(@isempty,ind_remembered_imm); indx_remembered_imm = tmp==0;
            
            for q = 1:length(stim_task)
%                 idx = find(strcmp([C{:}], 'a'))
                ind_forgotten_imm(q,1) = {find(strcmp(stim_task{q,1},stim_imm_for))};
            end            
            tmp = cellfun(@isempty,ind_forgotten_imm); indx_forgotten_imm = tmp==0;
            
            
            if eval(['d' num2str(d) 'm(id,2)==0'])
            disp('no task data')
            else
            for q = 1:length(stim_task)
%                 idx = find(strcmp([C{:}], 'a'))
                ind_remembered_del(q,1) = {find(strcmp(stim_task{q,1},stim_del_rem))};
            end            
            tmp = cellfun(@isempty,ind_remembered_del); indx_remembered_del = tmp==0;
            
            for q = 1:length(stim_task)
%                 idx = find(strcmp([C{:}], 'a'))
                ind_forgotten_del(q,1) = {find(strcmp(stim_task{q,1},stim_del_for))};
            end            
            tmp = cellfun(@isempty,ind_forgotten_del); indx_forgotten_del = tmp==0;
            end
            
            
            % sort onsets
            if eval(['d' num2str(d) 'm(id,2)==0']) % there is no delayed bit
            disp('no task data')
            onsets_temp{id,d}.stim.remembered  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(indx_remembered_imm==1 )-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.stim.forgotten   = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(indx_forgotten_imm==1 )-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.fb.remembered  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.cue(indx_remembered_imm==1 )-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.fb.forgotten   = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.cue(indx_forgotten_imm==1 )-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.stim.reward_rem  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.stim.reward_for  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.stim.neutral_rem = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.stim.neutral_for = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.fb.reward_rem    = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.fb.reward_for    = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);

            onsets_temp{id,d}.fb.neutral_rem   = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.fb.neutral_for   = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            
            
            else % there is delayed bit
                
            onsets_temp{id,d}.stim.remembered  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(indx_remembered_imm==1 | indx_remembered_del==1)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.stim.forgotten   = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(indx_forgotten_imm==1 | indx_forgotten_del==1)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.fb.remembered  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.cue(indx_remembered_imm==1 | indx_remembered_del==1)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.fb.forgotten   = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.cue(indx_forgotten_imm==1 | indx_forgotten_del==1)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.stim.reward_rem  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.stim.reward_for  = eval(['expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.stim.neutral_rem = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.stim.neutral_for = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.fb.reward_rem    = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.fb.reward_for    = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);

            onsets_temp{id,d}.fb.neutral_rem   = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            onsets_temp{id,d}.fb.neutral_for   = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            
            end
            
            onsets_temp{id,d}.fixX         = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.fix-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
            
            onsets_temp{id,d}.resp         = eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim+2.5-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
            
            
            if d==1
                if d1m(id,2)==0 % there is no delayed bit
                disp('no task data')
                onsets_temp{id,d}.stim.FreqReward_rem  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.FreqReward_for  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.stim.RareNeutral_rem = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.RareNeutral_for = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.FreqReward_rem    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.FreqReward_for    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.RareNeutral_rem   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.RareNeutral_for   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                else
                onsets_temp{id,d}.stim.FreqReward_rem  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.FreqReward_for  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.stim.RareNeutral_rem = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.RareNeutral_for = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.FreqReward_rem    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.FreqReward_for    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.RareNeutral_rem   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.RareNeutral_for   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                end
                
            elseif d==2
                if d2m(id,2)==0 % there is no delayed bit
                disp('no task data')
                onsets_temp{id,d}.stim.RareReward_rem  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.RareReward_for  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.stim.FreqNeutral_rem = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.FreqNeutral_for = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.RareReward_rem    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.RareReward_for    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.FreqNeutral_rem   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.FreqNeutral_for   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 ))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                else
                onsets_temp{id,d}.stim.RareReward_rem  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.RareReward_for  = eval(['expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.stim.FreqNeutral_rem = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.stim.FreqNeutral_for = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.RareReward_rem    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.RareReward_for    = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(1) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                
                onsets_temp{id,d}.fb.FreqNeutral_rem   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_remembered_imm==1 | indx_remembered_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                onsets_temp{id,d}.fb.FreqNeutral_for   = eval...
                    (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(cell2mat(expdat{id,d}.dat.day' num2str(d)...
                    '.maintask.results.trl(:,2))==contingency{id,d}(2) & (indx_forgotten_imm==1 | indx_forgotten_del==1))-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
                end
                
            end
        end
    end
    %% make physio + realignment parameters
    
    
    % ----- already made ------ %
    
%     % read physio parameters
%     physioFile  = [paths.physio num2str(IDs(id)) '/' num2str(IDs(id)) '.txt'];
%     phystemp    = importdata(physioFile);
%     physiodata  = phystemp.data; clear phystemp
%     
%     % read realignment parameters
%     realignFile = [paths.preproc num2str(IDs(id)) '/rp_all_ua' num2str(IDs(id)) '.txt'];
%     realigndata = importdata(realignFile);
%     
%     % concatenate
%     R = [physiodata realigndata];
%     
%     save([paths.analyses num2str(IDs(id)) '/multiple_regressors.mat'],'R');
%     clear R realigndata realignFile physiodata physioFile
    
    
    %             % physio+realignment regressors
    %             Rs{id,d}=load([paths.preproc num2str(IDs(id)) '_' num2str(d) '/multiple_regressors.mat']);
    %             if d==1
    %                 Rs{id,d}.R= [Rs{id,d}.R, ones(length(Rs{id,d}.R),1)];
    %             else
    %                 Rs{id,d}.R= [Rs{id,d}.R, zeros(length(Rs{id,d}.R),1)];
    %             end
    
    %% merge info together
    
    
    % merge onsets
    if days(id,1) == 0
        %Freq reward and rare neutral
        
        %remembered
        onsets{id,1}.stim.remembered      = [0];
        onsets{id,1}.fb.remembered        = [0];
        onsets{id,1}.stim.rewardFreq_rem  = [0];
        onsets{id,1}.stim.neutralRare_rem = [0];
        onsets{id,1}.fb.rewardFreq_rem    = [0];
        onsets{id,1}.fb.neutralRare_rem   = [0];
        
        %forgotten
        onsets{id,1}.stim.forgotten       = [0];
        onsets{id,1}.fb.forgotten         = [0];
        onsets{id,1}.stim.rewardFreq_for  = [0];
        onsets{id,1}.stim.neutralRare_for = [0];
        onsets{id,1}.fb.rewardFreq_for    = [0];
        onsets{id,1}.fb.neutralRare_for   = [0];
        
        %remembered
        onsets{id,1}.stim.rewardRare_rem  = [onsets_temp{id,2}.stim.reward_rem+length_scan1(id)];
        onsets{id,1}.stim.neutralFreq_rem = [onsets_temp{id,2}.stim.neutral_rem+length_scan1(id)];
        onsets{id,1}.fb.rewardRare_rem    = [onsets_temp{id,2}.fb.reward_rem+length_scan1(id)];
        onsets{id,1}.fb.neutralFreq_rem   = [onsets_temp{id,2}.fb.neutral_rem+length_scan1(id)];
        
        %forgotten
        onsets{id,1}.stim.rewardRare_for  = [onsets_temp{id,2}.stim.reward_for+length_scan1(id)];
        onsets{id,1}.stim.neutralFreq_for = [onsets_temp{id,2}.stim.neutral_for+length_scan1(id)];
        onsets{id,1}.fb.rewardRare_for    = [onsets_temp{id,2}.fb.reward_for+length_scan1(id)];
        onsets{id,1}.fb.neutralFreq_for   = [onsets_temp{id,2}.fb.neutral_for+length_scan1(id)];
        
        onsets{id,1}.fixX         = [onsets_temp{id,2}.fixX+length_scan1(id)];
        onsets{id,1}.resp         = [onsets_temp{id,2}.resp+length_scan1(id)];
        
    elseif days(id,2) == 0
        %Freq reward and rare neutral
        
        %remembered
        onsets{id,1}.stim.remembered      = [onsets_temp{id,1}.stim.remembered];
        onsets{id,1}.fb.remembered        = [onsets_temp{id,1}.fb.remembered];
        onsets{id,1}.stim.rewardFreq_rem  = [onsets_temp{id,1}.stim.reward_rem];
        onsets{id,1}.stim.neutralRare_rem = [onsets_temp{id,1}.stim.neutral_rem];
        onsets{id,1}.fb.rewardFreq_rem    = [onsets_temp{id,1}.fb.reward_rem];
        onsets{id,1}.fb.neutralRare_rem   = [onsets_temp{id,1}.fb.neutral_rem];
        %forgotten
        onsets{id,1}.stim.forgotten       = [onsets_temp{id,1}.stim.forgotten];
        onsets{id,1}.fb.forgotten         = [onsets_temp{id,1}.fb.forgotten];
        onsets{id,1}.stim.rewardFreq_for  = [onsets_temp{id,1}.stim.reward_for];
        onsets{id,1}.stim.neutralRare_for = [onsets_temp{id,1}.stim.neutral_for];
        onsets{id,1}.fb.rewardFreq_for    = [onsets_temp{id,1}.fb.reward_for];
        onsets{id,1}.fb.neutralRare_for   = [onsets_temp{id,1}.fb.neutral_for];
        
        %remembered
        onsets{id,1}.stim.rewardRare_rem  = [0];
        onsets{id,1}.stim.neutralFreq_rem = [0];
        onsets{id,1}.fb.rewardRare_rem    = [0];
        onsets{id,1}.fb.neutralFreq_rem   = [0];
        %forgotten
        onsets{id,1}.stim.rewardRare_for  = [0];
        onsets{id,1}.stim.neutralFreq_for = [0];
        onsets{id,1}.fb.rewardRare_for    = [0];
        onsets{id,1}.fb.neutralFreq_for   = [0];
        
        onsets{id,1}.fixX         = [onsets_temp{id,1}.fixX];
        onsets{id,1}.resp         = [onsets_temp{id,1}.resp];
        
    else
        %Freq reward and rare neutral
        
        onsets{id,1}.stim.remembered = [onsets_temp{id,1}.stim.remembered; onsets_temp{id,2}.stim.remembered+length_scan1(id)];
        onsets{id,1}.stim.forgotten  = [onsets_temp{id,1}.stim.forgotten; onsets_temp{id,2}.stim.forgotten+length_scan1(id)];
        onsets{id,1}.fb.remembered = [onsets_temp{id,1}.fb.remembered; onsets_temp{id,2}.fb.remembered+length_scan1(id)];
        onsets{id,1}.fb.forgotten  = [onsets_temp{id,1}.fb.forgotten; onsets_temp{id,2}.fb.forgotten+length_scan1(id)];

        %remembered
        onsets{id,1}.stim.rewardFreq_rem  = [onsets_temp{id,1}.stim.reward_rem];
        onsets{id,1}.stim.neutralRare_rem = [onsets_temp{id,1}.stim.neutral_rem];
        onsets{id,1}.fb.rewardFreq_rem    = [onsets_temp{id,1}.fb.reward_rem];
        onsets{id,1}.fb.neutralRare_rem   = [onsets_temp{id,1}.fb.neutral_rem];
        %forgotten
        onsets{id,1}.stim.rewardFreq_for  = [onsets_temp{id,1}.stim.reward_for];
        onsets{id,1}.stim.neutralRare_for = [onsets_temp{id,1}.stim.neutral_for];
        onsets{id,1}.fb.rewardFreq_for    = [onsets_temp{id,1}.fb.reward_for];
        onsets{id,1}.fb.neutralRare_for   = [onsets_temp{id,1}.fb.neutral_for];
        
        %remembered
        onsets{id,1}.stim.rewardRare_rem  = [onsets_temp{id,2}.stim.reward_rem+length_scan1(id)];
        onsets{id,1}.stim.neutralFreq_rem = [onsets_temp{id,2}.stim.neutral_rem+length_scan1(id)];
        onsets{id,1}.fb.rewardRare_rem    = [onsets_temp{id,2}.fb.reward_rem+length_scan1(id)];
        onsets{id,1}.fb.neutralFreq_rem   = [onsets_temp{id,2}.fb.neutral_rem+length_scan1(id)];
        %forgotten
        onsets{id,1}.stim.rewardRare_for  = [onsets_temp{id,2}.stim.reward_for+length_scan1(id)];
        onsets{id,1}.stim.neutralFreq_for = [onsets_temp{id,2}.stim.neutral_for+length_scan1(id)];
        onsets{id,1}.fb.rewardRare_for    = [onsets_temp{id,2}.fb.reward_for+length_scan1(id)];
        onsets{id,1}.fb.neutralFreq_for   = [onsets_temp{id,2}.fb.neutral_for+length_scan1(id)];
        
        onsets{id,1}.fixX         = [onsets_temp{id,1}.fixX; onsets_temp{id,2}.fixX+length_scan1(id)];
        onsets{id,1}.resp         = [onsets_temp{id,1}.resp; onsets_temp{id,2}.resp+length_scan1(id)];
        
        
    end
     
end

fprintf('\n preparation done \n')
%% start, stim

spm fmri  % open progress window

for id = 1%:length(IDs) % id 15 --> no forgotten items for rare reward : figure out something for this id
    %% build model
    
%     if IDs(id) == 2202 && d==1
%         expdat{id,d}.dat.day1.maintask.results.SOT.raw.cue=...
%             expdat{id,d}.dat.day1.maintask.results.SOT.raw.stim+4.5+((expdat{id,d}.dat.day1.maintask.config.timing.dot)/1000)'
%     else
%     end
    
    clear dir_model list_scan volnum
    dir_model = [paths.analyses num2str(IDs(id)) '/'];
    volnum = length(spm_vol([paths.preproc '/s3func' num2str(IDs(id)) '.nii']));
    
    for v1 = 1:volnum
        list_scan{v1,1} = [paths.preproc '/s3func' num2str(IDs(id)) '.nii,' num2str(v1)];
    end
    clear matlabbatch
    
    spm_jobman('initcfg');      % initiate job manager
    
    matlabbatch{1}.spm.stats.fmri_spec.dir               = cellstr(dir_model); % where will the model be saved?
    matlabbatch{1}.spm.stats.fmri_spec.timing.units      = 'secs'; % scans / secs
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT         = TR; % TR in seconds
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t     = 51; % volume size
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0    = 25; % microtime onset
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans        = list_scan;
    
    
    % remembered
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'stim_remembered';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset= onsets{id,1}.stim.remembered;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'fb_remembered';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.fb.remembered;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    
    % forgotten
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'stim_forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= onsets{id,1}.stim.forgotten;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'fb_forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.fb.forgotten;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
    
    % else
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'FixX';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= onsets{id,1}.fixX;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'Response';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= onsets{id,1}.resp;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
    
    %
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr([paths.preproc num2str(IDs(id+29)) '_multiple_regressors.mat']); % movement parameters
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    spm_jobman('run', matlabbatch) % run batch
    
    
    %% estimate
    
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([dir_model 'SPM.mat']);
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch) % run batch
    
    %% statistical inference
    
    firstlvl_dir        = [dir_model 'SPM.mat'];
    
    % stims
    c1_stim_rem_for          = [1 -1];
    c1_stim_rem              = [1];
    c1_stim_for              = [0 1];
    
    c1_fb_rem_for            = [1 -1];
    c1_fb_rem                = [1];
    c1_fb_for                = [0 1];
    
    % batch setup
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim: remembered vs forgotten';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = c1_stim_rem_for;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'stim: remembered';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = c1_stim_rem;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'stim: forgotten';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = c1_stim_for;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'fb: remembered vs forgotten';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = c1_fb_rem_for;
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'fb: remembered ';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = c1_fb_rem;
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'fb: forgotten';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = c1_fb_for;
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    
    matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
    
    spm_jobman('run', matlabbatch) % run batch
    
end


fprintf('\ndone\n')

%% fMRI 2nd-level analysis pipeline
%% set environmental variables

clear;
warning('off','all');

% paths

paths = [];
paths.parent  = '/Users/alex/Documents/angela/';
paths.funx    = ['/Users/alex/Dropbox/literatures_IKND/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent];
paths.analyses= [paths.parent 'model_memory_immdel_epi2t1/'];
paths.behav   = '/Users/alex/Documents/angela/behav/';
paths.temp    = [paths.parent 'model_memory_immdel_epi2t1/2nd/'];
paths.save2nd = [paths.parent 'model_memory_immdel_epi2t1/2nd/'];

% IDs
IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130 2131 2132 2233 2234 2235 2236]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

TR = 3.6;

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window


list_stim_remembered_forgotten= []; cnt1=0;
for id = 1:length(IDs)
            cnt1 = cnt1+1;
            list_stim_remembered_forgotten{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0001_mni.nii,1'];
end
clear i1 d cnt1

list_stim_Rem = []; cnt1=0;
for id = 1:length(IDs)
            cnt1 = cnt1+1;
            list_stim_Rem{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0002_mni.nii,1'];
end
clear i1 d cnt1

list_stim_for = []; cnt1=0;
for id = 1:length(IDs)
            cnt1 = cnt1+1;
            list_stim_for{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0003_mni.nii,1'];
end
clear i1 d cnt1

list_fb_rem_v_for = []; cnt1=0;
for id = 1:length(IDs)
            cnt1=cnt1+1;
            list_fb_rem_v_for{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0004_mni.nii,1'];
end
clear i1 d cnt1

list_fb_rem = []; cnt1=0;
for id = 1:length(IDs)
            cnt1=cnt1+1;
            list_fb_rem{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0005_mni.nii,1'];
end
clear i1 d cnt1

list_fb_for = []; cnt1=0;
for id = 1:length(IDs)
            cnt1=cnt1+1;
            list_fb_for{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0006_mni.nii,1'];
end
clear i1 d cnt1


%% compute models

% ------- compute: stim vs. fixation cross ------- %

cd(paths.save2nd); mkdir('stim_remembered_forgotten'); cd stim_remembered_forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_remembered_forgotten) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim: remembered > forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > neutral ------- %

cd(paths.save2nd); mkdir('stim_remembered'); cd stim_remembered; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_Rem) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim: remembered';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: neutral > reward ------- %
cd(paths.save2nd); mkdir('stim_forgotten'); cd stim_forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_for) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim: forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %
%


% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('fb_remembered_forgotten'); cd fb_remembered_forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_rem_v_for) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb: remembered > forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('fb_remembered'); cd fb_remembered; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_rem) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb: remembered ';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > fix ------- %

cd(paths.save2nd); mkdir('fb_forgotten'); cd fb_forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_for) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb: forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

%% set environmental variables

clear;
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.physio  = [paths.parent 'E_data/EA_raw/EAE_physio/physio/3Tpilot/'];
paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAB_MRI/EABB_preprocessed/MainTask_all/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.analyses= [paths.parent 'E_data/EB_cleaned/EBC_mri/pilot3T_modelBF/'];
paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/pilot_3T_2sess/tmp/';
paths.temp    = [paths.parent 'E_data/EB_cleaned/EBC_mri/pilot3T_modelBF/'];
paths.save2nd = [paths.parent 'E_data/EB_cleaned/EBC_mri/pilot3T_modelBF/2ndLvL/'];
paths.doc     = '/Users/yeojin/Desktop/C_writings/CB_figures/pilot_3T_2sess/MainTask/2ndLvl_uncorr/';

% add toolboxes and functions
addpath(paths.spm)
addpath(paths.funx)

IDs  = [2202 2204 2205 2206 2208 2109 2110 2112 2113 2114 2115 2116 2218 2219 2221 2222 2223 2224 2125 2126]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

TR = 3.6;

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window


list_fb_remembered_forgotten= []; cnt1=0;
for id = 1:length(IDs)
            cnt1 = cnt1+1;
            list_fb_remembered_forgotten{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0001_mni.nii,1'];
end
clear i1 d cnt1

list_fb_rewRem_neuRem = []; cnt1=0;
for id = 1:length(IDs)
            cnt1 = cnt1+1;
            list_fb_rewRem_neuRem{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0002_mni.nii,1'];
end
clear i1 d cnt1

list_fb_rewRem_rewFor = []; cnt1=0;
for id = 1:length(IDs)
            cnt1 = cnt1+1;
            list_fb_rewRem_rewFor{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0003_mni.nii,1'];
end
clear i1 d cnt1

list_fb_neuRem_neuFor = []; cnt1=0;
for id = 1:length(IDs)
            cnt1=cnt1+1;
            list_fb_neuRem_neuFor{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0004_mni.nii,1'];
end
clear i1 d cnt1

list_fb_rareRem_rareFor = []; cnt1=0;
for id = 1:length(IDs)
            cnt1=cnt1+1;
            list_fb_rareRem_rareFor{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0005_mni.nii,1'];
end
clear i1 d cnt1

list_fb_freqRem_freqFor = []; cnt1=0;
for id = 1:length(IDs)
            cnt1=cnt1+1;
            list_fb_freqRem_freqFor{cnt1,1} = [paths.analyses num2str(IDs(id)) '/con_0006_mni.nii,1'];
end
clear i1 d cnt1

% list_fb_forgotten_rew = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_fb_forgotten_rew{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0007_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_fb_forgotten_neu = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_fb_forgotten_neu{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0008_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_fb_rewards_remembered_forgotten = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_fb_rewards_remembered_forgotten{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0009_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% 
% list_fb_neutrals_remembered_forgotten = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_fb_neutrals_remembered_forgotten{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0010_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_remembered_rew_neu = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_remembered_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0011_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_remembered_neu_rew = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_remembered_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0012_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_remembered_rew = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_remembered_rew{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0013_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_remembered_neu = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_remembered_neu{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0014_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_forgotten_rew_neu = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_forgotten_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0015_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_forgotten_neu_rew = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_forgotten_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0016_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_forgotten_rew = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_forgotten_rew{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0017_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_forgotten_neu = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_forgotten_neu{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0018_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% list_feedback_rewards_remembered_forgotten = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_rewards_remembered_forgotten{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0019_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1
% 
% 
% list_feedback_neutrals_remembered_forgotten = []; cnt1=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%         else
%             cnt1=cnt1+1;
%             list_feedback_neutrals_remembered_forgotten{cnt1,1} = [paths.analyses num2str(IDs(id)) '_' num2str(days(id,d)) '/con_0020_mni.nii,1'];
%         end
%     end
% end
% clear i1 d cnt1


%% compute models

% ------- compute: fb vs. fixation cross ------- %

cd(paths.save2nd); mkdir('fb_remembered_forgotten'); cd fb_remembered_forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_remembered_forgotten) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'remembered > forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > neutral ------- %

cd(paths.save2nd); mkdir('fb_rewRem_neuRem'); cd fb_rewRem_neuRem; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_rewRem_neuRem) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'remembered: rew > neu';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: neutral > reward ------- %
cd(paths.save2nd); mkdir('fb_rewRem_rewFor'); cd fb_rewRem_rewFor; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_rewRem_rewFor) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'rewards: remembered vs forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('fb_neuRem_neuFor'); cd fb_neuRem_neuFor; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_neuRem_neuFor) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'neutrals: remembered > forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('fb_rareRem_rareFor'); cd fb_rareRem_rareFor; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_rareRem_rareFor) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'rare items: remembered > forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > fix ------- %

cd(paths.save2nd); mkdir('fb_freqRem_freqFor'); cd fb_freqRem_freqFor; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_freqRem_freqFor) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'frequent items: remembered > forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %






%% generate results report

clc;clear;

IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 0 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

path_2nd = '/Users/yeojin/Desktop/E_data/EB_cleaned/EBC_mri/pilot3T_remembered/2ndLvL/';
tmp = dir([path_2nd]); list_2nd = tmp(~startsWith({tmp.name},'.'));
clear tmp

SPMmat=[]; mc=0;
for l1 = 1:length(list_2nd)
    if strcmp(list_2nd(l1).name,'compare_days')
        mc=mc+1; SPMmat{mc,1} = load([path_2nd list_2nd(l1).name '/stim/Rew_Neu/SPM.mat']);
        mc=mc+1; SPMmat{mc,1} = load([path_2nd list_2nd(l1).name '/stim/Neu_Rew/SPM.mat']);
        mc=mc+1; SPMmat{mc,1} = load([path_2nd list_2nd(l1).name '/fb/Rew_Neu/SPM.mat']);
        mc=mc+1; SPMmat{mc,1} = load([path_2nd list_2nd(l1).name '/fb/Neu_Rew/SPM.mat']);
    else
        mc=mc+1;
        SPMmat{mc,1} = load([path_2nd list_2nd(l1).name '/SPM.mat']);
    end
end

%% non cluster-corrected

FDRmat=[]; fdrc=0;
for l2 = 1:length(SPMmat)
    if ~isempty(strfind(SPMmat{l2,1}.SPM.swd,'compare_days'))
        for c1 = 1:length(SPMmat{l2,1}.SPM.xCon)
            clear matlabbatch
            spm_jobman('initcfg');
            matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l2,1}.SPM.swd '/SPM.mat']);
            matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l2,1}.SPM.xCon(c1).name;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = c1;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
            matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/c1mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
            matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = [ SPMmat{l2,1}.SPM.xCon(c1).name '_p01_c0_grey'];
            spm_jobman('run', matlabbatch) % run batch
            
            fdrc=fdrc+1;
            FDRmat{fdrc,1}=TabDat;
            
        end
        
    elseif ~isempty(strfind(SPMmat{l2,1}.SPM.swd,'anova'))
        
        for c2 = 1:length(SPMmat{l2,1}.SPM.xCon)
            clear matlabbatch
            spm_jobman('initcfg');
            matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l2,1}.SPM.swd '/SPM.mat']);
            matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l2,1}.SPM.xCon(c2).name;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = c2;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
            matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/c1mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
            matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            matlabbatch{1}.spm.stats.results.units = 1;
            clear slashes basefname
            slashes = strfind([ SPMmat{l2,1}.SPM.xCon(c2).name '_p01_c0_grey'], '/');
            basefname = [ SPMmat{l2,1}.SPM.xCon(c2).name '_p01_c0_grey'];
            basefname(slashes)='_';
            matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = basefname;
            spm_jobman('run', matlabbatch) % run batch
            
            fdrc=fdrc+1;
            FDRmat{fdrc,1}=TabDat;
        end
    else
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l2,1}.SPM.swd '/SPM.mat']);
        matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l2,1}.SPM.xCon.name;
        matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
        matlabbatch{1}.spm.stats.results.conspec.extent = 0;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/c1mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
        matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = [ SPMmat{l2,1}.SPM.xCon.name '_p01_c0_grey'];
        spm_jobman('run', matlabbatch) % run batch
        
        fdrc=fdrc+1;
        FDRmat{fdrc,1}=TabDat;
        
    end
    
end

save('/Users/yeojin/Desktop/E_data/EB_cleaned/EBC_mri/pilot3T_remembered/results_nii/FDRmat.mat','FDRmat')


% cluster corrected by FDRc
load('/Users/yeojin/Desktop/E_data/EB_cleaned/EBC_mri/pilot3T_remembered/results_nii/FDRmat.mat')
cc=0;
for l3 = 1:length(SPMmat)
    if ~isempty(strfind(SPMmat{l3,1}.SPM.swd,'compare_days'))
        for c1 = 1:length(SPMmat{l3,1}.SPM.xCon)
            clear matlabbatch
            spm_jobman('initcfg');
            matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l3,1}.SPM.swd '/SPM.mat']);
            matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l3,1}.SPM.xCon(c1).name;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = c1;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
            cc=cc+1;
            matlabbatch{1}.spm.stats.results.conspec.extent = FDRmat{cc,1}.ftr{5,2}(4);
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/c1mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
            matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = [ SPMmat{l3,1}.SPM.xCon(c1).name '_p01_c' num2str(FDRmat{cc,1}.ftr{5,2}(4)) '_grey'];
            spm_jobman('run', matlabbatch) % run batch
            
        end
        
    elseif ~isempty(strfind(SPMmat{l3,1}.SPM.swd,'anova'))
        
        for c2 = 2:length(SPMmat{l3,1}.SPM.xCon)
            clear matlabbatch
            spm_jobman('initcfg');
            matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l3,1}.SPM.swd '/SPM.mat']);
            matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l3,1}.SPM.xCon(c2).name;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = c2;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
            cc=cc+1;
            matlabbatch{1}.spm.stats.results.conspec.extent = FDRmat{cc,1}.ftr{5,2}(4);
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/c1mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
            matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            matlabbatch{1}.spm.stats.results.units = 1;
            clear slashes basefname
            slashes = strfind([ SPMmat{l3,1}.SPM.xCon(c2).name '_p01_c' num2str(FDRmat{cc,1}.ftr{5,2}(4)) '_grey'], '/');
            basefname = [ SPMmat{l3,1}.SPM.xCon(c2).name '_p01_c' num2str(FDRmat{cc,1}.ftr{5,2}(4)) '_grey'];
            basefname(slashes)='_';
            matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = basefname;
            spm_jobman('run', matlabbatch) % run batch
            
        end
    else
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l3,1}.SPM.swd '/SPM.mat']);
        matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l3,1}.SPM.xCon.name;
        matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
        cc=cc+1;
        matlabbatch{1}.spm.stats.results.conspec.extent = FDRmat{cc,1}.ftr{5,2}(4);
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/c1mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
        matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = [ SPMmat{l3,1}.SPM.xCon.name '_p01_c' num2str(FDRmat{cc,1}.ftr{5,2}(4)) '_grey'];
        spm_jobman('run', matlabbatch) % run batch
        
    end
    
end


% brainstem mask
for l4 = 1:length(SPMmat)
    if ~isempty(strfind(SPMmat{l4,1}.SPM.swd,'compare_days'))
        for c1 = 1:length(SPMmat{l4,1}.SPM.xCon)
            clear matlabbatch
            spm_jobman('initcfg');
            matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l4,1}.SPM.swd '/SPM.mat']);
            matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l4,1}.SPM.xCon(c1).name;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = c1;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
            matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/brainstemMask_short.nii,1'};
            matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = [ SPMmat{l4,1}.SPM.xCon(c1).name '_p01_c0_Brainstem'];
            spm_jobman('run', matlabbatch) % run batch
            
        end
        
    elseif ~isempty(strfind(SPMmat{l4,1}.SPM.swd,'anova'))
        
        for c2 = 1:length(SPMmat{l4,1}.SPM.xCon)
            clear matlabbatch
            spm_jobman('initcfg');
            matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l4,1}.SPM.swd '/SPM.mat']);
            matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l4,1}.SPM.xCon(c2).name;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = c2;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
            matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/brainstemMask_short.nii,1'};
            matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            matlabbatch{1}.spm.stats.results.units = 1;
            clear slashes basefname
            slashes = strfind([ SPMmat{l4,1}.SPM.xCon(c2).name '_p01_c0_Braianstem'], '/');
            basefname = [ SPMmat{l4,1}.SPM.xCon(c2).name '_p01_c0_Brainstem'];
            basefname(slashes)='_';
            matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = basefname;
            spm_jobman('run', matlabbatch) % run batch
            
        end
    else
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMmat{l4,1}.SPM.swd '/SPM.mat']);
        matlabbatch{1}.spm.stats.results.conspec.titlestr = SPMmat{l4,1}.SPM.xCon.name;
        matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
        matlabbatch{1}.spm.stats.results.conspec.extent = 0;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/brainstemMask_short.nii,1'};
        matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = [ SPMmat{l4,1}.SPM.xCon.name '_p01_c0_Brainstem'];
        spm_jobman('run', matlabbatch) % run batch
        
    end
    
end



% move the result images

for f1 = 1:length(SPMmat)

    path_source = [SPMmat{f1,1}.SPM.swd];
    path_dest   = ['/Users/yeojin/Desktop/E_data/EB_cleaned/EBC_mri/pilot3T_remembered/results_nii/'];
    
    clear tmplist
    tmplist     = dir([path_source '/spmT_00*_*.nii']);
    
    for ff1 = 1:length(tmplist)
        movefile(fullfile(tmplist(ff1).folder,tmplist(ff1).name),...
            fullfile(path_dest,tmplist(ff1).name))
    end
    
    clear tmplist
    tmplist     = dir([path_source '/spmF_00*_*.nii']);
    
    for ff2 = 1:length(tmplist)
        movefile(fullfile(tmplist(ff2).folder,tmplist(ff2).name),...
            fullfile(path_dest,tmplist(ff2).name))
    end
    
end



%% calculate tSNR

for k = 1:length(IDs)
    
    file_end = ['/Users/yeojin/Desktop/E_data/EB_cleaned/EBC_mri/pilot_3T_2sess/MainTask/' num2str(IDs(k)) '_2/'];
    cd(file_end)
    
    disp('tSNr computation')
    FileNameSPMMat = {fullfile(file_end,'SPM.mat')};
    load(char(FileNameSPMMat))
    % S=size(SPM.xX.X);
    numbeta = 1;
    
    % matlabbatch{1}.spm.util.imcalc.input = {fullfile(file_end,strcat('beta_',sprintf('%04d',S(2)),'.nii'));fullfile(file_end,'ResMS.nii')};
    matlabbatch{1}.spm.util.imcalc.input = {fullfile(file_end,strcat('beta_',sprintf('%04d',numbeta(1)),'.nii'));fullfile(file_end,'ResMS.nii')};
    
    matlabbatch{1}.spm.util.imcalc.output = ['tSNR_' num2str(IDs(k)) '_2_' strcat('beta_',sprintf('%04d',numbeta(1)))];
    matlabbatch{1}.spm.util.imcalc.outdir = {file_end};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1./sqrt(i2))';
    
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -3;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
end