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
paths.analyses= [paths.parent 'model_TaskMemory_immdel_rarefreq_stim_epi2t1/'];
paths.behav   = '/Users/alex/Documents/angela/behav/';
load('/Users/alex/Dropbox/literatures_IKND/length_scan_allsubj.mat')
length_scan1=length_scan1(:,2)


% IDs
IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130 2131 2132 2233 2234 2235 2236]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

% IDs  = [2202 2204 2205 2206 2208 2109 2110 2112 2113 2114 2115 2116 2218 2219 2220 2221 2222 2223 2224 2125 2126]; %% 2203 2207 2217 delayed memory test lost
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
% d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

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

for id = 28:length(IDs) % id 15 --> no forgotten items for rare reward : figure out something for this id
    %% build model
    
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
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'stim_RewardRare_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset= onsets{id,1}.stim.rewardRare_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'stim_RewardFrequent_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.stim.rewardFreq_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'stim_NeutralRare_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= onsets{id,1}.stim.neutralRare_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'stim_NeutralFrequent_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.stim.neutralFreq_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;

    % forgotten
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'stim_RewardRare_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= onsets{id,1}.stim.rewardRare_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'stim_RewardFrequent_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= onsets{id,1}.stim.rewardFreq_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name = 'stim_NeutralRare_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset= onsets{id,1}.stim.neutralRare_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name = 'stim_NeutralFrequent_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset= onsets{id,1}.stim.neutralFreq_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth = 1;

    
    % else
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name = 'FixX';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset= onsets{id,1}.fixX;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).orth = 1;
   
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).name = 'Response';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).onset= onsets{id,1}.resp;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).orth = 1;
    
    %
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr([paths.preproc num2str(IDs(id)) '_multiple_regressors.mat']); % movement parameters
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
    stim_rew_Remembered_Forgotten = [1 1 0 0 -1 -1];
    stim_rew_Forgotten_Remembered = [-1 -1 0 0 1 1];
    stim_neu_Remembered_Forgotten = [0 0 1 1 0 0 -1 -1];
    stim_neu_Forgotten_Remembered = [0 0 -1 -1 0 0 1 1];
    
    stim_rare_Remembered_Forgotten = [1 0 1 0 -1 0 -1];
    stim_rare_Forgotten_Remembered = [-1 0 -1 0 1 0 1];
    stim_freq_Remembered_Forgotten = [0 1 0 1 0 -1 0 -1];
    stim_freq_Forgotten_Remembered = [0 -1 0 -1 0 1 0 1];
    
    
%     stim_rem_Reward_Neutral  = [1 1 -1 -1];
    stim_rem_Rare_Frequent   = [1 -1 1 -1];
%     stim_rem_Neutral_Reward  = [-1 -1 1 1];
    stim_rem_Frequent_Rare   = [-1 1 -1 1];
%     
%     stim_for_Reward_Neutral  = [zeros(1,4) 1 1 -1 -1];
    stim_for_Rare_Frequent   = [zeros(1,4) 1 -1 1 -1];
%     stim_for_Neutral_Reward  = [zeros(1,4) -1 -1 1 1];
    stim_for_Frequent_Rare   = [zeros(1,4) -1 1 -1 1];
%     
%     stim_rem_RewardRare_RewardFreq     = [1 -1];
%     stim_rem_NeutralRare_NeutralFreq   = [0 0 1 -1];
%     stim_rem_RewardFreq_RewardRare     = [-1 1];
%     stim_rem_NeutralFreq_NeutralRare   = [0 0 -1 1];
%     
%     stim_for_RewardRare_RewardFreq     = [zeros(1,4) 1 -1];
%     stim_for_NeutralRare_NeutralFreq   = [zeros(1,4) 0 0 1 -1];
%     stim_for_RewardFreq_RewardRare     = [zeros(1,4) -1 1];
%     stim_for_NeutralFreq_NeutralRare   = [zeros(1,4) 0 0 -1 1];

    
    % batch setup
    
    if IDs(id) == 2220
        
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'stim, Neutrals: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = stim_neu_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'stim, Neutrals: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = stim_neu_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'stim, Frequent: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = stim_freq_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'stim, Frequent: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = stim_freq_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
    
%     matlabbatch{1}.spm.stats.con.delete = 0;    % 1 means delete, 0 means append
    matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
    
    spm_jobman('run', matlabbatch) % run batch
    
    else
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim, Rewards: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = stim_rew_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'stim, Rewards: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = stim_rew_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'stim, Neutrals: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = stim_neu_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'stim, Neutrals: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = stim_neu_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';    
    
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'stim, Rare: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = stim_rare_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'stim, Rare: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = stim_rare_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'stim, Frequent: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = stim_freq_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'stim, Frequent: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = stim_freq_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'stim, Remembered: Rare > Frequent';
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = stim_rem_Rare_Frequent;
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'stim, Forgotten: Rare > Frequent';
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = stim_for_Rare_Frequent;
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
    
%     matlabbatch{1}.spm.stats.con.delete = 0;    % 1 means delete, 0 means append
    matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append    
    spm_jobman('run', matlabbatch) % run batch
    
    
    end
    
    
end


fprintf('\ndone\n')

%%
% for id=1:length(IDs)
%     mkdir(['/Users/alex/Documents/angela/model_TaskMemory_immdel_stim/coreg/' num2str(IDs(id)) '/data/'])
%     for ctr=1:8
%         copyfile([paths.analyses num2str(IDs(id)) '/' strcat('con_',sprintf('%04d',ctr),'.nii')],...
%             ['/Users/alex/Documents/angela/model_TaskMemory_immdel_stim/coreg/' num2str(IDs(id)) '/data/' strcat('con_',sprintf('%04d',ctr),'.nii')])
%     end
% end
% 
%%
for id=1:length(IDs)
    for ctr=9:10
        copyfile([paths.analyses num2str(IDs(id)) '/' strcat('con_',sprintf('%04d',ctr),'.nii')],...
            ['/Users/alex/Desktop/pilotmemory/' num2str(IDs(id)) '/data/' strcat('con_',sprintf('%04d',ctr),'_stim.nii')])
    end
end


%%
for id=15:length(IDs)
    for ctr=9:10
        copyfile(['/Users/alex/Desktop/pilotmemory/' num2str(IDs(id)) '/data/' strcat('con_',sprintf('%04d',ctr),'_stim_mni.nii')],...
            [paths.analyses num2str(IDs(id)) '/' strcat('con_',sprintf('%04d',ctr),'_mni.nii')])
    end
end

%% fMRI 2nd-level analysis pipeline
%% set environmental variables

clear;
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/alex/Documents/angela/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.physio  = [paths.parent 'E_data/EA_raw/EAE_physio/physio/3Tpilot/'];
paths.funx    = ['/Users/alex/Dropbox/literatures_IKND/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent];
paths.analyses= [paths.parent 'model_TaskMemory_immdel_rarefreq_stim_epi2t1/'];
paths.behav   = '/Users/alex/Documents/angela/behav/';
paths.group    = [paths.parent 'model_TaskMemory_immdel_rarefreq_stim_epi2t1/2nd/'];
paths.group = [paths.parent 'model_TaskMemory_immdel_rarefreq_stim_epi2t1/2nd/'];

% IDs
IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130 2131 2132 2233 2234 2235 2236]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

TR = 3.6;

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window

% 1
list_stim_rew_Remembered_Forgotten= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_stim_rew_Remembered_Forgotten{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0001_mni.nii,1'];
end

% 2
list_stim_rew_Forgotten_Remembered= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_stim_rew_Forgotten_Remembered{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0002_mni.nii,1'];
end

% 3
list_stim_neu_Remembered_Forgotten = [];
for id = 1:length(IDs)
    list_stim_neu_Remembered_Forgotten{id,1} = [paths.analyses num2str(IDs(id)) '/con_0003_mni.nii,1'];
end

% 4
list_stim_neu_Forgotten_Remembered = [];
for id = 1:length(IDs)
    list_stim_neu_Forgotten_Remembered{id,1} = [paths.analyses num2str(IDs(id)) '/con_0004_mni.nii,1'];
end

% 5
list_stim_rare_Remembered_Forgotten= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_stim_rare_Remembered_Forgotten{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0005_mni.nii,1'];
end

% 6
list_stim_rare_Forgotten_Remembered= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_stim_rare_Forgotten_Remembered{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0006_mni.nii,1'];
end

% 7
list_stim_freq_Remembered_Forgotten = [];
for id = 1:length(IDs)
    list_stim_freq_Remembered_Forgotten{id,1} = [paths.analyses num2str(IDs(id)) '/con_0007_mni.nii,1'];
end

% 8
list_stim_freq_Forgotten_Remembered = [];
for id = 1:length(IDs)
    list_stim_freq_Forgotten_Remembered{id,1} = [paths.analyses num2str(IDs(id)) '/con_0008_mni.nii,1'];
end

% 9
list_stim_rem_rare_freq = [];
for id = 1:length(IDs)
    list_stim_rem_rare_freq{id,1} = [paths.analyses num2str(IDs(id)) '/con_0009_mni.nii,1'];
end

% 10
list_stim_for_rare_freq = [];
for id = 1:length(IDs)
    list_stim_for_rare_freq{id,1} = [paths.analyses num2str(IDs(id)) '/con_0010_mni.nii,1'];
end

%% compute models

% ------- compute: stim, Rewards: Remembered > Forgotten ------- %

cd(paths.group); mkdir('stim_rew_Remembered_Forgotten'); cd stim_rew_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_rew_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Rewards: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Rewards: Forgotten > Remembered ------- %

cd(paths.group); mkdir('stim_rew_Forgotten_Remembered'); cd stim_rew_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_rew_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Rewards: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Neutrals: Remembered > Forgotten ------- %

cd(paths.group); mkdir('stim_neu_Remembered_Forgotten'); cd stim_neu_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_neu_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Neutrals: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Neutrals: Forgotten > Remembered ------- %

cd(paths.group); mkdir('stim_neu_Forgotten_Remembered'); cd stim_neu_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_neu_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Neutrals: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Rare: Remembered > Forgotten ------- %

cd(paths.group); mkdir('stim_rare_Remembered_Forgotten'); cd stim_rare_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_rare_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Rare: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Rare: Forgotten > Remembered ------- %

cd(paths.group); mkdir('stim_rare_Forgotten_Remembered'); cd stim_rare_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_rare_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Rare: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Frequent: Remembered > Forgotten ------- %

cd(paths.group); mkdir('stim_freq_Remembered_Forgotten'); cd stim_freq_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_freq_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Frequent: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: stim, Frequent: Forgotten > Remembered ------- %

cd(paths.group); mkdir('stim_freq_Forgotten_Remembered'); cd stim_freq_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_freq_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Frequent: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %




% ------- compute: stim, Remembered: Rare > Freq ------- %

cd(paths.group); mkdir('stim_remembered_rare_freq'); cd stim_remembered_rare_freq; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_rem_rare_freq;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Remembered: Rare > Freq';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %




% ------- compute: stim, Forgotten: Rare > Freq ------- %

cd(paths.group); mkdir('stim_forgotten_rare_freq'); cd stim_forgotten_rare_freq; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_stim_for_rare_freq;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'stim, Forgotten: Rare > Freq';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



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



%% feedback part 

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
paths.analyses= [paths.parent 'model_TaskMemory_immdel_rarefreq_fb_epi2t1/'];
paths.behav   = '/Users/alex/Documents/angela/behav/';
load('/Users/alex/Dropbox/literatures_IKND/length_scan_allsubj.mat')
length_scan1=length_scan1(:,2)

% IDs
IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130 2131 2132 2233 2234 2235 2236]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

% IDs  = [2202 2204 2205 2206 2208 2109 2110 2112 2113 2114 2115 2116 2218 2219 2220 2221 2222 2223 2224 2125 2126]; %% 2203 2207 2217 delayed memory test lost
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
% d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

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


%% start, fb

% spm fmri  % open progress window

for id = 28:length(IDs) % id 15 --> no forgotten items for rare reward : figure out something for this id
    %% build model
    
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
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'fb_RewardRare_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset= onsets{id,1}.fb.rewardRare_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'fb_RewardFrequent_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.fb.rewardFreq_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'fb_NeutralRare_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= onsets{id,1}.fb.neutralRare_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'fb_NeutralFrequent_Rem';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.fb.neutralFreq_rem;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;

    % forgotten
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'fb_RewardRare_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= onsets{id,1}.fb.rewardRare_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'fb_RewardFrequent_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= onsets{id,1}.fb.rewardFreq_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name = 'fb_NeutralRare_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset= onsets{id,1}.fb.neutralRare_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name = 'fb_NeutralFrequent_Forg';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset= onsets{id,1}.fb.neutralFreq_for;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth = 1;

    
    % else
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name = 'FixX';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset= onsets{id,1}.fixX;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).orth = 1;
   
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).name = 'Response';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).onset= onsets{id,1}.resp;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).orth = 1;
    
    %
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr([paths.preproc num2str(IDs(id)) '_multiple_regressors.mat']); % movement parameters
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    spm_jobman('run', matlabbatch) % run batch
    
    
    %% efbate
    
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([dir_model 'SPM.mat']);
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch) % run batch
    
    %% statistical inference
    
    firstlvl_dir        = [dir_model 'SPM.mat'];
    
    % fbs
    fb_rew_Remembered_Forgotten = [1 1 0 0 -1 -1];
    fb_rew_Forgotten_Remembered = [-1 -1 0 0 1 1];
    fb_neu_Remembered_Forgotten = [0 0 1 1 0 0 -1 -1];
    fb_neu_Forgotten_Remembered = [0 0 -1 -1 0 0 1 1];
    
    fb_rare_Remembered_Forgotten = [1 0 1 0 -1 0 -1];
    fb_rare_Forgotten_Remembered = [-1 0 -1 0 1 0 1];
    fb_freq_Remembered_Forgotten = [0 1 0 1 0 -1 0 -1];
    fb_freq_Forgotten_Remembered = [0 -1 0 -1 0 1 0 1];
    
    
%     fb_rem_Reward_Neutral  = [1 1 -1 -1];
%     fb_rem_Rare_Frequent   = [1 -1 1 -1];
%     fb_rem_Neutral_Reward  = [-1 -1 1 1];
%     fb_rem_Frequent_Rare   = [-1 1 -1 1];
%     
%     fb_for_Reward_Neutral  = [zeros(1,4) 1 1 -1 -1];
%     fb_for_Rare_Frequent   = [zeros(1,4) 1 -1 1 -1];
%     fb_for_Neutral_Reward  = [zeros(1,4) -1 -1 1 1];
%     fb_for_Frequent_Rare   = [zeros(1,4) -1 1 -1 1];
%     
%     fb_rem_RewardRare_RewardFreq     = [1 -1];
%     fb_rem_NeutralRare_NeutralFreq   = [0 0 1 -1];
%     fb_rem_RewardFreq_RewardRare     = [-1 1];
%     fb_rem_NeutralFreq_NeutralRare   = [0 0 -1 1];
%     
%     fb_for_RewardRare_RewardFreq     = [zeros(1,4) 1 -1];
%     fb_for_NeutralRare_NeutralFreq   = [zeros(1,4) 0 0 1 -1];
%     fb_for_RewardFreq_RewardRare     = [zeros(1,4) -1 1];
%     fb_for_NeutralFreq_NeutralRare   = [zeros(1,4) 0 0 -1 1];

    
    % batch setup
    
    if IDs(id) == 2220
        
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'fb, Neutrals: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = fb_neu_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'fb, Neutrals: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = fb_neu_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'nullcontrast';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'fb, Frequent: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = fb_freq_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'fb, Frequent: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = fb_freq_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
    
    spm_jobman('run', matlabbatch) % run batch
    
    else
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb, Rewards: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = fb_rew_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'fb, Rewards: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = fb_rew_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'fb, Neutrals: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = fb_neu_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'fb, Neutrals: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = fb_neu_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';    
    
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'fb, Rare: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = fb_rare_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'fb, Rare: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = fb_rare_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'fb, Frequent: Remembered > Forgotten';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = fb_freq_Remembered_Forgotten;
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
   
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'fb, Frequent: Forgotten > Remembered';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = fb_freq_Forgotten_Remembered;
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
    
    spm_jobman('run', matlabbatch) % run batch
    
    
    end
    
    
end

fprintf('\ndone\n')
%% set environmental variables

clear;
warning('off','all');

% paths

% paths
paths = [];
paths.parent  = '/Users/alex/Documents/angela/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.physio  = [paths.parent 'E_data/EA_raw/EAE_physio/physio/3Tpilot/'];
paths.funx    = ['/Users/alex/Dropbox/literatures_IKND/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent];
paths.analyses= [paths.parent 'model_TaskMemory_immdel_rarefreq_fb_epi2t1/'];
paths.behav   = '/Users/alex/Documents/angela/behav/';
paths.group    = [paths.parent 'model_TaskMemory_immdel_rarefreq_fb_epi2t1/2nd/'];
paths.group = [paths.parent 'model_TaskMemory_immdel_rarefreq_fb_epi2t1/2nd/'];

% IDs
IDs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2129 2130 2131 2132 2233 2234 2235 2236]; %% 2203 2207 2217 delayed memory test lost
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];


TR = 3.6;

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window

% 1
list_fb_rew_Remembered_Forgotten= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_fb_rew_Remembered_Forgotten{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0001_mni.nii,1'];
end

% 2
list_fb_rew_Forgotten_Remembered= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_fb_rew_Forgotten_Remembered{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0002_mni.nii,1'];
end

% 3
list_fb_neu_Remembered_Forgotten = [];
for id = 1:length(IDs)
    list_fb_neu_Remembered_Forgotten{id,1} = [paths.analyses num2str(IDs(id)) '/con_0003_mni.nii,1'];
end

% 4
list_fb_neu_Forgotten_Remembered = [];
for id = 1:length(IDs)
    list_fb_neu_Forgotten_Remembered{id,1} = [paths.analyses num2str(IDs(id)) '/con_0004_mni.nii,1'];
end

% 5
list_fb_rare_Remembered_Forgotten= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_fb_rare_Remembered_Forgotten{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0005_mni.nii,1'];
end

% 6
list_fb_rare_Forgotten_Remembered= []; cnt=0;
for id = [1:14 16:length(IDs)]
    cnt=cnt+1;
    list_fb_rare_Forgotten_Remembered{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0006_mni.nii,1'];
end

% 7
list_fb_freq_Remembered_Forgotten = [];
for id = 1:length(IDs)
    list_fb_freq_Remembered_Forgotten{id,1} = [paths.analyses num2str(IDs(id)) '/con_0007_mni.nii,1'];
end

% 8
list_fb_freq_Forgotten_Remembered = [];
for id = 1:length(IDs)
    list_fb_freq_Forgotten_Remembered{id,1} = [paths.analyses num2str(IDs(id)) '/con_0008_mni.nii,1'];
end

%% compute models

% ------- compute: fb, Rewards: Remembered > Forgotten ------- %

cd(paths.group); mkdir('fb_rew_Remembered_Forgotten'); cd fb_rew_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_rew_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Rewards: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Rewards: Forgotten > Remembered ------- %

cd(paths.group); mkdir('fb_rew_Forgotten_Remembered'); cd fb_rew_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_rew_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Rewards: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Neutrals: Remembered > Forgotten ------- %

cd(paths.group); mkdir('fb_neu_Remembered_Forgotten'); cd fb_neu_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_neu_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Neutrals: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Neutrals: Forgotten > Remembered ------- %

cd(paths.group); mkdir('fb_neu_Forgotten_Remembered'); cd fb_neu_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_neu_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Neutrals: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Rare: Remembered > Forgotten ------- %

cd(paths.group); mkdir('fb_rare_Remembered_Forgotten'); cd fb_rare_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_rare_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Rare: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Rare: Forgotten > Remembered ------- %

cd(paths.group); mkdir('fb_rare_Forgotten_Remembered'); cd fb_rare_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_rare_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Rare: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Frequent: Remembered > Forgotten ------- %

cd(paths.group); mkdir('fb_freq_Remembered_Forgotten'); cd fb_freq_Remembered_Forgotten; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_freq_Remembered_Forgotten;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Frequent: Remembered > Forgotten';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %



% ------- compute: fb, Frequent: Forgotten > Remembered ------- %

cd(paths.group); mkdir('fb_freq_Forgotten_Remembered'); cd fb_freq_Forgotten_Remembered; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_fb_freq_Forgotten_Remembered;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'fb, Frequent: Forgotten > Remembered';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% -------------------------------------------------------------- %
