%% PET 1st-level analysis pipeline

%% work log

%   12-12-2019    created the script

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/alex/Dropbox/paperwriting/MRPET/data/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = ['/Volumes/korokdorf/MRPET/img/'];
paths.analyses= [paths.parent 'fMRI/1stLevel/Eachsessions_memory/'];
paths.behav   = [paths.parent 'behav/'];
% paths.multiprm= [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/physio/'];
% paths.physio  = [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/physio/'];


% add toolboxes and functions
% addpath(paths.spm)
% addpath(paths.funx)

% IDs
% IDs  = [4001 4002 4004 4005 4008 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4028 4030 4031 4032 4033];
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; 
% IDs_half = [4003 4006 4007 4009 4027 4029];
% days_half= [1 0; 0 2; 1 0; 0 2; 0 2; 1 0];
% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 ...
    4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; ...
    1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; ...
    1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; ...
    1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/7.4.1/bin:/Applications/freesurfer/7.4.1/fsfast/bin:/Users/alex/fsl/bin:/Users/alex/fsl/share/fsl/bin:/Applications/freesurfer/7.4.1/mni/bin:/Applications/freesurfer/7.4.1:/usr/local/bin/python3/:/Library/Frameworks/Python.framework/Versions/3.12/bin:/Users/alex/fsl/share/fsl/bin:/Users/alex/fsl/share/fsl/bin:/opt/homebrew/bin:/opt/homebrew/sbin:/Users/alex/ants/bin/:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin:/opt/X11/bin']);
setenv('ANTSPATH','/Users/alex/ants/bin/')

% FSL Setup
setenv( 'FSLDIR', '/Users/alex/fsl' );
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

% scan specifics
TR = 3.6;
slices=51;
microtime=0;

%% load experimental details

expdat = []; onsets_temp=[]; onsets=[]; length_scan1=[];
for id = 1:length(IDs)

    % set up workspace
    mkdir([paths.analyses num2str(IDs(id))])

    for d = 1:2
        if days(id,d) == 0
            expdat{id,d} = {NaN};
            contingency{id,d} = [];
        else

            expdat{id,d}        = load([paths.behav num2str(IDs(id)) '_' num2str(days(id,d)) '.mat']);

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

            % clean null trials
            tmponsets = eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw']);
            tmptrigs  = eval(['(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st)']);
            tmptrials = cell2mat(eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,2)'])); % reward 1 neutral 2
            ind_null  = isnan(tmptrials); tmptrials(ind_null)=[];
            ind_rew   = tmptrials==contingency{id,d}(1); 
            ind_neu   = tmptrials==contingency{id,d}(2);

            % subsequent memory analyses
            clear stim_task stim_del stim_imm stim_del_rem stim_imm_rem accuracies_imm accuracies_del stim_rewards stim_neutrals
            stim_task = eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,1)']); stim_task(ind_null)=[];

            stim_imm  = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl(:,[1 4])']); % 1st column filename / second column old(1)/new(2)
            accuracies_imm= eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.accu']);
            stim_imm_rem  = stim_imm(cell2mat(stim_imm(:,2))==1 & accuracies_imm==1,1); % old items and remembered
            stim_imm_for  = stim_imm(cell2mat(stim_imm(:,2))==1 & accuracies_imm==0,1); % old items and forgotten
            stim_imm_unresp = stim_imm(cell2mat(stim_imm(:,2))==1 & isnan(accuracies_imm),1); % old items and unresponded

            if eval(['d' num2str(d) 'm(id,2)==0'])
                disp('no task data')
            else
                stim_del  = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.trl(:,[1 4])']);
                accuracies_del= eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.accu']);
                stim_del_rem  = stim_del(cell2mat(stim_del(:,2))==1 & accuracies_del==1,1); % old items and remembered
                stim_del_for  = stim_del(cell2mat(stim_del(:,2))==1 & accuracies_del==0,1); % old items and forgotten
                stim_del_unresp = stim_del(cell2mat(stim_del(:,2))==1 & isnan(accuracies_del),1); % old items and unresponded
            end

            clear tmpi tmpd
            tmpi = zeros(size(stim_task)); tmpd = zeros(size(stim_task));

            clear ind_remembered_imm isRemembered
            [isRemembered, pos] = ismember(stim_task, stim_imm_rem);
            ind_remembered_imm = find(isRemembered);
            indx_remembered_imm=tmpi; indx_remembered_imm(ind_remembered_imm)=1;

            clear ind_forgotten_imm isForgotten
            [isForgotten, pos] = ismember(stim_task, stim_imm_for);
            ind_forgotten_imm = find(isForgotten);
            indx_forgotten_imm=tmpi; indx_forgotten_imm(ind_forgotten_imm)=1;

            clear ind_unresponded_imm isUnresponded
            [isUnresponded, pos] = ismember(stim_task, stim_imm_unresp);
            ind_unresponded_imm = find(isUnresponded);
            indx_unresponded_imm=tmpi; indx_unresponded_imm(ind_unresponded_imm)=1;

            % clear ind_remembered_imm
            % for q = 1:length(stim_task)
            %     ind_remembered_imm(q,1) = {find(strcmp(stim_task{q,1},stim_imm_rem))};
            % end
            % tmp = cellfun(@isempty,ind_remembered_imm); indx_remembered_imm = tmp==0; clear tmp
            % indx_remembered_imm(ind_null)=[];

            % clear ind_forgotten_imm
            % for q = 1:length(stim_task)
            %     ind_forgotten_imm(q,1) = {find(strcmp(stim_task{q,1},stim_imm_for))};
            % end
            % tmp = cellfun(@isempty,ind_forgotten_imm); indx_forgotten_imm = tmp==0; clear tmp
            % indx_forgotten_imm(ind_null)=[];

            % clear ind_unresponded_imm
            % for q = 1:length(stim_task)
            %     ind_unresponded_imm(q,1) = {find(strcmp(stim_task{q,1},stim_imm_unresp))};
            % end
            % tmp = cellfun(@isempty,ind_unresponded_imm); indx_unresponded_imm = tmp==0; clear tmp
            % indx_unresponded_imm(ind_null)=[];

            if eval(['d' num2str(d) 'm(id,2)==0'])
                disp('no task data')
            else
                
                clear ind_remembered_del isRemembered
                [isRemembered, pos] = ismember(stim_task, stim_del_rem);
                ind_remembered_del = find(isRemembered); 
                indx_remembered_del=tmpd; indx_remembered_del(ind_remembered_del)=1;

                clear ind_forgotten_del isForgotten
                [isForgotten, pos] = ismember(stim_task, stim_del_for);
                ind_forgotten_del = find(isForgotten);
                indx_forgotten_del=tmpd; indx_forgotten_del(ind_forgotten_del)=1;

                clear ind_unresponded_del isUnresponded
                [isUnresponded, pos] = ismember(stim_task, stim_del_unresp);
                ind_unresponded_del = find(isUnresponded);
                indx_unresponded_del=tmpd; indx_unresponded_del(ind_unresponded_del)=1;

                % clear ind_remembered_del
                % for q = 1:length(stim_task)
                %     ind_remembered_del(q,1) = {find(strcmp(stim_task{q,1},stim_del_rem))};
                % end
                % tmp = cellfun(@isempty,ind_remembered_del); indx_remembered_del = tmp==0; clear tmp
                % indx_remembered_del(ind_null)=[];
                % 
                % clear ind_forgotten_del
                % for q = 1:length(stim_task)
                %     ind_forgotten_del(q,1) = {find(strcmp(stim_task{q,1},stim_del_for))};
                % end
                % tmp = cellfun(@isempty,ind_forgotten_del); indx_forgotten_del = tmp==0; clear tmp
                % indx_forgotten_del(ind_null)=[];

                % clear ind_unresponded_del
                % for q = 1:length(stim_task)
                %     ind_unresponded_del(q,1) = {find(strcmp(stim_task{q,1},stim_del_unresp))};
                % end
                % tmp = cellfun(@isempty,ind_unresponded_del); indx_unresponded_del = tmp==0; clear tmp
                % end
                % indx_unresponded_del(ind_null)=[];
            end

            % sort onsets
            onsets_temp{id,d}.stim.reward  = tmponsets.stim(ind_rew)-tmptrigs;
            onsets_temp{id,d}.stim.neutral = tmponsets.stim(ind_neu)-tmptrigs;
            onsets_temp{id,d}.null         = tmponsets.null-tmptrigs;
            onsets_temp{id,d}.resp         = tmponsets.resp-tmptrigs;
            onsets_temp{id,d}.fb.reward    = tmponsets.cue(ind_rew)-tmptrigs;
            onsets_temp{id,d}.fb.neutral   = tmponsets.cue(ind_neu)-tmptrigs;

            if eval(['d' num2str(d) 'm(id,2)==0'])
                onsets_temp{id,d}.stim.remembered  = tmponsets.stim(indx_remembered_imm==1)-tmptrigs;
                onsets_temp{id,d}.fb.remembered    = tmponsets.cue(indx_remembered_imm==1)-tmptrigs;

                onsets_temp{id,d}.stim.forgotten  = tmponsets.stim(indx_forgotten_imm==1)-tmptrigs;
                onsets_temp{id,d}.fb.forgotten    = tmponsets.cue(indx_forgotten_imm==1)-tmptrigs;

                onsets_temp{id,d}.stim.reward_rem  = tmponsets.stim(ind_rew==1 & indx_remembered_imm==1)-tmptrigs;
                onsets_temp{id,d}.stim.neutral_rem = tmponsets.stim(ind_neu==1 & indx_remembered_imm==1)-tmptrigs;
                onsets_temp{id,d}.fb.reward_rem    = tmponsets.cue(ind_rew==1 & indx_remembered_imm==1)-tmptrigs;
                onsets_temp{id,d}.fb.neutral_rem   = tmponsets.cue(ind_neu==1 & indx_remembered_imm==1)-tmptrigs;

                onsets_temp{id,d}.stim.reward_for  = tmponsets.stim(ind_rew==1 & indx_forgotten_imm==1)-tmptrigs;
                onsets_temp{id,d}.stim.neutral_for = tmponsets.stim(ind_neu==1 & indx_forgotten_imm==1)-tmptrigs;
                onsets_temp{id,d}.fb.reward_for    = tmponsets.cue(ind_rew==1 & indx_forgotten_imm==1)-tmptrigs;
                onsets_temp{id,d}.fb.neutral_for   = tmponsets.cue(ind_neu==1 & indx_forgotten_imm==1)-tmptrigs;


            else
                onsets_temp{id,d}.stim.remembered  = tmponsets.stim(indx_remembered_imm==1 | indx_remembered_del==1)-tmptrigs;
                onsets_temp{id,d}.fb.remembered    = tmponsets.cue(indx_remembered_imm==1 | indx_remembered_del==1)-tmptrigs;

                onsets_temp{id,d}.stim.forgotten  = tmponsets.stim(indx_forgotten_imm==1 | indx_forgotten_del==1)-tmptrigs;
                onsets_temp{id,d}.fb.forgotten    = tmponsets.cue(indx_forgotten_imm==1 | indx_forgotten_del==1)-tmptrigs;

                onsets_temp{id,d}.stim.reward_rem  = tmponsets.stim(ind_rew==1 & (indx_remembered_imm==1 | indx_remembered_del==1))-tmptrigs;
                onsets_temp{id,d}.stim.neutral_rem = tmponsets.stim(ind_neu==1 & (indx_remembered_imm==1 | indx_remembered_del==1))-tmptrigs;
                onsets_temp{id,d}.fb.reward_rem    = tmponsets.cue(ind_rew==1 & (indx_remembered_imm==1 | indx_remembered_del==1))-tmptrigs;
                onsets_temp{id,d}.fb.neutral_rem   = tmponsets.cue(ind_neu==1 & (indx_remembered_imm==1 | indx_remembered_del==1))-tmptrigs;

                onsets_temp{id,d}.stim.reward_for  = tmponsets.stim(ind_rew==1 & (indx_forgotten_imm==1 | indx_forgotten_del==1))-tmptrigs;
                onsets_temp{id,d}.stim.neutral_for = tmponsets.stim(ind_neu==1 & (indx_forgotten_imm==1 | indx_forgotten_del==1))-tmptrigs;
                onsets_temp{id,d}.fb.reward_for    = tmponsets.cue(ind_rew==1 & (indx_forgotten_imm==1 | indx_forgotten_del==1))-tmptrigs;
                onsets_temp{id,d}.fb.neutral_for   = tmponsets.cue(ind_neu==1 & (indx_forgotten_imm==1 | indx_forgotten_del==1))-tmptrigs;

            end
        end
    end

    % sort onsets

    if sum(days(id,:))==3

        clear numvolSess1
        numvolSess1  = length(spm_vol([paths.preproc num2str(IDs(id)) '_1/' num2str(IDs(id)) '_MRI_4D_MT1.nii']));
        length_scan1(id,1) = numvolSess1;

        onsets{id,1}.stim.highDA_reward  = onsets_temp{id,1}.stim.reward;
        onsets{id,1}.stim.lowDA_reward   = onsets_temp{id,2}.stim.reward;
        onsets{id,1}.stim.highDA_neutral = onsets_temp{id,1}.stim.neutral;
        onsets{id,1}.stim.lowDA_neutral  = onsets_temp{id,2}.stim.neutral;
        onsets{id,1}.stim.all_reward     = [onsets_temp{id,1}.stim.reward; onsets_temp{id,2}.stim.reward];
        onsets{id,1}.stim.all_neutral    = [onsets_temp{id,1}.stim.neutral; onsets_temp{id,2}.stim.neutral];

        onsets{id,1}.fb.highDA_reward    = onsets_temp{id,1}.fb.reward;
        onsets{id,1}.fb.lowDA_reward     = onsets_temp{id,2}.fb.reward;
        onsets{id,1}.fb.highDA_neutral   = onsets_temp{id,1}.fb.neutral;
        onsets{id,1}.fb.lowDA_neutral    = onsets_temp{id,2}.fb.neutral;
        onsets{id,1}.fb.all_reward       = [onsets_temp{id,1}.fb.reward; onsets_temp{id,2}.fb.reward];
        onsets{id,1}.fb.all_neutral      = [onsets_temp{id,1}.fb.neutral; onsets_temp{id,2}.fb.neutral];

        onsets{id,1}.highDA_null         = onsets_temp{id,1}.null;
        onsets{id,1}.highDA_resp         = onsets_temp{id,1}.resp;
        onsets{id,1}.lowDA_null          = onsets_temp{id,2}.null;
        onsets{id,1}.lowDA_resp          = onsets_temp{id,2}.resp;
        onsets{id,1}.all_null            = [onsets_temp{id,1}.null; onsets_temp{id,2}.null];
        onsets{id,1}.all_resp            = [onsets_temp{id,1}.resp; onsets_temp{id,2}.resp];

        % memory
        onsets{id,1}.stim.highDA_rem  = onsets_temp{id,1}.stim.remembered;
        onsets{id,1}.stim.lowDA_rem   = onsets_temp{id,2}.stim.remembered;
        onsets{id,1}.fb.highDA_rem    = onsets_temp{id,1}.fb.remembered;
        onsets{id,1}.fb.lowDA_rem     = onsets_temp{id,2}.fb.remembered;

        onsets{id,1}.stim.highDA_reward_rem  = onsets_temp{id,1}.stim.reward_rem;
        onsets{id,1}.stim.lowDA_reward_rem   = onsets_temp{id,2}.stim.reward_rem;
        onsets{id,1}.stim.highDA_neutral_rem = onsets_temp{id,1}.stim.neutral_rem;
        onsets{id,1}.stim.lowDA_neutral_rem  = onsets_temp{id,2}.stim.neutral_rem;
        onsets{id,1}.stim.all_reward_rem     = [onsets_temp{id,1}.stim.reward_rem; onsets_temp{id,2}.stim.reward_rem];
        onsets{id,1}.stim.all_neutral_rem    = [onsets_temp{id,1}.stim.neutral_rem; onsets_temp{id,2}.stim.neutral_rem];

        onsets{id,1}.fb.highDA_reward_rem    = onsets_temp{id,1}.fb.reward_rem;
        onsets{id,1}.fb.lowDA_reward_rem     = onsets_temp{id,2}.fb.reward_rem;
        onsets{id,1}.fb.highDA_neutral_rem   = onsets_temp{id,1}.fb.neutral_rem;
        onsets{id,1}.fb.lowDA_neutral_rem    = onsets_temp{id,2}.fb.neutral_rem;
        onsets{id,1}.fb.all_reward_rem       = [onsets_temp{id,1}.fb.reward_rem; onsets_temp{id,2}.fb.reward_rem];
        onsets{id,1}.fb.all_neutral_rem      = [onsets_temp{id,1}.fb.neutral_rem; onsets_temp{id,2}.fb.neutral_rem];

        onsets{id,1}.stim.highDA_for  = onsets_temp{id,1}.stim.forgotten;
        onsets{id,1}.stim.lowDA_for   = onsets_temp{id,2}.stim.forgotten;
        onsets{id,1}.fb.highDA_for    = onsets_temp{id,1}.fb.forgotten;
        onsets{id,1}.fb.lowDA_for     = onsets_temp{id,2}.fb.forgotten;

        onsets{id,1}.stim.highDA_reward_for  = onsets_temp{id,1}.stim.reward_for;
        onsets{id,1}.stim.lowDA_reward_for   = onsets_temp{id,2}.stim.reward_for;
        onsets{id,1}.stim.highDA_neutral_for = onsets_temp{id,1}.stim.neutral_for;
        onsets{id,1}.stim.lowDA_neutral_for  = onsets_temp{id,2}.stim.neutral_for;
        onsets{id,1}.stim.all_reward_for     = [onsets_temp{id,1}.stim.reward_for; onsets_temp{id,2}.stim.reward_for];
        onsets{id,1}.stim.all_neutral_for    = [onsets_temp{id,1}.stim.neutral_for; onsets_temp{id,2}.stim.neutral_for];

        onsets{id,1}.fb.highDA_reward_for    = onsets_temp{id,1}.fb.reward_for;
        onsets{id,1}.fb.lowDA_reward_for     = onsets_temp{id,2}.fb.reward_for;
        onsets{id,1}.fb.highDA_neutral_for   = onsets_temp{id,1}.fb.neutral_for;
        onsets{id,1}.fb.lowDA_neutral_for    = onsets_temp{id,2}.fb.neutral_for;
        onsets{id,1}.fb.all_reward_for       = [onsets_temp{id,1}.fb.reward_for; onsets_temp{id,2}.fb.reward_for];
        onsets{id,1}.fb.all_neutral_for      = [onsets_temp{id,1}.fb.neutral_for; onsets_temp{id,2}.fb.neutral_for];


    elseif sum(days(id,:))==1
        clear numvolSess1
        numvolSess1  = length(spm_vol([paths.preproc num2str(IDs(id)) '_1/' num2str(IDs(id)) '_MRI_4D_MT1.nii']));
        length_scan1(id,1) = numvolSess1;

        onsets{id,1}.stim.highDA_reward  = onsets_temp{id,1}.stim.reward;
        onsets{id,1}.stim.lowDA_reward   = numvolSess1*TR-30; % null regressor
        onsets{id,1}.stim.highDA_neutral = onsets_temp{id,1}.stim.neutral;
        onsets{id,1}.stim.lowDA_neutral  = numvolSess1*TR-15; % null regressor
        onsets{id,1}.stim.all_reward     = [onsets_temp{id,1}.stim.reward];
        onsets{id,1}.stim.all_neutral    = [onsets_temp{id,1}.stim.neutral];

        onsets{id,1}.fb.highDA_reward    = onsets_temp{id,1}.fb.reward;
        onsets{id,1}.fb.lowDA_reward     = numvolSess1*TR-30+4; % null regressor
        onsets{id,1}.fb.highDA_neutral   = onsets_temp{id,1}.fb.neutral;
        onsets{id,1}.fb.lowDA_neutral    = numvolSess1*TR-15+4; % null regressor
        onsets{id,1}.fb.all_reward       = [onsets_temp{id,1}.fb.reward];
        onsets{id,1}.fb.all_neutral      = [onsets_temp{id,1}.fb.neutral];

        onsets{id,1}.highDA_null         = onsets_temp{id,1}.null;
        onsets{id,1}.highDA_resp         = onsets_temp{id,1}.resp;
        onsets{id,1}.lowDA_null          = numvolSess1*TR-TR; % null regressor
        onsets{id,1}.lowDA_resp          = numvolSess1*-15+4+2.5; % null regressor
        onsets{id,1}.all_null            = [onsets_temp{id,1}.null];
        onsets{id,1}.all_resp            = [onsets_temp{id,1}.resp];

        % memory
        onsets{id,1}.stim.highDA_remembered  = onsets_temp{id,1}.stim.remembered;
        onsets{id,1}.fb.highDA_remembered    = onsets_temp{id,1}.fb.remembered;

        onsets{id,1}.stim.highDA_reward_rem  = onsets_temp{id,1}.stim.reward_rem;
        onsets{id,1}.stim.lowDA_reward_rem   = numvolSess1*TR-30; % null regressor
        onsets{id,1}.stim.highDA_neutral_rem = onsets_temp{id,1}.stim.neutral_rem;
        onsets{id,1}.stim.lowDA_neutral_rem  = numvolSess1*TR-15; % null regressor
        onsets{id,1}.stim.all_reward_rem     = [onsets_temp{id,1}.stim.reward_rem];
        onsets{id,1}.stim.all_neutral_rem    = [onsets_temp{id,1}.stim.neutral_rem];

        onsets{id,1}.fb.highDA_reward_rem    = onsets_temp{id,1}.fb.reward_rem;
        onsets{id,1}.fb.lowDA_reward_rem     = numvolSess1*TR-30+4; % null regressor
        onsets{id,1}.fb.highDA_neutral_rem   = onsets_temp{id,1}.fb.neutral_rem;
        onsets{id,1}.fb.lowDA_neutral_rem    = numvolSess1*TR-15+4; % null regressor
        onsets{id,1}.fb.all_reward_rem       = [onsets_temp{id,1}.fb.reward_rem];
        onsets{id,1}.fb.all_neutral_rem      = [onsets_temp{id,1}.fb.neutral_rem];

        onsets{id,1}.stim.highDA_forgotten  = onsets_temp{id,1}.stim.forgotten;
        onsets{id,1}.fb.highDA_forgotten    = onsets_temp{id,1}.fb.forgotten;

        onsets{id,1}.stim.highDA_reward_for  = onsets_temp{id,1}.stim.reward_for;
        onsets{id,1}.stim.lowDA_reward_for   = numvolSess1*TR-30; % null regressor
        onsets{id,1}.stim.highDA_neutral_for = onsets_temp{id,1}.stim.neutral_for;
        onsets{id,1}.stim.lowDA_neutral_for  = numvolSess1*TR-15; % null regressor
        onsets{id,1}.stim.all_reward_for     = [onsets_temp{id,1}.stim.reward_for];
        onsets{id,1}.stim.all_neutral_for    = [onsets_temp{id,1}.stim.neutral_for];

        onsets{id,1}.fb.highDA_reward_for    = onsets_temp{id,1}.fb.reward_for;
        onsets{id,1}.fb.lowDA_reward_for     = numvolSess1*TR-30+4; % null regressor
        onsets{id,1}.fb.highDA_neutral_for   = onsets_temp{id,1}.fb.neutral_for;
        onsets{id,1}.fb.lowDA_neutral_for    = numvolSess1*TR-15+4; % null regressor
        onsets{id,1}.fb.all_reward_for       = [onsets_temp{id,1}.fb.reward_for];
        onsets{id,1}.fb.all_neutral_for      = [onsets_temp{id,1}.fb.neutral_for];

    elseif sum(days(id,:))==2
        clear numvolSess1
        numvolSess1  = length(spm_vol([paths.preproc num2str(IDs(id)) '_2/' num2str(IDs(id)) '_MRI_4D_MT2.nii']));
        length_scan1(id,1) = numvolSess1;

        onsets{id,1}.stim.lowDA_reward  = onsets_temp{id,2}.stim.reward;
        onsets{id,1}.stim.highDA_reward   = numvolSess1*TR-30; % null regressor
        onsets{id,1}.stim.lowDA_neutral = onsets_temp{id,2}.stim.neutral;
        onsets{id,1}.stim.highDA_neutral  = numvolSess1*TR-15; % null regressor
        onsets{id,1}.stim.all_reward     = [onsets_temp{id,2}.stim.reward];
        onsets{id,1}.stim.all_neutral    = [onsets_temp{id,2}.stim.neutral];

        onsets{id,1}.fb.lowDA_reward    = onsets_temp{id,2}.fb.reward;
        onsets{id,1}.fb.highDA_reward     = numvolSess1*TR-30+4; % null regressor
        onsets{id,1}.fb.lowDA_neutral   = onsets_temp{id,2}.fb.neutral;
        onsets{id,1}.fb.highDA_neutral    = numvolSess1*TR-15+4; % null regressor
        onsets{id,1}.fb.all_reward       = [onsets_temp{id,2}.fb.reward];
        onsets{id,1}.fb.all_neutral      = [onsets_temp{id,2}.fb.neutral];

        onsets{id,1}.lowDA_null         = onsets_temp{id,2}.null;
        onsets{id,1}.lowDA_resp         = onsets_temp{id,2}.resp;
        onsets{id,1}.highDA_null          = numvolSess1*TR-TR; % null regressor
        onsets{id,1}.highDA_resp          = numvolSess1*TR-15+4+2.5; % null regressor
        onsets{id,1}.all_null            = [onsets_temp{id,2}.null];
        onsets{id,1}.all_resp            = [onsets_temp{id,2}.resp];

        % memory
        onsets{id,1}.stim.lowDA_remembered  = onsets_temp{id,2}.stim.remembered;
        onsets{id,1}.fb.lowDA_remembered    = onsets_temp{id,2}.fb.remembered;

        onsets{id,1}.stim.lowDA_reward_rem  = onsets_temp{id,2}.stim.reward_rem;
        onsets{id,1}.stim.highDA_reward_rem   = numvolSess1*TR-30; % null regressor
        onsets{id,1}.stim.lowDA_neutral_rem = onsets_temp{id,2}.stim.neutral_rem;
        onsets{id,1}.stim.highDA_neutral_rem  = numvolSess1*TR-15; % null regressor
        onsets{id,1}.stim.all_reward_rem     = [onsets_temp{id,2}.stim.reward_rem];
        onsets{id,1}.stim.all_neutral_rem    = [onsets_temp{id,2}.stim.neutral_rem];

        onsets{id,1}.fb.lowDA_reward_rem    = onsets_temp{id,2}.fb.reward_rem;
        onsets{id,1}.fb.highDA_reward_rem     = numvolSess1*TR-30+4; % null regressor
        onsets{id,1}.fb.lowDA_neutral_rem   = onsets_temp{id,2}.fb.neutral_rem;
        onsets{id,1}.fb.highDA_neutral_rem    = numvolSess1*TR-15+4; % null regressor
        onsets{id,1}.fb.all_reward_rem       = [onsets_temp{id,2}.fb.reward_rem];
        onsets{id,1}.fb.all_neutral_rem      = [onsets_temp{id,2}.fb.neutral_rem];

        onsets{id,1}.stim.lowDA_forgotten  = onsets_temp{id,2}.stim.forgotten;
        onsets{id,1}.fb.lowDA_forgotten    = onsets_temp{id,2}.fb.forgotten;

        onsets{id,1}.stim.lowDA_reward_for  = onsets_temp{id,2}.stim.reward_for;
        onsets{id,1}.stim.highDA_reward_for   = numvolSess1*TR-30; % null regressor
        onsets{id,1}.stim.lowDA_neutral_for = onsets_temp{id,2}.stim.neutral_for;
        onsets{id,1}.stim.highDA_neutral_for  = numvolSess1*TR-15; % null regressor
        onsets{id,1}.stim.all_reward_for     = [onsets_temp{id,2}.stim.reward_for];
        onsets{id,1}.stim.all_neutral_for    = [onsets_temp{id,2}.stim.neutral_for];

        onsets{id,1}.fb.lowDA_reward_for    = onsets_temp{id,2}.fb.reward_for;
        onsets{id,1}.fb.highDA_reward_for     = numvolSess1*TR-30+4; % null regressor
        onsets{id,1}.fb.lowDA_neutral_for   = onsets_temp{id,2}.fb.neutral_for;
        onsets{id,1}.fb.highDA_neutral_for    = numvolSess1*TR-15+4; % null regressor
        onsets{id,1}.fb.all_reward_for       = [onsets_temp{id,2}.fb.reward_for];
        onsets{id,1}.fb.all_neutral_for      = [onsets_temp{id,2}.fb.neutral_for];


    end

    % make EPI mask
    %     eval(['!bet2 ' paths.preproc num2str(IDs(id)) '/meanall_ua' num2str(IDs(id)) '.nii ' paths.preproc num2str(IDs(id)) '/EPI_brainonly -m'])
    %     gunzip([paths.preproc num2str(IDs(id)) '/EPI_brainonly_mask.nii.gz'])
    %     copyfile([paths.preproc num2str(IDs(id)) '/EPI_brainonly_mask.nii'],[paths.preproc num2str(IDs(id)) '/mask.nii'])
    %     delete([paths.preproc num2str(IDs(id)) '/EPI_brainonly_mask.nii.gz']); delete([paths.preproc num2str(IDs(id)) '/EPI_brainonly.nii.gz']);

end

%% make physio+realignment+session regressors (already done)

for id = [27 29]%1:length(IDs)
    
%     copyfile([paths.preproc num2str(IDs(id)) '/reg_all.mat'],...
%         ['/Users/alex/Dropbox/paperwriting/MRPET/scripts/kalina_firstlevel/data/multiparam/' num2str(IDs(id)) '_reg_all.mat'])
%     copyfile([paths.behav fname_beh{id,d}],['/Users/alex/Dropbox/paperwriting/MRPET/scripts/kalina_firstlevel/data/behav/' fname_beh{id,d}])
    
%     %% make physio + realignment parameters
%             
%     % read physio parameters
%     % day1
%     if IDs(id)==4008 || IDs(id)==4017 || IDs(id)==4018 
%         clear R
%         load(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_1/reg_all.mat'])
%         physiodata1  = R(:,1:18);
%     elseif IDs(id)==4012 || IDs(id)==4019 % no pulse data
%         physioFile1  = [paths.preproc num2str(IDs(id)) '_1/multiple_regressors.txt'];
%         physiodata1  = importdata(physioFile1);
%         physiodata1  = [ zeros(length(physiodata1),10) physiodata1];
%     else
%         physioFile1  = [paths.preproc num2str(IDs(id)) '_1/multiple_regressors.txt'];
%         physiodata1  = importdata(physioFile1);
%     end
%     % day2
%     if IDs(id)==4011 || IDs(id)==4017 || IDs(id)==4018 || IDs(id)==4019 
%         load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/4011_2/reg_all.mat')
%         physiodata2  = R(:,1:18);
%     elseif IDs(id)==4013 % no pulse data
%         physioFile2  = [paths.preproc num2str(IDs(id)) '_2/multiple_regressors.txt'];
%         physiodata2  = importdata(physioFile2);
%         physiodata2  = [ zeros(length(physiodata2),10) physiodata2];
%     else
%         physioFile2  = [paths.preproc num2str(IDs(id)) '_2/multiple_regressors.txt'];
%         physiodata2  = importdata(physioFile2);
%     end
%     % assemble
%     physiodata = [physiodata1; physiodata2];
% 
%     % read realignment parameters
%     realignFile1 = [paths.preproc num2str(IDs(id)) '/rp_all_ua' num2str(IDs(id)) '.txt'];
%     realigndata1 = importdata(realignFile1);
% 
%     % session marker
%     clear numvolSess1 numvolSess2 sessionMarker
%     numvolSess1  = spm_vol([paths.preproc num2str(IDs(id)) '_1/' num2str(IDs(id)) '_MRI_4D_MT1.nii']);
%     numvolSess2  = spm_vol([paths.preproc num2str(IDs(id)) '_2/' num2str(IDs(id)) '_MRI_4D_MT2.nii']);
%     sessionMarker= [ones(length(numvolSess1),1); zeros(length(numvolSess2),1)];
% 
%     % concatenate
%     R = [physiodata realigndata1 sessionMarker];
% 
%     % read already made realignment+physio regressor
% %     multiregFile = [paths.physio num2str(IDs(id)) '/reg_all.txt'];
% %     R = importdata(multiregFile);
%     save([paths.preproc num2str(IDs(id)) '/reg_all.mat'],'R');
% 
%     clear R realigndata1 realignFile1 physiodata1 physiodata2 physiodata physioFile1 physioFile2
% 
end


fprintf('\n preparation done \n')
%% start

spm fmri  % open progress window

onsets{13,1}.stim.highDA_neutral_for = 1;

for id = 1:length(IDs)

    for d = 1:2
        if days(id,d) == 0
        else

            fprintf('\n model estimation for: ID %d \n', IDs(id))

            % set up workspace
            cd(paths.analyses); mkdir([num2str(IDs(id)) '_' num2str(d)]);
            clear dir_model
            dir_model = [paths.analyses num2str(IDs(id)) '_' num2str(d) '/'];
            dir_multiparam = [paths.preproc num2str(IDs(id)) '_' num2str(d) '/reg_all.mat'];

            %% build model

            clear matlabbatch

            spm_jobman('initcfg');      % initiate job manager

            matlabbatch{1}.spm.stats.fmri_spec.dir               = cellstr(dir_model); % where will the model be saved?
            matlabbatch{1}.spm.stats.fmri_spec.timing.units      = 'secs'; % scans / secs
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT         = TR; % TR in seconds
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t     = slices; % volume size
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0    = microtime; % microtime onset

            matlabbatch{1}.spm.stats.fmri_spec.sess.scans        = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(d) '/' ...
                'sua' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii']);

            if days(id,d)==1
            % high reward condition
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Stim_RewRem';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset= onsets{id,1}.stim.highDA_reward_rem;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Stim_RewForg';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.stim.highDA_reward_for;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Stim_NeuRem';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= onsets{id,1}.stim.highDA_neutral_rem;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Stim_NeuForg';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.stim.highDA_neutral_for;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;

            elseif days(id,d)==2
            % low reward condition
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Stim_RewRem';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset= onsets{id,1}.stim.lowDA_reward_rem;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Stim_RewForg';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.stim.lowDA_reward_for;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Stim_NeuRem';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= onsets{id,1}.stim.lowDA_neutral_rem;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Stim_NeuForg';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.stim.lowDA_neutral_for;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;

            end

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'NullTrials';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= onsets{id,1}.all_null;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'ButtonPresses';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= onsets{id,1}.all_resp;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});

            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(dir_multiparam); % movement parameters
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {[paths.preproc num2str(IDs(id)) '/mask.nii']};
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

            % here you set up the contrast matrices

            stim_Remembered_vs_Forgotten = [1 -1 1 -1];
            stim_Forgotten_vs_Remembered = [-1 1 -1 1];

            stim_Rew_Remembered_Forgotten = [1 -1];
            stim_Rew_Forgotten_Remembered = [-1 1];
            stim_Neu_Remembered_Forgotten = [0 0 1 -1];
            stim_Neu_Forgotten_Remembered = [0 0 -1 1];

            stim_Rem_Reward_Neutral  = [1 0 -1];
            stim_Rem_Neutral_Reward  = [-1 0 1];

            stim_Forg_Reward_Neutral  = [0 1 0 -1];
            stim_Forg_Neutral_Reward  = [0 -1 0 1];

            Nulls  = [0 0 0 0 1];
            Button = [0 0 0 0 0 1];
            stim_vs_null  = [ones(1,4)./4 -1];

            % % ANOVA contrasts
            stim_rew_rem   = [1];
            stim_neu_rem   = [0 0 1];

            stim_rew_forg  = [0 1];
            stim_neu_forg  = [0 0 0 1];

            % batch setup
            clear matlabbatch
            spm_jobman('initcfg');      % initiate job manager

            matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);

            % if IDs(id)==4007 | IDs(id)==4009
            %     paramlength=24-4;
            % else
            %     paramlength=34-4;
            % end

            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_Remembered_vs_Forgotten';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = stim_Remembered_vs_Forgotten;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'stim_Forgotten_vs_Remembered';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = stim_Forgotten_vs_Remembered;
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'stim_Rew_Remembered_Forgotten';
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = stim_Rew_Remembered_Forgotten;
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'stim_Rew_Forgotten_Remembered';
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = stim_Rew_Forgotten_Remembered;
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'stim_Neu_Remembered_Forgotten';
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = stim_Neu_Remembered_Forgotten;
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'stim_Neu_Forgotten_Remembered';
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = stim_Neu_Forgotten_Remembered;
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'stim_Rem_Reward_Neutral';
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = stim_Rem_Reward_Neutral;
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'stim_Rem_Neutral_Reward';
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = stim_Rem_Neutral_Reward;
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'stim_Forg_Reward_Neutral';
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = stim_Forg_Reward_Neutral;
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'stim_Forg_Neutral_Reward';
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = stim_Forg_Neutral_Reward;
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'Nulls';
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = Nulls;
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'Button';
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = Button;
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'Scenes_vs_Null';
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = stim_vs_null;% 29
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'stim: reward remembered';
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = stim_rew_rem;
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{15}.tcon.name = 'stim: neutral remembered';
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.weights = stim_neu_rem;
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{16}.tcon.name = 'stim: reward forgotten';
            matlabbatch{1}.spm.stats.con.consess{16}.tcon.weights = stim_rew_forg;
            matlabbatch{1}.spm.stats.con.consess{16}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.consess{17}.tcon.name = 'stim: neutral forgotten';
            matlabbatch{1}.spm.stats.con.consess{17}.tcon.weights = stim_neu_forg;
            matlabbatch{1}.spm.stats.con.consess{17}.tcon.sessrep = 'none';

            matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append

            spm_jobman('run', matlabbatch) % run batch
        end

    end
end

fprintf('\ndone\n')

%% move contrasts for registration

path_source      =[paths.analyses];
path_destination ='/Volumes/korokdorf/MRPET/coreg_mri/eachsession/';

mrpetID=[]; c1=0;
for id=1:length(IDs)

            c1=c1+1;
            mrpetID{c1,1}=[num2str(IDs(id))];
                        
            for cnt=1:37
            copyfile([path_source num2str(IDs(id)) '/con_00' sprintf('%02i',cnt) '.nii'],...
                [path_destination num2str(IDs(id)) '/data/con_00' sprintf('%02i',cnt) '.nii'])
            end

            % for cnt=1:28
            % copyfile([path_destination num2str(IDs(id)) '/data/con_00' sprintf('%02i',cnt) '_mem_mni.nii'],...
            %     [path_source num2str(IDs(id)) '/con_00' sprintf('%02i',cnt) '_mni.nii'])
            % end

            % copyfile([path_source num2str(IDs(id)) '/mask.nii'],...
            %     [path_destination num2str(IDs(id)) '/data/mask.nii'])
            
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/T1WB.nii'])
            
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/meanEPI.nii'])
            
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_GRE3D' num2str(d) '_2.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/MTw.nii'])

end


%% fMRI 2nd-level analysis pipeline
%% set environmental variables

clear;
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
% paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.analyses= [paths.parent '1stLevel/Bothsessions_memory/'];
paths.analyses_half = [paths.parent '1stLevel/Bothsessions_memory/'];
paths.behav   = ['/Users/alex/Dropbox/paperwriting/MRPET/data/behav/'];
% paths.temp    = [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/MRI/'];
paths.save2nd = [paths.parent '2ndLevel/2nd_template_CompleteDatasetsOnly_memory/'];


% add toolboxes and functions
% addpath(paths.spm)
% addpath(paths.funx)

% IDs
% IDs  = [4001 4002 4004 4005 4008 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4028 4030 4031 4032 4033];
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; 

IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2]; 

IDs_both  = [4001 4002 4004 4005 4008 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4028 4030 4031 4032 4033];
% days      = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; 

IDs_high = [4001 4002 4003 4004 4005 4007 4008 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4030 4031 4032 4033];
% days     = [1 2; 1 2; 1 0; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2]; 

IDs_low  = [4001 4002 4004 4005 4006 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4028 4029 4030 4031 4032 4033];
% days     = [1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2]; 


IDs_half = [4003 4006 4007 4009 4027 4029];
days_half= [1 0; 0 2; 1 0; 0 2; 1 0; 1 0];


% load experimental details
expdat = [];
for c1 = 1:length(IDs)
    for d = 1:2
        if days(c1,d) == 0
            fname_beh{c1,d} = {NaN};
            expdat{c1,d} = {NaN};
        else
            fname_beh{c1,d}     = [ num2str(IDs(c1)) '_' num2str(days(c1,d)) '.mat' ];
            expdat{c1,d} = load([paths.behav fname_beh{c1,d}]);
        end
    end
end

TR = 3.6;

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window

list_stim_bothsess_Remembered_vs_Forgotten_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_bothsess_Remembered_vs_Forgotten_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0001_template.nii,1'];
end

list_stim_bothsess_Forgotten_vs_Remembered_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_bothsess_Forgotten_vs_Remembered_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0002_template.nii,1'];
end

list_stim_RewAll_Remembered_Forgotten_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_RewAll_Remembered_Forgotten_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0003_template.nii,1'];
end

list_stim_RewAll_Forgotten_Remembered_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_RewAll_Forgotten_Remembered_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0004_template.nii,1'];
end

list_stim_NeuAll_Remembered_Forgotten_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_NeuAll_Remembered_Forgotten_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0005_template.nii,1'];
end

list_stim_NeuAll_Forgotten_Remembered_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_NeuAll_Forgotten_Remembered_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0006_template.nii,1'];
end

list_stim_highSess_Rew_Remembered_Forgotten = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Rew_Remembered_Forgotten{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0007_template.nii,1'];
end

list_stim_highSess_Rew_Forgotten_Remembered = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Rew_Forgotten_Remembered{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0008_template.nii,1'];
end

list_stim_highSess_Neu_Remembered_Forgotten = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Neu_Remembered_Forgotten{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0009_template.nii,1'];
end

list_stim_highSess_Neu_Forgotten_Remembered = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Neu_Forgotten_Remembered{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0010_template.nii,1'];
end

list_stim_lowSess_Rew_Remembered_Forgotten = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Rew_Remembered_Forgotten{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0011_template.nii,1'];
end

list_stim_lowSess_Rew_Forgotten_Remembered = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Rew_Forgotten_Remembered{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0012_template.nii,1'];
end

list_stim_lowSess_Neu_Remembered_Forgotten = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Neu_Remembered_Forgotten{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0013_template.nii,1'];
end

list_stim_lowSess_Neu_Forgotten_Remembered = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Neu_Forgotten_Remembered{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0014_template.nii,1'];
end

list_stim_RemAll_Reward_Neutral_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_RemAll_Reward_Neutral_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0015_template.nii,1'];
end

list_stim_RemAll_Neutral_Reward_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_RemAll_Neutral_Reward_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0016_template.nii,1'];
end

list_stim_highSess_Rem_Reward_Neutral = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Rem_Reward_Neutral{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0017_template.nii,1'];
end

list_stim_highSess_Rem_Neutral_Reward=[];
for c1 = 1:length(IDs_high)
    list_stim_highSess_Rem_Neutral_Reward{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0018_template.nii,1'];
end

list_stim_lowSess_Rem_Reward_Neutral = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Rem_Reward_Neutral{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0019_template.nii,1'];
end

list_stim_lowSess_Rem_Neutral_Reward = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Rem_Neutral_Reward{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0020_template.nii,1'];
end

list_stim_ForgAll_Reward_Neutral_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_ForgAll_Reward_Neutral_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0021_template.nii,1'];
end

list_stim_ForgAll_Neutral_Reward_bothonly = [];
for c1 = 1:length(IDs_both)
    list_stim_ForgAll_Neutral_Reward_bothonly{c1,1} = [paths.analyses num2str(IDs_both(c1)) '/con_0022_template.nii,1'];
end

list_stim_highSess_Forg_Reward_Neutral = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Forg_Reward_Neutral{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0023_template.nii,1'];
end

list_stim_highSess_Forg_Neutral_Reward = []; 
for c1 = 1:length(IDs_high)
    list_stim_highSess_Forg_Neutral_Reward{c1,1} = [paths.analyses num2str(IDs_high(c1)) '/con_0024_template.nii,1'];
end

list_stim_lowSess_Forg_Reward_Neutral = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Forg_Reward_Neutral{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0025_template.nii,1'];
end

list_stim_lowSess_Forg_Neutral_Reward = []; 
for c1 = 1:length(IDs_low)
    list_stim_lowSess_Forg_Neutral_Reward{c1,1} = [paths.analyses num2str(IDs_low(c1)) '/con_0026_template.nii,1'];
end

list_Nulls = [];
for c1 = 1:length(IDs)
    list_Nulls{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0027_template.nii,1'];
end

list_Button = [];
for c1 = 1:length(IDs)
    list_Button{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0028_template.nii,1'];
end

list_scene_vs_null = [];
for c1 = 1:length(IDs)
    list_scene_vs_null{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0029_template.nii,1'];
end


%% compute models


for CollapseAnalysis=1
% ------- compute: stim_bothsess_Remembered_vs_Forgotten_bothonly ------- %

cd(paths.save2nd); mkdir('stim_bothsess_Remembered_vs_Forgotten_bothonly'); cd stim_bothsess_Remembered_vs_Forgotten_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_bothsess_Remembered_vs_Forgotten_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_bothsess_Remembered_vs_Forgotten_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_bothsess_Forgotten_vs_Remembered_bothonly ------- %

cd(paths.save2nd); mkdir('stim_bothsess_Forgotten_vs_Remembered_bothonly'); cd stim_bothsess_Forgotten_vs_Remembered_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_bothsess_Forgotten_vs_Remembered_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_bothsess_Forgotten_vs_Remembered_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_RewAll_Remembered_Forgotten_bothonly ------- %
cd(paths.save2nd); mkdir('stim_RewAll_Remembered_Forgotten_bothonly'); cd stim_RewAll_Remembered_Forgotten_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_RewAll_Remembered_Forgotten_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_RewAll_Remembered_Forgotten_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_RewAll_Forgotten_Remembered_bothonly ------- %
cd(paths.save2nd); mkdir('stim_RewAll_Forgotten_Remembered_bothonly'); cd stim_RewAll_Forgotten_Remembered_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_RewAll_Forgotten_Remembered_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_RewAll_Forgotten_Remembered_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: fb_lowDA_reward_vs_neutral ------- %
cd(paths.save2nd); mkdir('stim_NeuAll_Remembered_Forgotten_bothonly'); cd stim_NeuAll_Remembered_Forgotten_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_NeuAll_Remembered_Forgotten_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_NeuAll_Remembered_Forgotten_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: fb_rewardAll_vs_neutralAll ------- %

cd(paths.save2nd); mkdir('stim_NeuAll_Forgotten_Remembered_bothonly'); cd stim_NeuAll_Forgotten_Remembered_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_NeuAll_Forgotten_Remembered_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_NeuAll_Forgotten_Remembered_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_highSess_Rew_Remembered_Forgotten ------- %
cd(paths.save2nd); mkdir('stim_highSess_Rew_Remembered_Forgotten'); cd stim_highSess_Rew_Remembered_Forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Rew_Remembered_Forgotten) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Rew_Remembered_Forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_highSess_Rew_Forgotten_Remembered ------- %
cd(paths.save2nd); mkdir('stim_highSess_Rew_Forgotten_Remembered'); cd stim_highSess_Rew_Forgotten_Remembered; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Rew_Forgotten_Remembered) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Rew_Forgotten_Remembered';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_highSess_Neu_Remembered_Forgotten ------- %
cd(paths.save2nd); mkdir('stim_highSess_Neu_Remembered_Forgotten'); cd stim_highSess_Neu_Remembered_Forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Neu_Remembered_Forgotten) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Neu_Remembered_Forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_highSess_Neu_Forgotten_Remembered ------- %
cd(paths.save2nd); mkdir('stim_highSess_Neu_Forgotten_Remembered'); cd stim_highSess_Neu_Forgotten_Remembered; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Neu_Forgotten_Remembered) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Neu_Forgotten_Remembered';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Rew_Remembered_Forgotten ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Rew_Remembered_Forgotten'); cd stim_lowSess_Rew_Remembered_Forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Rew_Remembered_Forgotten) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Rew_Remembered_Forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Rew_Forgotten_Remembered ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Rew_Forgotten_Remembered'); cd stim_lowSess_Rew_Forgotten_Remembered; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Rew_Forgotten_Remembered) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Rew_Forgotten_Remembered';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Neu_Remembered_Forgotten ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Neu_Remembered_Forgotten'); cd stim_lowSess_Neu_Remembered_Forgotten; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Neu_Remembered_Forgotten) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Neu_Remembered_Forgotten';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Neu_Forgotten_Remembered ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Neu_Forgotten_Remembered'); cd stim_lowSess_Neu_Forgotten_Remembered; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Neu_Forgotten_Remembered) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Neu_Forgotten_Remembered';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_RemAll_Reward_Neutral_bothonly ------- %
cd(paths.save2nd); mkdir('stim_RemAll_Reward_Neutral_bothonly'); cd stim_RemAll_Reward_Neutral_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_RemAll_Reward_Neutral_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_RemAll_Reward_Neutral_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_RemAll_Neutral_Reward_bothonly ------- %
cd(paths.save2nd); mkdir('stim_RemAll_Neutral_Reward_bothonly'); cd stim_RemAll_Neutral_Reward_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_RemAll_Neutral_Reward_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_RemAll_Neutral_Reward_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_highSess_Rem_Reward_Neutral ------- %
cd(paths.save2nd); mkdir('stim_highSess_Rem_Reward_Neutral'); cd stim_highSess_Rem_Reward_Neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Rem_Reward_Neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Rem_Reward_Neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_highSess_Rem_Neutral_Reward ------- %
cd(paths.save2nd); mkdir('stim_highSess_Rem_Neutral_Reward'); cd stim_highSess_Rem_Neutral_Reward; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Rem_Neutral_Reward) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Rem_Neutral_Reward';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Rem_Reward_Neutral ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Rem_Reward_Neutral'); cd stim_lowSess_Rem_Reward_Neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Rem_Reward_Neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Rem_Reward_Neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Rem_Neutral_Reward ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Rem_Neutral_Reward'); cd stim_lowSess_Rem_Neutral_Reward; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Rem_Neutral_Reward) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Rem_Neutral_Reward';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_ForgAll_Reward_Neutral_bothonly ------- %
cd(paths.save2nd); mkdir('stim_ForgAll_Reward_Neutral_bothonly'); cd stim_ForgAll_Reward_Neutral_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_ForgAll_Reward_Neutral_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_ForgAll_Reward_Neutral_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_ForgAll_Neutral_Reward_bothonly ------- %
cd(paths.save2nd); mkdir('stim_ForgAll_Neutral_Reward_bothonly'); cd stim_ForgAll_Neutral_Reward_bothonly; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_ForgAll_Neutral_Reward_bothonly) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_ForgAll_Neutral_Reward_bothonly';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_highSess_Forg_Reward_Neutral ------- %
cd(paths.save2nd); mkdir('stim_highSess_Forg_Reward_Neutral'); cd stim_highSess_Forg_Reward_Neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Forg_Reward_Neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Forg_Reward_Neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_highSess_Forg_Neutral_Reward ------- %
cd(paths.save2nd); mkdir('stim_highSess_Forg_Neutral_Reward'); cd stim_highSess_Forg_Neutral_Reward; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highSess_Forg_Neutral_Reward) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highSess_Forg_Neutral_Reward';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_lowSess_Forg_Reward_Neutral ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Forg_Reward_Neutral'); cd stim_lowSess_Forg_Reward_Neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Forg_Reward_Neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Forg_Reward_Neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowSess_Forg_Neutral_Reward ------- %
cd(paths.save2nd); mkdir('stim_lowSess_Forg_Neutral_Reward'); cd stim_lowSess_Forg_Neutral_Reward; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowSess_Forg_Neutral_Reward) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowSess_Forg_Neutral_Reward';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: Nulls ------- %
cd(paths.save2nd); mkdir('Nulls'); cd Nulls; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_Nulls) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Nulls';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: Button ------- %
cd(paths.save2nd); mkdir('Button'); cd Button; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_Button) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Button';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: scene_vs_null ------- %
cd(paths.save2nd); mkdir('scene_vs_null'); cd scene_vs_null; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_scene_vs_null) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'scene_vs_null';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

end

fprintf('\ndone\n')


%% anova contrast lists

% 1
list_stim_rew_high_rem= []; cnt=0;
for id = 1:length(IDs_high)
    cnt=cnt+1;
    list_stim_rew_high_rem{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0030_template.nii,1'];
end

% 2
list_stim_rew_low_rem= []; cnt=0;
for id = 1:length(IDs_low)
    cnt=cnt+1;
    list_stim_rew_low_rem{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0031_template.nii,1'];
end

% 3
list_stim_neu_high_rem = [];
for id = 1:length(IDs_high)
    list_stim_neu_high_rem{id,1} = [paths.analyses num2str(IDs(id)) '/con_0032_template.nii,1'];
end

% 4
list_stim_neu_low_rem = [];
for id = 1:length(IDs_low)
    list_stim_neu_low_rem{id,1} = [paths.analyses num2str(IDs(id)) '/con_0033_template.nii,1'];
end

% 5
list_stim_rew_high_forg= []; cnt=0;
for id = 1:length(IDs_high)
    cnt=cnt+1;
    list_stim_rew_high_forg{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0034_template.nii,1'];
end

% 6
list_stim_rew_low_forg= []; cnt=0;
for id = 1:length(IDs_low)
    cnt=cnt+1;
    list_stim_rew_low_forg{cnt,1} = [paths.analyses num2str(IDs(id)) '/con_0035_template.nii,1'];
end

% 7
list_stim_neu_high_forg = [];
for id = 1:length(IDs_high)
    list_stim_neu_high_forg{id,1} = [paths.analyses num2str(IDs(id)) '/con_0036_template.nii,1'];
end

% 8
list_stim_neu_low_forg = [];
for id = 1:length(IDs_low)
    list_stim_neu_low_forg{id,1} = [paths.analyses num2str(IDs(id)) '/con_0037_template.nii,1'];
end


%% full factorial ANOVA 3 way

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {[paths.parent '2ndLevel/fullfactorial_template_3way_HalfIncluded/']};

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'HighLow';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % they are dependent because they're repeated measures
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'RewNeu';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'RemForg';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1 1 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = list_stim_rew_high_rem;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [2 1 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = list_stim_rew_low_rem;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [1 2 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = list_stim_neu_high_rem;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2 2 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = list_stim_neu_low_rem;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).levels = [1 1 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).scans = list_stim_rew_high_forg;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).levels = [2 1 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).scans = list_stim_rew_low_forg;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).levels = [1 2 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).scans = list_stim_neu_high_forg;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).levels = [2 2 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).scans = list_stim_neu_low_forg;
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch) % run batch

% estimate

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([paths.parent '2ndLevel/fullfactorial_mni_3way_HalfIncluded/SPM.mat']);
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch) % run batch


%%

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager

matlabbatch{1}.spm.stats.factorial_design.dir = {'/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/2ndLevel/fullfactorial_template_3way_HalfIncluded'};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'HighLow';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'RewNeu';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'RemForg';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1
                                                                    1];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4003/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4007/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4027/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0030_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0030_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    1
                                                                    1];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4006/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4009/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4029/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0031_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0031_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [1
                                                                    2
                                                                    1];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4003/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4007/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4027/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0032_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0032_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                    2
                                                                    1];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4006/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4009/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4029/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0033_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0033_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).levels = [1
                                                                    1
                                                                    2];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4003/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4007/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4027/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0034_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0034_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).levels = [2
                                                                    1
                                                                    2];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4006/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4009/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4029/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0035_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0035_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).levels = [1
                                                                    2
                                                                    2];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4003/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4007/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4027/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0036_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0036_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).levels = [2
                                                                    2
                                                                    2];

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).scans = {
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4001/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4002/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4004/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4005/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4006/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4008/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4009/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4010/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4011/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4012/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4013/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4014/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4015/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4016/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4017/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4018/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4019/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4020/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4021/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4022/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4023/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4024/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4025/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4026/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4028/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4029/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4030/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4031/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4032/con_0037_template.nii,1'
                                                                   '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions_memory/4033/con_0037_template.nii,1'
                                                                   };

matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
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

spm_jobman('run', matlabbatch) % run batch

%%

for id=1:length(IDs)

    copyfile(['/Volumes/korokdorf/MRPET/coreg_mri/bothsessions/' num2str(IDs(id)) '/data/NLantsreg_meanEPI_to_template_t12epi.nii.gz'],...
        ['/Users/alex/Dropbox/paperwriting/MRPET/data/meanEPI/' num2str(IDs(id)) '_meanfunc_template2.nii.gz'])

end


%% calculate tSNR
% 
% for k = 1:length(IDs)
%     
%     file_end = ['/Users/yeojin/Desktop/E_data/EB_cleaned/EBC_mri/pilot_3T_2sess/MainTask/' num2str(IDs(k)) '_2/'];
%     cd(file_end)
%     
%     disp('tSNr computation')
%     FileNameSPMMat = {fullfile(file_end,'SPM.mat')};
%     load(char(FileNameSPMMat))
%     % S=size(SPM.xX.X);
%     numbeta = 1;
%     
%     % matlabbatch{1}.spm.util.imcalc.input = {fullfile(file_end,strcat('beta_',sprintf('%04d',S(2)),'.nii'));fullfile(file_end,'ResMS.nii')};
%     matlabbatch{1}.spm.util.imcalc.input = {fullfile(file_end,strcat('beta_',sprintf('%04d',numbeta(1)),'.nii'));fullfile(file_end,'ResMS.nii')};
%     
%     matlabbatch{1}.spm.util.imcalc.output = ['tSNR_' num2str(IDs(k)) '_2_' strcat('beta_',sprintf('%04d',numbeta(1)))];
%     matlabbatch{1}.spm.util.imcalc.outdir = {file_end};
%     matlabbatch{1}.spm.util.imcalc.expression = '(i1./sqrt(i2))';
%     
%     matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%     matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
%     matlabbatch{1}.spm.util.imcalc.options.mask = 0;
%     matlabbatch{1}.spm.util.imcalc.options.interp = -3;
%     matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
%     spm_jobman('run', matlabbatch);
%     clear matlabbatch
%     
% end