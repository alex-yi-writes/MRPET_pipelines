%% PET 1st-level analysis pipeline

%% work log

%   12-12-2019    created the script

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.analyses= [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/MRI/'];
paths.behav   = [paths.parent 'E_data/EA_raw/EAC_behav/MRPET/'];

% add toolboxes and functions
addpath(paths.spm)
addpath(paths.funx)

% IDs
IDs  = [4001];
days = [0 2];
d1m  = [1 2]; % 1=immediate 2=delayed
d2m  = [1 2];

% load experimental details
expdat = [];
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            fname_beh{i1,d} = {NaN};
            expdat{i1,d} = {NaN};
        else
            fname_beh{i1,d}     = [ num2str(IDs(i1)) '_' num2str(days(i1,d)) '.mat' ];
            expdat{i1,d} = load([paths.behav fname_beh{i1,d}]);
        end
    end
end

TR = 3.6;

fprintf('\n preparation done \n')
%% start

spm fmri  % open progress window

for id = 1:length(IDs)
    
    for d = 1:2
        if days(id,d) == 0
        else
            
            fprintf('\n model estimation for: ID %d \n', IDs(id))
            
            % define contingency
            clear rewcond neucond
            eval(['rewcond = str2num(expdat{id,d}.dat.day' num2str(d) '.RewardCategory(9))']);
            if rewcond == 1
                neucond = 2;
            elseif rewcond == 3
                neucond = 4;
            elseif rewcond == 2
                neucond = 1;
            elseif rewcond == 4
                neucond = 3;
            end
            trlinfo = eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,2:3))']);
            trlinfo(trlinfo(:,2)==0,:) = [];            
            
            % set up workspace
            cd(paths.analyses); mkdir([num2str(IDs(id)) '_' num2str(d)]);
            clear dir_model
            dir_model = [paths.analyses num2str(IDs(id)) '_' num2str(d) '/'];
            dir_movep = [paths.preproc num2str(IDs(id)) '_' num2str(d) '/rp_a' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.txt'];
            
            %% build model
            
            clear matlabbatch
            
            spm_jobman('initcfg');      % initiate job manager
            
            matlabbatch{1}.spm.stats.fmri_spec.dir               = cellstr(dir_model); % where will the model be saved?
            matlabbatch{1}.spm.stats.fmri_spec.timing.units      = 'secs'; % scans / secs
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT         = TR; % TR in seconds
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t     = 51; % volume size
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0    = 25; % microtime onset
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans        = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(d) '/' ...
                'swra' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii']);
            
            
            % reward condition
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Rewards';
            eval...
                (['matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset=expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(trlinfo(:,1)==rewcond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Neutrals';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(trlinfo(:,1)==neucond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Nulls';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.null-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Response';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim+2.5-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'Rewards_feedback';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(trlinfo(:,1)==rewcond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'neutrals_feedback';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(trlinfo(:,1)==neucond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(dir_movep); % movement parameters
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
            
            % here you set up the contrast matrices: Rewards - Neutrals - Fix
            c_t_1   = [1 1 -2];        % stim vs fix
            c_t_2   = [1 -1];          % rew > neu
            c_t_3   = [-1 1];          % neu > rew
            c_t_4   = [1 1 1];         % when there's something on the screen
            c_t_5   = [0 0 0 1];
            
            c_t_r   = [1 0 -1];        % rew > fix
            c_t_p   = [0 1 -1];        % neu > fix
            
            c_t_6   = [0 0 0 0 1 -1];  % feedbacks rew > neu
            c_t_7   = [0 0 0 0 -1 1];  % feedbacks neu > rew
            c_t_8   = [0 0 0 0 1 0];   % feedbacks rew
            c_t_9   = [0 0 0 0 0 1];   % feedbacks neu
            
            c_t_10  = [-1 0 0 0 1];    % reward: feedback > stim
            c_t_11  = [0 -1 0 0 0 1];  % neutral: feedback > stim
            
            c_t_12  = [1 0 0 0 -1];    % reward: stim > feedback
            c_t_13  = [0 1 0 0 0 -1];  % neutral: stim > feedback
            
            %     c_F_1 	  = [2/3 -1/3 -1/3; ...
            %                 -1/3 2/3 -1/3; ...
            %                 -1/3 -1/3 2/3 ];  % ess_0005: rew vs neu vs fix
            
            % batch setup
            clear matlabbatch
            spm_jobman('initcfg');      % initiate job manager
            
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
            
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stim-Null';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = c_t_1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'rew > neu';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = c_t_2;
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'neu > rew';
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = c_t_3;
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'manipulation v intercept';
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = c_t_4;
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'response v intercept';
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = c_t_5;
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'reward v null';
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = c_t_r;
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'neutral v null';
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = c_t_p;
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'feedbacks: rew > neu';
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = c_t_6;
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'feedbacks: neu > rew';
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = c_t_7;
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'feedbacks: rew';
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = c_t_8;
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'feedbacks: neu';
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = c_t_9;
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'reward: feedback > stim';
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = c_t_10;
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'neutral: feedback > stim';
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = c_t_11;
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'reward: stim > feedback';
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = c_t_12;
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.name = 'neutral: stim > feedback';
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.weights = c_t_13;
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
            
            %     matlabbatch{1}.spm.stats.con.consess{6}.fcon.name = 'ANOVA (rew vs neu vs fix)';
            %     matlabbatch{1}.spm.stats.con.consess{6}.fcon.weights = c_F_1;
            %     matlabbatch{1}.spm.stats.con.consess{6}.fcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
            
            spm_jobman('run', matlabbatch) % run batch
            
        end
    end
    
end

fprintf('\ndone\n')

%% fMRI 2nd-level analysis pipeline
%% set environmental variables

clear;
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.analyses= [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/MRI/'];
paths.behav   = [paths.parent 'E_data/EA_raw/EAC_behav/MRPET/'];

paths.temp    = [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/MRI/'];
paths.save2nd = [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/MRI/2ndLvL/'];
paths.doc     = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/MRI/'];

% add toolboxes and functions
addpath(paths.spm)
addpath(paths.funx)

% IDs
IDs  = [4001];
days = [0 2];
d1m  = [1 2]; % 1=immediate 2=delayed
d2m  = [1 2];

% load experimental details
expdat = [];
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            fname_beh{i1,d} = {NaN};
            expdat{i1,d} = {NaN};
        else
            fname_beh{i1,d}     = [ num2str(IDs(i1)) '_' num2str(days(i1,d)) '.mat' ];
            expdat{i1,d} = load([paths.behav fname_beh{i1,d}]);
        end
    end
end

TR = 3.6;

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window


list_stim_fix = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_stim_fix{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0001_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_rew_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0002_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_neu_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0003_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_manipulation = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_manipulation{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0004_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_resp = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_resp{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0005_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0006_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0007_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_feedback_rew_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_feedback_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0008_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_feedback_neu_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_feedback_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0009_mni.nii,1'];
        end
    end
end
clear i1 d cnt1


list_feedback_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_feedback_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0010_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_feedback_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1=cnt1+1;
            list_feedback_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_' num2str(days(i1,d)) '/con_0011_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

%% compute models

% ------- compute: stim vs. fixation cross ------- %

cd(paths.save2nd); mkdir('Stim_v_Fix'); cd Stim_v_Fix; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_fix) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stim-Fixation';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Stim_v_Fix_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Stim_v_Fix_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > neutral ------- %

cd(paths.save2nd); mkdir('Rew_Neu'); cd Rew_Neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_rew_neu) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Rew>Neu';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Rew_Neu_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Rew_Neu_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: neutral > reward ------- %
cd(paths.save2nd); mkdir('Neu_Rew'); cd Neu_Rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_neu_rew) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Neu>Rew';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Neu_Rew_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_Rew_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('manipulation'); cd manipulation; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_manipulation) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'manipulation_intercept';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('manipulation_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'manipulation_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('response'); cd response; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_resp) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'response_intercept';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('response_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'response_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > fix ------- %

cd(paths.save2nd); mkdir('Rew'); cd Rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_rew_neu) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Rew>fix';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Rew_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Rew_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: neutral > fix ------- %
cd(paths.save2nd); mkdir('Neu'); cd Neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_neu_rew) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Neu>fix';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Neu_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: neutral > fix ------- %
cd(paths.save2nd); mkdir('Neu'); cd Neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_neu_rew) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Neu>fix';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Neu_uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: feedback: rewards > neutral ------- %
cd(paths.save2nd); mkdir('fb_rew_neu'); cd fb_rew_neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_feedback_rew_neu) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Feedback: Rew>Neu';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Feedback: Rew>Neu uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Feedback: Rew>Neu uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: feedback: neutral > rewards ------- %
cd(paths.save2nd); mkdir('fb_neu_rew'); cd fb_neu_rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_feedback_rew_neu) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Feedback: Neu>Rew';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% document results
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('Feedback: Neu>Rew uncorr');
matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'pdf';
matlabbatch{1}.spm.stats.results.write.none = 1;
spm_jobman('run', matlabbatch) % run batch
clear matlabbatch

doc_list_uncorr = dir('*.pdf'); % list docs
tmpdir = pwd;
for i3= 1:length(doc_list_uncorr) % move docs
    movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Feedback: Neu>Rew uncorr.pdf']);
end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

fprintf('\ndone\n')


%% compare days in 2nd-lvl

% stim
% make lists
list_D1_stim_rew_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D1_stim_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0002_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_D2_stim_rew_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D2_stim_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0002_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_D1_stim_neu_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D1_stim_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0003_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_D2_stim_neu_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D2_stim_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0003_mni.nii,1'];
        end
    end
end
clear i1 d cnt1


% compute
% ------- compute: reward > neutral ------- %

cd(paths.save2nd); mkdir('compare_days'); cd('compare_days'); mkdir('stim'); cd('stim'); mkdir('Rew_Neu'); cd Rew_Neu; dir_spm = pwd;

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_spm);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = cellstr(list_D1_stim_rew_neu);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = cellstr(list_D2_stim_rew_neu);
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch) % run batch

fMRI_estimate([dir_spm '/SPM.mat'])

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Rew>Neu: D1>D2';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Rew>Neu: D2>D1';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

% ------- compute: neutral > reward ------- %

cd(paths.save2nd); cd('compare_days');cd('stim'); mkdir('Neu_Rew'); cd Neu_Rew; dir_spm = pwd;

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_spm);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = cellstr(list_D1_stim_neu_rew);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = cellstr(list_D2_stim_neu_rew);
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch) % run batch

fMRI_estimate([dir_spm '/SPM.mat'])

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Neu>Rew: D1>D2';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Neu>Rew: D2>D1';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% fb
% make lists
list_D1_fb_rew_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D1_fb_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0008_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_D2_fb_rew_neu = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D2_fb_rew_neu{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0008_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_D1_fb_neu_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D1_fb_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0009_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_D2_fb_neu_rew = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_D2_fb_neu_rew{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0009_mni.nii,1'];
        end
    end
end
clear i1 d cnt1


% compute
% ------- compute: reward > neutral ------- %

cd(paths.save2nd); cd('compare_days'); mkdir('fb'); cd('fb'); mkdir('Rew_Neu'); cd Rew_Neu; dir_spm = pwd;

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_spm);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = cellstr(list_D1_fb_rew_neu);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = cellstr(list_D2_fb_rew_neu);
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch) % run batch

fMRI_estimate([dir_spm '/SPM.mat'])

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Rew>Neu: D1>D2';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Rew>Neu: D2>D1';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

% ------- compute: neutral > reward ------- %

cd(paths.save2nd); cd('compare_days'); cd('fb'); mkdir('Neu_Rew'); cd Neu_Rew; dir_spm = pwd;

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_spm);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = cellstr(list_D1_fb_neu_rew);
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = cellstr(list_D2_fb_neu_rew);
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch) % run batch

fMRI_estimate([dir_spm '/SPM.mat'])

clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Neu>Rew: D1>D2';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Neu>Rew: D2>D1';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


%% stim vs. feedback

list_fb_stim_rewards = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_fb_stim_rewards{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0012_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_fb_stim_neutrals = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_fb_stim_neutrals{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0013_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_stim_fb_rewards = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_stim_fb_rewards{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0014_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_stim_fb_neutrals = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_stim_fb_neutrals{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0015_mni.nii,1'];
        end
    end
end
clear i1 d cnt1


% ------- compute: rewards: feedback > stimulus ------- %
cd(paths.save2nd); mkdir('fb_vs_stim'); cd fb_vs_stim; mkdir('rew'); cd rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_stim_rewards) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'rewards: feedback > stimulus';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% % document results
% clear matlabbatch
% spm_jobman('initcfg');      % initiate job manager
% matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
% matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('rewards: feedback > stimulus');
% matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
% matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
% matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
% matlabbatch{1}.spm.stats.results.conspec.extent = 0;
% matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
% matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
% matlabbatch{1}.spm.stats.results.units = 1;
% matlabbatch{1}.spm.stats.results.print = 'pdf';
% matlabbatch{1}.spm.stats.results.write.none = 1;
% spm_jobman('run', matlabbatch) % run batch
% clear matlabbatch

% doc_list_uncorr = dir('*.pdf'); % list docs
% tmpdir = pwd;
% for i3= 1:length(doc_list_uncorr) % move docs
%     movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_Rew_uncorr.pdf']);
% end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

% ------- compute: neutrals: feedback > stimulus ------- %
cd(paths.save2nd); cd fb_vs_stim; mkdir('neu'); cd neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_stim_neutrals) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'neutrals: feedback > stimulus';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% % document results
% clear matlabbatch
% spm_jobman('initcfg');      % initiate job manager
% matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
% matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('rewards: feedback > stimulus');
% matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
% matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
% matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
% matlabbatch{1}.spm.stats.results.conspec.extent = 0;
% matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
% matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
% matlabbatch{1}.spm.stats.results.units = 1;
% matlabbatch{1}.spm.stats.results.print = 'pdf';
% matlabbatch{1}.spm.stats.results.write.none = 1;
% spm_jobman('run', matlabbatch) % run batch
% clear matlabbatch

% doc_list_uncorr = dir('*.pdf'); % list docs
% tmpdir = pwd;
% for i3= 1:length(doc_list_uncorr) % move docs
%     movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_Rew_uncorr.pdf']);
% end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: rewards: stim > fb ------- %
cd(paths.save2nd); cd stim_vs_fb; cd rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_fb_rewards) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'rewards: stim > fb';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% % document results
% clear matlabbatch
% spm_jobman('initcfg');      % initiate job manager
% matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
% matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('rewards: feedback > stimulus');
% matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
% matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
% matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
% matlabbatch{1}.spm.stats.results.conspec.extent = 0;
% matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
% matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
% matlabbatch{1}.spm.stats.results.units = 1;
% matlabbatch{1}.spm.stats.results.print = 'pdf';
% matlabbatch{1}.spm.stats.results.write.none = 1;
% spm_jobman('run', matlabbatch) % run batch
% clear matlabbatch

% doc_list_uncorr = dir('*.pdf'); % list docs
% tmpdir = pwd;
% for i3= 1:length(doc_list_uncorr) % move docs
%     movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_Rew_uncorr.pdf']);
% end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %

% ------- compute: neutrals: stim > fb ------- %
cd(paths.save2nd); cd stim_vs_fb; cd neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_fb_neutrals) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'neutrals: stim > fb';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

% % document results
% clear matlabbatch
% spm_jobman('initcfg');      % initiate job manager
% matlabbatch{1}.spm.stats.results.spmmat = cellstr([dir_spm '/SPM.mat']);
% matlabbatch{1}.spm.stats.results.conspec.titlestr = cellstr('rewards: feedback > stimulus');
% matlabbatch{1}.spm.stats.results.conspec.contrasts = inf;
% matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
% matlabbatch{1}.spm.stats.results.conspec.thresh = 0.01;
% matlabbatch{1}.spm.stats.results.conspec.extent = 0;
% matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
% matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
% matlabbatch{1}.spm.stats.results.units = 1;
% matlabbatch{1}.spm.stats.results.print = 'pdf';
% matlabbatch{1}.spm.stats.results.write.none = 1;
% spm_jobman('run', matlabbatch) % run batch
% clear matlabbatch

% doc_list_uncorr = dir('*.pdf'); % list docs
% tmpdir = pwd;
% for i3= 1:length(doc_list_uncorr) % move docs
%     movefile([tmpdir '/' doc_list_uncorr(i3,1).name],[paths.doc 'Neu_Rew_uncorr.pdf']);
% end
clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


%% Two-Way ANOVA

% collect the contrast images
% stim
list_less_neu_stim = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_less_neu_stim{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0007_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_less_rew_stim = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_less_rew_stim{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0006_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_more_neu_stim = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_more_neu_stim{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0007_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_more_rew_stim = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_more_rew_stim{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0006_mni.nii,1'];
        end
    end
end
clear i1 d cnt1


% fb
list_less_neu_fb = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_less_neu_fb{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0011_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_less_rew_fb = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_less_rew_fb{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_2/con_0010_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_more_neu_fb = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_more_neu_fb{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0011_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

list_more_rew_fb = []; cnt1=0;
for i1 = 1:length(IDs)
    for d = 1
        if days(i1,d) == 0
        else
            cnt1 = cnt1+1;
            list_more_rew_fb{cnt1,1} = [paths.analyses num2str(IDs(i1)) '_1/con_0010_mni.nii,1'];
        end
    end
end
clear i1 d cnt1

%% VERSION 1

%%%%%%%%%%%%%%%%%%   factors table   %%%%%%%%%%%%%%%
%
%                  less(f1)          more(f1)
%
%
% neutral(f2)         1  1            2  1  
%
% reward(f2)          1  2            2  2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%   design matrix   %%%%%%%%%%%%%%%
%
% less & 
%    neutral
%
%
%             less & 
%                 reward
%
%
%                         more &       
%                            neutral  
%  
%
%                                      more &
%                                          reward
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr([paths.analyses '2ndLvL/two-way_anova1/']);
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'less/more'; % factor 1, 2 levels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'neutral/reward'; % factor 2, 2 levels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1]; % 1st level of 1st factor and 1st level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = cellstr(list_less_neu_stim); % frequent neutral

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1
                                                                    2]; % 1st level of 1st factor and 2nd level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = cellstr(list_less_rew_stim); % rare rewards

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2
                                                                    1]; % 2nd level of 1st factor and 1st level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = cellstr(list_more_neu_stim); % rare neutral

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                    2]; % 2nd level of 1st factor and 2nd level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = cellstr(list_more_rew_stim); % frequent neutral
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

fMRI_estimate([paths.analyses '2ndLvL/two-way_anova1/SPM.mat'])



%% VERSION 2
% variables dependent


%%%%%%%%%%%%%%%%%%   factors table   %%%%%%%%%%%%%%%
%
%                  less(f1)          more(f1)
%
%
% reward (f2)         1  1            2  1  
%
% neutral(f2)         1  2            2  2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%   design matrix   %%%%%%%%%%%%%%%
%
% less & 
%    reward
%
%
%             less & 
%                 neutral
%
%
%                         more &       
%                            reward  
%  
%
%                                      more &
%                                          neutral
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear matlabbatch
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr([paths.analyses '2ndLvL/two-way_anova2']);

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'less/more'; % factor 1, 2 levels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'neutral/reward'; % factor 2, 2 levels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1]; % 1st level of 1st factor and 1st level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = cellstr(list_less_rew_stim);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1
                                                                    2]; % 1st level of 1st factor and 2nd level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = cellstr(list_less_neu_stim);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2
                                                                    1]; % 2nd level of 1st factor and 1st level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = cellstr(list_more_rew_stim);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                    2]; % 2nd level of 1st factor and 2nd level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = cellstr(list_more_neu_stim);
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

fMRI_estimate([paths.analyses '2ndLvL/two-way_anova2/SPM.mat'])


%% version 2-2
% variables independent


clear matlabbatch
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr([paths.analyses '2ndLvL/two-way_anova2-2/']);

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'less/more'; % factor 1, 2 levels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'neutral/reward'; % factor 2, 2 levels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1]; % 1st level of 1st factor and 1st level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = cellstr(list_less_rew_stim);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1
                                                                    2]; % 1st level of 1st factor and 2nd level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = cellstr(list_less_neu_stim);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2
                                                                    1]; % 2nd level of 1st factor and 1st level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = cellstr(list_more_rew_stim);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                    2]; % 2nd level of 1st factor and 2nd level of 2nd factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = cellstr(list_more_neu_stim);
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

fMRI_estimate([paths.analyses '2ndLvL/two-way_anova2-2/SPM.mat'])

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