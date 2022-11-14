%% PET 1st-level analysis pipeline

%% work log

%   12-12-2019    created the script

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Volumes/ALEX3/MRPET/';
paths.funx    = '/Users/alex/Dropbox/literatures_IKND/BB_analyses/BBC_MRI/analyses_functions/';
paths.preproc = [paths.parent 'img/'];
paths.analyses= [paths.parent 'analysis/'];
paths.behav   = [paths.parent 'behav/'];

% add toolboxes and functions
% addpath(paths.spm)
addpath(paths.funx)

% IDs
IDs = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2]; 
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 2; 1 2];

%% load experimental details

expdat = []; onsets_temp=[]; onsets=[]; length_scan1=[];
for id = 1:length(IDs)
    
    for d = 1:2
        if days(id,d) == 0
            fname_beh{id,d} = {NaN};
            expdat{id,d} = {NaN};
            contingency{id,d} = [];
        else
            
             % set up workspace
            mkdir([paths.analyses num2str(IDs(id)) '_' num2str(d)])
            
            
            fname_beh{id,d}     = [ num2str(IDs(id)) '_' num2str(days(id,d)) '.mat' ];
            expdat{id,d}        = load([paths.behav fname_beh{id,d}]);
            
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
            tmptrials = cell2mat(eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,3)']));
            ind_null  = tmptrials==0;
            ind_rew   = tmptrials==contingency{id,d}(1); ind_rew(ind_null)=[];
            ind_neu   = tmptrials==contingency{id,d}(2); ind_neu(ind_null)=[];
            
            
            % sort onsets
            onsets{id,d}.stim.reward  = tmponsets.stim(ind_rew);
            onsets{id,d}.stim.neutral = tmponsets.stim(ind_neu);
            onsets{id,d}.null         = tmponsets.null;
            onsets{id,d}.resp         = tmponsets.resp;
            onsets{id,d}.fb.reward    = tmponsets.cue(ind_rew);
            onsets{id,d}.fb.neutral   = tmponsets.cue(ind_neu);
            
            %% make physio + realignment parameters
            
            % read physio parameters
            physioFile  = [paths.preproc num2str(IDs(id)) '_' num2str(d) '/multiple_regressors.txt'];
            physiodata    = importdata(physioFile);   
            
            % read realignment parameters
            realignFile = [paths.preproc num2str(IDs(id)) '_' num2str(d) '/rp_a' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.txt'];
            realigndata = importdata(realignFile);
            
            % concatenate
            R = [physiodata realigndata];
            
            save([paths.preproc num2str(IDs(id)) '_' num2str(d) '/reg_all.mat'],'R');
            clear R realigndata realignFile physiodata physioFile
            
        end
    end
end

TR = 3.6;

fprintf('\n preparation done \n')

%% treat each as single session data?

IDs_singleSession=[];
cnt=0;
for id=1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            disp('missing')
        else
            
            cnt=cnt+1;
            IDs_singleSession{cnt,1}=[num2str(IDs(id)) num2str(d)];
            
        end
    end
    
end


%% start

spm fmri  % open progress window

for id = 9:length(IDs)
    
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
            cd(paths.analyses); mkdir([num2str(IDs(id)) num2str(d)]);
            clear dir_model
            dir_model = [paths.analyses num2str(IDs(id)) num2str(d) '/'];
            dir_multiparam = [paths.preproc num2str(IDs(id)) '_' num2str(d) '/reg_all.mat'];
            
            %% build model
            
            clear matlabbatch
            
            spm_jobman('initcfg');      % initiate job manager
            
            matlabbatch{1}.spm.stats.fmri_spec.dir               = cellstr(dir_model); % where will the model be saved?
            matlabbatch{1}.spm.stats.fmri_spec.timing.units      = 'secs'; % scans / secs
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT         = TR; % TR in seconds
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t     = 51; % volume size
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0    = 25; % microtime onset
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans        = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(d) '/' ...
                'sua' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii']);
            
            
            % reward condition
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'StimRewards';
            eval...
                (['matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset=expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.stim(trlinfo(:,1)==rewcond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'StimNeutrals';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(trlinfo(:,1)==neucond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'FeedbackRewards';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(trlinfo(:,1)==rewcond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'FeedbackNeutrals';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.cue(trlinfo(:,1)==neucond)-(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'Nulls';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.null-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
            
            
            % --------- left and right response differences -------- %
            
            clear respLind respRind tmpresp
            tmpresp=eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.keypress']);
            tmpresp(tmpresp==0)=[];
            respLind=tmpresp==28;
            respRind=tmpresp==29;
            
            % ------------------------------------------------------ %
            
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'Response_L';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(respLind)+2.5-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'Response_R';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= eval...
                (['expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.stim(respRind)+2.5-(expdat{id,d}.dat.day' num2str(d)...
                '.maintask.results.SOT.raw.trig_1st);']);
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
            
            %% statistical inference: model E
            
            firstlvl_dir        = [dir_model 'SPM.mat'];
            
            % here you set up the contrast matrices
            
            c_stim_v_null   = [1 1 1 1 -4];        % stim vs null
            
            %%%---
            c_stim_rew_v_neu   = [1 -1];          % stim: rew > neu (PET DA-nergic, 2nd)
            c_stim_neu_v_rew   = [-1 1];          % stim: neu > rew
            %%%---
            
            %%%---
            c_fb_rew_v_neu   = [0 0 1 -1];  % feedback: rew > neu (PET DA-nergic, 1st)
            c_fb_neu_v_rew   = [0 0 1 -1];  % feedback: neu > rew
            %%%---

            %%%---
            c_reward_fb_v_stim  = [-1 0 1];    % reward: feedback > stim
            c_reward_stim_v_fb  = [1 0 -1];    % reward: stim > feedback
            
            c_neutral_fb_v_stim  = [0 -1 0 1];    % neutral: feedback > stim
            c_neutral_stim_v_fb  = [0 1 0 -1];    % neutral: stim > feedback
            %%%---
            
            c_resp_L    = [0 0 0 0 0 1];
            c_resp_R    = [0 0 0 0 0 0 1];
            
            c_stim_rew_v_null = [1 0 0 0 -1];
            c_stim_neu_v_null = [0 1 0 0 -1];
            c_fb_rew_v_null   = [0 0 1 0 -1];
            c_fb_neu_v_null   = [0 0 0 1 -1];
            

            % batch setup
            clear matlabbatch
            spm_jobman('initcfg');      % initiate job manager
            
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);
            
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stim-Null';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = c_stim_v_null;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Stimuli: rew > neu';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = c_stim_rew_v_neu;
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Stimuli: neu > rew';
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = c_stim_neu_v_rew;
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Feedback: rew > neu';
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = c_fb_rew_v_neu;
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Feedback: neu > rew';
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = c_fb_neu_v_rew;
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'reward: feedback > stim';
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = c_reward_fb_v_stim;
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'reward: stim > feedback';
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = c_reward_stim_v_fb;
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'neutral: feedback > stim';
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = c_neutral_fb_v_stim;
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'neutral: stim > feedback';
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = c_neutral_stim_v_fb;
            matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'response Left';
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = c_resp_L;
            matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'response Right';
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = c_resp_R;
            matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'StimReward > Null';
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = c_stim_rew_v_null;
            matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'StimNeutral > Null';
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = c_stim_neu_v_null;
            matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'FeedbackReward > Null';
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = c_fb_rew_v_null;
            matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.name = 'FeedbackNeutral > Null';
            matlabbatch{1}.spm.stats.con.consess{15}.tcon.weights = c_fb_neu_v_null;
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


%% move contrasts for registration

path_source      ='/Volumes/ALEX3/MRPET/analysis/';
path_destination ='/Volumes/ALEX3/MRPET/coreg/';

mrpetID=[]; c1=0;
for id=1:length(IDs)
    
    for d = 1:2
        if days(id,d) == 0
        else
            
            c1=c1+1;
            mrpetID{c1,1}=[num2str(IDs(id)) num2str(d)];
            
%             mkdir([path_destination num2str(IDs(id)) num2str(d) '/data/'])
            
            for cnt=1:15
            copyfile([path_source num2str(IDs(id)) num2str(d) '/con_00' sprintf('%02i',cnt) '.nii'],...
                [path_destination num2str(IDs(id)) num2str(d) '/data/con_00' sprintf('%02i',cnt) '.nii'])
            end
            copyfile([path_source num2str(IDs(id)) num2str(d) '/mask.nii'],...
                [path_destination num2str(IDs(id)) num2str(d) '/data/mask.nii'])
            
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/T1WB.nii'])
%             
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/meanEPI.nii'])
%             
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_GRE3D' num2str(d) '_2.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/MTw.nii'])
%             
    
        end
    end
    

end

%%


path_source      ='/Volumes/ALEX3/MRPET/analysis/';
path_destination ='/Volumes/ALEX3/MRPET/coreg/';

for id=1:length(IDs)
    
    for d = 1:2
        if days(id,d) == 0
        else
            
%             c1=c1+1;
%             mrpetID{c1,1}=[num2str(IDs(id)) num2str(d)];
            
%             mkdir([path_destination num2str(IDs(id)) num2str(d) '/data/'])
            
            for cnt=1:15
            copyfile([path_destination num2str(IDs(id)) num2str(d) '/data/con_00' sprintf('%02i',cnt) '_mni.nii'],...
                [path_source num2str(IDs(id)) num2str(d) '/con_00' sprintf('%02i',cnt) '_mni.nii'])
            end
            copyfile([path_destination num2str(IDs(id)) num2str(d) '/data/mask_mni.nii'],...
                [path_source num2str(IDs(id)) num2str(d) '/mask_mni.nii'])
            
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/T1WB.nii'])
%             
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/meanEPI.nii'])
%             
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_GRE3D' num2str(d) '_2.nii'],...
%                 [path_destination num2str(IDs(id)) num2str(d) '/data/MTw.nii'])
%             
    
        end
    end
    

end


%% fMRI 2nd-level analysis pipeline
%% set environmental variables

clear;
warning('off','all');

% paths
paths = [];
paths.parent  = '/Volumes/ALEX3/MRPET/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
paths.funx    = ['/Users/alex/Dropbox/literatures_IKND/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent 'img/'];
paths.analyses= [paths.parent 'analysis/'];
paths.behav   = [paths.parent 'behav/'];

paths.temp    = [paths.parent 'img/'];
paths.save2nd = [paths.parent 'analysis/2nd/'];
% paths.doc     = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/MRI/'];

% add toolboxes and functions
% addpath(paths.spm)
addpath(paths.funx)

% IDs
IDs = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2]; 
% d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
% d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 2];
load('/Volumes/ALEX3/MRPET/IDs_singleSession.mat')

% load experimental details
expdat = [];
for i1 = 1:length(IDs_singleSession)
    expdat{i1,1} = load([paths.behav IDs_singleSession{i1}(1:4) '_' IDs_singleSession{i1}(5) '.mat']);
end

TR = 3.6;
IDhighDA=[]; IDlowDA=[]; h1=0;l1=0;
for i1 = 1:length(IDs_singleSession)
    
    if str2num(IDs_singleSession{i1}(5))==1
        h1=h1+1;
        IDhighDA{h1,1}=IDs_singleSession{i1};
    elseif str2num(IDs_singleSession{i1}(5))==2
        l1=l1+1;
        IDlowDA{l1,1}=IDs_singleSession{i1};
    end 
end

fprintf('\n preparation done \n')
%% start

% spm fmri  % open progress window

list_c_stim_v_null = [];
for i1 = 1:length(IDs_singleSession)
    list_c_stim_v_null{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0001_mni.nii,1'];
end


list_c_stim_rew_v_neu = [];
for i1 = 1:length(IDs_singleSession)
    list_c_stim_rew_v_neu{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0002_mni.nii,1'];
end


list_c_stim_neu_v_rew = [];
for i1 = 1:length(IDs_singleSession)
    list_c_stim_neu_v_rew{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0003_mni.nii,1'];
end


list_c_fb_rew_v_neu = [];
for i1 = 1:length(IDs_singleSession)
    list_c_fb_rew_v_neu{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0004_mni.nii,1'];
end


list_c_fb_neu_v_rew = [];
for i1 = 1:length(IDs_singleSession)
    list_c_fb_neu_v_rew{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0005_mni.nii,1'];
end


list_c_reward_fb_v_stim = [];
for i1 = 1:length(IDs_singleSession)
    list_c_reward_fb_v_stim{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0006_mni.nii,1'];
end


list_c_reward_stim_v_fb = [];
for i1 = 1:length(IDs_singleSession)
    list_c_reward_stim_v_fb{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0007_mni.nii,1'];
end


list_c_neutral_fb_v_stim = [];
for i1 = 1:length(IDs_singleSession)
    list_c_neutral_fb_v_stim{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0008_mni.nii,1'];
end


list_c_neutral_stim_v_fb = [];
for i1 = 1:length(IDs_singleSession)
    list_c_neutral_stim_v_fb{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0009_mni.nii,1'];
end



list_c_resp_L = [];
for i1 = 1:length(IDs_singleSession)
    list_c_resp_L{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0010_mni.nii,1'];
end


list_c_resp_R = [];
for i1 = 1:length(IDs_singleSession)
    list_c_resp_R{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0011_mni.nii,1'];
end

list_c_stim_rew_v_null = [];
for i1 = 1:length(IDs_singleSession)
    list_c_stim_rew_v_null{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0012_mni.nii,1'];
end

list_c_stim_neu_v_null = [];
for i1 = 1:length(IDs_singleSession)
    list_c_stim_neu_v_null{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0013_mni.nii,1'];
end

list_c_fb_rew_v_null = [];
for i1 = 1:length(IDs_singleSession)
    list_c_fb_rew_v_null{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0014_mni.nii,1'];
end

list_c_fb_neu_v_null = [];
for i1 = 1:length(IDs_singleSession)
    list_c_fb_neu_v_null{i1,1} = [paths.analyses IDs_singleSession{i1} '/con_0015_mni.nii,1'];
end

%% compute models

% ------- compute: stim vs. fixation cross ------- %

cd(paths.save2nd); mkdir('Stim-Null'); cd Stim-Null; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_stim_v_null) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stim-Null';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > neutral ------- %

cd(paths.save2nd); mkdir('stim_rew_v_neu'); cd stim_rew_v_neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_stim_rew_v_neu) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stimuli: rew > neu';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: neutral > reward ------- %
cd(paths.save2nd); mkdir('stim_neu_v_rew'); cd stim_neu_v_rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_stim_neu_v_rew) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stimuli: neu > rew';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('fb_rew_v_neu'); cd fb_rew_v_neu; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_fb_rew_v_neu) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Feedback: rew > neu';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: sanity check ------- %
cd(paths.save2nd); mkdir('fb_neu_v_rew'); cd fb_neu_v_rew; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_fb_neu_v_rew) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Feedback: neu > rew';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: reward > fix ------- %

cd(paths.save2nd); mkdir('reward_fb_v_stim'); cd reward_fb_v_stim; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_reward_fb_v_stim) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'reward: feedback > stim';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: neutral > fix ------- %
cd(paths.save2nd); mkdir('reward_stim_v_fb'); cd reward_stim_v_fb; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_reward_stim_v_fb) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'reward: stim > feedback';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: feedback: rewards > neutral ------- %
cd(paths.save2nd); mkdir('neutral_fb_v_stim'); cd neutral_fb_v_stim; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_neutral_fb_v_stim) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'neutral: feedback > stim';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: feedback: neutral > rewards ------- %
cd(paths.save2nd); mkdir('neutral_stim_v_fb'); cd neutral_stim_v_fb; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_neutral_stim_v_fb) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'neutral: stim > feedback';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: feedback: rewards > null ------- %
cd(paths.save2nd); mkdir('resp_R'); cd resp_R; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_resp_R) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'response Right';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: feedback: rewards > null ------- %
cd(paths.save2nd); mkdir('stim_rew_v_null'); cd stim_rew_v_null; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_stim_rew_v_null) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'StimReward > Null';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: feedback: rewards > null ------- %
cd(paths.save2nd); mkdir('stim_neu_v_null'); cd stim_neu_v_null; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_stim_neu_v_null) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'StimNeutral > Null';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: feedback: rewards > null ------- %
cd(paths.save2nd); mkdir('fb_rew_v_null'); cd fb_rew_v_null; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_fb_rew_v_null) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'FeedbackReward > Null';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %




% ------- compute: feedback: rewards > null ------- %
cd(paths.save2nd); mkdir('fb_neu_v_null'); cd fb_neu_v_null; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_c_fb_neu_v_null) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'FeedbackNeutral > Null';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



fprintf('\ndone\n')


%% high and low reward conditions

% for paired t-test, we measure the high/low DA conditions in the same
% subject
IDs = [4001 4002 4004 4005 4008 4010 4012 4013 4015 4016 4017];
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; 

mkdir('/Volumes/ALEX3/MRPET/analysis/2nd/PairedTtest_Feedback_Rew_v_Neu/')

clear matlabbatch
spm_jobman('initcfg')
matlabbatch{1}.spm.stats.factorial_design.dir = {'/Volumes/ALEX3/MRPET/analysis/2nd/PairedTtest_Feedback_Rew_v_Neu/'};
for id=1:length(IDs)    
matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(id).scans = {
                                                                  [paths.analyses num2str(IDs(id)) '1/con_0004_mni.nii,1']
                                                                  [paths.analyses num2str(IDs(id)) '2/con_0004_mni.nii,1']
                                                                  };
end
matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/Volumes/ALEX3/MRPET/analysis/MNI_brainOnlyMask.nii,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',matlabbatch)
fMRI_estimate(['/Volumes/ALEX3/MRPET/analysis/2nd/PairedTtest_Feedback_Rew_v_Neu/SPM.mat'])

%% for two-sample t-test, we treat high/low DA conditions as two separate
% groups

mkdir('/Volumes/ALEX3/MRPET/analysis/2nd/TwoSampleT_Feedback_Rew_v_Neu/')


clear matlabbatch
spm_jobman('initcfg')
matlabbatch{1}.spm.stats.factorial_design.dir = {'/Volumes/ALEX3/MRPET/analysis/2nd/TwoSampleT_Feedback_Rew_v_Neu/'};
listHighDA=[];
for id=1:length(IDhighDA)
    listHighDA{id,1}=['/Volumes/ALEX3/MRPET/analysis/' IDhighDA{id} '/con_0004_mni.nii,1'];
end
listLowDA=[];
for id=1:length(IDlowDA)
    listLowDA{id,1}=['/Volumes/ALEX3/MRPET/analysis/' IDlowDA{id} '/con_0004_mni.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = listHighDA;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = listLowDA;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;

%%% no cov %%%
% matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

%%% yes cov %%%

% make cov matrix
concatenatedIDs_num=str2num(cell2mat(cellfun( @ (x) x(1:4),[IDhighDA;IDlowDA], 'UniformOutput',false)));
IDhighDA_num=str2num(cell2mat(cellfun( @ (x) x(1:4),IDhighDA, 'UniformOutput',false)));
IDlowDA_num=str2num(cell2mat(cellfun( @ (x) x(1:4),IDlowDA, 'UniformOutput',false)));
covmat = IDs==concatenatedIDs_num;
%%%%%%

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov.files = {'/Volumes/ALEX3/MRPET/analysis/2nd/Cov_Subjects.mat'};
matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run',matlabbatch)
fMRI_estimate(['/Volumes/ALEX3/MRPET/analysis/2nd/TwoSampleT_Feedback_Rew_v_Neu/SPM.mat'])

