%% PET 1st-level analysis pipeline

%% work log

%   12-12-2019    created the script

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Volumes/ALEX3/MRPET/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
paths.preproc = [paths.parent 'img/'];
paths.analyses= [paths.parent 'analysis/bothSessions/1st/'];
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
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2]; 


% IDs = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4026];
% days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0; 0 2]; 
% d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 0 0; 0 0; 1 2; 0 0]; % 1=immediate 2=delayed
% d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 2; 1 2;1 2; 0 0; 1 2; 1 2; 1 2; 0 0; 0 0; 1 2];

setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Applications/freesurfer/7.2.0/bin:/Users/ikndadmin/abin:/usr/local/fsl']);
setenv('ANTSPATH','/usr/local/bin')

% FSL Setup
setenv( 'FSLDIR', '/usr/local/fsl' );
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
            fname_beh{id,d} = {NaN};
            expdat{id,d} = {NaN};
            contingency{id,d} = [];
        else
            
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
            tmptrigs  = eval(['(expdat{id,d}.dat.day' num2str(d) '.maintask.results.SOT.raw.trig_1st)']);
            tmptrials = cell2mat(eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,2)']));
            ind_null  = isnan(tmptrials);
            ind_rew   = tmptrials==contingency{id,d}(1); ind_rew(ind_null)=[];
            ind_neu   = tmptrials==contingency{id,d}(2); ind_neu(ind_null)=[];
            
            
            % sort onsets
            onsets_temp{id,d}.stim.reward  = tmponsets.stim(ind_rew)-tmptrigs;
            onsets_temp{id,d}.stim.neutral = tmponsets.stim(ind_neu)-tmptrigs;
            onsets_temp{id,d}.null         = tmponsets.null-tmptrigs;
            onsets_temp{id,d}.resp         = tmponsets.resp-tmptrigs;
            onsets_temp{id,d}.fb.reward    = tmponsets.cue(ind_rew)-tmptrigs;
            onsets_temp{id,d}.fb.neutral   = tmponsets.cue(ind_neu)-tmptrigs;

        end 
    end

    % sort onsets
    
    if sum(days(id,:))==3
    
    clear numvolSess1 
    numvolSess1  = length(spm_vol([paths.preproc num2str(IDs(id)) '_1/' num2str(IDs(id)) '_MRI_4D_MT1.nii']));
    length_scan1(id,1) = numvolSess1;

    onsets{id,1}.stim.highDA_reward  = onsets_temp{id,1}.stim.reward;
    onsets{id,1}.stim.lowDA_reward   = onsets_temp{id,2}.stim.reward+numvolSess1*TR;
    onsets{id,1}.stim.highDA_neutral = onsets_temp{id,1}.stim.neutral;
    onsets{id,1}.stim.lowDA_neutral  = onsets_temp{id,2}.stim.neutral+numvolSess1*TR;
    onsets{id,1}.stim.all_reward     = [onsets_temp{id,1}.stim.reward; onsets_temp{id,2}.stim.reward+numvolSess1*TR];
    onsets{id,1}.stim.all_neutral    = [onsets_temp{id,1}.stim.neutral; onsets_temp{id,2}.stim.neutral+numvolSess1*TR];

    onsets{id,1}.fb.highDA_reward    = onsets_temp{id,1}.fb.reward;
    onsets{id,1}.fb.lowDA_reward     = onsets_temp{id,2}.fb.reward+numvolSess1*TR;
    onsets{id,1}.fb.highDA_neutral   = onsets_temp{id,1}.fb.neutral;
    onsets{id,1}.fb.lowDA_neutral    = onsets_temp{id,2}.fb.neutral+numvolSess1*TR;
    onsets{id,1}.fb.all_reward       = [onsets_temp{id,1}.fb.reward; onsets_temp{id,2}.fb.reward+numvolSess1*TR];
    onsets{id,1}.fb.all_neutral      = [onsets_temp{id,1}.fb.neutral; onsets_temp{id,2}.fb.neutral+numvolSess1*TR];

    onsets{id,1}.highDA_null         = onsets_temp{id,1}.null;
    onsets{id,1}.highDA_resp         = onsets_temp{id,1}.resp;
    onsets{id,1}.lowDA_null          = onsets_temp{id,2}.null+numvolSess1*TR;
    onsets{id,1}.lowDA_resp          = onsets_temp{id,2}.resp+numvolSess1*TR;
    onsets{id,1}.all_null            = [onsets_temp{id,1}.null; onsets_temp{id,2}.null+numvolSess1*TR];
    onsets{id,1}.all_resp            = [onsets_temp{id,1}.resp; onsets_temp{id,2}.resp+numvolSess1*TR];
    
    
    elseif sum(days(id,:))==1
    clear numvolSess1 
    numvolSess1  = length(spm_vol([paths.preproc num2str(IDs(id)) '_1/' num2str(IDs(id)) '_MRI_4D_MT1.nii']));
    length_scan1(id,1) = numvolSess1;

    onsets{id,1}.stim.highDA_reward  = onsets_temp{id,1}.stim.reward;
    onsets{id,1}.stim.lowDA_reward   = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.stim.highDA_neutral = onsets_temp{id,1}.stim.neutral;
    onsets{id,1}.stim.lowDA_neutral  = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.stim.all_reward     = [onsets_temp{id,1}.stim.reward];
    onsets{id,1}.stim.all_neutral    = [onsets_temp{id,1}.stim.neutral];

    onsets{id,1}.fb.highDA_reward    = onsets_temp{id,1}.fb.reward;
    onsets{id,1}.fb.lowDA_reward     = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.fb.highDA_neutral   = onsets_temp{id,1}.fb.neutral;
    onsets{id,1}.fb.lowDA_neutral    = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.fb.all_reward       = [onsets_temp{id,1}.fb.reward];
    onsets{id,1}.fb.all_neutral      = [onsets_temp{id,1}.fb.neutral];

    onsets{id,1}.highDA_null         = onsets_temp{id,1}.null;
    onsets{id,1}.highDA_resp         = onsets_temp{id,1}.resp;
    onsets{id,1}.lowDA_null          = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.lowDA_resp          = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.all_null            = [onsets_temp{id,1}.null];
    onsets{id,1}.all_resp            = [onsets_temp{id,1}.resp];
    
    elseif sum(days(id,:))==2
    clear numvolSess1 
    numvolSess1  = length(spm_vol([paths.preproc num2str(IDs(id)) '_2/' num2str(IDs(id)) '_MRI_4D_MT2.nii']));
    length_scan1(id,1) = numvolSess1;

    onsets{id,1}.stim.lowDA_reward  = onsets_temp{id,2}.stim.reward;
    onsets{id,1}.stim.highDA_reward   = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.stim.lowDA_neutral = onsets_temp{id,2}.stim.neutral;
    onsets{id,1}.stim.highDA_neutral  = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.stim.all_reward     = [onsets_temp{id,2}.stim.reward];
    onsets{id,1}.stim.all_neutral    = [onsets_temp{id,2}.stim.neutral];

    onsets{id,1}.fb.lowDA_reward    = onsets_temp{id,2}.fb.reward;
    onsets{id,1}.fb.highDA_reward     = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.fb.lowDA_neutral   = onsets_temp{id,2}.fb.neutral;
    onsets{id,1}.fb.highDA_neutral    = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.fb.all_reward       = [onsets_temp{id,2}.fb.reward];
    onsets{id,1}.fb.all_neutral      = [onsets_temp{id,2}.fb.neutral];

    onsets{id,1}.lowDA_null         = onsets_temp{id,2}.null;
    onsets{id,1}.lowDA_resp         = onsets_temp{id,2}.resp;
    onsets{id,1}.highDA_null          = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.highDA_resp          = numvolSess1*TR+TR; % null regressor
    onsets{id,1}.all_null            = [onsets_temp{id,2}.null];
    onsets{id,1}.all_resp            = [onsets_temp{id,2}.resp];
    
    end

    % make EPI mask
%     eval(['!bet2 ' paths.preproc num2str(IDs(id)) '/meanall_ua' num2str(IDs(id)) '.nii ' paths.preproc num2str(IDs(id)) '/EPI_brainonly -m'])
%     gunzip([paths.preproc num2str(IDs(id)) '/EPI_brainonly_mask.nii.gz'])
%     copyfile([paths.preproc num2str(IDs(id)) '/EPI_brainonly_mask.nii'],[paths.preproc num2str(IDs(id)) '/mask.nii'])
%     delete([paths.preproc num2str(IDs(id)) '/EPI_brainonly_mask.nii.gz']); delete([paths.preproc num2str(IDs(id)) '/EPI_brainonly.nii.gz']); 

end

%% make physio+realignment+session regressors (already done)

for id = 1:length(IDs)
    
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

for id = 1%:length(IDs)

    fprintf('\n model estimation for: ID %d \n', IDs(id))   
    
    % set up workspace
    cd(paths.analyses); mkdir([num2str(IDs(id))]);
    clear dir_model
    dir_model = [paths.analyses num2str(IDs(id)) '/'];
    dir_multiparam = [paths.preproc num2str(IDs(id)) '/reg_all.mat'];

    %% build model

    clear matlabbatch

    spm_jobman('initcfg');      % initiate job manager

    matlabbatch{1}.spm.stats.fmri_spec.dir               = cellstr(dir_model); % where will the model be saved?
    matlabbatch{1}.spm.stats.fmri_spec.timing.units      = 'secs'; % scans / secs
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT         = TR; % TR in seconds
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t     = slices; % volume size
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0    = microtime; % microtime onset

    matlabbatch{1}.spm.stats.fmri_spec.sess.scans        = cellstr([paths.preproc num2str(IDs(id)) '/' ...
        'srall_ua' num2str(IDs(id)) '.nii']);


    % reward condition
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Stim_highDA_Rewards';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset= onsets{id,1}.stim.highDA_reward;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Stim_highDA_Neutrals';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.stim.highDA_neutral;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Stim_lowDA_Rewards';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset= onsets{id,1}.stim.lowDA_reward;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Stim_lowDA_Neutrals';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.stim.lowDA_neutral;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'FB_highDA_Rewards';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset= onsets{id,1}.fb.highDA_reward;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'FB_highDA_Neutrals';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset= onsets{id,1}.fb.highDA_neutral;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name = 'FB_lowDA_Rewards';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset= onsets{id,1}.fb.lowDA_reward;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name = 'FB_lowDA_Neutrals';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset= onsets{id,1}.fb.lowDA_neutral;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name = 'NullTrials';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset= onsets{id,1}.all_null;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).name = 'ButtonPresses';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).onset= onsets{id,1}.all_resp;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(10).orth = 1;

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

    % here you set up the contrast matrices: Rewards - Neutrals - Fix
    stim_highDA_reward_vs_neutral  = [1 -1];
    stim_lowDA_reward_vs_neutral   = [0 0 1 -1];
    stim_rewardAll_vs_neutralAll   = [1 -1 1 -1];

    fb_highDA_reward_vs_neutral    = [0 0 0 0 1 -1];
    fb_lowDA_reward_vs_neutral     = [0 0 0 0 0 0 1 -1];
    fb_rewardAll_vs_neutralAll     = [0 0 0 0 1 -1 1 -1];

    stim_Rewards_highDA_vs_lowDA   = [1 0 -1];
    fb_Rewards_highDA_vs_lowDA     = [0 0 0 0 1 0 -1];
    stim_Neutrals_highDA_vs_lowDA  = [0 1 0 -1];
    fb_Neutrals_highDA_vs_lowDA    = [0 0 0 0 0 1 0 -1];

    stim_highDA_vs_lowDA           = [1 1 -1 -1]; 
    fb_highDA_vs_lowDA             = [0 0 0 0 1 1 -1 -1];

    Stimuli_vs_Feedbacks           = [1 1 1 1 -1 -1 -1 -1];   
    Feedbacks_vs_Stimuli           = [-1 -1 -1 -1 1 1 1 1];

    Nulls  = [0 0 0 0 0 0 0 0 1]; 
    Button = [0 0 0 0 0 0 0 0 0 1]; 

    % batch setup
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager

    matlabbatch{1}.spm.stats.con.spmmat = cellstr(firstlvl_dir);

    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highDA_reward_vs_neutral';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = stim_highDA_reward_vs_neutral;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'stim_lowDA_reward_vs_neutral';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = stim_lowDA_reward_vs_neutral;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'stim_rewardAll_vs_neutralAll';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = stim_rewardAll_vs_neutralAll;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'fb_highDA_reward_vs_neutral';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = fb_highDA_reward_vs_neutral;
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'fb_lowDA_reward_vs_neutral';
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = fb_lowDA_reward_vs_neutral;
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'fb_rewardAll_vs_neutralAll';
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = fb_rewardAll_vs_neutralAll;
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'stim_Rewards_highDA_vs_lowDA';
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = stim_Rewards_highDA_vs_lowDA;
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'fb_Rewards_highDA_vs_lowDA';
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = fb_Rewards_highDA_vs_lowDA;
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'stim_Neutrals_highDA_vs_lowDA';
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = stim_Neutrals_highDA_vs_lowDA;
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'fb_Neutrals_highDA_vs_lowDA';
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = fb_Neutrals_highDA_vs_lowDA;
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'stim_highDA_vs_lowDA';
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = stim_highDA_vs_lowDA;
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'fb_highDA_vs_lowDA';
    matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = fb_highDA_vs_lowDA;
    matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'Stimuli_vs_Feedbacks';
    matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = Stimuli_vs_Feedbacks;
    matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'Feedbacks_vs_Stimuli';
    matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = Feedbacks_vs_Stimuli;
    matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{15}.tcon.name = 'Nulls';
    matlabbatch{1}.spm.stats.con.consess{15}.tcon.weights = Nulls;
    matlabbatch{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.consess{16}.tcon.name = 'Button';
    matlabbatch{1}.spm.stats.con.consess{16}.tcon.weights = Button;
    matlabbatch{1}.spm.stats.con.consess{16}.tcon.sessrep = 'none';

    %     matlabbatch{1}.spm.stats.con.consess{6}.fcon.name = 'ANOVA (rew vs neu vs fix)';
    %     matlabbatch{1}.spm.stats.con.consess{6}.fcon.weights = c_F_1;
    %     matlabbatch{1}.spm.stats.con.consess{6}.fcon.sessrep = 'none';

    matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append

    spm_jobman('run', matlabbatch) % run batch

end

fprintf('\ndone\n')

%% move contrasts for registration

path_source      =[paths.analyses];
path_destination ='/Volumes/ALEX3/MRPET/coreg_mri/bothsessions/';

mrpetID=[]; c1=0;
for id=1:length(IDs)

            c1=c1+1;
            mrpetID{c1,1}=[num2str(IDs(id))];
                        
            for cnt=1:16
            copyfile([path_source num2str(IDs(id)) '/con_00' sprintf('%02i',cnt) '.nii'],...
                [path_destination num2str(IDs(id)) '/data/con_00' sprintf('%02i',cnt) '.nii'])
            end
            copyfile([path_source num2str(IDs(id)) '/mask.nii'],...
                [path_destination num2str(IDs(id)) '/data/mask.nii'])
            
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
paths.parent  = '/Volumes/ALEX3/MRPET/analysis/bothSessions/';
% paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
% paths.funx    = [paths.parent 'B_scripts/BB_analyses/BBC_MRI/analyses_functions/'];
% paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.analyses= [paths.parent '1st/'];
paths.analyses_half = '/Volumes/ALEX3/MRPET/analysis/eachSession/1st/';
paths.behav   = ['/Volumes/ALEX3/MRPET/behav/'];
% paths.temp    = [paths.parent 'E_data/EB_cleaned/EBD_mrpet/RewardTask/MRI/'];
paths.save2nd = [paths.parent '2nd_mni/'];


% add toolboxes and functions
% addpath(paths.spm)
% addpath(paths.funx)

% IDs
% IDs  = [4001 4002 4004 4005 4008 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4028 4030 4031 4032 4033];
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; 

IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2]; 

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

list_stim_highDA_reward_vs_neutral = []; cnt=0;
for c1 = 1:length(IDs)
    if sum(days(c1,:))==3
        cnt=cnt+1;
    list_stim_highDA_reward_vs_neutral{cnt,1} = [paths.analyses num2str(IDs(c1)) '/con_0001_mni.nii,1'];
    elseif sum(days(c1,:))==1
        cnt=cnt+1;
    list_stim_highDA_reward_vs_neutral{cnt,1} = [paths.analyses_half num2str(IDs(c1)) '_1/con_0002_mni.nii,1'];
    else
    disp('skipped')
    end
end

list_stim_lowDA_reward_vs_neutral = []; cnt=0;
for c1 = 1:length(IDs)
    if sum(days(c1,:))==3
        cnt=cnt+1;
    list_stim_lowDA_reward_vs_neutral{cnt,1} = [paths.analyses num2str(IDs(c1)) '/con_0002_mni.nii,1'];
    elseif sum(days(c1,:))==2
        cnt=cnt+1;
    list_stim_lowDA_reward_vs_neutral{cnt,1} = [paths.analyses_half num2str(IDs(c1)) '_2/con_0002_mni.nii,1'];
    else
    disp('skipped')
    end
end

list_stim_rewardAll_vs_neutralAll = [];
for c1 = 1:length(IDs)
    if sum(days(c1,:))==3
    list_stim_rewardAll_vs_neutralAll{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0003_mni.nii,1'];
    elseif sum(days(c1,:))==1
    list_stim_rewardAll_vs_neutralAll{c1,1} = [paths.analyses_half num2str(IDs(c1)) '_1/con_0002_mni.nii,1'];
    elseif sum(days(c1,:))==2
    list_stim_rewardAll_vs_neutralAll{c1,1} = [paths.analyses_half num2str(IDs(c1)) '_2/con_0002_mni.nii,1'];
    end
end

list_fb_highDA_reward_vs_neutral = []; cnt=0;
for c1 = 1:length(IDs)
    if sum(days(c1,:))==3
        cnt=cnt+1;
    list_fb_highDA_reward_vs_neutral{cnt,1} = [paths.analyses num2str(IDs(c1)) '/con_0004_mni.nii,1'];
    elseif sum(days(c1,:))==2
        cnt=cnt+1;
    list_fb_highDA_reward_vs_neutral{cnt,1} = [paths.analyses num2str(IDs(c1)) '_1/con_0008_mni.nii,1'];
    else
    disp('skipped')
    end
end

list_fb_lowDA_reward_vs_neutral = []; cnt=0;
for c1 = 1:length(IDs)
    if sum(days(c1,:))==3
        cnt=cnt+1;
    list_fb_lowDA_reward_vs_neutral{cnt,1} = [paths.analyses num2str(IDs(c1)) '/con_0005_mni.nii,1'];
    elseif sum(days(c1,:))==2
        cnt=cnt+1;
    list_fb_lowDA_reward_vs_neutral{cnt,1} = [paths.analyses num2str(IDs(c1)) '_2/con_0008_mni.nii,1'];
    else
    disp('skipped')
    end
end

list_fb_rewardAll_vs_neutralAll = [];
for c1 = 1:length(IDs)
    if sum(days(c1,:))==3
    list_fb_rewardAll_vs_neutralAll{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0006_mni.nii,1'];
    elseif sum(days(c1,:))==1
    list_fb_rewardAll_vs_neutralAll{c1,1} = [paths.analyses num2str(IDs(c1)) '_1/con_0008_mni.nii,1'];
    elseif sum(days(c1,:))==2
    list_fb_rewardAll_vs_neutralAll{c1,1} = [paths.analyses num2str(IDs(c1)) '_2/con_0008_mni.nii,1'];
    end
end

list_stim_Rewards_highDA_vs_lowDA = [];
for c1 = 1:length(IDs)
    list_stim_Rewards_highDA_vs_lowDA{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0007_mni.nii,1'];
end

list_fb_Rewards_highDA_vs_lowDA = [];
for c1 = 1:length(IDs)
    list_fb_Rewards_highDA_vs_lowDA{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0008_mni.nii,1'];
end

list_stim_Neutrals_highDA_vs_lowDA = [];
for c1 = 1:length(IDs)
    list_stim_Neutrals_highDA_vs_lowDA{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0009_mni.nii,1'];
end

list_fb_Neutrals_highDA_vs_lowDA = [];
for c1 = 1:length(IDs)
    list_fb_Neutrals_highDA_vs_lowDA{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0010_mni.nii,1'];
end

list_stim_highDA_vs_lowDA = [];
for c1 = 1:length(IDs)
    list_stim_highDA_vs_lowDA{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0011_mni.nii,1'];
end

list_fb_highDA_vs_lowDA = [];
for c1 = 1:length(IDs)
    list_fb_highDA_vs_lowDA{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0012_mni.nii,1'];
end

list_Stimuli_vs_Feedbacks = [];
for c1 = 1:length(IDs)
    list_Stimuli_vs_Feedbacks{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0013_mni.nii,1'];
end

list_Feedbacks_vs_Stimuli = [];
for c1 = 1:length(IDs)
    list_Feedbacks_vs_Stimuli{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0014_mni.nii,1'];
end

list_Nulls = [];
for c1 = 1:length(IDs)
    list_Nulls{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0015_mni.nii,1'];
end

list_Button = [];
for c1 = 1:length(IDs)
    list_Button{c1,1} = [paths.analyses num2str(IDs(c1)) '/con_0016_mni.nii,1'];
end

%% compute models

% ------- compute: stim_highDA_reward_vs_neutral ------- %

cd(paths.save2nd); mkdir('stim_highDA_reward_vs_neutral'); cd stim_highDA_reward_vs_neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highDA_reward_vs_neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highDA_reward_vs_neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_lowDA_reward_vs_neutral ------- %

cd(paths.save2nd); mkdir('stim_lowDA_reward_vs_neutral'); cd stim_lowDA_reward_vs_neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_lowDA_reward_vs_neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_lowDA_reward_vs_neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_rewardAll_vs_neutralAll ------- %
cd(paths.save2nd); mkdir('stim_rewardAll_vs_neutralAll'); cd stim_rewardAll_vs_neutralAll; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_rewardAll_vs_neutralAll) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_rewardAll_vs_neutralAll';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: fb_highDA_reward_vs_neutral ------- %
cd(paths.save2nd); mkdir('fb_highDA_reward_vs_neutral'); cd fb_highDA_reward_vs_neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_highDA_reward_vs_neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb_highDA_reward_vs_neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: fb_lowDA_reward_vs_neutral ------- %
cd(paths.save2nd); mkdir('fb_lowDA_reward_vs_neutral'); cd fb_lowDA_reward_vs_neutral; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_lowDA_reward_vs_neutral) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb_lowDA_reward_vs_neutral';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: fb_rewardAll_vs_neutralAll ------- %

cd(paths.save2nd); mkdir('fb_rewardAll_vs_neutralAll'); cd fb_rewardAll_vs_neutralAll; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_rewardAll_vs_neutralAll) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb_rewardAll_vs_neutralAll';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: stim_Rewards_highDA_vs_lowDA ------- %
cd(paths.save2nd); mkdir('stim_Rewards_highDA_vs_lowDA'); cd stim_Rewards_highDA_vs_lowDA; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_Rewards_highDA_vs_lowDA) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_Rewards_highDA_vs_lowDA';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %



% ------- compute: fb_Rewards_highDA_vs_lowDA ------- %
cd(paths.save2nd); mkdir('fb_Rewards_highDA_vs_lowDA'); cd fb_Rewards_highDA_vs_lowDA; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_Rewards_highDA_vs_lowDA) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb_Rewards_highDA_vs_lowDA';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_Neutrals_highDA_vs_lowDA ------- %
cd(paths.save2nd); mkdir('stim_Neutrals_highDA_vs_lowDA'); cd stim_Neutrals_highDA_vs_lowDA; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_Neutrals_highDA_vs_lowDA) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_Neutrals_highDA_vs_lowDA';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: fb_Neutrals_highDA_vs_lowDA ------- %
cd(paths.save2nd); mkdir('fb_Neutrals_highDA_vs_lowDA'); cd fb_Neutrals_highDA_vs_lowDA; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_Neutrals_highDA_vs_lowDA) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb_Neutrals_highDA_vs_lowDA';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: stim_highDA_vs_lowDA ------- %
cd(paths.save2nd); mkdir('stim_highDA_vs_lowDA'); cd stim_highDA_vs_lowDA; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_stim_highDA_vs_lowDA) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'stim_highDA_vs_lowDA';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: fb_highDA_vs_lowDA ------- %
cd(paths.save2nd); mkdir('fb_highDA_vs_lowDA'); cd fb_highDA_vs_lowDA; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_fb_highDA_vs_lowDA) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'fb_highDA_vs_lowDA';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: Stimuli_vs_Feedbacks ------- %
cd(paths.save2nd); mkdir('Stimuli_vs_Feedbacks'); cd Stimuli_vs_Feedbacks; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_Stimuli_vs_Feedbacks) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stimuli_vs_Feedbacks';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
spm_jobman('run', matlabbatch) % run batch

clear dir_spm tmpdir i3 doc_list_uncorr

% ------------------------------------------------ %


% ------- compute: Feedbacks_vs_Stimuli ------- %
cd(paths.save2nd); mkdir('Feedbacks_vs_Stimuli'); cd Feedbacks_vs_Stimuli; dir_spm = pwd;
secondlvl_onesampleT(dir_spm,list_Feedbacks_vs_Stimuli) % run
fMRI_estimate([dir_spm '/SPM.mat'])
clear matlabbatch
spm_jobman('initcfg');      % initiate job manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr([dir_spm '/SPM.mat']);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Feedbacks_vs_Stimuli';
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

fprintf('\ndone\n')


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