%% extract T-values from the con_00##.nii

%   ------ both sessions
%   consess{1}.tcon.name = 'stim_highDA_reward_vs_neutral';
%   consess{2}.tcon.name = 'stim_lowDA_reward_vs_neutral';
%   consess{3}.tcon.name = 'stim_rewardAll_vs_neutralAll';
%   consess{4}.tcon.name = 'fb_highDA_reward_vs_neutral';
%   consess{5}.tcon.name = 'fb_lowDA_reward_vs_neutral';
%   consess{6}.tcon.name = 'fb_rewardAll_vs_neutralAll';
%   consess{7}.tcon.name = 'stim_Rewards_highDA_vs_lowDA';
%   consess{8}.tcon.name = 'fb_Rewards_highDA_vs_lowDA';
%   consess{9}.tcon.name = 'stim_Neutrals_highDA_vs_lowDA';
%   consess{10}.tcon.name = 'fb_Neutrals_highDA_vs_lowDA';
%   consess{11}.tcon.name = 'stim_highDA_vs_lowDA';;
%   consess{12}.tcon.name = 'fb_highDA_vs_lowDA';
%   consess{13}.tcon.name = 'Stimuli_vs_Feedbacks';
%   consess{14}.tcon.name = 'Feedbacks_vs_Stimuli';
%   consess{15}.tcon.name = 'Nulls';
%   consess{16}.tcon.name = 'Button';

%   ------ each session: task
%   consess{1}.tcon.name = 'Stim-Null';
%   consess{2}.tcon.name = 'rew > neu';
%   consess{3}.tcon.name = 'neu > rew';
%   consess{4}.tcon.name = 'manipulation v intercept';
%   consess{5}.tcon.name = 'response v intercept';
%   consess{6}.tcon.name = 'reward v null';
%   consess{7}.tcon.name = 'neutral v null';
%   consess{8}.tcon.name = 'feedbacks: rew > neu';
%   consess{9}.tcon.name = 'feedbacks: neu > rew';
%   consess{10}.tcon.name = 'feedbacks: rew > null';
%   consess{11}.tcon.name = 'feedbacks: neu > null';
%   consess{12}.tcon.name = 'reward: feedback > stim';
%   consess{13}.tcon.name = 'neutral: feedback > stim';
%   consess{14}.tcon.name = 'reward: stim > feedback';
%   consess{15}.tcon.name = 'neutral: stim > feedback';
%   consess{16}.tcon.name = 'Stim_v_FB';
%   consess{17}.tcon.name = 'FB_v_Stim';

%   ----- each session: memory
%   consess{1}.tcon.name = 'stim_Remembered_vs_Forgotten';
%   consess{2}.tcon.name = 'stim_Forgotten_vs_Remembered';
%   consess{3}.tcon.name = 'stim_Rew_Remembered_Forgotten';
%   consess{4}.tcon.name = 'stim_Rew_Forgotten_Remembered';
%   consess{5}.tcon.name = 'stim_Neu_Remembered_Forgotten';
%   consess{6}.tcon.name = 'stim_Neu_Forgotten_Remembered';
%   consess{7}.tcon.name = 'stim_Rem_Reward_Neutral';
%   consess{8}.tcon.name = 'stim_Rem_Neutral_Reward';
%   consess{9}.tcon.name = 'stim_Forg_Reward_Neutral';
%   consess{10}.tcon.name = 'stim_Forg_Neutral_Reward';
%   consess{11}.tcon.name = 'Nulls';
%   consess{12}.tcon.name = 'Button';
%   consess{13}.tcon.name = 'Scenes_vs_Null';
%   consess{14}.tcon.name = 'stim: reward remembered';
%   consess{15}.tcon.name = 'stim: neutral remembered';
%   consess{16}.tcon.name = 'stim: reward forgotten';
%   consess{17}.tcon.name = 'stim: neutral forgotten';

clc;clear;

path_1st_both    = '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Bothsessions/';
path_1st_each1    = '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Eachsession_additional/';
path_1st_each2    = '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Eachsession/';
path_1st_each_mem = '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Eachsessions_memory/';

path_mnimask     = '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/masks/fMRI_clustermasks/';

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 ...
    4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2;...
    1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

%% pipeline: full factorial anova

% 2way, context * memory: 
% (1) bilateral SN

clear tmp_masklist masklist
tmp_masklist = dir([path_mnimask 'FullFactANOVA/PosInteraction_2way/*.nii']);
masklist = {tmp_masklist.name};

for m=1:length(masklist)

    clear mask_img mask_ind
    mask_img        = spm_read_vols(spm_vol([ path_mnimask 'FullFactANOVA/PosInteraction_2way/' masklist{m}]));
    mask_mni_ind    = mask_img~=0;

    for id=1:length(IDs)
        for d=1:2
            if days(id,d)==1
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each_mem num2str(IDs(id)) '_1/con_0001_mni.nii']));
                eval(['Tval_' masklist{m}(16:end-4) '_HighRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])

            elseif days(id,d)==2
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each_mem num2str(IDs(id)) '_2/con_0001_mni.nii']));
                eval(['Tval_' masklist{m}(16:end-4) '_LowRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])
            else
                if d==1
                eval(['Tval_' masklist{m}(16:end-4) '_HighRewSess(id,1)=NaN;'])
                elseif d==2
                eval(['Tval_' masklist{m}(16:end-4) '_LowRewSess(id,1)=NaN;'])
                end
            end
        end
    end
end



% 3way, context * memory * reward: 
% (1) left HPC
% (2) left VTA
% (3) bilateral SN

clear tmp_masklist masklist
tmp_masklist = dir([path_mnimask 'FullFactANOVA/PosInteraction_3way/*.nii']);
masklist = {tmp_masklist.name};

for m=1:length(masklist)

    clear mask_img mask_ind
    mask_img        = spm_read_vols(spm_vol([ path_mnimask 'FullFactANOVA/PosInteraction_3way/' masklist{m}]));
    mask_mni_ind    = mask_img~=0;

    for id=1:length(IDs)
        for d=1:2
            if days(id,d)==1
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each_mem num2str(IDs(id)) '_1/con_0014_mni.nii']));
                eval(['Tval_' masklist{m}(16:end-4) '_HighRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])

            elseif days(id,d)==2
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each_mem num2str(IDs(id)) '_2/con_0014_mni.nii']));
                eval(['Tval_' masklist{m}(16:end-4) '_LowRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])
            else
                if d==1
                eval(['Tval_' masklist{m}(16:end-4) '_HighRewSess(id,1)=NaN;'])
                elseif d==2
                eval(['Tval_' masklist{m}(16:end-4) '_LowRewSess(id,1)=NaN;'])
                end
            end
        end
    end
end


%% pipeline: Main Task Only

% Scenes & Rewards, highDA vs LowDA: 
% (1) bilateral LC

clear tmp_masklist masklist
tmp_masklist = dir([path_mnimask 'MainTask/SceneRewards_highRew_v_lowRew_LC_bi_cluster.nii']);
masklist = {tmp_masklist.name};

for m=1:length(masklist)

    clear mask_img mask_ind
    mask_img        = spm_read_vols(spm_vol([ path_mnimask 'MainTask/' masklist{m}]));
    mask_mni_ind    = mask_img~=0;

    for id=1:length(IDs)
        for d=1:2
            if days(id,d)==1
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0002_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])

            elseif days(id,d)==2
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0002_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])
            else
                if d==1
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=NaN;'])
                elseif d==2
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=NaN;'])    
                end
            end
        end
    end
end


% Scenes & neutrals, highDA vs LowDA: 
% (1) bilateral LC
% (2) left amygdala

clear tmp_masklist masklist
tmp_masklist = dir([path_mnimask 'MainTask/SceneNeutrals_highRew_v_lowRew_*.nii']);
masklist = {tmp_masklist.name};

for m=1:length(masklist)

    clear mask_img mask_ind
    mask_img        = spm_read_vols(spm_vol([ path_mnimask 'MainTask/' masklist{m}]));
    mask_mni_ind    = mask_img~=0;

    for id=1:length(IDs)
        for d=1:2
            if days(id,d)==1
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0002_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])

            elseif days(id,d)==2
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0002_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])
            else
                if d==1
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=NaN;'])
                elseif d==2
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=NaN;'])    
                end
            end
        end
    end
end



% FB, Rew vs Neu: 
% (1) right VTA
% (2) right thalamus

clear tmp_masklist masklist
tmp_masklist = dir([path_mnimask 'MainTask/FB_RewAll_v_NeuAll_*.nii']);
masklist = {tmp_masklist.name};

for m=1:length(masklist)

    clear mask_img mask_ind
    mask_img        = spm_read_vols(spm_vol([ path_mnimask 'MainTask/' masklist{m}]));
    mask_mni_ind    = mask_img~=0;

    for id=1:length(IDs)
        for d=1:2
            if days(id,d)==1
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0008_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])

            elseif days(id,d)==2
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0008_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])
            else
                if d==1
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=NaN;'])
                elseif d==2
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=NaN;'])    
                end
            end
        end
    end
end


% FB, highRew: Rew vs Neu
% (1) bilateral VTA
% (2) bilateral LC
% (3) right thalamus

clear tmp_masklist masklist
tmp_masklist = dir([path_mnimask 'MainTask/FB_HighRew_Rew_v_Neu_*.nii']);
masklist = {tmp_masklist.name};

for m=1:length(masklist)

    clear mask_img mask_ind
    mask_img        = spm_read_vols(spm_vol([ path_mnimask 'MainTask/' masklist{m}]));
    mask_mni_ind    = mask_img~=0;

    for id=1:length(IDs)
        for d=1:2
            if days(id,d)==1
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0008_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])

            elseif days(id,d)==2
                clear tmp_img
                tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0008_mni.nii']));
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=nanmean(tmp_img(mask_mni_ind));'])
            else
                if d==1
                eval(['Tval_' masklist{m}(1:end-4) '_HighRewSess(id,1)=NaN;'])
                elseif d==2
                eval(['Tval_' masklist{m}(1:end-4) '_LowRewSess(id,1)=NaN;'])    
                end
            end
        end
    end
end

%% now make a table

vars = who('Tval_*');
T = table();
for i = 1:length(vars)
    varName = vars{i}; 
    varData = eval(varName);
    T.(varName) = varData; 
end
writetable(T, ['/Users/alex/Dropbox/paperwriting/MRPET/data/Tvals_fMRIcluters_' date '.xls']);

%% both session table

% Step 1: Identify variables
highRewVars = who('Tval_*_HighRewSess');
lowRewVars = who('Tval_*_LowRewSess');

% Step 2: Loop through HighRewSess and find corresponding LowRewSess variables
for i = 1:length(highRewVars)
    % Extract the variable name for HighRewSess
    highRewVarName = highRewVars{i};
    
    % Derive the corresponding LowRewSess variable name
    baseName = strrep(highRewVarName, '_HighRewSess', ''); % Get wildcard part
    lowRewVarName = [baseName, '_LowRewSess'];
    
    % Check if the LowRewSess variable exists
    if ismember(lowRewVarName, lowRewVars)
        % Step 3: Concatenate vertically
        highRewData = eval(highRewVarName);
        lowRewData = eval(lowRewVarName);
        bothSessData = [highRewData; lowRewData]; % Concatenate vertically
        
        % Create a new variable name with '_BothSess'
        newVarName = [baseName, '_BothSess'];
        
        % Assign the concatenated data to the new variable
        assignin('base', newVarName, bothSessData);
    else
        % If the corresponding LowRewSess variable doesn't exist, print a warning
        warning('LowRewSess variable %s not found for %s', lowRewVarName, highRewVarName);
    end
end

%%
varsboth = who('Tval_*_BothSess');
Tboth = table();
for i = 1:length(varsboth)
    varNameboth = varsboth{i}; 
    varDataboth = eval(varNameboth);
    Tboth.(varNameboth) = varDataboth; 
end
writetable(Tboth, ['/Users/alex/Dropbox/paperwriting/MRPET/data/Tvals_fMRIcluters_' date '_BothSession.xls']);


%% extract: both sessions

fMRIclusters1_mask = spm_read_vols(spm_vol([ path_mnimask 'mni_icbm152_Stim_vs_FB_p001_c255_HPConly.nii']));
fMRIclusters1_ind  = fMRIclusters1_mask~=0;

Tval_fMRIclusters1=[];
for id=1:length(IDs)

    if sum(days(id,:))==3
        clear tmp_img
        tmp_img = spm_read_vols(spm_vol([ path_1st_both num2str(IDs(id)) '/con_0013_mni.nii']));
        Tval_fMRIclusters1(id,1)=nanmean(tmp_img(fMRIclusters1_ind));

    else
        disp('skipped')
        Tval_fMRIclusters1(id,1)= NaN;
        
    end
end

%% extract: each session (STIM vs FB)

% but used both session mask
fMRIclusters1_mask = spm_read_vols(spm_vol([ path_mnimask 'mni_icbm152_Stim_vs_FB_p001_c255_HPConly.nii']));
fMRIclusters1_ind  = fMRIclusters1_mask~=0;

Tval_fMRIclusters1=[]; Tval_fMRIclusters1_highDA=[]; Tval_fMRIclusters1_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each1 num2str(IDs(id)) '_1/con_0016_mni.nii']));
            Tval_fMRIclusters1(id,d)=nanmean(tmp_img(fMRIclusters1_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each1 num2str(IDs(id)) '_2/con_0016_mni.nii']));
            Tval_fMRIclusters1(id,d)=nanmean(tmp_img(fMRIclusters1_ind));
        else
            Tval_fMRIclusters1(id,d)=NaN;
        end
    end
end
Tval_fMRIclusters1_highDA=Tval_fMRIclusters1(:,1);
Tval_fMRIclusters1_lowDA=Tval_fMRIclusters1(:,2);


%% extract: each session (STIM: Rew vs Neu)

fMRIclusters2_mask = spm_read_vols(spm_vol([ path_mnimask 'EachSession_STIM_Rew_v_Neu_p005_c370_HPConly.nii']));
fMRIclusters2_ind  = fMRIclusters2_mask~=0;

Tval_fMRIclusters2=[]; Tval_fMRIclusters2_highDA=[]; Tval_fMRIclusters2_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0002_mni.nii']));
            Tval_fMRIclusters2(id,d)=nanmean(tmp_img(fMRIclusters2_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0002_mni.nii']));
            Tval_fMRIclusters2(id,d)=nanmean(tmp_img(fMRIclusters2_ind));
        else
            Tval_fMRIclusters2(id,d)=NaN;
        end
    end
end
Tval_fMRIclusters2_highDA=Tval_fMRIclusters2(:,1);
Tval_fMRIclusters2_lowDA=Tval_fMRIclusters2(:,2);

%% each session, hippocampaus mask (1)

fMRIclusters1_mask = spm_read_vols(spm_vol([ path_mnimask 'mni_icbm152_AAL3_hippocampus_bi.nii']));
fMRIclusters1_ind  = fMRIclusters1_mask~=0;

Tval_fMRI_stim_v_fb_HPC=[]; Tval_fMRI_stim_v_fb_HPC_highDA=[]; Tval_fMRI_stim_v_fb_HPC_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each1 num2str(IDs(id)) '_1/con_0016_mni.nii']));
            Tval_fMRI_stim_v_fb_HPC(id,d)=nanmean(tmp_img(fMRIclusters1_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each1 num2str(IDs(id)) '_2/con_0016_mni.nii']));
            Tval_fMRI_stim_v_fb_HPC(id,d)=nanmean(tmp_img(fMRIclusters1_ind));
        else
            Tval_fMRI_stim_v_fb_HPC(id,d)=NaN;
        end
    end
end
Tval_fMRI_stim_v_fb_HPC_highDA=Tval_fMRI_stim_v_fb_HPC(:,1);
Tval_fMRI_stim_v_fb_HPC_lowDA=Tval_fMRI_stim_v_fb_HPC(:,2);

%% each session, hippocampaus mask (2)

fMRIclusters1_mask = spm_read_vols(spm_vol([ path_mnimask 'mni_icbm152_AAL3_hippocampus_bi.nii']));
fMRIclusters1_ind  = fMRIclusters1_mask~=0;

Tval_fMRI_stim_NEUvsREW_HPC=[]; Tval_fMRI_stim_NEUvsREW_highDA=[]; Tval_fMRI_stim_NEUvsREW_HPC_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0003_mni.nii']));
            Tval_fMRI_stim_NEUvsREW_HPC(id,d)=nanmean(tmp_img(fMRIclusters1_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0003_mni.nii']));
            Tval_fMRI_stim_NEUvsREW_HPC(id,d)=nanmean(tmp_img(fMRIclusters1_ind));
        else
            Tval_fMRI_stim_NEUvsREW_HPC(id,d)=NaN;
        end
    end
end
Tval_fMRI_stim_NEUvsREW_HPC_highDA=Tval_fMRI_stim_NEUvsREW_HPC(:,1);
Tval_fMRI_stim_NEUvsREW_HPC_lowDA=Tval_fMRI_stim_NEUvsREW_HPC(:,2);

Tval_fMRI_fb_NEUvsREW_HPC=[]; Tval_fMRI_fb_NEUvsREW_highDA=[]; Tval_fMRI_fb_NEUvsREW_HPC_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0009_mni.nii']));
            Tval_fMRI_fb_NEUvsREW_HPC(id,d)=nanmean(tmp_img(fMRIclusters1_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0009_mni.nii']));
            Tval_fMRI_fb_NEUvsREW_HPC(id,d)=nanmean(tmp_img(fMRIclusters1_ind));
        else
            Tval_fMRI_fb_NEUvsREW_HPC(id,d)=NaN;
        end
    end
end
Tval_fMRI_fb_NEUvsREW_HPC_highDA=Tval_fMRI_fb_NEUvsREW_HPC(:,1);
Tval_fMRI_fb_NEUvsREW_HPC_lowDA=Tval_fMRI_fb_NEUvsREW_HPC(:,2);

disp('done')


%% each session, AMG cluster

fMRIclusters1_mask = spm_read_vols(spm_vol([ '/Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/2ndLevel/2nd_mni_CompleteDatasetsOnly/stim_Neutrals_highDA_vs_lowDA_bothonly/StimNeutrals_highRew_v_lowRew_p005_c510_AMG.nii']));
fMRIclusters1_ind  = fMRIclusters1_mask~=0;

Tval_fMRI_stim_NEUvsNULL_AMGc=[]; Tval_fMRI_stim_NEUvsNULL_highDA=[]; Tval_fMRI_stim_NEUvsNULL_AMGc_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0007_mni.nii']));
            Tval_fMRI_stim_NEUvsNULL_AMGc(id,d)=nanmean(tmp_img(fMRIclusters1_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0007_mni.nii']));
            Tval_fMRI_stim_NEUvsNULL_AMGc(id,d)=nanmean(tmp_img(fMRIclusters1_ind));
        else
            Tval_fMRI_stim_NEUvsNULL_AMGc(id,d)=NaN;
        end
    end
end
Tval_fMRI_stim_NEUvsNULL_AMGc_highDA=Tval_fMRI_stim_NEUvsNULL_AMGc(:,1);
Tval_fMRI_stim_NEUvsNULL_AMGc_lowDA=Tval_fMRI_stim_NEUvsNULL_AMGc(:,2);

Tval_fMRI_fb_NEUvsNULL_AMGc=[]; Tval_fMRI_fb_NEUvsNULL_highDA=[]; Tval_fMRI_fb_NEUvsNULL_AMGc_lowDA=[];
for id=1:length(IDs)

    for d=1:2
        if days(id,d)==1
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_1/con_0011_mni.nii']));
            Tval_fMRI_fb_NEUvsNULL_AMGc(id,d)=nanmean(tmp_img(fMRIclusters1_ind));

        elseif days(id,d)==2
    clear tmp_img
            tmp_img = spm_read_vols(spm_vol([ path_1st_each2 num2str(IDs(id)) '_2/con_0011_mni.nii']));
            Tval_fMRI_fb_NEUvsNULL_AMGc(id,d)=nanmean(tmp_img(fMRIclusters1_ind));
        else
            Tval_fMRI_fb_NEUvsNULL_AMGc(id,d)=NaN;
        end
    end
end
Tval_fMRI_fb_NEUvsNULL_AMGc_highDA=Tval_fMRI_fb_NEUvsNULL_AMGc(:,1);
Tval_fMRI_fb_NEUvsNULL_AMGc_lowDA=Tval_fMRI_fb_NEUvsNULL_AMGc(:,2);


disp('done')

