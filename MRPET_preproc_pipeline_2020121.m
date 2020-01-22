%% PET and fMRI preprocessing pipeline main script
%  order of operation:
%   1. reslice - MRI
%   2. generate fieldmap - MRI
%   3. realign / unwarp - MRI, PET
%   4. coregister (estimate & write) - MRI, PET
%   5. coregister individual functional data to coregistered EPI (estimate) - MRI, PET
%   6. normalise - MRI
%   7. smoothe - MRI

%% work log

%   03-12-2019    created the script
%   18-01-2020    modified PET coregistration bit: now it incorporates T1
%       segmentation process as well
%   19-01-2020    added In-Flow and baseline preprocessing as well

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.spm     = '/Users/yeojin/Documents/MATLAB/spm12';
paths.funx_MRI= [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
paths.funx_PET= [paths.parent 'B_scripts/BA_preprocessing/BAC_PET/preproc_functions/'];
paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';
paths.seg     = [paths.parent 'E_data/EA_raw/EAB_MRI/EABD_segmented/'];
paths.TACs    = [paths.parent 'E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask/'];
paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/']

% add toolboxes and functions
addpath(genpath('/Users/yeojin/Documents/MATLAB/spm12'))
addpath(paths.funx_MRI)
addpath(paths.funx_PET)

% IDs
IDs = [4001];
days = [0 2];


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

% record history
% flags: 0 not done / 1 done
fg_1_resliced      = 0;
fg_2_realigned     = 0;
fg_3_coregistered  = 0;
fg_4_segmented     = 0;
fg_5_normalised    = 0;
fg_6_smoothed      = 0;

config = [];

fprintf('\n preparation done \n')
%% start preprocessing

spm fmri  % open progress window

%% MRI part

%% run

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            fprintf('\n*** ID %d being preprocessed: %2.0f out of %2.0f ***\n', IDs(id), id, length(IDs))
            
            %% manage configuration
            
            config.ID   = IDs(id);
            config.exp  = expdat{id,d}.dat;
            
            %% reslice
            
            if fg_1_resliced == 0
                
                clear nvols
                cd([paths.raw num2str(IDs(id)) '_' num2str(days(id,d))])
                tmp = dir('*study*'); cd(tmp(1).name); clear tmp;
                tmp = dir('*fMRI*'); cd(tmp(1).name);
                nvols = length(dir('MR*')); clear tmp
                
                for cc = 1:nvols
                    flist_PETtask{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end
%                 flist    = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                nslices  = 51;
                TR       = 3.6;
                refslice = 25;
                
                [config]        = preproc_reslice(IDs(id),flist_PETtask,nslices,TR,refslice,[1:2:51,2:2:50]);
                
                fg_1_resliced = 1;
                config.flags.resliced = 1;
                fprintf('\n *** resliced ***\n')
            else
                fprintf('\n Already resliced batch \n')
            end
            
            %% realign
            
            %     %% calculate VDM, then realign and unwarp
            %
            %     if fg_2_realigned == 0
            %
            %         shortTE         = 10;
            %         longTE          = 11.02;
            %         Total_TReadouot = 52.4698;
            %
            %         fMRI            = cellstr([paths.preproc num2str(IDs(id)) '/' num2str(IDs(id)) '_BL_4D_EmoMem2.nii']);
            %         Phase_FMap      = cellstr([paths.preproc num2str(IDs(id)) '/' num2str(IDs(id)) '_BL_4D_b1map2.nii']);
            %         Magnitude_FMap  = cellstr([paths.preproc num2str(IDs(id)) '/' num2str(IDs(id)) '_BL_4D_b1map1.nii']);
            %
            %         [config]        = preproc_realign_unwarp(IDs(id),Phase_FMap,Magnitude_FMap,fMRI,shortTE,longTE,Total_TReadout)
            %
            %         fg_2_realigned  = 1;
            %         config.flags.realigned = 1;
            %         fprintf('\n *** realigned and unwarped ***\n')
            %
            %     else
            %         fprintf('\n Already realigned batch \n')
            %     end
            %
            
            if fg_2_realigned == 0
                
                flist_PETtask = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/a' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                [config]        = preproc_realign(IDs(id),flist_PETtask);
                
                fg_2_realigned  = 1;
                config.flags.realigned = 1;
                fprintf('\n *** realigned ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            %% coregister
            
%             if fg_3_coregistered == 0
%                 
%                 fMRI_mean= cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
%                 try
%                     sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
%                     fMRI     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
%                 
%                     [config] = preproc_coregister(IDs(id),fMRI_mean,sMRI,fMRI);
%                 catch
%                     sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
%                     fMRI     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
%                 
%                     [config] = preproc_coregister(IDs(id),fMRI_mean,sMRI,fMRI);
%                 end
%                 fg_3_coregistered = 1;
%                 config.flags.coregistered = 1;
%                 fprintf('\n *** coregistered ***\n')
%             else
%                 fprintf('\n Already realigned batch \n')
%             end
            
            
            %% segment
            
%             if fg_4_segmented == 0
%                 
%                 [config] = preproc_segment(IDs(id),sMRI); % header information has changed from the previous step
%                 
%                 fg_4_segmented = 1;
%                 config.flags.segmented = 1;
%                 fprintf('\n *** segmented ***\n')
%             else
%                 fprintf('\n Already segmented batch \n')
%             end
            
            %% normalise
            
%             if fg_5_normalised == 0
%                 
%                 try
%                     Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
%                     [config] = preproc_normalise(IDs(id),Def_field,fMRI);
%                 catch
%                     Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
%                     [config] = preproc_normalise(IDs(id),Def_field,fMRI);
%                 end
%                 
%                 fg_5_normalised = 1;
%                 config.flags.normalised = 1;
%                 fprintf('\n *** normalised ***\n')
%             else
%                 fprintf('\n Already normalised batch \n')
%             end
            
            %% smoothe
            
            if fg_6_smoothed == 0
                
                fMRI      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
%                 fMRI      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/wra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
                
                [config] = preproc_smoothe(IDs(id),fMRI,4);
                
                fg_6_smoothed = 1;
                config.flags.smoothed = 1;
                fprintf('\n *** smoothed ***\n')
            else
                fprintf('\n Already smoothed batch \n')
            end
            
            
            %% wrap up
            
            % save configuration
            save([paths.history num2str(IDs(id)) '_' num2str(days(id,d)) '_config_MRI.mat'],'config')
            
            % flags hoisted again
            fg_1_resliced      = 0;
            fg_2_realigned     = 0;
            fg_3_coregistered  = 0;
            fg_4_segmented     = 0;
            fg_5_normalised    = 0;
            fg_6_smoothed      = 0;
            
        end
    end
    fprintf('\n***********\nID %d MRI preprocessing done\n***********\n',IDs(id))
    
end

%% PET part

%% run

clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            fprintf('\n*** ID %d being preprocessed: %2.0f out of %2.0f ***\n', IDs(id), id, length(IDs))
            
            %% manage configuration
            
            config.ID   = IDs(id);
            config.exp  = expdat{id,d}.dat;
            
            % prepare
            cd([paths.preproc num2str(IDs(id)) '_' num2str(d)]);
            seriesnum     = numel(dir([num2str(IDs(id)) '*PET_3D_T*MT*.nii'])); % see how many binned PET data there are
            binsize       = 300;
            
            %% reslice
            
%             if fg_1_resliced == 0
%                 
%                 clear nvols
%                 cd([paths.raw num2str(IDs(id)) '_' num2str(days(id,d))])
%                 tmp = dir('*study*'); cd(tmp(1).name); clear tmp;
%                 tmp = dir('*fMRI*'); cd(tmp(1).name);
%                 nvols = length(dir('MR*')); clear tmp
%                 
%                 for cc = 1:nvols
%                     flist{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
%                 end
% %                 flist    = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
%                 nslices  = 51;
%                 TR       = 3.6;
%                 refslice = 25;
%                 
%                 [config]        = preproc_reslice(IDs(id),flist,nslices,TR,refslice,[1:2:51,2:2:50]);
%                 
%                 fg_1_resliced = 1;
%                 config.flags.resliced = 1;
%                 fprintf('\n *** resliced ***\n')
%             else
%                 fprintf('\n Already resliced batch \n')
%             end
            
            %% realign
            
            if fg_2_realigned == 0
                
                % PET task
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']))
                    flist_PETtask{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end; clear cc
                [config]        = preproc_realign_PET(IDs(id),flist_PETtask);
                [config]        = preproc_realign_estwrt_PET(IDs(id),flist_PETtask);
                
                % PET inflow
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii']))
                    flist_PETflow{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii,' num2str(cc)];
                end;clear cc
                [config]        = preproc_realign_PET(IDs(id),flist_PETflow);
                [config]        = preproc_realign_estwrt_PET(IDs(id),flist_PETflow);
                
                % PET baseline
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii']))
                    flist_PETbsl{cc,1}     = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii,' num2str(cc)];
                end;clear cc
                [config]        = preproc_realign_PET(IDs(id),flist_PETbsl);
                [config]        = preproc_realign_estwrt_PET(IDs(id),flist_PETbsl);
                
                fg_2_realigned  = 1;
                config.flags.realigned = 1;
                fprintf('\n *** realigned ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            %% coregister
            
            % have you segmented your T1 image? this needs to be done as
            % well. how we go about this process is:
            %   (1) segment cerebellum with FSL, using our T1
            %   (2) reslice the T1 to fit in the segmentation generated
            %       from (1)
            %   (3) reslice and coregister 4D PET to this T1
            
            if fg_3_coregistered == 0
                
                PETtask_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']);
                PETflow_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii']);
                PETbsl_mean  = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii']);
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']))
                    PETtask{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)]; % estimated
                end; clear cc
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii']))
                    PETflow{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii,' num2str(cc)];
                end;clear cc
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii']))
                    PETbsl{cc,1}     = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii,' num2str(cc)];
                end;clear cc
                
                try
                    sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
                    
                    % step (2)
                    clear matlabbatch
                    spm_jobman('initcfg')
                    matlabbatch{1}.spm.spatial.coreg.write.ref = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/aparc+aseg.nii,1'}; % example FSL segmentation
                    matlabbatch{1}.spm.spatial.coreg.write.source = sMRI;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
                    spm_jobman('run', matlabbatch);
                    clear matlabbatch
                    sMRI_resliced = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
                    
                    % step (3)
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PETtask_mean,sMRI_resliced,PETtask);
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PETflow_mean,sMRI_resliced,PETflow);
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PETbsl_mean,sMRI_resliced,PETbsl);
                    
                catch
                    sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
                    
                    % step (2)
                    clear matlabbatch
                    spm_jobman('initcfg')
                    matlabbatch{1}.spm.spatial.coreg.write.ref = {'/Users/yeojin/Desktop/E_data/EE_atlases_templates/aparc+aseg.nii,1'};
                    matlabbatch{1}.spm.spatial.coreg.write.source = sMRI;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
                    spm_jobman('run', matlabbatch);
                    clear matlabbatch
                    sMRI_resliced = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
                    
                    % step (3)
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PETtask_mean,sMRI_resliced,PETtask);
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PETflow_mean,sMRI_resliced,PETflow);
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PETbsl_mean,sMRI_resliced,PETbsl);
                    
                end
                fg_3_coregistered = 1;
                config.flags.coregistered = 1;
                fprintf('\n *** coregistered ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            
            %% segment
            
%             if fg_4_segmented == 0
%                 
%                 [config] = preproc_segment_PET(IDs(id),sMRI); % header information has changed from the previous step
%                 
%                 fg_4_segmented = 1;
%                 config.flags.segmented = 1;
%                 fprintf('\n *** segmented ***\n')
%             else
%                 fprintf('\n Already segmented batch \n')
%             end
            
            %% normalise
            
%             if fg_5_normalised == 0
%                 
%                 
%                 try
%                     Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
%                     [config] = preproc_normalise(IDs(id),Def_field,PET);
%                 catch
%                     Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
%                     [config] = preproc_normalise(IDs(id),Def_field,PET);
%                 end
%                 
%                 fg_5_normalised = 1;
%                 config.flags.normalised = 1;
%                 fprintf('\n *** normalised ***\n')
%             else
%                 fprintf('\n Already normalised batch \n')
%             end
            
            %% smoothe -> do we have to?
            
%             if fg_6_smoothed == 0
%                 
% %                 for series = 1:seriesnum
% %                     fprintf('\n series %d\n', series)
%                     PET      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/cr' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
%                     [config] = preproc_smoothe(IDs(id),PET,4);
% %                 end
%                 
%                 fg_6_smoothed = 1;
%                 config.flags.smoothed = 1;
%                 fprintf('\n *** smoothed ***\n')
%             else
%                 fprintf('\n Already smoothed batch \n')
%             end
            
            
            %% wrap up
            
            % save configuration
            save([paths.history num2str(IDs(id)) '_' num2str(days(id,d)) '_config_PET.mat'],'config')
            
            % flags hoisted again
            fg_1_resliced      = 0;
            fg_2_realigned     = 0;
            fg_3_coregistered  = 0;
            fg_4_segmented     = 0;
            fg_5_normalised    = 0;
            fg_6_smoothed      = 0;
            
        end
    end
    fprintf('\n***********\nID %d preprocessing done\n***********\n',IDs(id))
    
end

%% Extract TACs

%% Define input parameters

addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')

% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
TrimPerc=15;

clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask   =[ paths.seg num2str(IDs(id)) num2str(d) '/mri/aparc+aseg.nii'];
            
            % 4-D PET file
            PETtask = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/c' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii'];
            PETflow = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/c' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii'];
            PETbsl  = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/c' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
            maskidx=[];
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,2}).LongName=strtrim(ROIDef{ROITabidx,3});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1}+1000; % Add 1000 to the index
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,1}+2000; % Add 2000 to the index
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}+1000, ROIDef{ROITabidx,1}+2000];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,2}).(Hemi{1})=IndList;
                end
            end
            
            
            [A B ROIDef]=xlsread('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/ExtractPETTACs/Subcortical_Dictionary.xls');
            
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,3}).LongName=strtrim(ROIDef{ROITabidx,4});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1};
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,2};
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}, ROIDef{ROITabidx,2}];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,3}).(Hemi{1})=IndList;
                end
            end
            
            %% Read in 4D-PET data and extract ROI-averages in each frame
            
            % task
            DynPET=load_nii(PETtask);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'RewardTask';
            save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
            % inFlow
            DynPET=load_nii(PETflow);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'inFlow';
            save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
            % baseline
            DynPET=load_nii(PETbsl);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'Baseline';
            save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
            %% plot TACs and save
            
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Times=[tt1; tt2+95*60; tt3+115*60]
            
            Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
            Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
            Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
            tmid=mean(Times,2)/60;
            plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
            xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
            legend('Cerebellum','Putamen','Caudate');
            cd(paths.figures)
            print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '.pdf']);
            
        end
    end
end

