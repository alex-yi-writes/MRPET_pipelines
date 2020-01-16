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

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.spm     = '/Users/yeojin/Documents/MATLAB/spm12';
paths.funx    = [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';

% add toolboxes and functions
addpath(paths.spm)
addpath(paths.funx)

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
                    flist{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end
%                 flist    = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                nslices  = 51;
                TR       = 3.6;
                refslice = 25;
                
                [config]        = preproc_reslice(IDs(id),flist,nslices,TR,refslice,[1:2:51,2:2:50]);
                
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
                
                flist = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/a' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                [config]        = preproc_realign(IDs(id),flist);
                
                fg_2_realigned  = 1;
                config.flags.realigned = 1;
                fprintf('\n *** realigned ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            %% coregister
            
            if fg_3_coregistered == 0
                
                fMRI_mean= cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                try
                    sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
                    fMRI     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                
                    [config] = preproc_coregister(IDs(id),fMRI_mean,sMRI,fMRI);
                catch
                    sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
                    fMRI     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                
                    [config] = preproc_coregister(IDs(id),fMRI_mean,sMRI,fMRI);
                end
                fg_3_coregistered = 1;
                config.flags.coregistered = 1;
                fprintf('\n *** coregistered ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            
            %% segment
            
            if fg_4_segmented == 0
                
                [config] = preproc_segment(IDs(id),sMRI); % header information has changed from the previous step
                
                fg_4_segmented = 1;
                config.flags.segmented = 1;
                fprintf('\n *** segmented ***\n')
            else
                fprintf('\n Already segmented batch \n')
            end
            
            %% normalise
            
            if fg_5_normalised == 0
                
                try
                    Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
                    [config] = preproc_normalise(IDs(id),Def_field,fMRI);
                catch
                    Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
                    [config] = preproc_normalise(IDs(id),Def_field,fMRI);
                end
                
                fg_5_normalised = 1;
                config.flags.normalised = 1;
                fprintf('\n *** normalised ***\n')
            else
                fprintf('\n Already normalised batch \n')
            end
            
            %% smoothe
            
            if fg_6_smoothed == 0
                
%                 fMRI      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
                fMRI      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/wra' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
                
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
    fprintf('\n***********\nID %d preprocessing done\n***********\n',IDs(id))
    
end

%% PET part

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
                
                for cc = 1:seriesnum
                    flist{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end
                
%                 flist = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']);
                [config]        = preproc_realign_PET(IDs(id),flist);
                
                fg_2_realigned  = 1;
                config.flags.realigned = 1;
                fprintf('\n *** realigned ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            %% coregister
            
            if fg_3_coregistered == 0
                
                PET_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']);
                try
                    sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '.nii']);
                    for cc = 1:seriesnum
                        PET{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                    end
%                     sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/NLreg_T1WB_to_MNI.nii']);
%                     PET  = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']);
                    
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PET_mean,sMRI,PET);
                catch
                    sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
                    for cc = 1:seriesnum
                        PET{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                    end
%                     PET  = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']);
%                     PET     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/r' num2str(IDs(id)) '_PET_4D_T' num2str((series-1)*binsize) '_MT' num2str(days(id,d)) '.nii']);
                    
                    [config] = preproc_coregister_estwrt_PET(IDs(id),PET_mean,sMRI,PET);
                end
                fg_3_coregistered = 1;
                config.flags.coregistered = 1;
                fprintf('\n *** coregistered ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            
            %% segment
            
            if fg_4_segmented == 0
                
                [config] = preproc_segment_PET(IDs(id),sMRI); % header information has changed from the previous step
                
                fg_4_segmented = 1;
                config.flags.segmented = 1;
                fprintf('\n *** segmented ***\n')
            else
                fprintf('\n Already segmented batch \n')
            end
            
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