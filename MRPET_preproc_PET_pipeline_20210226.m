%% PET data preprocessing pipeline
%  
%% work log

%   18-02-2020      separated this part of the script from the previous
%               pipeline that was both fMRI and PET

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
paths.seg     = [paths.parent 'E_data/EA_raw/EAD_PET/EADD_segmented/'];
paths.TACs    = [paths.parent 'E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask/'];
paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/']

% add toolboxes and functions
% addpath(genpath('/Users/yeojin/Documents/MATLAB/spm12'))
addpath(paths.funx_MRI)
addpath(paths.funx_PET)

% IDs
IDs     = [4016];%[4008 4009 4010 4011 4012 4013 4014 4015];
days    = [1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2]; 

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
                flist_PETtask = [];
                    for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']))
                        flist_PETtask{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                    end; clear cc
%                 end
                [config]        = preproc_realign_PET(IDs(id),flist_PETtask);
                [config]        = preproc_realign_estwrt_PET(IDs(id),flist_PETtask);
                disp('MT realigned')
                
                % PET inflow
                %%%%%%%%%%%%%%%%%
                % somehow inflow images look funny and throw errors so i'm
                % treating it differently
                %%%%%%%%%%%%%%%%%
                
                % how many volumes are there?
                flist_PETflow=[];
                for cc = 10:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii']))
                    flist_PETflow{cc-9,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii,' num2str(cc)];
                end;clear cc
                [config]        = preproc_realign_PET(IDs(id),flist_PETflow);
                [config]        = preproc_realign_estwrt_PET(IDs(id),flist_PETflow);
                disp('inflow realigned')
                
                
                
                % PET baseline
                flist_PETbsl=[];
                try
                    for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii']))
                        flist_PETbsl{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii,' num2str(cc)];
                    end;clear cc
                catch
                    for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '_1.nii']))
                        flist_PETbsl{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '_1.nii,' num2str(cc)];
                    end;clear cc
                end
                [config]        = preproc_realign_PET(IDs(id),flist_PETbsl);
                [config]        = preproc_realign_estwrt_PET(IDs(id),flist_PETbsl);
                disp('baseline realigned')
                
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
                
                % InFlow
                PETflow_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii']);
                PETflow=[];
                for cc = 1:(length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii']))-9)
                    PETflow{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii,' num2str(cc)];
                end;clear cc
                
                % baseline
                try
                    PETbsl_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii']);
                    PETbsl=[];
                    for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii']))
                        PETbsl{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii,' num2str(cc)];
                    end;clear cc
                catch
                    PETbsl_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '_1.nii']);
                    PETbsl=[];
                    for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '_1.nii']))
                        PETbsl{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '_1.nii,' num2str(cc)];
                    end;clear cc
                end
                
                % task
                PETtask_mean = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/mean' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']);
                for cc = 1:length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii']))
                    PETtask{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)]; % estimated
                end; clear cc
                                               
                
                % now coregister really
                sMRI_pt1      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_pt1.nii']);
                sMRI_pt2      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_pt2.nii']);
                seg_image_pt1 = {[paths.seg num2str(IDs(id)) num2str(d) '1/aparc+aseg.nii']};
                seg_image_pt2 = {[paths.seg num2str(IDs(id)) num2str(d) '2/aparc+aseg.nii']};
                
                % step (2) - 1: reslice T1 part 1 (InFlow part) with the corresponding segmentation.
                clear matlabbatch
                spm_jobman('initcfg')
                matlabbatch{1}.spm.spatial.coreg.write.ref = seg_image_pt1; % example FSL segmentation
                matlabbatch{1}.spm.spatial.coreg.write.source = sMRI_pt1;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'seg_';
                spm_jobman('run', matlabbatch);
                clear matlabbatch
                sMRI_resliced_pt1 = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/seg_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_pt1.nii']);
                
                % step (2) - 3: reslice T1 part 2 (Baseline and Task part) with the corresponding segmentation.
                clear matlabbatch
                spm_jobman('initcfg')
                matlabbatch{1}.spm.spatial.coreg.write.ref = seg_image_pt2; % example FSL segmentation
                matlabbatch{1}.spm.spatial.coreg.write.source = sMRI_pt2;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'seg_';
                spm_jobman('run', matlabbatch);
                clear matlabbatch
                sMRI_resliced_pt2 = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/seg_' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(days(id,d)) '_pt2.nii']);
                
                % step (3): coregister the dynamic images to the resliced T1
                [config] = preproc_coregister_estwrt_PET(IDs(id),PETflow_mean,sMRI_resliced_pt1,PETflow);
                [config] = preproc_coregister_estwrt_PET(IDs(id),PETbsl_mean,sMRI_resliced_pt2,PETbsl);
                [config] = preproc_coregister_estwrt_PET(IDs(id),PETtask_mean,sMRI_resliced_pt2,PETtask);
%                 [config] = preproc_coregister_wrt_PET(IDs(id),sMRI_resliced_pt1,PETflow);
%                 [config] = preproc_coregister_wrt_PET(IDs(id),sMRI_resliced_pt2,PETtask);
%                 [config] = preproc_coregister_wrt_PET(IDs(id),sMRI_resliced_pt2,PETbsl);
                
                
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
% this part is very heavy for some computers


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
paths.seg     = [paths.parent 'E_data/EA_raw/EAD_PET/EADD_segmented/'];
paths.TACs    = [paths.parent 'E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask_deccorr_2/'];
paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/']

% add toolboxes and functions
% addpath(genpath('/Users/yeojin/Documents/MATLAB/spm12'))
addpath(paths.funx_MRI)
addpath(paths.funx_PET)

% IDs
IDs     = [4008 4009 4010 4011 4012 4013 4014 4015];
days    = [1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2]; 

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
%%

addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')
opengl hardwarebasic


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
            Mask1   =[ paths.seg num2str(IDs(id)) num2str(d) '1/aparc+aseg.nii'];
            
            % 4-D PET file
            PETflow = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/c' num2str(IDs(id)) '_PET_4D_InFlow' num2str(days(id,d)) '.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask1);
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
        end
    end
end


clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask2   =[ paths.seg num2str(IDs(id)) num2str(d) '2/aparc+aseg.nii']; % task and the baseline
            
            % 4-D PET file
            PETtask = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/c' num2str(IDs(id)) '_PET_4D_MT' num2str(days(id,d)) '.nii'];
            PETbsl  = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/c' num2str(IDs(id)) '_PET_4D_Baseline' num2str(days(id,d)) '.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask2);
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
            
%             % inFlow
%             DynPET=load_nii(PETflow);
%             temp=size(DynPET.img);
%             ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
%             for ROI=fieldnames(maskidx)'
%                 TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
%                 for Hemi={'Left','Right','Bilateral'}
%                     maskidxcur = maskidx.(ROI{1}).(Hemi{1});
%                     TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
%                     TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
%                 end
%                 
%             end
%             TACDATA.info = 'inFlow';
%             save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
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
        end
    end
end


%% retro correct the decay

Appointments = [1 2 2 1 1 2 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    % bilaterals
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^(t0_frame_bsl/109));
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^(t0_frame_task/109));
                    
                    % lefts
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^(t0_frame_bsl/109));
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^(t0_frame_task/109));
                    
                    % rights
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^(t0_frame_bsl/109));
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^(t0_frame_task/109));
                end
            end
            
            save(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask_deccorr_2/' num2str(IDs(id)) num2str(d) '_TACs_deccorr.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
            try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '.pdf']);
                
            catch
                
                Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '.pdf']);

                
            end
            
            
        end
    end
end

%% problematic cases

Appointments = [1 2 2 1 1 2 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 6%1:length(IDs)
    
    for d = 1%1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    % bilaterals
%                     TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^(t0_frame_bsl/109));
%                     TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^(t0_frame_task/109));
                    
                    % lefts
%                     TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^(t0_frame_bsl/109));
%                     TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^(t0_frame_task/109));
                    
                    % rights
%                     TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^(t0_frame_bsl/109));
%                     TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^(t0_frame_task/109));
                end
            end
            
            save(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask_deccorr_2/' num2str(IDs(id)) num2str(d) '_TACs_deccorr.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
            try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '.pdf']);
                
            catch
                
                Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_noDecayCorrection.pdf']);

                
            end
            
            
        end
    end
end


%% modelling


clear; clc; close all
warning('off','all');

% set paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
% paths.spm     = '/Users/yeojin/Documents/MATLAB/spm12';
paths.funx_MRI= [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
paths.funx_PET= [paths.parent 'B_scripts/BA_preprocessing/BAC_PET/preproc_functions/'];
paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';
paths.seg     = [paths.parent 'E_data/EA_raw/EAB_MRI/EABD_segmented/'];
paths.TACs    = [paths.parent 'E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask_deccorr_2/'];
paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/'];

% set envs
% IDs
IDs     = [4008 4009 4010 4011 4012 4013 4014 4015];
days    = [1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2]; 
addpath('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/modelling')

% load TACs
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            TACs{i1,d} = [];
        else
            TACs{i1,d} = load([paths.TACs num2str(IDs(i1)) num2str(d) '_TACs_deccorr.mat']);
        end
    end
end

disp('prep done')

%% run the model

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg={'CerC','Put','Caud','Nac','Hipp','ThalProp'}
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
            end
            
            Lengths=[10*ones(30,1); 60*ones(55,1)]; % time windows
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
            
            tmid=mean(times,2);
            tmidMin=tmid/60;
            t_points    = length(tmid);
            dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
            break_point=find(times(:,1)>=115*60,1,'first'); %% TIme of activation start
            
            PlotStrFit=1;
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg={'Striatum','Caudate','Putamen','Accumbens-area','Thalamus','Hippocampus'}
                    if PlotStrFit
                        figure('Position',[100 100 800 1200],'Visible','on'); hold on;
                        spidx=1;
                    end
                    reftac=TACDATA.CerC;
                    mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
                    %% Set the SRTM part of A
                    ASRTM = zeros(t_points ,3);
                    ASRTM(:,1)  = reftac;
                    for k = 1:t_points
                        ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
                    end
                    refauc=sum(mreftac.*dt);
                    
                    if PlotStrFit
                        subplot(3,2,[1 2]);
                        plot(tmidMin,reftac,'k-'); hold on;
                        legendtext={'Cerebellum'};
                    end
                    switch(reg{1})
                        case 'Striatum'
                            tempTac=[];
                            tempVol=[];
                            for subReg={'Caud','Put'}
                                tempTac=[tempTac, TACDATA.(subReg{1})];
                                tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol];
                            end
                            roitac=sum(tempTac.*tempVol,2)./sum(tempVol);
                        case 'Caudate'
                            roitac=TACDATA.Caud;
                        case 'Putamen'
                            roitac=TACDATA.Put;
                        case 'Accumbens-area'
                            roitac=TACDATA.Nac;
                        case 'Thalamus'
                            roitac=TACDATA.ThalProp;
                        case 'Hippocampus'
                            roitac=TACDATA.Hipp;
                            
                    end
                    
                    mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2];
                    ASTRM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mroitac(1:k).*dt(1:k));
                    end
                    %LSQ-estimation using lscov
                    [parest se_srtm mse_srtm]   = lscov(ASRTM,roitac);
                    fittac=ASRTM*parest;
                    BP=parest(2)/parest(3)-1;
                    k2p=parest(2)/parest(1);
                    
                    %%%%% DO real SRTM
                    options = optimset('MaxFunEvals',1000);
                    weighs=[0.25*ones(30,1); ones(t_points-30,1)];
                    fobj = @(x) norm((simESRTM_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-roitac).*weighs);
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    fobj = @(x) norm((simESRTM_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-roitac(1:end-11)).*weighs(1:end-11));
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1_bl=parest_srtm(1);
                    k2_bl=parest_srtm(2);
                    BP_bl=parest_srtm(3);
                    modfit_esrtm_bl=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    %     weighs=[0.25*ones(30,1); ones(15+46,1); 100*ones(11,1)];
                    %     fobj = @(x) norm((simESRTM_1_0_0(tmidMin([1:76 92:t_points]),reftac([1:76 92:t_points]),t_points-15,x(1),x(2),x(3)*ones(t_points-15,1))-roitac([1:76 92:t_points])).*weighs([1:76 92:t_points]));
                    %     [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    %     R1_task=parest_srtm(1);
                    %     k2_task=parest_srtm(2);
                    %     BP_task=parest_srtm(3);
                    %     modfit_esrtm_task=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    
                    %%% Do lp-ntPET
                    Alpntpet=zeros(t_points,4);
                    Alpntpet(:,1:3)=ASRTM;
                    best_mse=10^20;
                    best_parest=[];
                    best_se=[];
                    % for estimating the optimal peak and alpha parameters
                    for point_rise=break_point:length(tmid)-1
                        t2p_index = find(tmid>tmid(point_rise)); %times(point_rise,1)
                        if length(t2p_index)>1
                            for t_ind=t2p_index(1):t2p_index(end-1)  %t=1:0.1:20 %
                                for alpha=[0.25 1 4]
                                    t_peak=tmid(t_ind)-tmid(point_rise);   %  times(point_rise,1)
                                    p = [1 alpha tmid(point_rise)+t_peak tmid(point_rise)];
                                    actfun = zeros(size(tmid));
                                    actfun(point_rise:t_points) = gamma_variate_madsen(p,tmid(point_rise:t_points));
                                    roitac_gamma = roitac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(roitac))+roitac_gamma(1:length(roitac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,roitac);
                                    
                                    %estimated model TAC
                                    modelfit = Alpntpet*parest;
                                    
                                    if (best_mse > mse)
                                        best_mse = mse;
                                        best_parest=parest;
                                        best_modelfit=modelfit;
                                        best_actfun=actfun;
                                        best_se=se;
                                    end
                                end
                            end                                        %                        breakpoint2=breakpoint;
                        end
                    end
                    BP_lp=best_parest(2)/best_parest(3)-1;
                    
                    fprintf(1,'Here\n')
                    
                    if PlotStrFit
                        legendtext{end+1}=[reg{1} ' raw'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{Baseline}'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{All}'];
                        legendtext{end+1}=[reg{1} ' fit lpnt-PET_{All}'];
                        legendtext{end+1}=['Activation function'];
                        SE=abs(BP)*sqrt((se_srtm(2)/parest(2))^2+(se_srtm(3)/parest(3))^2);
                        subplot(3,2,[1 2]);
                        yyaxis left;
                        plot(tmidMin,roitac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
                        %xlim([0 60]);
                        xlabel('Time (min)');
                        ylabel('Radioactivity concentration');
                        legend(legendtext,'Location','south');
                        yyaxis right;
                        plot(tmidMin,best_parest(4)*best_actfun);
                        ylabel('Compensatory function');
                        %ylim([0 4*10^(-4)]);
                        title([Subj{r} ' compartmental fits: BP_{Baseline}=' num2str(BP_bl,'%1.2f') ', BP_{All}=' num2str(BP__,'%1.2f') ', BP_{lp-nt}=' num2str(BP_lp,'%1.2f')]);
                        subplot(3,2,3);
                        %         plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'go',tmidMin,roitac-modfit_esrtm_bl,'co',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'bo',tmidMin,roitac-modfit_esrtm_bl,'ro',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        [h p]=runstest(roitac-fittac);
                        [h1 p1]=runstest(roitac-modfit_esrtm);
                        [h2 p2]=runstest(roitac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
                        title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
                        xlabel('Time (min)');
                        
                        subplot(3,2,4);
                        y=-ASRTM(:,3)./roitac;
                        x=ASRTM(:,2)./roitac;
                        [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                        [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                        plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
                        title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        xlabel(['\int REF/ROI']);
                        ylabel('\int ROI/ROI')
                        subplot(3,2,[5 6]);
                        plot(tmidMin,roitac./reftac,'ko',[0 180],[0 0],'k--');
                        %ylim([0.5 5]);
                        title('Target to reference ratio');
                        xlabel('Time (min)');
                        
                        print('-dpsc2','-append','-bestfit',fullfile(paths.figures, [ num2str(IDs(id)) num2str(d) '_TAC_Fit_lpntpet_logan_framesDropped.ps']));
%                                 pause
                        close(gcf)
                        BPdata.BP_mrtm(Subj{r})=BP;
                        BPdata.BP_srtm(Subj{r})=BP__;
                        BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                        BPdata.BP_lpnt(Subj{r})=BP_lp;
                        BPdata.BP_logan(Subj{r})=pp(1)-1;
                        %return
                        
                        

                        
%                         continue;
                    end
                end
                
                
            end
            close all
            keep IDs days paths id d TACs
        end
    end
end


disp('done')



