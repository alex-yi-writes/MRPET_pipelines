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


%% PET part

%% run

clear id d

for id = 3%:length(IDs)
    
    for d = 2%1:2
        
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

addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')
opengl hardwarebasic

% IDs
IDs     = [4008 4009 4010 4011 4012 4013 4014 4015];
days    = [1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2]; 

% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
TrimPerc=15;

clear id d

for id = 4:length(IDs)
    
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
            
%             % task
%             DynPET=load_nii(PETtask);
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
%             TACDATA.info = 'RewardTask';
%             save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
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
%             
%             % baseline
%             DynPET=load_nii(PETbsl);
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
%             TACDATA.info = 'Baseline';
%             save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end

%%
clear id d

for id = 3:length(IDs)
    
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


% IDs     = [4008 4009 4010 4011 4012 4013 4014 4015];
% days    = [1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2]; 
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
            
            save(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask_decaycorrected/' num2str(IDs(id)) num2str(d) '_TACs_deccorr.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
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
                
            end
            
            
        end
    end
end


%% plot TACs and save

close all

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
            
            try
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; (TACDATA_Baseline.CerC.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.CerC.Bilateral.tac).*(2^(t0_frame_task/109))];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; (TACDATA_Baseline.Put.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Put.Bilateral.tac).*(2^(t0_frame_task/109))];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; (TACDATA_Baseline.Caud.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Caud.Bilateral.tac).*(2^(t0_frame_task/109))];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                ylim([0 30000]);
                clear titlestring
                if Appointments(id) == 1
                    titlestring = ['ID ' num2str(IDs(id)) ', morning appointment, Session type ' num2str(d)];
                else
                    titlestring = ['ID ' num2str(IDs(id)) ', afternoon appointment, Session type ' num2str(d)];
                end
                title(titlestring,'Fontsize',20,'Fontweight','bold')
                print('-dpdf','-bestfit', fullfile(paths.figures, [ num2str(IDs(id)) num2str(d) '_TAC.pdf']));
                
            catch
                Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); (TACDATA_Baseline.CerC.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.CerC.Bilateral.tac).*(2^(t0_frame_task/109))];
                Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); (TACDATA_Baseline.Put.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Put.Bilateral.tac).*(2^(t0_frame_task/109))];
                Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); (TACDATA_Baseline.Caud.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Caud.Bilateral.tac).*(2^(t0_frame_task/109))];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                ylim([0 30000]);
                clear titlestring
                if Appointments(id) == 1
                    titlestring = ['ID ' num2str(IDs(id)) ', morning appointment, Session type ' num2str(d)];
                else
                    titlestring = ['ID ' num2str(IDs(id)) ', afternoon appointment, Session type ' num2str(d)];
                end
                title(titlestring,'Fontsize',20,'Fontweight','bold')
                print('-dpdf','-bestfit', fullfile(paths.figures, [ num2str(IDs(id)) num2str(d) '_TAC.pdf']));
                
            end

        end
    end
end

