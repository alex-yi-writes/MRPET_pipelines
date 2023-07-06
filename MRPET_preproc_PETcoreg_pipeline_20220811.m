%% new PET registrations

clc;clear;

% IDs
% IDs         = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4026];
% days        = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0; 0 2];

% IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
% days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];
% d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
% d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; 
IDs  = [4014];
days = [1 0];

% set env
path_parent = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/';
path_segs   = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/';

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin:/Applications/freesurfer/7.2.0/bin:/usr/local/bin']);
setenv('ANTSPATH','/usr/local/bin')

%% convert seg files

for id=1:length(IDs)

    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            try                
            gunzip([path_segs num2str(IDs(id)) num2str(d) '1/T1pt1.nii.gz'])
            gunzip([path_segs num2str(IDs(id)) num2str(d) '2/T1pt2.nii.gz'])
            catch
            end

            try
                
                copyfile([path_segs num2str(IDs(id)) num2str(d) '1/T1pt1.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'])
                copyfile([path_segs num2str(IDs(id)) num2str(d) '2/T1pt2.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'])
                copyfile([path_segs num2str(IDs(id)) num2str(d) '1/mri/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'])
                copyfile([path_segs num2str(IDs(id)) num2str(d) '2/mri/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'])
            catch
%                 copyfile([path_segs num2str(IDs(id)) num2str(d) '1/' num2str(IDs(id)) num2str(d) '1/T1WB.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'])
%                 copyfile([path_segs num2str(IDs(id)) num2str(d) '2/' num2str(IDs(id)) num2str(d) '2/T1WB.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'])
                try
                    copyfile([path_segs num2str(IDs(id)) num2str(d) '1/' num2str(IDs(id)) num2str(d) '1/mri/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'])
                    copyfile([path_segs num2str(IDs(id)) num2str(d) '2/' num2str(IDs(id)) num2str(d) '2/mri/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'])
                catch
                    eval(['!mri_convert -it mgz -ot nii ' path_segs num2str(IDs(id)) num2str(d) '1/' num2str(IDs(id)) num2str(d) '1/mri/aparc+aseg.mgz ' path_segs num2str(IDs(id)) num2str(d) '1/' num2str(IDs(id)) num2str(d) '1/mri/aparc+aseg.nii'])
                    eval(['!mri_convert -it mgz -ot nii ' path_segs num2str(IDs(id)) num2str(d) '2/' num2str(IDs(id)) num2str(d) '2/mri/aparc+aseg.mgz ' path_segs num2str(IDs(id)) num2str(d) '2/' num2str(IDs(id)) num2str(d) '2/mri/aparc+aseg.nii'])
                    copyfile([path_segs num2str(IDs(id)) num2str(d) '1/' num2str(IDs(id)) num2str(d) '1/mri/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'])
                    copyfile([path_segs num2str(IDs(id)) num2str(d) '2/' num2str(IDs(id)) num2str(d) '2/mri/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'])
                end
            end

        end
    end
end
%% extra

% for id=1:length(IDs)
%     
%     for d=1:2
%         if days(id,d)==0
%             warning('skipped')
%         else
%             
%             mkdir(['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d)])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
%                 ['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d) '/T1pt2.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2_nat.nii'],...
%                 ['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat.nii'])
%         end
%         
%     end
% end

%% move files for ANTs coreg

% for id=1:length(IDs)
%     
%     for d=1:2
%         if days(id,d)==0
%             warning('skipped')
%         else
%             
%             mkdir(['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/'])
% 
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/T1pt1.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/T1pt2.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/aparc+aseg_pt1.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/aparc+aseg_pt2.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/meanEPI.nii'])
%         end
%         
%     end
% end


%% convert seg files

% for id=1:length(IDs)
%     
%     for d=1:2
%         if days(id,d)==0
%             warning('skipped')
%         else
%             
% 
%             movefile([path_segs num2str(IDs(id)) num2str(d) '1/T1.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'])
%             movefile([path_segs num2str(IDs(id)) num2str(d) '2/T1.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'])
%             movefile([path_segs num2str(IDs(id)) num2str(d) '1/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'])
%             movefile([path_segs num2str(IDs(id)) num2str(d) '2/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'])
%             
%         end
%         
%     end
% end



%% first, realign and create
% id 23 d 2 --> do it again

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            disp('***************************')
            disp('***************************')
            disp(['******  ' num2str(IDs(id)) '_' num2str(d)  '  *******'])
            disp('***************************')
            disp('***************************')

            %% register FSL T1 and the segmentation to the native T1 
            eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii -f '...
                path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_FSL1toNative_' ])
            eval(['!antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_FSL1toNative_0GenericAffine.mat -i '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii -r '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/aparc+aseg_pt1_nat.nii'])
                
            eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii -f '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_FSL2toNative_' ])
            eval(['!antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_FSL2toNative_0GenericAffine.mat -i '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii -r '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/aparc+aseg_pt2_nat.nii'])
            
            %% inflow
            
            list_inflow=[]; clear numvol
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
            for v1=10:numvol
                list_inflow{v1-9,1} = [path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii,' num2str(v1)];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_inflow};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            spm_jobman('run',matlabbatch)
            
            
            mkdir([path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/mean' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii -f '...
                path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlowtoT1pt1_' ])
            % apply transformations to the images
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlowtoT1pt1_0GenericAffine.mat -i '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '_' sprintf('%05i',v2) '.nii -r '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_inflow=[];
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1_' sprintf('%05i',v1) '.nii,1'];
            end
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_inflow;
            matlabbatch{1}.spm.util.cat.name = ['coreg_InFlow' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1.nii'],...
                [path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlow' num2str(d) '_on_T1.nii'])
            rmdir([path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D'],'s')
            
            disp('===========')
            disp('inflow done')
            disp('===========')

            %% baseline
            list_bsl=[];
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(d) '.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_Baseline' num2str(d) '.nii,' num2str(v1)];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_bsl};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            spm_jobman('run',matlabbatch)
            
            mkdir([path_parent num2str(IDs(id)) '_' num2str(d) '/bsl3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_Baseline' num2str(d) '.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent num2str(IDs(id)) '_' num2str(d) '/bsl3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/mean' num2str(IDs(id)) '_PET_4D_Baseline' num2str(d) '.nii -f '...
                path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_BaselinetoT1pt2_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_BaselinetoT1pt2_0GenericAffine.mat -i '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/bsl3D/r' num2str(IDs(id)) '_PET_4D_Baseline' num2str(d) '_' sprintf('%05i',v2) '.nii -r '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/bsl3D/coreg_Baseline' num2str(d) '_on_T1_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_bsl=[];
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d)  '/bsl3D/coreg_Baseline' num2str(d) '_on_T1_' sprintf('%05i',v1) '.nii,1'];
            end
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_bsl;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Baseline' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent num2str(IDs(id)) '_' num2str(d)  '/bsl3D/coreg_Baseline' num2str(d) '_on_T1.nii'],...
                [path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_Baseline' num2str(d) '_on_T1.nii'])
            rmdir([path_parent num2str(IDs(id)) '_' num2str(d) '/bsl3D'],'s')

            disp('=============')
            disp('baseline done')
            disp('=============')
            
            %% task
            list_task=[];
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(d) '.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_MT' num2str(d) '.nii,' num2str(v1)];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_task};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            spm_jobman('run',matlabbatch)
            
            mkdir([path_parent num2str(IDs(id)) '_' num2str(d) '/task3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_MT' num2str(d) '.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent num2str(IDs(id)) '_' num2str(d) '/task3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/mean' num2str(IDs(id)) '_PET_4D_MT' num2str(d) '.nii -f '...
                path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_MTtoT1pt2_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_MTtoT1pt2_0GenericAffine.mat -i '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/task3D/r' num2str(IDs(id)) '_PET_4D_MT' num2str(d) '_' sprintf('%05i',v2) '.nii -r '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/task3D/coreg_MT' num2str(d) '_on_T1_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_task=[];
            for v1=1:numvol
                list_task{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d)  '/task3D/coreg_MT' num2str(d) '_on_T1_' sprintf('%05i',v1) '.nii,1'];
            end
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_task;
            matlabbatch{1}.spm.util.cat.name = ['coreg_MT' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent num2str(IDs(id)) '_' num2str(d)  '/task3D/coreg_MT' num2str(d) '_on_T1.nii'],...
                [path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_MT' num2str(d) '_on_T1.nii'])
            rmdir([path_parent num2str(IDs(id)) '_' num2str(d) '/task3D'],'s')

            disp('=========')
            disp('task done')
            disp('=========')
            
        end

        disp('***************************')
        disp('***************************')
        disp(['***  ' num2str(IDs(id)) '_' num2str(d)  ' done  ****'])
        disp('***************************')
        disp('***************************')

    end
    
end

%% write out TACs

addpath('/Users/yeojin/Desktop/B_scripts/BA_preprocessing/BAC_PET/preproc_functions/')


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
            Mask1   =[ path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1_nat.nii'];
            
            % 4-D PET file
            PETflow = [path_parent num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_InFlow' num2str(d) '_on_T1.nii'];
            
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
                disp(maskidx.(ROI{1}).LongName)
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'inFlow';
            mkdir(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d)])
            save([ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData            
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
            Mask2   =[ path_parent num2str(IDs(id)) '_' num2str(days(id,d)) '/aparc+aseg_pt2_nat.nii']; % task and the baseline
            
            % 4-D PET file
            PETtask = [path_parent num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_MT' num2str(d) '_on_T1.nii'];
            PETbsl  = [path_parent num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_Baseline' num2str(d) '_on_T1.nii'];
            
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
            save([ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
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
            save([ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end


%% retro correct the decay
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            if IDs(id)==4020 && d==2
            Lengths=[10*ones(30-9,1); 60*ones(55,1)];
            else
            Lengths=[10*ones(30,1); 60*ones(55,1)];    
            end
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
            
            save(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/TACs/RewardTask_ANTscoreg/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_decorr.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            if IDs(id)==4020 && d==2
            Lengths=[10*ones(30-9,1); 60*ones(55,1)];
            else
            Lengths=[10*ones(30,1); 60*ones(55,1)];    
            end
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
                print('-dpdf','-bestfit',[ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/' num2str(IDs(id)) num2str(d) '_' date '.pdf']);
                
            catch
                
                if IDs(id)==2031 && d==1
                    Cer=[vertcat(TACDATA_InFlow.CerC.Bilateral.tac,nan(14,1)); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                    Put=[vertcat(TACDATA_InFlow.Put.Bilateral.tac,nan(14,1)); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                    Caud=[vertcat(TACDATA_InFlow.Caud.Bilateral.tac,nan(14,1)); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                    tmid=mean(Times,2)/60;

                else
                    Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                    Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                    Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                    tmid=mean(Times,2)/60;

                end
                
                % now draw
                figure('Renderer', 'painters')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/' num2str(IDs(id)) num2str(d) '_' date '.pdf']);
                
            end
            
            
        end
    end
end
