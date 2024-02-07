%% new PET registrations

clc;clear;

path_parent = '/Users/yyi/Desktop/new_recon/';
path_segs   = '/Users/yyi/Desktop/new_recon/segs/';

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
% days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

%% copy files to the workspace

for id=12:length(IDs)
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else

            mkdir([path_parent 'segs/' num2str(IDs(id)) num2str(d)])
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d)])

            copyfile(['/Volumes/korokdorf/MRPET/coreg_roi/' num2str(IDs(id)) num2str(d) '/data/T1pt1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii'])
            copyfile(['/Volumes/korokdorf/MRPET/coreg_roi/' num2str(IDs(id)) num2str(d) '/data/T1pt2.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii'])

            copyfile(['/Volumes/korokdorf/MRPET/coreg_roi/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii'])
            copyfile(['/Volumes/korokdorf/MRPET/coreg_roi/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii'])

        end
    end
end

%% convert dicoms to niftiis

for id=29:length(IDs)
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else

            % convert
            clear tmp kuerzel studyfolder
            tmp=dir([path_parent 'raw/' num2str(IDs(id)) num2str(d) '/*.zip']);
            kuerzel = tmp.name(1:4);
            clear tmp
            tmp=dir([path_parent 'raw/' num2str(IDs(id)) num2str(d) '/' kuerzel '/study*'])
            studyfolder=tmp.name;

            eval(['!/Applications/MRIcron.app/Contents/Resources/dcm2niix ' ...
                '-f "%f_%p_%t_%s" -p y -z n -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) ' ' ...
                path_parent 'raw/' num2str(IDs(id)) num2str(d) '/' kuerzel '/' studyfolder])

            % rename
            clear tmp
            tmp=dir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/*InFlow*.nii']);
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/' tmp.name ], ['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ])
            movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ], [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow.nii' ] )

            clear tmp
            tmp=dir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/*Baseline*.nii']);
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/' tmp.name ], ['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ])
            movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ], [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl.nii' ] )

            clear tmp
            tmp=dir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/*Task*.nii']);
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/' tmp.name ], ['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ])
            movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ], [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task.nii' ] )



        end
    end

end


%% set env

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/bin']);
setenv('ANTSPATH','/usr/local/bin')

%% rename

% for id=5:12%length(IDs)
%     
%     for d=1:2
%         if days(id,d)==0
%             warning('skipped')
%         else
%             clear tmp
%             tmp=dir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/*InFlow*.nii']);
%             movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/' tmp.name ], ['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ])
%             movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ], [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow.nii' ] )
% 
%             clear tmp
%             tmp=dir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/*Baseline*.nii']);
%             movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/' tmp.name ], ['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ])
%             movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ], [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl.nii' ] )
% 
%             clear tmp
%             tmp=dir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/*Task*.nii']);
%             movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/' tmp.name ], ['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ])
%             movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/' tmp.name ], [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task.nii' ] )
% 
%         end
%     end
% end

%% extra

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            
            mkdir(['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d)])
            copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
                ['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d) '/T1pt2.nii'])
            copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2_nat_labelled.nii'],...
                ['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii'])
        end
        
    end
end


%% first, realign and create

for id=25:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            %% start

            disp('**************************')
            disp(['Starting ID ' num2str(IDs(id)) num2str(d)])
            disp('**************************')


            %% baseline

            list_bsl=[]; clear numvol
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl.nii,' num2str(v1)];
            end

            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D'])

            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D']};
            spm_jobman('run',matlabbatch)

            % register each PET frame to T1: linear interpolation by
            % default
            for v2=1:numvol
                eval(['!antsRegistrationSyNQuick.sh -d 3 -n 4 -t r -f ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -m '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/bsl_' sprintf('%05i',v2) '.nii -o '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/rbsl_' sprintf('%05i',v2) '_'])
            end

            % cleanup
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/waste'])
            for v2=1:numvol
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/rbsl_' sprintf('%05i',v2) '_0GenericAffine.mat'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/waste/rbsl_' sprintf('%05i',v2) '_0GenericAffine.mat'])
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/rbsl_' sprintf('%05i',v2) '_InverseWarped.nii.gz'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/waste/rbsl_' sprintf('%05i',v2) '_InverseWarped.nii.gz'])
            end
            eval(['!gzip ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/bsl_*.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/waste'],'s')

            for v2=1:numvol
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/rbsl_' sprintf('%05i',v2) '_Warped.nii.gz' ], ...
                    ['/Users/yyi/Desktop/new_recon/preproc/tmp/rbsl_' sprintf('%05i',v2) '_Warped.nii.gz' ])
                movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/rbsl_' sprintf('%05i',v2) '_Warped.nii.gz'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/bsl_on_t1pt2_' sprintf('%05i',v2) '.nii.gz' ] )
            end


            %% task

            list_task=[]; clear numvol
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task.nii,' num2str(v1)];
            end

            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D'])

            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D']};
            spm_jobman('run',matlabbatch)

            % register each PET frame to T1
            for v2=1:numvol
                eval(['!antsRegistrationSyNQuick.sh -d 3 -n 4 -t r -f ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -m '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/task_' sprintf('%05i',v2) '.nii -o '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/rtask_' sprintf('%05i',v2) '_'])
            end

            % cleanup
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/waste'])
            for v2=1:numvol
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/rtask_' sprintf('%05i',v2) '_0GenericAffine.mat'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/waste/rtask_' sprintf('%05i',v2) '_0GenericAffine.mat'])
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/rtask_' sprintf('%05i',v2) '_InverseWarped.nii.gz'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/waste/rtask_' sprintf('%05i',v2) '_InverseWarped.nii.gz'])
            end
            eval(['!gzip ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/task_*.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/waste'],'s')

            for v2=1:numvol
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/rtask_' sprintf('%05i',v2) '_Warped.nii.gz' ], ...
                    ['/Users/yyi/Desktop/new_recon/preproc/tmp/rtask_' sprintf('%05i',v2) '_Warped.nii.gz' ])
                movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/rtask_' sprintf('%05i',v2) '_Warped.nii.gz'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/task_on_t1pt2_' sprintf('%05i',v2) '.nii.gz' ] )
            end

            %% inflow
            
            list_inflow=[]; clear numvol
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow.nii,' num2str(v1)];
            end

            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D'])

            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D']};
            spm_jobman('run',matlabbatch)

            % register each PET frame to T1
            for v2=1:numvol
                eval(['!antsRegistrationSyNQuick.sh -d 3 -n 4 -t r -f ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii -m '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/inflow_' sprintf('%05i',v2) '.nii -o '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/rinflow_' sprintf('%05i',v2) '_'])
            end

            % cleanup
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/waste'])
            for v2=1:numvol
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/rinflow_' sprintf('%05i',v2) '_0GenericAffine.mat'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/waste/rinflow_' sprintf('%05i',v2) '_0GenericAffine.mat'])
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/rinflow_' sprintf('%05i',v2) '_InverseWarped.nii.gz'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/waste/rinflow_' sprintf('%05i',v2) '_InverseWarped.nii.gz'])
            end
            eval(['!gzip ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/inflow_*.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/waste'],'s')

            for v2=1:numvol
                movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/rinflow_' sprintf('%05i',v2) '_Warped.nii.gz' ], ...
                    ['/Users/yyi/Desktop/new_recon/preproc/tmp/rinflow_' sprintf('%05i',v2) '_Warped.nii.gz' ])
                movefile(['/Users/yyi/Desktop/new_recon/preproc/tmp/rinflow_' sprintf('%05i',v2) '_Warped.nii.gz'],...
                    [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/inflow_on_t1pt1_' sprintf('%05i',v2) '.nii.gz' ] )
            end
            
        end
    end
    
end

%% smoothe - absolute

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            list_inflow=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Inflow' num2str(d) '_on_T1.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Inflow' num2str(d) '_on_T1.nii,'  num2str(v1)];
            end
            list_bsl=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Baseline' num2str(d) '_on_T1.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Baseline' num2str(d) '_on_T1.nii,'  num2str(v1)];
            end
            list_task=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Task' num2str(d) '_on_T1.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Task' num2str(d) '_on_T1.nii,'  num2str(v1)];
            end
            
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = list_inflow;
            matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's3_';
            matlabbatch{2}.spm.spatial.smooth.data = list_bsl;
            matlabbatch{2}.spm.spatial.smooth.fwhm = [3 3 3];
            matlabbatch{2}.spm.spatial.smooth.dtype = 0;
            matlabbatch{2}.spm.spatial.smooth.im = 0;
            matlabbatch{2}.spm.spatial.smooth.prefix = 's3_';
            matlabbatch{3}.spm.spatial.smooth.data = list_task;
            matlabbatch{3}.spm.spatial.smooth.fwhm = [3 3 3];
            matlabbatch{3}.spm.spatial.smooth.dtype = 0;
            matlabbatch{3}.spm.spatial.smooth.im = 0;
            matlabbatch{3}.spm.spatial.smooth.prefix = 's3_';
            spm_jobman('run',matlabbatch)
            
        end
    end
end


%% write out TACs

% addpath('/Users/yeojin/Desktop/B_scripts/BA_preprocessing/BAC_PET/preproc_functions/')


% addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')
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
            Mask1   =[ path_parent 'segs/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii'];
            
            % 4-D PET file
%             PETflow = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/coreg_InFlow' num2str(d) '_on_T1.nii'];
            PETflow = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/s3_coreg_InFlow' num2str(d) '_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask1);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
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
            
            
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
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
%             mkdir(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d)])
%             save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_abs.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData            
            save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_abs_s3.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData            
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
            Mask2   =[ path_parent 'segs/' num2str(IDs(id)) num2str(days(id,d)) '/aparc+aseg_pt2_nat_labelled.nii']; % task and the baseline
            
            % 4-D PET file
%             PETtask = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/coreg_Task' num2str(d) '_on_T1.nii'];
%             PETbsl  = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/coreg_Baseline' num2str(d) '_on_T1.nii'];
            PETtask = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/s3_coreg_Task' num2str(d) '_on_T1.nii'];
            PETbsl  = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/s3_coreg_Baseline' num2str(d) '_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask2);
            ROIMask=round(ROI.img); 
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
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
            
            
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
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
%             save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_abs.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_abs_s3.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
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
%             save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_abs.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_abs_s3.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end

%% originals

close all

Appointments = [1 2 2 1 1 2 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_abs.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_abs.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_abs.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    % bilaterals
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac);
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac);
                    
                    % lefts
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac);
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac);
                    
                    % rights
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac);
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac);
                end
            end
            
            save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_uncorr_abs.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
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
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_uncorrected_abs.pdf']);
                
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
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_uncorrected_abs.pdf']);

                
            end
            
            
        end
    end
end


%% retro correct the decay

close all

Appointments = [1 2 2 1 1 2 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
%             load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_abs.mat']);
%             TACDATA_InFlow=TACDATA; clear TACDATA
%             load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_abs.mat']);
%             TACDATA_Baseline=TACDATA; clear TACDATA
%             load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_abs.mat']);
%             TACDATA_Task=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_abs_s3.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_abs_s3.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_abs_s3.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
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
            
%             save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_abs.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_abs_s3.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
%             try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_abs_s3.pdf']);
                
%             catch
%                 
%                 Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
%                 Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
%                 Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
%                 tmid=mean(Times,2)/60;
%                 
%                 % now draw
%                 figure('Renderer', 'painters ')
%                 plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
%                 xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
%                 legend('Cerebellum','Putamen','Caudate');
%                 title([ num2str(IDs(id)) '_' num2str(d) ]);
%                 ax = gca; ax.YAxis.Exponent = 0;
%                 print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_abs.pdf']);
% 
%                 
%             ends
            
            
        end
    end
end


%% relative

%% first, realign and create

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            %% inflow
%             
%             list_inflow=[]; clear numvol
%             numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
%             for v1=10:numvol
%                 list_inflow{v1-9,1} = [path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii,' num2str(v1)];
%             end
%             
%             clear matlabbatch
%             spm_jobman('initcfg')
%             matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_inflow};
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
%             matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
%             matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
%             matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
%             matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
%             matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
%             matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
%             spm_jobman('run',matlabbatch)
%             
%             
%             mkdir([path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D'])
%             
%             % unpack 4d for registration
%             clear matlabbatch
%             spm_jobman('initcfg')
%             matlabbatch{1}.spm.util.split.vol = {[path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']};
%             matlabbatch{1}.spm.util.split.outdir = {[path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D']};
%             spm_jobman('run',matlabbatch)
%             
%             % register meanPET to T1
%             eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/mean' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii -f '...
%                 path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlowtoT1pt1_' ])
%             % apply transformations to the images
%             numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
%             for v2=1:numvol
%                 eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlowtoT1pt1_0GenericAffine.mat -i '...
%                     path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '_' sprintf('%05i',v2) '.nii -r '...
%                     path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1_' sprintf('%05i',v2) '.nii'])
%             end
%             
%             % assemble again
%             list_inflow=[];
%             numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
%             for v1=1:numvol
%                 list_inflow{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1_' sprintf('%05i',v1) '.nii,1'];
%             end
%             clear matlabbatch
%             spm_jobman('initcfg')
%             matlabbatch{1}.spm.util.cat.vols = list_inflow;
%             matlabbatch{1}.spm.util.cat.name = ['coreg_InFlow' num2str(d) '_on_T1.nii'];
%             matlabbatch{1}.spm.util.cat.dtype = 4;
%             spm_jobman('run',matlabbatch)
%             
%             movefile([path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1.nii'],...
%                 [path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlow' num2str(d) '_on_T1.nii'])
%             rmdir([path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D'],'s')
            
            %% realign
            
            list_inflow=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow_relative.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow_relative.nii,'  num2str(v1)];
            end
            
            list_bsl=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl_relative.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl_relative.nii,'  num2str(v1)];
            end
            
            list_task=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task_relative.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task_relative.nii,'  num2str(v1)];
            end
            
            clear matlabbatch
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_inflow};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            matlabbatch{2}.spm.spatial.realign.estwrite.data = {list_task};
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 1;
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            matlabbatch{3}.spm.spatial.realign.estwrite.data = {list_bsl};
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 1;
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

            
            spm_jobman('run',matlabbatch)

            %% inflow
            
            list_inflow=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rinflow_relative.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rinflow_relative.nii,'  num2str(v1)];
            end
                        
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rinflow_relative.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/meaninflow_relative.nii -f '...
                path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Inflow_relative_toT1pt1_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Inflow_relative_toT1pt1_0GenericAffine.mat -i '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/rinflow_relative_' sprintf('%05i',v2) '.nii -r '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/inflow3D/crinflow_relative_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_inflow=[];
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/inflow3D/crinflow_relative_' sprintf('%05i',v1) '.nii,1'];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_inflow;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Inflow' num2str(d) '_relative_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/inflow3D/coreg_Inflow' num2str(d) '_relative_on_T1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Inflow' num2str(d) '_relative_on_T1.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D'],'s')

            %% baseline
            list_bsl=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rbsl_relative.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rbsl_relative.nii,'  num2str(v1)];
            end
                        
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rbsl_relative.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/meanbsl_relative.nii -f '...
                path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Baseline_relative_toT1pt2_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Baseline_relative_toT1pt2_0GenericAffine.mat -i '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/rbsl_relative_' sprintf('%05i',v2) '.nii -r '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/bsl3D/crbsl_relative_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_bsl=[];
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/bsl3D/crbsl_relative_' sprintf('%05i',v1) '.nii,1'];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_bsl;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Baseline' num2str(d) '_relative_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/bsl3D/coreg_Baseline' num2str(d) '_relative_on_T1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Baseline' num2str(d) '_relative_on_T1.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D'],'s')
            
            
            %% task
            list_task=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rtask_relative.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rtask_absolute.nii,'  num2str(v1)];
            end
                        
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rtask_relative.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/meantask_relative.nii -f '...
                path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Task_relative_toT1pt2_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Task_relative_toT1pt2_0GenericAffine.mat -i '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/rtask_relative_' sprintf('%05i',v2) '.nii -r '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/task3D/crtask_relative_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_task=[];
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/task3D/crtask_relative_' sprintf('%05i',v1) '.nii,1'];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_task;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Task' num2str(d) '_relative_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/task3D/coreg_Task' num2str(d) '_relative_on_T1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Task' num2str(d) '_relative_on_T1.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D'],'s')
            
            
        end
    end
    
end


%% write out TACs

clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask1   =[ path_parent 'segs/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii'];
            
            % 4-D PET file
%             PETflow = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/coreg_InFlow' num2str(d) '_relative_on_T1.nii'];
            PETflow = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/s3_coreg_InFlow' num2str(d) '_relative_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask1);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
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
            
            
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
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
%             mkdir(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d)])
            save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_rel_s3.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData            
%             save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_rel.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData            
        end
    end
end
%

clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask2   =[ path_parent 'segs/' num2str(IDs(id)) num2str(days(id,d)) '/aparc+aseg_pt2_nat_labelled.nii']; % task and the baseline
            
            % 4-D PET file
%             PETtask = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/coreg_Task' num2str(d) '_relative_on_T1.nii'];
%             PETbsl  = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/coreg_Baseline' num2str(d) '_relative_on_T1.nii'];
            PETtask = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/s3_coreg_Task' num2str(d) '_relative_on_T1.nii'];
            PETbsl  = [path_parent 'preproc/' num2str(IDs(id)) num2str(days(id,d)) '/s3_coreg_Baseline' num2str(d) '_relative_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask2);
            ROIMask=round(ROI.img); 
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
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
            
            
            [A B ROIDef]=xlsread('/Users/alex/Documents/MATLAB/BBD_PET/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
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
%             save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_rel.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_rel_s3.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData

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
%             save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_rel.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            save([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_rel_s3.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData

        
        end
    end
end

%% originals

close all

Appointments = [1 2 2 1 1 2 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_rel.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_rel.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_rel.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    % bilaterals
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac);
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac);
                    
                    % lefts
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac);
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac);
                    
                    % rights
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac);
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac);
                end
            end
            
            save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_uncorr_rel.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
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
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_uncorrected_rel.pdf']);
                
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
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_uncorrected_rel.pdf']);

                
            end
            
            
        end
    end
end


%% retro correct the decay

close all

Appointments = [1 2 2 1 1 2 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
%             load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_rel.mat']);
%             TACDATA_InFlow=TACDATA; clear TACDATA
%             load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_rel.mat']);
%             TACDATA_Baseline=TACDATA; clear TACDATA
%             load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_rel.mat']);
%             TACDATA_Task=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_rel_s3.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_rel_s3.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_rel_s3.mat']);
            TACDATA_Task=TACDATA; clear TACDATA

            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
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
            
%             save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_rel.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_rel_s3.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
%             try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_rel_s3.pdf']);
                
%             catch
%                 
%                 Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
%                 Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
%                 Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
%                 tmid=mean(Times,2)/60;
%                 
%                 % now draw
%                 figure('Renderer', 'painters ')
%                 plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
%                 xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
%                 legend('Cerebellum','Putamen','Caudate');
%                 title([ num2str(IDs(id)) '_' num2str(d) ]);
%                 ax = gca; ax.YAxis.Exponent = 0;
%                 print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_abs.pdf']);
% 
%                 
%             ends
            
            
        end
    end
end



%% modelling

clear; clc; close all
warning('off','all');

% set paths
paths = [];
paths.parent  = '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/';
paths.TACs    = '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/';
paths.TACs_new= '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/tacs/';
paths.ROImask = '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/segs/';
paths.figure  = '/Users/alex/Dropbox/paperwriting/MRPET/figures/modellingFigures/';
paths.rawImg  = '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/';
paths.exports = '/Users/alex/Dropbox/paperwriting/MRPET/data/';

addpath(genpath('/Users/alex/Dropbox/paperwriting/MRPET/scripts/modelling'))

% IDs
IDs         = [4001 4002 4004];
days        = [1 2; 1 2; 1 2];

%% absolute
%% load TACs for modelling

for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            TACs{i1,d} = [];
        else
            TACs{i1,d} = load([paths.TACs_new num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr_abs_s3.mat']);
        end
    end
end


%% run the model: condensed ROIs

set(0, 'DefaultLineLineWidth', 1);

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg={'CerC','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA',...
                        'caudalanteriorcingulate','entorhinal','inferiortemporal','parahippocampal','precuneus',...
                        'rostralanteriorcingulate','superiortemporal','temporalpole','middletemporal','insula'}
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
            end
            
            Lengths=[60*ones(60,1)]; % time windows
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
            if IDs(id)==4032 && d==2
            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3(2:11,:)+115*60];
            else
            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
            end
            
            tmid=mean(times,2);
            tmidMin=tmid/60;
            t_points    = length(tmid);
            dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
            break_point=find(times(:,1)>=115*60,1,'first'); %% Time of activation start
            
            %%%%%%%%%%%
            PlotStrFit=1;
            %%%%%%%%%%%
            
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg={'Striatum','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA',...
                        'caudalanteriorcingulate','entorhinal','inferiortemporal','parahippocampal','precuneus',...
                        'rostralanteriorcingulate','superiortemporal','temporalpole','middletemporal','insula','CerC'}
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
                            ROItac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                        case 'Caud'
                            ROItac=TACDATA.Caud;
                            savestr = 'Caud';
                        case 'Put'
                            ROItac=TACDATA.Put;
                            savestr = 'Put';
                        case 'Nacc'
                            ROItac=TACDATA.Nac;
                            savestr = 'Nac';
                        case 'ThalProp'
                            ROItac=TACDATA.ThalProp;
                            savestr = 'ThalProp';
                        case 'Hipp'
                            ROItac=TACDATA.Hipp;
                            savestr = 'Hipp';
                        case 'Amy'
                            ROItac=TACDATA.Amy;
                            savestr = 'Amy';
                        case 'Pall'
                            ROItac=TACDATA.Pall;
                        case 'SN'
                            ROItac=TACDATA.SN;
                            savestr = 'SN';
                        case 'LC'
                            ROItac=TACDATA.LC;
                            savestr = 'LC';
                        case 'VTA'
                            ROItac=TACDATA.VTA;
                            savestr = 'VTA';    
                        case 'caudalanteriorcingulate'
                            ROItac=TACDATA.caudalanteriorcingulate;
                            savestr = 'caudalanteriorcingulate';
                        case 'entorhinal'
                            ROItac=TACDATA.entorhinal;
                            savestr = 'entorhinal';
                        case 'inferiortemporal'
                            ROItac=TACDATA.inferiortemporal;
                            savestr = 'inferiortemporal';
                        case 'parahippocampal'
                            ROItac=TACDATA.parahippocampal;
                            savestr = 'parahippocampal';
                        case 'precuneus'
                            ROItac=TACDATA.precuneus;
                            savestr = 'precuneus';
                        case 'rostralanteriorcingulate'
                            ROItac=TACDATA.rostralanteriorcingulate;
                            savestr = 'rostralanteriorcingulate';
                        case 'superiortemporal'
                            ROItac=TACDATA.superiortemporal;
                            savestr = 'superiortemporal';
                        case 'temporalpole'
                            ROItac=TACDATA.temporalpole;
                            savestr = 'temporalpole';
                        case 'middletemporal'
                            ROItac=TACDATA.middletemporal;
                            savestr = 'middletemporal';
                        case 'insula'
                            ROItac=TACDATA.insula;
                            savestr = 'insula';
                        case 'CerC'
                            ROItac=TACDATA.CerC;
                            savestr = 'CerC';    
                    end
                    
                    mROItac  = [ROItac(1)/2; (ROItac(2:end)+ROItac(1:end-1))/2];
                    ASRTM(:,3)=zeros(t_points,1);
%                     ASTRM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mROItac(1:k).*dt(1:k));
                    end
                    %LSQ-estimation using lscov
                    [parest se_srtm mse_srtm]   = lscov(ASRTM,ROItac); % parest(1)->R1=K1/K1', parest(2)->V_T, parest(3)->V_ND
                    fittac=ASRTM*parest;
                    BP=parest(2)/parest(3)-1;
                    k2p=parest(2)/parest(1); % k2 of the reference region
                    
                    %%%%% DO real SRTM
                    options = optimset('MaxFunEvals',1000);
                    weighs=[0.25*ones(30,1); ones(t_points-30,1)];
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-ROItac).*weighs);
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-ROItac(1:end-11)).*weighs(1:end-11));
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1_bl=parest_srtm(1);
                    k2_bl=parest_srtm(2);
                    BP_bl=parest_srtm(3);
                    modfit_esrtm_bl=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
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
                                    roitac_gamma = ROItac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(ROItac))+roitac_gamma(1:length(ROItac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,ROItac);
                                    
                                    %estimated model TAC
                                    modelfit = Alpntpet*parest;
                                    
                                    %%%%%% magic %%%%%
                                    if (best_mse > mse)
                                        best_mse = mse;
                                        best_parest=parest; % R1 - k2 - k2a - gamma
                                        best_modelfit=modelfit;
                                        best_actfun=actfun;
                                        best_se=se;
                                    end
                                end
                            end
        %                        breakpoint2=breakpoint;
                        end
                    end
                    
                    y=-ASRTM(:,3)./ROItac;
                    x=ASRTM(:,2)./ROItac;
                    [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                    [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                    
                    % save these and work out receptor occupancy later in
                    % the script

                    k2=best_parest(2);
                    k2a=best_parest(3);
                    BP_lp=k2/k2a-1; % baseline binding potential of lpntPET
                    BP_srtm=BP__; % binding potential until the end of the task
                    BP_srtm_bsl=BP_bl; % binding potential until the end of baseline
                    G=best_parest(4); % gamma, obviously...
                    BPND = ((R1__*k2p)/k2a)-1; % correct? - not sure
%                     DBP =  k2./(k2a + best_actfun) - 1; % dynamic binding potential
                    DBP =  k2./(k2a + G*best_actfun) - 1; % is it this way...? confused
                    OCC = 100*(1-DBP/BP_lp);% receptor occupancy
%                     OCC2 = (1-DBP2/BP_lp);% receptor occupancy
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.' savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.' savestr '=OCC;' ])
%                     eval(['Occupancy2{id,d}.' savestr '=OCC2;' ])
                    eval(['BPND_save{id,d}.' savestr '=BPND;' ])
                    eval(['BP_srtm_save{id,d}.' savestr '=BP_srtm;' ])
                    eval(['BP_srtm_Bsl_save{id,d}.' savestr '=BP_srtm_bsl;' ])
                    eval(['BestActivationFunction{id,d}.' savestr '=best_actfun;'])
                    eval(['BestModelFit{id,d}.' savestr '=best_modelfit;'])
                    eval(['gamma_lp_save{id,d}.' savestr '=best_parest(4);' ])
                    eval(['Residuals{id,d}.' savestr '.BaselineFit=ROItac-fittac;' ])
                    eval(['Residuals{id,d}.' savestr '.AllFit_ESRTM=ROItac-modfit_esrtm;' ])
                    eval(['Residuals{id,d}.' savestr '.lpntPETFit=ROItac-best_modelfit;' ])
                    eval(['CompensatoryFunction{id,d}.' savestr '=best_parest(4)*best_actfun;' ])
                    
                    
                    % Compute residuals
%                     residuals_lpnt = ROItac-best_modelfit;
%                     
%                     % Apply weights
%                     Weight = [ones(1,100) ones(1,11).*2];
%                     weights = 1./sqrt(Weight);
%                     weighted_resid_lpnt = residuals_lpnt .* weights;
%                     
%                     % Square and sum weighted residuals to obtain WRSS
%                     WRSS_lpnt = sum(weighted_resid_lpnt.^2);
                    

%                     BP_lp=best_parest(2)/best_parest(3)-1; % is this the BP from lpntPET?
%                     BP_lpntPET=(best_parest(2)/(best_parest(3)+(best_parest(4)).*best_actfun))-1; % isn't this the way?
%                     BP_bsl=BP_bl; % baseline from SRTM
                    % save BP for later analysis
%                     eval(['BP_lp_orig{id,d}.' savestr '=BP_lp;' ])
%                     eval(['BP_lp_save{id,d}.' savestr '=BP_lpntPET;' ])
%                     eval(['BP_bsl_save{id,d}.' savestr '=BP_bsl;' ])
                    
                    
                    BPdata.BP_mrtm(Subj{r})=BP;
                    BPdata.BP_srtm(Subj{r})=BP__;
                    BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                    BPdata.BP_lpnt(Subj{r})=BP_lp;
                    BPdata.BP_logan(Subj{r})=pp(1)-1;
                    
                    % save
                    BPdataSave{id,d}.(reg{1}).BP_mrtm=BP;
                    BPdataSave{id,d}.(reg{1}).BP_srtm=BP__;
                    BPdataSave{id,d}.(reg{1}).BP_srtm_bl=BP_bl;
                    BPdataSave{id,d}.(reg{1}).BP_lpnt=BP_lp;
                    BPdataSave{id,d}.(reg{1}).BP_logan=pp(1)-1;
                    BPdataSave{id,d}.(reg{1}).BPND=BPND;
                    BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
%                     BPdataSave{id,d}.(reg{1}).Occupancy2=OCC2;
                    BPdataSave{id,d}.(reg{1}).DBP=DBP;
                    BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);


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
                        plot(tmidMin,ROItac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
%                         xlim([0 60]);
                        xlabel('Time (min)');
                        ylabel('Radioactivity concentration');
                        legend(legendtext,'Location','best');
                        yyaxis right;
                        plot(tmidMin,best_parest(4)*best_actfun,'c-'); % plot gamma
                        ylabel('Compensatory function','Color','c');
                        ylim([0 4*10^(-4)]);
                        title([Subj{r} ' compartmental fits: BP_{Baseline}=' num2str(BP_bl,'%1.2f') ', BP_{All}=' num2str(BP__,'%1.2f') ', BP_{lp-nt}=' num2str(BP_lp,'%1.2f')]);
                        subplot(3,2,3);
                        plot(tmidMin,ROItac-fittac,'bo',tmidMin,ROItac-modfit_esrtm,'go',tmidMin,ROItac-modfit_esrtm_bl,'co',tmidMin,ROItac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        plot(tmidMin,ROItac-fittac,'bo',tmidMin,ROItac-modfit_esrtm,'bo',tmidMin,ROItac-modfit_esrtm_bl,'ro',tmidMin,ROItac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        [h p]=runstest(ROItac-fittac);
                        [h1 p1]=runstest(ROItac-modfit_esrtm);
                        [h2 p2]=runstest(ROItac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
                        title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
                        xlabel('Time (min)');

                        subplot(3,2,4);
                        y=-ASRTM(:,3)./ROItac;
                        x=ASRTM(:,2)./ROItac;
                        [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                        [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                        plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
                        title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        xlabel(['\int REF/ROI']);
                        ylabel('\int ROI/ROI')
                        subplot(3,2,[5 6]);
                        plot(tmidMin,ROItac./reftac,'ko',[0 180],[0 0],'k--');
                        ylim([0.5 30]); xlim([0 190])
                        title('Target to reference ratio');
                        xlabel('Time (min)');

                        print('-dpsc2','-append','-bestfit',fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_TAC_AbsoluteData_Fit_lpntpet_logan_' date '_s3.ps']));
                        close(gcf)
                        BPdata.BP_mrtm(Subj{r})=BP;
                        BPdata.BP_srtm(Subj{r})=BP__;
                        BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                        BPdata.BP_lpnt(Subj{r})=BP_lp;
                        BPdata.BP_logan(Subj{r})=pp(1)-1;
                        
                        % save
                        BPdataSave{id,d}.(reg{1}).BP_mrtm=BP;
                        BPdataSave{id,d}.(reg{1}).BP_srtm=BP__;
                        BPdataSave{id,d}.(reg{1}).BP_srtm_bl=BP_bl;
                        BPdataSave{id,d}.(reg{1}).BP_lpnt=BP_lp;
                        BPdataSave{id,d}.(reg{1}).BP_logan=pp(1)-1;
                        BPdataSave{id,d}.(reg{1}).BPND=BPND;
                        BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
%                         BPdataSave{id,d}.(reg{1}).Occupancy2=OCC2;
                        BPdataSave{id,d}.(reg{1}).DBP=DBP;
%                         BPdataSave{id,d}.(reg{1}).DBP2=DBP2;
                        BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);
                        
                        
                        continue;
                    end
                end
                
                
            end
            close all
            keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI BPdataSave BP_lp_save DBP_save Occupancy BPND_save BP_srtm_save BP_srtm_Bsl_save BestActivationFunction...
    BestModelFit gamma_lp_save Residuals CompensatoryFunction
        end
    end
end

disp('done')

save(['/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/TACs/MRPET_BPpackage_condensed_' date '_absolute_s3.mat'],...
    'BPdataSave', 'BP_lp_save', 'DBP_save', 'Occupancy', 'BPND_save', 'BP_srtm_save', 'BP_srtm_Bsl_save', 'BestActivationFunction',...
    'BestModelFit', 'gamma_lp_save', 'Residuals', 'CompensatoryFunction')








