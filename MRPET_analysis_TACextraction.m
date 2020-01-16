%% Define input parameters
clc;clear;
addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')

% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
TrimPerc=15;
% Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
Mask='/Users/yeojin/Desktop/E_data/EA_raw/EAB_MRI/EABD_segmented/4001_2/40012/mri/aparc+aseg.nii';
% Mask=['/Users/yeojin/Desktop/E_data/EA_raw/EAB_MRI/EABD_segmented/4001_2/40012/mri/aparc+aseg.nii'];
% 4-D PET file
CurPET='/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/4001_2/cr4001_PET_4D_MT2.nii'; 

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
DynPET=load_nii(CurPET);
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
save('TACDATA_40012.mat','TACDATA');
