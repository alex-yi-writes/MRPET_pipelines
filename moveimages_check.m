%% moving PET files for checks!

clc;clear;

% IDs
IDs = [4001 4002 4003 4004 4005 4006 4007];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0]; 

% paths
destpath    = '/Users/yeojin/Dropbox/PET_checks/';
imgpath     = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/';
segpath     = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/';

%% loops

for i = 1:length(IDs)
    
    for d = 1:2
        
        if days(i,d) == 0
        else
            
            mkdir([destpath num2str(IDs(i)) num2str(d)])
            
            % move PET images
            copyfile([imgpath num2str(IDs(i)) '_' num2str(d) '/cmean' num2str(IDs(i)) '_PET_4D_InFlow' num2str(d) '.nii'],...
                [destpath num2str(IDs(i)) num2str(d) '/cmean' num2str(IDs(i)) '_PET_4D_InFlow' num2str(d) '.nii'])
            copyfile([imgpath num2str(IDs(i)) '_' num2str(d) '/cmean' num2str(IDs(i)) '_PET_4D_Baseline' num2str(d) '.nii'],...
                [destpath num2str(IDs(i)) num2str(d) '/cmean' num2str(IDs(i)) '_PET_4D_Baseline' num2str(d) '.nii'])
            copyfile([imgpath num2str(IDs(i)) '_' num2str(d) '/cmean' num2str(IDs(i)) '_PET_4D_MT' num2str(d) '.nii'],...
                [destpath num2str(IDs(i)) num2str(d) '/cmean' num2str(IDs(i)) '_PET_4D_MT' num2str(d) '.nii'])
            
            % move T1 images
            copyfile([imgpath num2str(IDs(i)) '_' num2str(d) '/seg_' num2str(IDs(i)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii'],...
                [destpath num2str(IDs(i)) num2str(d) '/seg_' num2str(IDs(i)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii'])
            copyfile([imgpath num2str(IDs(i)) '_' num2str(d) '/seg_' num2str(IDs(i)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
                [destpath num2str(IDs(i)) num2str(d) '/seg_' num2str(IDs(i)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'])
            
            % move segmentation images
            for seg = 1:2
                copyfile([segpath num2str(IDs(i)) num2str(d) num2str(seg) '/mri/aparc+aseg.nii'],...
                    [destpath num2str(IDs(i)) num2str(d) '/segmentation_pt' num2str(seg) '.nii'])
            end
            
        end
        
    end
    
end