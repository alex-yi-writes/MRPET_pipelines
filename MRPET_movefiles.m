%% move around files for the MRPET dataset

%% preparation

close all;clc;clear;
warning('off','all');

% paths
% paths
paths = [];
paths.parent      = '/Users/yeojin/Desktop/';
paths.raw         = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.source      = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.destination = [paths.parent 'E_data/EA_raw/EAB_MRI/EABD_segmented/transmat/'];

% IDs
IDs = [4001 4002 4003 4004 4005 4006 4007];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0]; 

%% move T1 whole-brain images for FSL segmentation
% info: ID # # # #   #       #
%          <-ID->   day      pt
%          4=OAS,   1=big    1=inflow 
%          3=YAs,   2=small  2=bsl,task


for id = 1:length(IDs)
    
     for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            for pts = 1:2
                
                clear src destination
                % from
                src = [paths.source num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) ...
                    '_MRI_4D_MPRAGE' num2str(d) '_pt' num2str(pts) '.nii'];
                % to
                mkdir([paths.destination num2str(IDs(id)) num2str(d) num2str(pts)])
                destination = [paths.destination num2str(IDs(id)) num2str(d) num2str(pts) '/T1WB.nii'];
                
                copyfile(src,destination)
                
            end
            
        end
     end
end


%% write bash scripts

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            for pts = 1:2
                
                fileID = fopen([num2str(IDs(id)) num2str(d) num2str(pts) '.sh'],'w');
                fprintf(fileID,'#!/bin/bash \n\n');
                fprintf(fileID,['recon-all -i /mnt/work/yyi/temp/parcellation/MRPET/' num2str(IDs(id)) num2str(d) num2str(pts)...
                    '/T1WB.nii -s ' num2str(IDs(id)) num2str(d) num2str(pts) ' -sd /mnt/work/yyi/temp/parcellation/MRPET/'...
                    num2str(IDs(id)) num2str(d) num2str(pts) ' -all \n\n']);
                fprintf(fileID,['mri_convert --in_type mgz --out_type nii --out_orientation RAS /mnt/work/yyi/temp/parcellation/MRPET/' ...
                    num2str(IDs(id)) num2str(d) num2str(pts) '/mri/aparc+aseg.mgz /mnt/work/yyi/temp/parcellation/MRPET/' ...
                    num2str(IDs(id)) num2str(d) num2str(pts) '/mri/aparc+aseg.nii \n']);
                fclose(fileID);
                
            end
            
        end
    end
end
