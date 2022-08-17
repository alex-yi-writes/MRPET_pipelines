%% extract LC contrast ratio

clc;clear

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin']);
setenv('ANTSPATH','/usr/local/antsbin/bin')

load('/Users/alex/Documents/LCaverage/IDs.mat')

path_parent='/Users/alex/Documents/LCaverage/';
pons='/Users/alex/Dropbox/Masks/Masks/mni_icbm152/ponsmask_25vox.nii';

%% registration T1->LCslab averaged

for id=1:length(ID)
    
    clear FixedImage MovingImage OutputImage
    MovingImage     = [path_parent ID{id} '/T1WB_corrected.nii'];
    FixedImage      = [path_parent ID{id} '/LCslab_averaged.nii'];
    OutputPrefix    = [path_parent ID{id} '/T1toLCslab_'];
    
    % register
%      eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -f ' FixedImage ' -m ' MovingImage ' -o ' OutputPrefix])
    % transform
    eval(['!antsApplyTransforms -d 3 -v 0 -n NearestNeighbor  -t ' path_parent ID{id}...
        '/T1toLCslab_0GenericAffine.mat -t [' path_parent ID{id} '/NLreg_T1WB_to_template_0GenericAffine.mat,1] -t ' ...
        path_parent 'NLreg_template_to_MNI_1InverseWarp.nii.gz -t [' path_parent '/NLreg_template_to_MNI_0GenericAffine.mat,1] -i '...
        pons ' -r ' FixedImage ' -o ' path_parent ID{id} '/pons_native.nii'])
    
end

%% extract LC CR

for id=1:length(ID)
    
    clear LC S_LC_mean S_pons_mean S_LC_max S_pons_max
    
    % load images
    LC          = spm_read_vols(spm_vol([path_parent ID{id} '/LCmask.nii']));
    pons_nat    = spm_read_vols(spm_vol([path_parent ID{id} '/pons_native.nii']));
    LCslab      = spm_read_vols(spm_vol([path_parent ID{id} '/LCslab_averaged.nii']));
    
    % define masks
    [x_lc,y_lc,z_lc]= ind2sub(size(LC),find(LC~=0));
    [x_p,y_p,z_p]   = ind2sub(size(pons_nat),find(pons_nat~=0));
    coord_LC        = [x_lc,y_lc,z_lc];
    coord_pons      = [x_p,y_p,z_p];
    slices          = unique(z_lc); % identify slices in LC masks
    
    % extract voxel values from each mask, in each slice
    for sl=1:length(slices)
        
        clear tmp_lc voxvals_lc tmp_p voxvals_p
        
        % LC
        tmp_lc  = coord_LC(z_lc==slices(sl),:);
        tmp_p   = coord_pons(z_p==slices(sl),:);
        for vox=1:size(tmp_lc,1)
            voxvals_lc(vox,1)=LCslab(tmp_lc(vox,1),tmp_lc(vox,2),tmp_lc(vox,3));
        end
        S_LC_mean(sl,1) = mean(voxvals_lc(:));
        S_LC_max(sl,1)  = max(voxvals_lc(:));
        
        % pons
        for vox=1:size(tmp_p,1)
            voxvals_pons(vox,1)=LCslab(tmp_p(vox,1),tmp_p(vox,2),tmp_p(vox,3));
        end
        S_pons_mean(sl,1) = mean(voxvals_pons(:));
        S_pons_max(sl,1)  = max(voxvals_pons(:));
        
    end
    
    CR_mean(id,1) = mean((S_LC_mean - S_pons_mean)./ S_pons_mean);
    CR_max(id,1)  = max((S_LC_max - S_pons_max)./ S_pons_max);
    
end
