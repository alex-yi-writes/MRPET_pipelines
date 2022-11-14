#!/usr/local/bin/bash


# actually it's MNI to T1 now
MNI=/Volumes/ALEX3/mni_icbm152_t1_tal_nlin_asym_09c.nii
template=/Users/alex/Documents/pilot_template.nii.gz

for K in 40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40121 40122 40131 40132 40142 40151 40152 40161 40162 #40171 40172
do

# subject ID
ID=${K}
# set up folders and group space images
folder_trx=/Volumes/ALEX3/MRPET/coreg_mri/"${ID}"/
folder_pet=/Volumes/ALEX3/MRPET/coreg_roi/"${ID}"/

# -------------- start the transformation -------------- #

# example T1
#antsApplyTransforms -d 3 -v 0 -n Linear -t ["${folder_trx}"data/coreg_meanEPI_to_T1WB_0GenericAffine.mat,1] -i "${folder_trx}"data/T1WB_corrected.nii -r "${folder_trx}"data/meanEPI_corrected.nii -o "${folder_pet}"3D/InFlow/T1_on_EPI.nii

# segmentations
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t "${folder_pet}"data/NLreg_MNI2T1pt2_1Warp.nii.gz -t "${folder_pet}"data/NLreg_MNI2T1pt2_0GenericAffine.mat -r "${folder_pet}"data/T1pt2.nii -o "${folder_pet}"fMRIclusters_on_t1pt2.nii -i /Volumes/ALEX3/MRPET/coreg_roi/FB_rew_v_neu_p005_c492_noMask.nii

antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t "${folder_pet}"data/NLreg_MNI2T1pt1_1Warp.nii.gz -t "${folder_pet}"data/NLreg_MNI2T1pt1_0GenericAffine.mat -r "${folder_pet}"data/T1pt1.nii -o "${folder_pet}"fMRIclusters_on_t1pt1.nii -i /Volumes/ALEX3/MRPET/coreg_roi/FB_rew_v_neu_p005_c492_noMask.nii


#antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t "${folder_pet}"data/NLreg_MNI2T1pt1_1InverseWarp.nii.gz -t ["${folder_pet}"data/NLreg_MNI2T1pt1_0GenericAffine.mat,1] -i "${folder_pet}"aparc+aseg_pt1_nat_labelled.nii -r /Volumes/ALEX3/mni_icbm152_t1_tal_nlin_asym_09c.nii -o "${folder_pet}"aparc+aseg_pt1_MNI_labelled.nii
#antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t "${folder_pet}"data/NLreg_MNI2T1pt2_1InverseWarp.nii.gz -t ["${folder_pet}"data/NLreg_MNI2T1pt2_0GenericAffine.mat,1] -i "${folder_pet}"aparc+aseg_pt2_nat_labelled.nii -r /Volumes/ALEX3/mni_icbm152_t1_tal_nlin_asym_09c.nii -o "${folder_pet}"aparc+aseg_pt2_MNI_labelled.nii


######### example
#antsApplyTransforms -d 3 -v 1 -n BSpline[4] -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder}"data/T1WB_corrected.nii -r "${MNI}" -o "${folder}"data/NLreg_T1WB_to_MNI.nii
########

# ------------------------------------------------------ #

done
