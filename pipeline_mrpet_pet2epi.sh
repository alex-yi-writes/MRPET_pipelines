#!/usr/local/bin/bash

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
#antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ["${folder_trx}"data/coreg_meanEPI_to_T1WB_0GenericAffine.mat,1] -i "${folder_pet}"aparc+aseg_pt1_nat_labelled.nii -r "${folder_trx}"data/meanEPI_corrected.nii -o "${folder_pet}"aparc+aseg_pt1_EPI_labelled.nii
#antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ["${folder_trx}"data/coreg_meanEPI_to_T1WB_0GenericAffine.mat,1] -i "${folder_pet}"aparc+aseg_pt2_nat_labelled.nii -r "${folder_trx}"data/meanEPI_corrected.nii -o "${folder_pet}"aparc+aseg_pt2_EPI_labelled.nii


######### example
#antsApplyTransforms -d 3 -v 1 -n BSpline[4] -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder}"data/T1WB_corrected.nii -r "${MNI}" -o "${folder}"data/NLreg_T1WB_to_MNI.nii
########


# first, inflow
for I in {01..85}
do
  antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t /Volumes/ALEX3/MRPET/coreg_mri/NLreg_template_to_MNI_1Warp.nii.gz -t /Volumes/ALEX3/MRPET/coreg_mri/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder_trx}"data/NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder_trx}"data/NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder_pet}"3D/InFlow/coreg_InFlowonT1_000${I}.nii -r "${MNI}" -o "${folder_pet}"3D/InFlow/coreg_InFlowonEPI_000${I}.nii
done

# second, baseline
for I in {01..15}
do
  antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t /Volumes/ALEX3/MRPET/coreg_mri/NLreg_template_to_MNI_1Warp.nii.gz -t /Volumes/ALEX3/MRPET/coreg_mri/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder_trx}"data/NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder_trx}"data/NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder_pet}"3D/BSL/coreg_BaselineonT1_000${I}.nii -r "${MNI}" -o "${folder_pet}"3D/BSL/coreg_BaselineonEPI_000${I}.nii
done

# last, task
for I in {01..11}
do
  antsApplyTransforms -d 3 -v 0 -n Linear -t /Volumes/ALEX3/MRPET/coreg_mri/NLreg_template_to_MNI_1Warp.nii.gz -t /Volumes/ALEX3/MRPET/coreg_mri/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder_trx}"data/NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder_trx}"data/NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder_pet}"3D/Task/coreg_MTonT1_000${I}.nii -r "${MNI}" -o "${folder_pet}"3D/Task/coreg_MTonEPI_000${I}.nii
done

# ------------------------------------------------------ #

done
