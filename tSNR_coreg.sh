#!/bin/bash
# coregister tSNR betas

# the subject folders
path=/mnt/work/yyi/temp/ED_coreg/pilot_RT

# MNI template with filename
MNI=/mnt/work/yyi/temp/ED_coreg/mni_icbm152_t1_tal_nlin_asym_09c.nii

# template with filename
template=/mnt/work/yyi/temp/ED_coreg/hc_template.nii.gz

# assign number of CPUs used
CPUs=3

# subject folders in ${path}
ls -d "${path}"/*/ > "${path}"/participants_3T.txt

while read folder; do

	# subject ID
	ID=$(echo "${folder}" | grep -o -E '[0-9][0-9][0-9][0-9][0-9]')

	# move contrast mask
	antsApplyTransforms -d 3 -v 1 -n Linear -t "${path}"/NLreg_template_to_MNI_1Warp.nii.gz -t "${path}"/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1mean_to_template_0GenericAffine.mat -t ["${folder}"data/coreg_meanEPI3D1_to_T1mean_0GenericAffine.mat, 1] -i "${folder}"data/tSNR_${ID}_beta_0013.nii -r "${MNI}" -o "${folder}"data/ctSNR_${ID}_beta_0013.nii

	# move contrast mask
	antsApplyTransforms -d 3 -v 1 -n Linear -t "${path}"/NLreg_template_to_MNI_1Warp.nii.gz -t "${path}"/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1mean_to_template_0GenericAffine.mat -t ["${folder}"data/coreg_meanEPI3D1_to_T1mean_0GenericAffine.mat, 1] -i "${folder}"data/con_0016.nii -r "${MNI}" -o "${folder}"data/con_0016_mni.nii

	cp "${folder}"data/ctSNR_${ID}_beta_0013.nii /mnt/work/yyi/temp/ED_coreg/tSNR/ctSNR_${ID}_beta_0013.nii

  echo "ID ${ID} DONE!"

done < "${path}"/participants_3T.txt
