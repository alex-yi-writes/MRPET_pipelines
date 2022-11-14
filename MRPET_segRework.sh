#!/bin/bash

CPUs=1
for K in 40112 #40171 40172 40181 40182 40191 40201 40202 40211 40212 40222 40231 40241 40262 # 40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092

do

ID=${K}

folder=/Users/yeojin/Desktop/E_data/ED_coreg/mrpet/${ID}/data/

MNI=/Users/yeojin/Desktop/E_data/EE_atlases_templates/mni_icbm152_t1_tal_nlin_asym_09c.nii
LC=/Users/yeojin/Desktop/E_data/EE_atlases_templates/mni_icbm152_LCmetaMask_final.nii
SN=/Users/yeojin/Desktop/E_data/EE_atlases_templates/mni_icbm152_SN.nii
VTA=/Users/yeojin/Desktop/E_data/EE_atlases_templates/mni_icbm152_VTA.nii

T1pt1=${folder}T1pt1.nii
T1pt2=${folder}T1pt2.nii

# register MNI to T1s
antsRegistrationSyNQuick.sh -d 3 -t s -f $T1pt1 -m $MNI -o ${folder}NLreg_MNI2T1pt1_
antsRegistrationSyNQuick.sh -d 3 -t s -f $T1pt2 -m $MNI -o ${folder}NLreg_MNI2T1pt2_

# LC, pt1
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder}NLreg_MNI2T1pt1_1Warp.nii.gz -t ${folder}NLreg_MNI2T1pt1_0GenericAffine.mat -i "${LC}" -r "${T1pt1}" -o "${folder}"LC_onT1pt1.nii
# LC, pt2
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder}NLreg_MNI2T1pt2_1Warp.nii.gz -t ${folder}NLreg_MNI2T1pt2_0GenericAffine.mat -i "${LC}" -r "${T1pt2}" -o "${folder}"LC_onT1pt2.nii

# SN, pt1
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder}NLreg_MNI2T1pt1_1Warp.nii.gz -t ${folder}NLreg_MNI2T1pt1_0GenericAffine.mat -i "${SN}" -r "${T1pt1}" -o "${folder}"SN_onT1pt1.nii
# SN, pt2
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder}NLreg_MNI2T1pt2_1Warp.nii.gz -t ${folder}NLreg_MNI2T1pt2_0GenericAffine.mat -i "${SN}" -r "${T1pt2}" -o "${folder}"SN_onT1pt2.nii

# VTA, pt1
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder}NLreg_MNI2T1pt1_1Warp.nii.gz -t ${folder}NLreg_MNI2T1pt1_0GenericAffine.mat -i "${VTA}" -r "${T1pt1}" -o "${folder}"VTA_onT1pt1.nii
# VTA, pt2
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder}NLreg_MNI2T1pt2_1Warp.nii.gz -t ${folder}NLreg_MNI2T1pt2_0GenericAffine.mat -i "${VTA}" -r "${T1pt2}" -o "${folder}"VTA_onT1pt2.nii

echo "${ID} done"

done
