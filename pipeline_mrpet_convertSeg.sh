#!/usr/local/bin/bash

MNI=/Volumes/ALEX3/mni_icbm152_t1_tal_nlin_asym_09c.nii
template=/Volumes/ALEX3/MRPET/mrpet-template.nii.gz

for K in 40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40121 40122 40131 40132 40142 40151 40152 40161 40162 40171 40172
do

# subject ID
ID=${K}
# set up folders and group space images
folder_origMNI=/Volumes/ALEX3/MRPET/MRPETseg/MNIspace/"${ID}"/mri
folder_origtemplate1=/Volumes/ALEX3/MRPET/MRPETseg/templateSpace/"${ID}"pt1/mri
folder_origtemplate2=/Volumes/ALEX3/MRPET/MRPETseg/templateSpace/"${ID}"pt2/mri
folder_dest=/Volumes/ALEX3/MRPET/coreg_roi/"${ID}"

#mri_convert --in_type mgz --out_type nii ${folder_origMNI}/aparc+aseg.mgz ${folder_origMNI}/aparc+aseg.nii
#mri_convert --in_type mgz --out_type nii ${folder_origMNI}/T1.mgz ${folder_origMNI}/T1.nii

#antsRegistrationSyN.sh -d 3 -n 1 -t r -m ${folder_origMNI}/T1.nii -f ${MNI} -o ${folder_dest}/FS_to_MNI_


#mri_convert --in_type mgz --out_type nii ${folder_origtemplate1}/aparc+aseg.mgz ${folder_origtemplate1}/aparc+aseg.nii
#mri_convert --in_type mgz --out_type nii ${folder_origtemplate1}/T1.mgz ${folder_origtemplate1}/T1.nii
#mri_convert --in_type mgz --out_type nii ${folder_origtemplate2}/aparc+aseg.mgz ${folder_origtemplate2}/aparc+aseg.nii
#mri_convert --in_type mgz --out_type nii ${folder_origtemplate2}/T1.mgz ${folder_origtemplate2}/T1.nii

#antsRegistrationSyN.sh -d 3 -n 1 -t r -m ${folder_origtemplate1}/T1.nii -f ${template} -o ${folder_dest}/FSpt1_to_Template_
#antsRegistrationSyN.sh -d 3 -n 1 -t r -m ${folder_origtemplate2}/T1.nii -f ${template} -o ${folder_dest}/FSpt2_to_Template_

antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder_dest}/FS_to_MNI_0GenericAffine.mat -r ${MNI} -i ${folder_origMNI}/aparc+aseg.nii -o ${folder_dest}/aparc+aseg_on_MNI.nii

antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder_dest}/FSpt1_to_Template_0GenericAffine.mat -r ${template} -i ${folder_origtemplate1}/aparc+aseg.nii -o ${folder_dest}/aparc+aseg_pt1_on_template.nii

antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ${folder_dest}/FSpt2_to_Template_0GenericAffine.mat -r ${template} -i ${folder_origtemplate2}/aparc+aseg.nii -o ${folder_dest}/aparc+aseg_pt2_on_template.nii

done
