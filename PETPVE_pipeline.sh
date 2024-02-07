#!/bin/bash	


export FSLDIR=/usr/local/fsl
export PATH=$PATH:/usr/local/fsl/bin
source $FSLDIR/etc/fslconf/fsl.sh

for K in  40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40121
do

	# subject ID
	ID=${K}

	# set up folders and group space images
	folder_pet=/Users/yyi/Desktop/new_recon/preproc/"${ID}"
	folder_pvc=/Users/yyi/Desktop/PVEc/"${ID}"

	#mkdir "${folder_pet}"/inflow3D
	#mkdir "${folder_pet}"/bsl3D
	#mkdir "${folder_pet}"/task3D

	#pvc_make4d -i "${folder_pet}"/aparc+aseg_pt2_nat_labelled.nii -o "${folder_pet}"/aparc+aseg_pt2_nat_labelled_4D.nii.gz
	#pvc_make4d -i "${folder_pet}"/aparc+aseg_pt1_nat_labelled.nii -o "${folder_pet}"/aparc+aseg_pt1_nat_labelled_4D.nii.gz

	#fslsplit "${folder_pet}"/task_on_t1pt2.nii "${folder_pet}"/task/task_on_t1pt2_
	#fslsplit "${folder_pet}"/baseline_on_t1pt2.nii "${folder_pet}"/baseline/baseline_on_t1pt2_
	#fslsplit "${folder_pet}"/inflow_on_t1pt1.nii "${folder_pet}"/inflow/inflow_on_t1pt1_


	for I in {1..15}
	do
		wtf=$(printf "%05d" "$I")
		petpvc -i "${folder_pet}"/bsl3D/bsl_on_t1pt2_"${wtf}".nii.gz -m "${folder_pvc}"/aparc+aseg_pt2_nat_labelled_4D.nii.gz -o "${folder_pet}"/bsl3D/PVEc_RBV_baseline_on_t1pt2_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
	done

	fslmerge -t "${folder_pet}"/PVEc_RBV_baseline_on_t1pt2.nii.gz "${folder_pet}"/bsl3D/PVEc_RBV_baseline_on_t1pt2*.nii.gz
	rm "${folder_pet}"/bsl3D/PVEc_RBV_baseline_on_t1pt2*.nii.gz

	for I in {1..55}
	do
		wtf=$(printf "%05d" "$I")
		petpvc -i "${folder_pet}"/task3D/task_on_t1pt2_"${wtf}".nii.gz -m "${folder_pvc}"/aparc+aseg_pt2_nat_labelled_4D.nii.gz -o "${folder_pet}"/task3D/PVEc_RBV_task_on_t1pt2_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
	done

	fslmerge -t "${folder_pet}"/PVEc_RBV_task_on_t1pt2.nii.gz "${folder_pet}"/task3D/PVEc_RBV_task_on_t1pt2*.nii.gz
	rm "${folder_pet}"/task3D/PVEc_RBV_task_on_t1pt2*.nii.gz



	for I in {1..60}
	do
		wtf=$(printf "%05d" "$I")
		echo $wtf
		petpvc -i "${folder_pet}"/inflow3D/inflow_on_t1pt1_"${wtf}".nii.gz -m "${folder_pvc}"/aparc+aseg_pt1_nat_labelled_4D.nii.gz -o "${folder_pet}"/inflow3D/PVEc_RBV_inflow_on_t1pt1_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
	done

	fslmerge -t "${folder_pet}"/PVEc_RBV_inflow_on_t1pt1.nii.gz "${folder_pet}"/inflow3D/PVEc_RBV_inflow_on_t1pt1*.nii.gz
	rm "${folder_pet}"/inflow3D/PVEc_RBV_inflow_on_t1pt1*.nii.gz


done