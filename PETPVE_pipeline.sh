#!/bin/bash	


export FSLDIR=/usr/local/fsl
export PATH=$PATH:/usr/local/fsl/bin
source $FSLDIR/etc/fslconf/fsl.sh

for K in  40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40112 40121 40122 40131 40132 40141 40142 40151 40152 40161 40162 40171 40172 40181 40182 40191 40192 40201 40202 40211 40212 40221 40222 40231 40232 40241 40242 40251 40252 40261 40262 40271 40281 40282 40292 40301 40302 40311 40312 40321 40322 40331 40332 40011 40012
do

	# subject ID
	ID=${K}

	# set up folders and group space images
	folder_pet=/Users/yyi/Desktop/PVEc/"${ID}"

	mkdir "${folder_pet}"/inflow
	#mkdir "${folder_pet}"/baseline
	#mkdir "${folder_pet}"/task

	#pvc_make4d -i "${folder_pet}"/aparc+aseg_pt2_nat_labelled.nii -o "${folder_pet}"/aparc+aseg_pt2_nat_labelled_4D.nii.gz
	pvc_make4d -i "${folder_pet}"/aparc+aseg_pt1_nat_labelled.nii -o "${folder_pet}"/aparc+aseg_pt1_nat_labelled_4D.nii.gz

	#fslsplit "${folder_pet}"/task_on_t1pt2.nii "${folder_pet}"/task/task_on_t1pt2_
	#fslsplit "${folder_pet}"/baseline_on_t1pt2.nii "${folder_pet}"/baseline/baseline_on_t1pt2_
	fslsplit "${folder_pet}"/inflow_on_t1pt1.nii "${folder_pet}"/inflow/inflow_on_t1pt1_


	#for I in {0..10}
	#do
	#	wtf=$(printf "%04d" "$I")
	#	petpvc -i "${folder_pet}"/task/task_on_t1pt2_"${wtf}".nii.gz -m "${folder_pet}"/aparc+aseg_pt2_nat_labelled_4D.nii.gz -o "${folder_pet}"/task/PVEc_RBV_task_on_t1pt2_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
	#done

	#fslmerge -t "${folder_pet}"/PVEc_RBV_task_on_t1pt2.nii.gz "${folder_pet}"/task/PVEc_RBV_task_on_t1pt2*.nii.gz
	#rm "${folder_pet}"/task/task_on_t1pt2_*.nii.gz



	#for I in {0..14}
	#do
	#	wtf=$(printf "%04d" "$I")
	#	petpvc -i "${folder_pet}"/baseline/baseline_on_t1pt2_"${wtf}".nii.gz -m "${folder_pet}"/aparc+aseg_pt2_nat_labelled_4D.nii.gz -o "${folder_pet}"/baseline/PVEc_RBV_baseline_on_t1pt2_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
	#done

	#fslmerge -t "${folder_pet}"/PVEc_RBV_baseline_on_t1pt2.nii.gz "${folder_pet}"/baseline/PVEc_RBV_baseline_on_t1pt2*.nii.gz
	#rm "${folder_pet}"/baseline/baseline_on_t1pt2_*.nii.gz


	for I in {0..84}
	do
		wtf=$(printf "%04d" "$I")
		echo $wtf
		petpvc -i "${folder_pet}"/inflow/inflow_on_t1pt1_"${wtf}".nii.gz -m "${folder_pet}"/aparc+aseg_pt1_nat_labelled_4D.nii.gz -o "${folder_pet}"/inflow/PVEc_RBV_inflow_on_t1pt1_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
	done

	fslmerge -t "${folder_pet}"/PVEc_RBV_inflow_on_t1pt1.nii.gz "${folder_pet}"/inflow/PVEc_RBV_inflow_on_t1pt1*.nii.gz
	rm "${folder_pet}"/inflow/inflow_on_t1pt1_*.nii.gz


done