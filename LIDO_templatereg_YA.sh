#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=168:00:00,h_vmem=8G,mem_free=8G
#$ -q work.q

# make transformation files for LIDO data, in each corresponding templates

template=/mnt/work/yyi/LIDO_YA_template.nii.gz
STUDYtemplte=/mnt/work/yyi/LIDO_STUDY_template.nii.gz
HCtemplate=/mnt/work/yyi/LIDO_HC_template.nii.gz

for K in 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 127 128 129 130 131 132 133 134 135
do

	# subject ID
	ID=${K}

	# set up folders
	folder_img=/mnt/work/yyi/LIDO_reg/"${ID}"
	T1=${folder_img}/T1_corrected.nii
	EPI=${folder_img}/meanEPI_EM_biascorr.nii


	antsRegistrationSyN.sh -d 3 -n 2 -t s -f ${template} -m ${T1} -o ${folder_img}/NLreg_T1mean_to_YAtemplate_
	antsRegistrationSyN.sh -d 3 -n 2 -t s -f ${HCtemplate} -m ${T1} -o ${folder_img}/NLreg_T1mean_to_HCtemplate_
	antsRegistrationSyN.sh -d 3 -n 2 -t r -f ${T1} -m ${EPI} -o ${folder_img}/coreg_meanEPI_EM_to_T1mean_
	antsRegistrationSyN.sh -d 3 -n 2 -t r -m ${T1} -f ${EPI} -o ${folder_img}/coreg_T1mean_meanEPI_EM_

	mkdir ${folder_img}/YAtemplate
	mkdir ${folder_img}/STUDYtemplate


	# register images to the group template
	antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t ${folder_img}/NLreg_T1mean_to_YAtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_YAtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/${EPI} -r "${template}" -o "${folder_img}"/YAtemplate/NLreg_EPI_EM_on_YAtemplate.nii.gz
	for I in {01..08}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_YAtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_YAtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_memory.nii -r "${template}" -o "${folder_img}"/YAtemplate/con_00${I}_memory_YAtemplate.nii.gz
	done

	for I in {01..12}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_YAtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_YAtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutral.nii -r "${template}" -o "${folder_img}"/YAtemplate/con_00${I}_emoneutral_YAtemplate.nii.gz
	done

	for I in {01..07}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_YAtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_YAtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutralnofix.nii -r "${template}" -o "${folder_img}"/YAtemplate/con_00${I}_emoneutralnofix_YAtemplate.nii.gz
	done



	# register images to the healthy control template
	antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t ${folder_img}/NLreg_T1mean_to_HCtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_HCtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/${EPI} -r "${HCtemplate}" -o "${folder_img}"/HCtemplate/NLreg_EPI_EM_on_HCtemplate.nii.gz
	for I in {01..08}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_HCtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_HCtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_memory.nii -r "${HCtemplate}" -o "${folder_img}"/HCtemplate/con_00${I}_memory_HCtemplate.nii.gz
	done

	for I in {01..12}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_HCtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_HCtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutral.nii -r "${HCtemplate}" -o "${folder_img}"/HCtemplate/con_00${I}_emoneutral_HCtemplate.nii.gz
	done

	for I in {01..07}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_HCtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_HCtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutralnofix.nii -r "${HCtemplate}" -o "${folder_img}"/HCtemplate/con_00${I}_emoneutralnofix_HCtemplate.nii.gz
	done



	# register images to the study template
	antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t ${folder_img}/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_template_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/${EPI} -r "${STUDYtemplte}" -o "${folder_img}"/STUDYtemplate/NLreg_EPI_EM_on_STUDYtemplate.nii.gz
	for I in {01..08}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_template_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_memory.nii -r "${STUDYtemplte}" -o "${folder_img}"/STUDYtemplate/con_00${I}_memory_STUDYtemplate.nii.gz
	done

	for I in {01..12}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_template_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutral.nii -r "${STUDYtemplte}" -o "${folder_img}"/STUDYtemplate/con_00${I}_emoneutral_STUDYtemplate.nii.gz
	done

	for I in {01..07}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_template_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutralnofix.nii -r "${STUDYtemplte}" -o "${folder_img}"/STUDYtemplate/con_00${I}_emoneutralnofix_STUDYtemplate.nii.gz
	done

done