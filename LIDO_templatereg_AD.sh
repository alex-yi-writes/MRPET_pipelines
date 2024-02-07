#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=168:00:00,h_vmem=8G,mem_free=8G
#$ -q work.q

# make transformation files for LIDO data, in each corresponding templates

template=/mnt/work/yyi/LIDO_AD_template.nii.gz
STUDYtemplte=/mnt/work/yyi/LIDO_STUDY_template.nii.gz


for K in 302 303 309 310 312 314 317 319 320 326 327 329 330 331 301 304 305 306 307 308 311 313 315 316 318 321 323 324 325 328 332 333
do

	# subject ID
	ID=${K}

	# set up folders
	folder_img=/mnt/work/yyi/LIDO_reg/"${ID}"
	T1=${folder_img}/T1_corrected.nii
	EPI=${folder_img}/meanEPI_EM_biascorr.nii


	antsRegistrationSyN.sh -d 3 -n 2 -t s -f ${template} -m ${T1} -o ${folder_img}/NLreg_T1mean_to_ADtemplate_
	antsRegistrationSyN.sh -d 3 -n 2 -t r -f ${T1} -m ${EPI} -o ${folder_img}/coreg_meanEPI_EM_to_T1mean_
	antsRegistrationSyN.sh -d 3 -n 2 -t r -m ${T1} -f ${EPI} -o ${folder_img}/coreg_T1mean_meanEPI_EM_

	mkdir ${folder_img}/ADtemplate
	mkdir ${folder_img}/STUDYtemplate


	# register images to the group template
	antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t ${folder_img}/NLreg_T1mean_to_ADtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_ADtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/${EPI} -r "${template}" -o "${folder_img}"/ADtemplate/NLreg_EPI_EM_on_ADtemplate.nii.gz
	for I in {01..08}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_ADtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_ADtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_memory.nii -r "${template}" -o "${folder_img}"/ADtemplate/con_00${I}_memory_ADtemplate.nii.gz
	done

	for I in {01..12}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_ADtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_ADtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutral.nii -r "${template}" -o "${folder_img}"/ADtemplate/con_00${I}_emoneutral_ADtemplate.nii.gz
	done

	for I in {01..07}
	do
	antsApplyTransforms -d 3 -v 0 -n Linear -t ${folder_img}/NLreg_T1mean_to_ADtemplate_1Warp.nii.gz -t "${folder_img}"/NLreg_T1mean_to_ADtemplate_0GenericAffine.mat -t "${folder_img}"/coreg_meanEPI_EM_to_T1mean_0GenericAffine.mat -i "${folder_img}"/con_00${I}_emoneutralnofix.nii -r "${template}" -o "${folder_img}"/ADtemplate/con_00${I}_emoneutralnofix_ADtemplate.nii.gz
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