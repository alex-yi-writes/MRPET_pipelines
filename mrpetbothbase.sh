#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=168:00:00,h_vmem=8G,mem_free=8G
#$ -q work.q

# subject ID





# set up folders and group space images
folder=/mnt/work/yyi/mrpet_both/"${ID}"/data/
MNI=/mnt/work/yyi/mni_icbm152_t1_tal_nlin_asym_09c.nii
template=/mnt/work/yyi/mrpet_template.nii.gz


# -------------- start the transformation -------------- #

# T1WB -> study template
#antsRegistrationSyN.sh -d 3 -t s -f "${template}" -m "${folder}"T1WB.nii -o "${folder}"NLreg_T1WB_to_template_

# T1 slab -> T1WB
#antsRegistrationSyN.sh -d 3 -t r -m "${folder}"t1slab.nii -f "${folder}"T1WB.nii -o "${folder}"coreg_t1slab_to_T1WB_

# EPI -> T1WB
#antsRegistrationSyN.sh -d 3 -t r -m "${folder}"meanEPI.nii -f "${folder}"T1WB.nii -o "${folder}"coreg_meanEPI_to_T1WB_

# T1WB -> EPI
#antsRegistrationSyN.sh -d 3 -t r -f "${folder}"meanEPI.nii -m "${folder}"T1WB.nii -o "${folder}"coreg_T1WB_to_meanEPI_

# T1WB -> MNI
#antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_1Warp.nii.gz -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder}"T1WB.nii -r "${MNI}" -o "${folder}"NLreg_T1WB_to_MNI.nii

# EPI -> MNI
#antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_1Warp.nii.gz -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folder}"meanEPI.nii -r "${MNI}" -o "${folder}"NLreg_meanEPI_to_MNI.nii

#for I in {01..16}
#do
  #antsApplyTransforms -d 3 -v 0 -n Linear -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_1Warp.nii.gz -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folder}"con_00${I}.nii -r "${MNI}" -o "${folder}"con_00${I}_mni.nii
#done

for I in {01..16}
do
antsApplyTransforms -d 3 -v 0 -n Linear -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folder}"con_00${I}.nii -r "${template}" -o "${folder}"con_00${I}_mni.nii
done


# ------------------------------------------------------ #
