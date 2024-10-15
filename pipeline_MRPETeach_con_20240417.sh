#!/bin/bash

  # set up folders and group space images
  folder=/Volumes/korokdorf/MRPET/coreg_mri/eachsession/${ID}/data/
  MNI=/Users/alex/Documents/mni_icbm152_t1_tal_nlin_asym_09c.nii
  template=/Users/alex/Documents/mrpet_template.nii.gz

for K in 40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40112 40121 40122 40131 40132 40141 40142 40151 40152 40161 40162 40171 40172 40181 40182 40191 40192 40201 40202 40211 40212 40221 40222 40231 40232 40241 40242 40251 40252 40261 40262 40271 40281 40282 40292 40301 40302 40311 40312 40321 40322 40331 40332

do

  ID=${K:0:4}   # Extract first 4 digits of K
  DAY=${K:4:1}  # Extract last digit of K
  IDALL=$K

  for I in {1..17}
  do
    K=$(printf %04d ${I})
    antsApplyTransforms -d 3 -v 0 -n Linear -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/NLreg_template_to_MNI_1Warp.nii.gz -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/NLreg_template_to_MNI_0GenericAffine.mat -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/${IDALL}/data/NLreg_T1WB_to_template_1Warp.nii.gz -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/${IDALL}/data/NLreg_T1WB_to_template_0GenericAffine.mat -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/${IDALL}/data/coreg_meanEPI_to_T1WB_0GenericAffine.mat -i /Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Eachsessions_memory/${ID}_${DAY}/con_${K}.nii -r "${MNI}" -o /Users/alex/Dropbox/paperwriting/MRPET/data/fMRI/1stLevel/Eachsessions_memory/${ID}_${DAY}/con_${K}_mni.nii
  done

done