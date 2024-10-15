#!/bin/bash

folder=/Volumes/korokdorf/MRPET/coreg_roi
template=/Volumes/korokdorf/MRPET/mrpet_template.nii.gz
maskdir=/Users/alex/Dropbox/paperwriting/MRPET/data/BPdat_fMRIclusters/clustermasks

for mask in ${maskdir}/*.nii.gz
do
  maskname=$(basename "$mask" .nii.gz)

  for K in 40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40112 40121 40122 40131 40132 40141 40142 40151 40152 40161 40162 40171 40172 40181 40182 40191 40192 40201 40202 40211 40212 40221 40222 40231 40232 40241 40242 40251 40252 40261 40262 40271 40281 40282 40292 40301 40302 40311 40312 40321 40322 40331 40332

  do

    ID=${K}
    path=${folder}/${ID}


    #pt 1
    #antsRegistrationSyNQuick.sh -n 4 -d 3 -f ${template} -m "${path}"/data/T1pt1.nii -t s -o "${path}"/NLreg_T1pt1_to_template_
    antsApplyTransforms -d 3 -v 1 -n NearestNeighbor -r "${path}"/data/T1pt1.nii -o "${path}"/data/"${maskname}".nii.gz -i "${mask}" -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/${ID}/data/NLreg_T1pt1_to_template_1InverseWarp.nii.gz -t [/Volumes/korokdorf/MRPET/coreg_mri/eachsession/${ID}/data/NLreg_T1pt1_to_template_0GenericAffine.mat,1] -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/NLreg_MNI_to_template_1Warp.nii.gz -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/NLreg_MNI_to_template_0GenericAffine.mat

    #pt 2
    #antsRegistrationSyNQuick.sh -n 4 -d 3 -f ${template} -m "${path}"/data/T1pt2.nii -t s -o "${path}"/NLreg_T1pt2_to_template_
    antsApplyTransforms -d 3 -v 1 -n NearestNeighbor -r "${path}"/data/T1pt2.nii -o "${path}"/data/"${maskname}".nii.gz -i "${mask}" -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/${ID}/data/NLreg_T1pt2_to_template_1InverseWarp.nii.gz -t [/Volumes/korokdorf/MRPET/coreg_mri/eachsession/${ID}/data/NLreg_T1pt2_to_template_0GenericAffine.mat,1] -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/NLreg_MNI_to_template_1Warp.nii.gz -t /Volumes/korokdorf/MRPET/coreg_mri/eachsession/NLreg_MNI_to_template_0GenericAffine.mat
  done

done
