#!/bin/bash

# setup paths
path_parent=/mnt/work/yyi/temp/parcellation/MRPET
path_subjlist=/mnt/work/yyi/temp/parcellation/subjlist

# assign number of CPUs used
CPUs=2

# read in subjects
ls -d ${path_parent}/* > ${path_subjlist}/MRPET_subj1.txt

while read folder; do

  ID=$(echo ${folder} | grep -o -E '[0-9][0-9][0-9][0-9][0-9]')
  recon-all \
  -i  ${folder}/T1WB.nii \
  -s  ${ID} \
  -sd  ${folder} \
  -all

  echo "ID ${ID} DONE!"

  mri_convert --in_type mgz --out_type nii --out_orientation RAS ${folder}/40012/mri/aparc+aseg.mgz ${folder}/40012/mri/aparc+aseg.nii


done < ${path_subjlist}/MRPET_subj1.txt
