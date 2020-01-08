#!/bin/bash

# setup paths
path_parent=/storage/users/yyi/segmentation/

# assign number of CPUs used
CPUs=3

# read in subjects
ls -d "${path_parent}"/*/ > "${path_parent}"/subjs.txt

while read folder; do

  ID=$(echo "${folder}" | grep -o -E '[0-9][0-9][0-9][0-9][0-9][0-9]')
  recon-all \
  -i  ${path_parent}T1w_clusterIKND/${ID}_T1w.nii.gz \
  -s  ${ID} \
  -sd ${folder} \
  -all

  echo "ID ${ID} DONE!"

done < "${path}"/subjs.txt
