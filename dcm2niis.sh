#!/bin/bash

flist=`ls ./PI*`
for eachfile in $flist
do
   echo $eachfile
   mri_convert --in_type dicom --out_type nii
   /Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADY_originals/DOPE/4001_2/study_2_20191126/series_40004__PET-Baseline_HD_I3_256_Z3_PRR_AC_Images/$eachfile /Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/4001_2/$eachfile.nii
done
