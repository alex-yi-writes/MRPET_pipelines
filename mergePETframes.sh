#!/bin/bash	

export FSLDIR=/usr/local/fsl
export PATH=$PATH:/usr/local/fsl/bin
source $FSLDIR/etc/fslconf/fsl.sh

for K in 40011 40012 40021 40022 40031 40041 40042 40051 40052 40062 40071 40081 40082 40092 40101 40102 40111 40112 40121 40122 40131 40132 40141 40142 40151 40152 40161 40162 40171 40172 40181 40182 40191 40192 40201 40202 40211 40212 40221 40222 40231 40232 40241 40242 #40251 40252 40261 40262 40271 40281 40282 40292 40301 40302 40311 40312 40321 40322 40331 40332 
do
	mkdir /Volumes/ALEX5/MRPET/newRecon/preproc/${K}
	fslmerge -t /Volumes/ALEX5/MRPET/newRecon/preproc/${K}/inflow_on_t1.nii.gz /Users/yyi/Desktop/new_recon/preproc/${K}/inflow3D/inflow_on_t1pt1_*.nii.gz
	fslmerge -t /Volumes/ALEX5/MRPET/newRecon/preproc/${K}/bsl_on_t1.nii.gz /Users/yyi/Desktop/new_recon/preproc/${K}/bsl3D/bsl_on_t1pt2_*.nii.gz
	fslmerge -t /Volumes/ALEX5/MRPET/newRecon/preproc/${K}/task_on_t1.nii.gz /Users/yyi/Desktop/new_recon/preproc/${K}/task3D/task_on_t1pt2_*.nii.gz
done