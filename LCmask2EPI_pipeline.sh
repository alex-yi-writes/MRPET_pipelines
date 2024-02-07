#!/bin/bash

for K in 4282 #3109 3110 3113 3115 3121 3123 3125 3165 3171 3172 3173 3174 3208 3211 3212 3214 3216 3217 3218 3220 3222 3252 3253 3264 3269 3270 3275 3276 4103 4104 4108 4156 4173 4176 4177 4183 4184 4185 4202 4205 4210 4211 4268 4269 4274 4275 4279 4280 4281

do
	ID=${K}

	mkdir /Users/yyi/Desktop/LCmask_on_nativeEPI/${ID}/

	antsApplyTransforms -d 3 -n NearestNeighbor -i /Volumes/ALEX_DATA4/mareike_coreg/LCmask.nii -r /Volumes/ALEX_DATA4/mareike_coreg/${ID}/data/meanEPI3D1.nii -o /Users/yyi/Desktop/LCmask_on_nativeEPI/${ID}/LCmask_on_nativeEPI.nii -t /Volumes/ALEX_DATA4/mareike_coreg/${ID}/data/coreg_T1mean_to_meanEPI3D1_0GenericAffine.mat -t /Volumes/ALEX_DATA4/mareike_coreg/${ID}/data/NLreg_T1mean_to_template_1InverseWarp.nii.gz -t [/Volumes/ALEX_DATA4/mareike_coreg/${ID}/data/NLreg_T1mean_to_template_0GenericAffine.mat,1] -t /Volumes/ALEX_DATA4/mareike_coreg/NLreg_template_to_MNI_1InverseWarp.nii.gz -t [/Volumes/ALEX_DATA4/mareike_coreg/NLreg_template_to_MNI_0GenericAffine.mat,1]

done