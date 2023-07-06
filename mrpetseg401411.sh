#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=168:00:00,h_vmem=4G,mem_free=4G
#$ -q work.q

# subject ID



ID=401411






folder=/mnt/work/yyi/MRPETseg/$ID

recon-all \
-i  ${folder}/T1pt1.nii.gz \
-s  ${ID} \
-sd  ${folder} \
-all
