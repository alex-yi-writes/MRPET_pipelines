#!/usr/local/bin/bash

# work in progress

# actually its MNI to T1 now
MNI=/Volumes/ALEX3/mni_icbm152_t1_tal_nlin_asym_09c.nii
template=/Users/alex/Documents/pilot_template.nii.gz

for K in 401811 401812 401821 401822 401911 401912 402011 402012 402021 402022 402111 402112 402121 402122 402221 402222

do
  # subject ID
  ID=${K}
  # set up folders and group space images
  folder_trx=/Volumes/ALEX3/MRPET/coreg_mri/"${ID}"/
  folder_pet=/Volumes/ALEX3/MRPET/coreg_roi/"${ID}"/


done
