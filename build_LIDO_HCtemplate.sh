#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=24:00:00,h_vmem=6G,mem_free=6G
#$ -q test.q

# if you submit to cluster, use this
# cd /mnt/work/

# if only in hxdra, use whatever path
# cd /DataTempVolatile/lancini4template/
cd /mnt/work/yyi/LIDO_HCtemplate/

# write a list of images
find "$PWD" -name '*nii' -o -name '*.nii' > images.txt

mv images.txt /mnt/work/yyi/LIDO_HCtemplate/output/

#Run from own dir
cd /mnt/work/yyi/LIDO_HCtemplate/output/

# generate a group template
/mnt/work/yyi/antsMultivariateTemplateConstruction2_v3.sh -d 3 -v 6gb -a 1 -i 8 -k 1 -r 1 -n 1 -f 6x4x2x1 -s 4x2x1x0vox -q 100x100x70x20 -t SyN -m CC -c 1 -o zz_ images.txt
