#!/bin/bash

CPUs=2

# write a list of images
find . -name '*nii' -o -name '*.nii.gz' > images.txt

# generate a group template
antsMultivariateTemplateConstruction2.sh -d 3 -a 1 -i 8 -k 1 -r 1 -f 6x4x2x1 -s 4x2x1x0vox -q 100x100x70x20 -t SyN -m CC -c 1 -o images.txt
