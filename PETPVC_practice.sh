#!/bin/bash	

#pvc_make4d -i /Users/yyi/Desktop/testdata/aparc+aseg_pt2_nat_labelled.nii -o /Users/yyi/Desktop/testdata/aparc+aseg_pt2_nat_labelled_4D.nii.gz
#pvc_make4d -i /Users/yyi/Desktop/testdata/aparc+aseg_pt2_nat_labelled.nii -o /Users/yyi/Desktop/testdata/aparc+aseg_pt1_nat_labelled_4D.nii.gz

#fslsplit /Users/yyi/Desktop/testdata/task_on_t1pt2.nii /Users/yyi/Desktop/testdata/task/task_on_t1pt2_
#fslsplit /Users/yyi/Desktop/testdata/baseline_on_t1pt2.nii /Users/yyi/Desktop/testdata/baseline/baseline_on_t1pt2_
#fslsplit /Users/yyi/Desktop/testdata/inflow_on_t1pt1.nii /Users/yyi/Desktop/testdata/inflow/inflow_on_t1pt1_


for I in {0..10}
do
	wtf=$(printf "%04d" "$I")
	petpvc -i /Users/yyi/Desktop/testdata/task/task_on_t1pt2_"${wtf}".nii.gz -m /Users/yyi/Desktop/testdata/aparc+aseg_pt2_nat_labelled_4D.nii.gz -o /Users/yyi/Desktop/testdata/task/PVEc_RBV_task_on_t1pt2_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
done

fslmerge -t /Users/yyi/Desktop/testdata/PVEc_RBV_task_on_t1pt2.nii.gz /Users/yyi/Desktop/testdata/task/PVEc_RBV_task_on_t1pt2*.nii.gz
rm /Users/yyi/Desktop/testdata/task/task_on_t1pt2_*.nii.gz



for I in {0..14}
do
	wtf=$(printf "%04d" "$I")
	petpvc -i /Users/yyi/Desktop/testdata/baseline/baseline_on_t1pt2_"${wtf}".nii.gz -m /Users/yyi/Desktop/testdata/aparc+aseg_pt2_nat_labelled_4D.nii.gz -o /Users/yyi/Desktop/testdata/baseline/PVEc_RBV_baseline_on_t1pt2_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
done

fslmerge -t /Users/yyi/Desktop/testdata/PVEc_RBV_baseline_on_t1pt2.nii.gz /Users/yyi/Desktop/testdata/baseline/PVEc_RBV_baseline_on_t1pt2*.nii.gz
rm /Users/yyi/Desktop/testdata/baseline/baseline_on_t1pt2_*.nii.gz



for I in {0..84}
do
	wtf=$(printf "%04d" "I")
	petpvc -i /Users/yyi/Desktop/testdata/inflow/inflow_on_t1pt1_"${wtf}".nii.gz -m /Users/yyi/Desktop/testdata/aparc+aseg_pt1_nat_labelled_4D.nii.gz -o /Users/yyi/Desktop/testdata/inflow/PVEc_RBV_inflow_on_t1pt1_"${wtf}".nii.gz --pvc RBV -x 3.0 -y 3.0 -z 3.0
done

fslmerge -t /Users/yyi/Desktop/testdata/PVEc_RBV_inflow_on_t1pt1.nii.gz /Users/yyi/Desktop/testdata/inflow/PVEc_RBV_inflow_on_t1pt1*.nii.gz
rm /Users/yyi/Desktop/testdata/inflow/inflow_on_t1pt1_*.nii.gz

