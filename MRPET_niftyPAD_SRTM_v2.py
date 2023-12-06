##### why doesn't it work?!?!?!?!?!?!?!?!?!?!?!?!?!

import numpy as np
import nibabel as nib
import niftypad as npd

from niftypad.models import feng_srtm
from niftypad.models import srtmb_basis
from niftypad.models import mrtm

# constants for decay correction: siemens scanner cannot handle breaks inbetween the scans
t0_frame_bsl = 95 # minutes after injection
t0_frame_task = 115 # minutes after injection
half_life = 109 # fallypride half-life
decay_factor_bsl = 2 ** (t0_frame_bsl / half_life)
decay_factor_task = 2 ** (t0_frame_task / half_life)

# time and frame duration (in seconds)
time = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,360,420,480,540,600,660,720,780,840,900,960,1020,1080,1140,1200,1260,1320,1380,1440,1500,1560,1620,1680,1740,1800,1860,1920,1980,2040,2100,2160,2220,2280,2340,2400,2460,2520,2580,2640,2700,2760,2820,2880,2940,3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,5700,5760,5820,5880,5940,6000,6060,6120,6180,6240,6300,6360,6420,6480,6540,6900,7200,7500,7800,8100,8400,8700,9000,9300,9600,9900])
frame_duration = np.array([10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,300,300,300,300,300,300,300,300,300,300,300])

# frame weights for modelling purpose
weights = np.concatenate([0.25 * np.ones(30), np.ones(15 + 46), 100 * np.ones(11)])

# ROI labels dictionary - now streamlined: it's long but i don't know any other way
roi_labels = {
#    "bankssts": [1001, 2001],
    "caudalanteriorcingulate": [1002, 2002],
#    "caudalmiddlefrontal": [1003, 2003],
#    "cuneus": [1005, 2005],
    "entorhinal": [1006, 2006],
#    "fusiform": [1007, 2007],
#    "inferiorparietal": [1008, 2008],
    "inferiortemporal": [1009, 2009],
#    "isthmuscingulate": [1010, 2010],
#    "lateraloccipital": [1011, 2011],
#    "lateralorbitofrontal": [1012, 2012],
#    "lingual": [1013, 2013],
#    "medialorbitofrontal": [1014, 2014],
    "middletemporal": [1015, 2015],
    "parahippocampal": [1016, 2016],
#    "paracentral": [1017, 2017],
#    "parsopercularis": [1018, 2018],
#    "parsorbitalis": [1019, 2019],
#    "parstriangularis": [1020, 2020],
#    "pericalcarine": [1021, 2021],
#    "postcentral": [1022, 2022],
#    "posteriorcingulate": [1023, 2023],
#    "precentral": [1024, 2024],
#    "precuneus": [1025, 2025],
    "rostralanteriorcingulate": [1026, 2026],
#    "rostralmiddlefrontal": [1027, 2027],
#    "superiorfrontal": [1028, 2028],
#    "superiorparietal": [1029, 2029],
    "superiortemporal": [1030, 2030],
#    "supramarginal": [1031, 2031],
#    "frontalpole": [1032, 2032],
    "temporalpole": [1033, 2033],
#    "transversetemporal": [1034, 2034],
    "insula": [1035, 2035],
    "CerC": [8, 47],
    "ThalProp": [10, 49],
    "Caud": [11, 50],
    "Put": [12, 51],
    "Pall": [13, 52],
    "Hipp": [17, 53],
    "Amy": [18, 54],
    "Nac": [26, 58],
    "SN": [991, 992],
    "LC": [993, 994],
    "VTA": [995, 996],
 #   "SubcorticalROIs": [997, 998],
 #   "SC": [888, 887]
}

# let's try with one subject first...
ID=["4001"]#,"4002","4003","4004","4005","4006","4007","4008","4009","4010","4011","4012","4013","4014","4015","4016","4017","4018","4019","4020","4021","4022","4023","4024","4025","4026","4027","4028","4029","4030","4031","4032","4033"]

# loop over subjects
for subject_id in ID: 
    for acquisition in [1, 2]:  # two acquisitions for each subject
        for roi_name, roi_idx in roi_labels.items():
            try:

                # load all three scans: InFlow, Baseline, Task(MT),
                img_inflow = nib.load(f'/Volumes/ALEX5/MRPET/img/{subject_id}_{acquisition}/coreg_InFlow{acquisition}_on_T1.nii').get_fdata()
                img_baseline = nib.load(f'/Volumes/ALEX5/MRPET/img/{subject_id}_{acquisition}/coreg_Baseline{acquisition}_on_T1.nii').get_fdata()
                img_task = nib.load(f'/Volumes/ALEX5/MRPET/img/{subject_id}_{acquisition}/coreg_MT{acquisition}_on_T1.nii').get_fdata()

                # load ROI mask unique to this subject and acquisition
                roi_mask_inflow = nib.load(f'/Volumes/ALEX5/MRPET/coreg_roi/{subject_id}{acquisition}/aparc+aseg_pt1_nat_labelled.nii').get_fdata()
                roi_mask_baseline_n_task = nib.load(f'/Volumes/ALEX5/MRPET/coreg_roi/{subject_id}{acquisition}/aparc+aseg_pt2_nat_labelled.nii').get_fdata()

                # create boolean masks based on label indices for the ROI
                roi_mask_inflow_bool = np.isin(roi_mask_inflow, roi_idx)
                roi_mask_baseline_n_task_bool = np.isin(roi_mask_baseline_n_task, roi_idx)

                # extract TAC from ROI for InFlow scan
                inflow_tac = np.mean(img_inflow[roi_mask_inflow == roi_idx], axis=0)  # No decay correction for InFlow

                # decay correction for baseline and task
                baseline_tac = np.mean(img_baseline[roi_mask_baseline_n_task == roi_idx], axis=0) * decay_factor_bsl
                task_tac = np.mean(img_task[roi_mask_baseline_n_task == roi_idx], axis=0) * decay_factor_task

                # ref region TAC
                CerC_indices = roi_labels['CerC']
                CerCtac_baseline = np.mean(img_baseline[np.isin(roi_mask_baseline_n_task, CerC_indices)], axis=0) * decay_factor_bsl
                CerCtac_task = np.mean(img_task[np.isin(roi_mask_baseline_n_task, CerC_indices)], axis=0) * decay_factor_task
                CerCtac_inflow = np.mean(img_inflow[np.isin(roi_mask_inflow, CerC_indices)], axis=0)
                CerCtac = np.concatenate([CerCtac_inflow, CerCtac_baseline, CerCtac_task])

                ############## SRTM ##############
                srtm_result = npd.models.srtm(
                    tac=[inflow_tac, baseline_tac, task_tac],
                    w=weights,
                    dt=time,
                    inputf1=CerCtac
                )
                # extract and save SRTM modelling results
                R1, k2, BP_ND = srtm_result['R1'], srtm_result['k2'], srtm_result['BP_ND']
                np.save(f'/Volumes/ALEX5/MRPET/img/{subject_id}_{acquisition}/srtm_results_{subject_id}_{acquisition}_{roi_name}.npy', {'R1': R1, 'k2': k2, 'BP_ND': BP_ND})
                ##################################

                ############## MRTM: doesn't work and i'm too tired to try more
#                mrtm_result = npd.models.mrtm(
#                    tac=[inflow_tac, baseline_tac, task_tac],
#                    dt=time,
#                )
                # extract and save MRTM modelling results
#                k2a, R1, BP_ND, BP = mrtm_result['k2a'], mrtm_result['R1'], mrtm_result['BP_ND'], mrtm_result['BP']
#                np.save(f'mrtm_results_{subject_id}_{acquisition}_{roi_name}.npy', {'k2a': k2a, 'R1': R1, 'BP_ND': BP_ND, 'BP': BP})
                ##################################

            except FileNotFoundError:
                print(f"dataset missing for subject {subject_id}, day {acquisition}, ROI {roi_name}. skipping...")
                continue
