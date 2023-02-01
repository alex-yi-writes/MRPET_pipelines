%% write out TACs

clc;clear

% 4001_2 --> inflow bugged
% IDs     = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4026]; % 4009, visit 1 has no task data
% days    = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0; 0 2];
IDs  = [4028 4029 4030 4031];
days = [0 2; 0 2; 1 2; 0 2; 1 0]; 

% % move files
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%             warning('Skipped')
%         else
%             mkdir(['/Users/yeojin/Desktop/E_data/ED_coreg/mrpet/' num2str(IDs(id)) num2str(d) '/data/'])
%             copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'],...
%                 ['/Users/yeojin/Desktop/E_data/ED_coreg/mrpet/' num2str(IDs(id)) num2str(d) '/data/T1pt2.nii'])
% 
%         end
%     end
% end
% 

addpath('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/modelling/')
addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')
opengl hardwarebasic

path_parent = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/';
path_rawImg   = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/';


% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
addpath('/Users/yeojin/Documents/MATLAB/NIfTI_20140122')
opengl hardwarebasic

%%
% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
TrimPerc=15;

clear id d

for id = 2:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask1   =[ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii'];
            
            % 4-D PET file
            PETflow = [path_rawImg num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_InFlow' num2str(d) '_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask1);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
            maskidx=[];
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,2}).LongName=strtrim(ROIDef{ROITabidx,3});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1}+1000; % Add 1000 to the index
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,1}+2000; % Add 2000 to the index
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}+1000, ROIDef{ROITabidx,1}+2000];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,2}).(Hemi{1})=IndList;
                end
            end
            
            
            [A B ROIDef]=xlsread('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,3}).LongName=strtrim(ROIDef{ROITabidx,4});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1};
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,2};
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}, ROIDef{ROITabidx,2}];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,3}).(Hemi{1})=IndList;
                end
            end
            
            %% Read in 4D-PET data and extract ROI-averages in each frame
            
            % inFlow
            DynPET=load_nii(PETflow);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'inFlow';
            save([ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end

%
clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask2   =[ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(days(id,d)) '/aparc+aseg_pt2_nat_labelled.nii']; % task and the baseline
            
            % 4-D PET file
            PETtask = [path_rawImg num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_MT' num2str(d) '_on_T1.nii'];
            PETbsl  = [path_rawImg num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_Baseline' num2str(d) '_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask2);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
            maskidx=[];
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,2}).LongName=strtrim(ROIDef{ROITabidx,3});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1}+1000; % Add 1000 to the index
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,1}+2000; % Add 2000 to the index
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}+1000, ROIDef{ROITabidx,1}+2000];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,2}).(Hemi{1})=IndList;
                end
            end
            
            
            [A B ROIDef]=xlsread('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBD_PET/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,3}).LongName=strtrim(ROIDef{ROITabidx,4});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1};
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,2};
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}, ROIDef{ROITabidx,2}];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,3}).(Hemi{1})=IndList;
                end
            end
            
            %% Read in 4D-PET data and extract ROI-averages in each frame
            
            % task
            DynPET=load_nii(PETtask);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'RewardTask';
            save([ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
            %             % inFlow
            %             DynPET=load_nii(PETflow);
            %             temp=size(DynPET.img);
            %             ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            %             for ROI=fieldnames(maskidx)'
            %                 TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
            %                 for Hemi={'Left','Right','Bilateral'}
            %                     maskidxcur = maskidx.(ROI{1}).(Hemi{1});
            %                     TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
            %                     TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
            %                 end
            %
            %             end
            %             TACDATA.info = 'inFlow';
            %             save([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
            % baseline
            DynPET=load_nii(PETbsl);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'Baseline';
            save([ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end


%% retro correct the decay

cd(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned'])

t0_frame_bsl    = 95;
t0_frame_task   = 115;

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    % bilaterals
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^(t0_frame_bsl/109));
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^(t0_frame_task/109));
                    
                    % lefts
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^(t0_frame_bsl/109));
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^(t0_frame_task/109));
                    
                    % rights
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^(t0_frame_bsl/109));
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^(t0_frame_task/109));
                end
            end
            
            save([path_parent 'coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
            
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
            try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_' date '.pdf']);
                
            catch
                
                Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_' date '.pdf']);
                
                
            end
            
            
        end
    end
end


%% modelling

clear; clc; close all
warning('off','all');

% set paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/';
% paths.spm     = '/Users/yeojin/Documents/MATLAB/spm12';
% paths.funx_MRI= [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
% paths.funx_PET= [paths.parent 'B_scripts/BA_preprocessing/BAC_PET/preproc_functions/'];
% paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
% paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
% paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
% paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
% paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
% paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';
% paths.seg     = [paths.parent 'E_data/EA_raw/EAB_MRI/EABD_segmented/'];
paths.TACs    = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/'];
% paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/'];


% paths.TACs    = ['/Users/Alex/Dropbox/literatures_IKND/EADC_TACs/modelling/RewardTask_deccorr_2/'];
% paths.figures = ['/Users/Alex/Dropbox/literatures_IKND/EADC_TACs/modelling/RewardTask_deccorr_2/'];
% addpath(genpath('/Users/Alex/Dropbox/literatures_IKND/EADC_TACs/modelling'))
addpath(genpath('/Users/alex/Documents/MATLAB/BBD_PET/modelling'))

% set envs
% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

% load TACs
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            TACs{i1,d} = [];
        else
            TACs{i1,d} = load([paths.TACs num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr.mat']);
        end
    end
end

load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/ROIs_LCSNVTA.mat')
ROI=ROI';

disp('prep done')

%% run the model: condensed ROIs

set(0, 'DefaultLineLineWidth', 1);

for id = [19:length(IDs)]%:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg={'CerC','Put','Caud','Nac','Hipp','ThalProp'}
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
            end
            
            Lengths=[10*ones(30,1); 60*ones(55,1)]; % time windows
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
            
            tmid=mean(times,2);
            tmidMin=tmid/60;
            t_points    = length(tmid);
            dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
            break_point=find(times(:,1)>=115*60,1,'first'); %% TIme of activation start
            
            %%%%%%%%%%%
            PlotStrFit=1;
            %%%%%%%%%%%
            
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg={'Striatum','Caudate','Putamen','Accumbens-area','Thalamus','Hippocampus'}
                    if PlotStrFit
                        figure('Position',[100 100 800 1200],'Visible','on'); hold on;
                        spidx=1;
                    end
                    reftac=TACDATA.CerC;
                    mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
                    %% Set the SRTM part of A
                    ASRTM = zeros(t_points ,3);
                    ASRTM(:,1)  = reftac;
                    for k = 1:t_points
                        ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
                    end
                    refauc=sum(mreftac.*dt);
                    
                    if PlotStrFit
                        subplot(3,2,[1 2]);
                        plot(tmidMin,reftac,'k-'); hold on;
                        legendtext={'Cerebellum'};
                    end
                    switch(reg{1})
                        case 'Striatum'
                            tempTac=[];
                            tempVol=[];
                            for subReg={'Caud','Put'}
                                tempTac=[tempTac, TACDATA.(subReg{1})];
                                tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol];
                            end
                            roitac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                        case 'Caudate'
                            roitac=TACDATA.Caud;
                            savestr = 'Caudate';
                        case 'Putamen'
                            roitac=TACDATA.Put;
                            savestr = 'Putamen';
                        case 'Accumbens-area'
                            roitac=TACDATA.Nac;
                            savestr = 'Accumbensarea';
                        case 'Thalamus'
                            roitac=TACDATA.ThalProp;
                            savestr = 'Thalamus';
                        case 'Hippocampus'
                            roitac=TACDATA.Hipp;
                            savestr = 'Hippocampus';
                    end
                    
                    mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2];
                    ASRTM(:,3)=zeros(t_points,1);
                    %                     ASTRM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mroitac(1:k).*dt(1:k));
                    end
                    %LSQ-estimation using lscov
                    [parest se_srtm mse_srtm]   = lscov(ASRTM,roitac);
                    fittac=ASRTM*parest;
                    BP=parest(2)/parest(3)-1;
                    k2p=parest(2)/parest(1);
                    
                    %%%%% DO real SRTM
                    options = optimset('MaxFunEvals',1000);
                    weighs=[0.25*ones(30,1); ones(t_points-30,1)];
                    fobj = @(x) norm((simESRTM_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-roitac).*weighs);
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    fobj = @(x) norm((simESRTM_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-roitac(1:end-11)).*weighs(1:end-11));
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1_bl=parest_srtm(1);
                    k2_bl=parest_srtm(2);
                    BP_bl=parest_srtm(3);
                    modfit_esrtm_bl=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    %     weighs=[0.25*ones(30,1); ones(15+46,1); 100*ones(11,1)];
                    %     fobj = @(x) norm((simESRTM_1_0_0(tmidMin([1:76 92:t_points]),reftac([1:76 92:t_points]),t_points-15,x(1),x(2),x(3)*ones(t_points-15,1))-roitac([1:76 92:t_points])).*weighs([1:76 92:t_points]));
                    %     [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    %     R1_task=parest_srtm(1);
                    %     k2_task=parest_srtm(2);
                    %     BP_task=parest_srtm(3);
                    %     modfit_esrtm_task=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    
                    %%% Do lp-ntPET
                    Alpntpet=zeros(t_points,4);
                    Alpntpet(:,1:3)=ASRTM;
                    best_mse=10^20;
                    best_parest=[];
                    best_se=[];
                    
                    % for estimating the optimal peak and alpha parameters
                    for point_rise=break_point:length(tmid)-1
                        t2p_index = find(tmid>tmid(point_rise)); %times(point_rise,1)
                        if length(t2p_index)>1
                            for t_ind=t2p_index(1):t2p_index(end-1)  %t=1:0.1:20 %
                                for alpha=[0.25 1 4]
                                    t_peak=tmid(t_ind)-tmid(point_rise);   %  times(point_rise,1)
                                    p = [1 alpha tmid(point_rise)+t_peak tmid(point_rise)];
                                    actfun = zeros(size(tmid));
                                    actfun(point_rise:t_points) = gamma_variate_madsen(p,tmid(point_rise:t_points));
                                    roitac_gamma = roitac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(roitac))+roitac_gamma(1:length(roitac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,roitac);
                                    
                                    %estimated model TAC
                                    modelfit = Alpntpet*parest;
                                    
                                    %%%%%% magic %%%%%
                                    if (best_mse > mse)
                                        best_mse = mse;
                                        best_parest=parest; % R1 - k2 - k2a - gamma
                                        best_modelfit=modelfit;
                                        best_actfun=actfun;
                                        best_se=se;
                                    end
                                end
                            end                                        %                        breakpoint2=breakpoint;
                        end
                    end
                    
                    % save these and work out receptor occupancy later in
                    % the script

                    k2=best_parest(2);
                    k2a=best_parest(3);
                    BP_lp=k2/k2a-1;
                    G=best_parest(4);
                    DBP =  k2./(k2a + G*best_actfun) - 1; % dynamic binding potential
                    OCC = 100*(1-DBP/BP_lp); % receptor occupancy
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.' savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.' savestr '=OCC;' ])

                    %                     BP_lp=best_parest(2)/best_parest(3)-1; % is this the BP from lpntPET?
                    %                     BP_lpntPET=(best_parest(2)/(best_parest(3)+(best_parest(4)).*best_actfun))-1; % isn't this the way?
                    %                     BP_bsl=BP_bl; % baseline from SRTM
                    % save BP for later analysis
                    %                     eval(['BP_lp_orig{id,d}.' savestr '=BP_lp;' ])
                    %                     eval(['BP_lp_save{id,d}.' savestr '=BP_lpntPET;' ])
                    %                     eval(['BP_bsl_save{id,d}.' savestr '=BP_bsl;' ])
                    eval(['gamma_lp_save{id,d}.' savestr '=best_parest(4);' ])


                    fprintf(1,'Here\n')

                    if PlotStrFit
                        legendtext{end+1}=[reg{1} ' raw'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{Baseline}'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{All}'];
                        legendtext{end+1}=[reg{1} ' fit lpnt-PET_{All}'];
                        legendtext{end+1}=['Activation function'];
                        SE=abs(BP)*sqrt((se_srtm(2)/parest(2))^2+(se_srtm(3)/parest(3))^2);
                        subplot(3,2,[1 2]);
                        yyaxis left;
                        plot(tmidMin,roitac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
%                         xlim([0 60]);
                        xlabel('Time (min)');
                        ylabel('Radioactivity concentration');
                        legend(legendtext,'Location','south');
                        yyaxis right;
                        plot(tmidMin,best_parest(4)*best_actfun,'c-'); % plot gamma
                        ylabel('Compensatory function','Color','c');
                        ylim([0 4*10^(-4)]);
                        title([Subj{r} ' compartmental fits: BP_{Baseline}=' num2str(BP_bl,'%1.2f') ', BP_{All}=' num2str(BP__,'%1.2f') ', BP_{lp-nt}=' num2str(BP_lp,'%1.2f')]);
                        subplot(3,2,3);
                        plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'go',tmidMin,roitac-modfit_esrtm_bl,'co',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'bo',tmidMin,roitac-modfit_esrtm_bl,'ro',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        [h p]=runstest(roitac-fittac);
                        [h1 p1]=runstest(roitac-modfit_esrtm);
                        [h2 p2]=runstest(roitac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
                        title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
                        xlabel('Time (min)');

                        subplot(3,2,4);
                        y=-ASRTM(:,3)./roitac;
                        x=ASRTM(:,2)./roitac;
                        [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                        [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                        plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
                        title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        xlabel(['\int REF/ROI']);
                        ylabel('\int ROI/ROI')
                        subplot(3,2,[5 6]);
                        plot(tmidMin,roitac./reftac,'ko',[0 180],[0 0],'k--');
                        ylim([0.5 5]);
                        title('Target to reference ratio');
                        xlabel('Time (min)');

                        print('-dpsc2','-append','-bestfit',fullfile('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADC_cleaned/modelling/', [ num2str(IDs(id)) num2str(d) '_TAC_Fit_lpntpet_logan_20221109.ps']));
                        close(gcf)
                        BPdata.BP_mrtm(Subj{r})=BP;
                        BPdata.BP_srtm(Subj{r})=BP__;
                        BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                        BPdata.BP_lpnt(Subj{r})=BP_lp;
                        BPdata.BP_logan(Subj{r})=pp(1)-1;
                        continue;
                    end
                end
                
                
            end
            close all
%             keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI
        end
    end
end


disp('done')

%% calculate occupancy per ROI

% occupancies = [];
reg={'Striatum','Caudate','Putamen','Accumbensarea','Thalamus','Hippocampus'};
for id = 1:length(IDs)
    for d=1:2
        if days(id,d) == 0
            disp('no measurement for this date')
            %             occupancies{id,d} = [];
        else
            figure;
            for r=[1 4 5]%:length(reg)
                eval(['plot(Occupancy{id,d}.' reg{r} ',''LineWidth'',3);hold on'])
            end
            title(['ID: ' num2str(IDs(id)) ' Session ' num2str(d)],'FontSize',30); ylim([-5 100]); %xlim([70 102]);
            xline(102-15,'--','Linewidth',2) % task onset
            xlabel('Frame','Fontsize',20);
            ylabel('Percentage (%)','Fontsize',20)
            legend({'Striatum','Accumbensarea','Hippocampus','Task Onset'},'Location','northwest','Fontsize',15);
            set(gca,'fontsize',20)
        end
    end
end

%% all ROIs

load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')

for id = 19:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg=ROI%{'CerC','Put','Caud','Nac','Hipp','ThalProp'}
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
            end
            
            Lengths=[10*ones(30,1); 60*ones(55,1)]; % time windows
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
            %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
            
            tmid=mean(times,2);
            tmidMin=tmid/60;
            t_points    = length(tmid);
            dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
            break_point=find(times(:,1)>=115*60,1,'first'); %% TIme of activation start
            
            %%%%%%%%%%%
            PlotStrFit=1;
            %%%%%%%%%%%
            
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg=ROI
                    if PlotStrFit
                        figure('Position',[100 100 800 1200],'Visible','on'); hold on;
                        spidx=1;
                    end
                    reftac=TACDATA.CerC;
                    mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
                    %% Set the SRTM part of A
                    ASRTM = zeros(t_points ,3);
                    ASRTM(:,1)  = reftac;
                    for k = 1:t_points
                        ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
                    end
                    refauc=sum(mreftac.*dt);
                    
                    if PlotStrFit
                        subplot(3,2,[1 2]);
                        plot(tmidMin,reftac,'k-'); hold on;
                        legendtext={'Cerebellum'};
                    end
                    switch(reg{1})
                        %                         case 'Striatum'
                        
                        case 'Caud'
                            roitac=TACDATA.Caud;
                            savestr = 'Caudate';
                        case 'Put'
                            roitac=TACDATA.Put;
                            savestr = 'Putamen';
                            
                            tempTac=[];
                            tempVol=[];
                            for subReg={'Caud','Put'}
                                tempTac=[tempTac, TACDATA.(subReg{1})];
                                tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol];
                            end
                            roitac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                            
                        case 'Nac'
                            roitac=TACDATA.Nac;
                            savestr = 'Nac';
                        case 'ThalProp'
                            roitac=TACDATA.ThalProp;
                            savestr = 'ThalProp';
                        case 'Hipp'
                            roitac=TACDATA.Hipp;
                            savestr = 'Hipp';
                        case 'bankssts'
                            roitac=TACDATA.bankssts;
                            savestr = 'bankssts';
                        case 'caudalanteriorcingulate'
                            roitac=TACDATA.caudalanteriorcingulate;
                            savestr = 'caudalanteriorcingulate';
                        case 'caudalmiddlefrontal'
                            roitac=TACDATA.caudalmiddlefrontal;
                            savestr = 'caudalmiddlefrontal';
                        case 'cuneus'
                            roitac=TACDATA.cuneus;
                            savestr = 'cuneus';
                        case 'entorhinal'
                            roitac=TACDATA.entorhinal;
                            savestr = 'entorhinal';
                        case 'fusiform'
                            roitac=TACDATA.fusiform;
                            savestr = 'fusiform';
                        case 'inferiorparietal'
                            roitac=TACDATA.inferiorparietal;
                            savestr = 'inferiorparietal';
                        case 'inferiortemporal'
                            roitac=TACDATA.inferiortemporal;
                            savestr = 'inferiortemporal';
                        case 'isthmuscingulate'
                            roitac=TACDATA.isthmuscingulate;
                            savestr = 'isthmuscingulate';
                        case 'lateraloccipital'
                            roitac=TACDATA.lateraloccipital;
                            savestr = 'lateraloccipital';
                        case 'lateralorbitofrontal'
                            roitac=TACDATA.lateralorbitofrontal;
                            savestr = 'lateralorbitofrontal';
                        case 'lingual'
                            roitac=TACDATA.lingual;
                            savestr = 'lingual';
                        case 'medialorbitofrontal'
                            roitac=TACDATA.medialorbitofrontal;
                            savestr = 'medialorbitofrontal';
                        case 'parahippocampal'
                            roitac=TACDATA.parahippocampal;
                            savestr = 'parahippocampal';
                        case 'paracentral'
                            roitac=TACDATA.paracentral;
                            savestr = 'paracentral';
                        case 'precuneus'
                            roitac=TACDATA.precuneus;
                            savestr = 'precuneus';
                        case 'parsopercularis'
                            roitac=TACDATA.parsopercularis;
                            savestr = 'parsopercularis';
                        case 'parsorbitalis'
                            roitac=TACDATA.parsorbitalis;
                            savestr = 'parsorbitalis';
                        case 'parstriangularis'
                            roitac=TACDATA.parstriangularis;
                            savestr = 'parstriangularis';
                        case 'pericalcarine'
                            roitac=TACDATA.pericalcarine;
                            savestr = 'pericalcarine';
                        case 'postcentral'
                            roitac=TACDATA.postcentral;
                            savestr = 'postcentral';
                        case 'posteriorcingulate'
                            roitac=TACDATA.posteriorcingulate;
                            savestr = 'posteriorcingulate';
                        case 'precentral'
                            roitac=TACDATA.precentral;
                            savestr = 'precentral';
                        case 'rostralanteriorcingulate'
                            roitac=TACDATA.rostralanteriorcingulate;
                            savestr = 'rostralanteriorcingulate';
                        case 'rostralmiddlefrontal'
                            roitac=TACDATA.rostralmiddlefrontal;
                            savestr = 'rostralmiddlefrontal';
                        case 'superiorfrontal'
                            roitac=TACDATA.superiorfrontal;
                            savestr = 'superiorfrontal';
                        case 'superiorparietal'
                            roitac=TACDATA.superiorparietal;
                            savestr = 'superiorparietal';
                        case 'superiortemporal'
                            roitac=TACDATA.superiortemporal;
                            savestr = 'superiortemporal';
                        case 'supramarginal'
                            roitac=TACDATA.supramarginal;
                            savestr = 'supramarginal';
                        case 'frontalpole'
                            roitac=TACDATA.frontalpole;
                            savestr = 'frontalpole';
                        case 'temporalpole'
                            roitac=TACDATA.temporalpole;
                            savestr = 'temporalpole';
                        case 'transversetemporal'
                            roitac=TACDATA.transversetemporal;
                            savestr = 'transversetemporal';
                        case 'insula'
                            roitac=TACDATA.insula;
                            savestr = 'insula';
                        case 'CerC'
                            roitac=TACDATA.CerC;
                            savestr = 'CerC';
                        case 'Pall'
                            roitac=TACDATA.Pall;
                            savestr = 'Pall';
                        case 'Amy'
                            roitac=TACDATA.Amy;
                            savestr = 'Amy';
                        case 'SN'
                            roitac=TACDATA.SN;
                            savestr = 'SubstantiaNigra';
                        case 'LC'
                            roitac=TACDATA.LC;
                            savestr = 'LocusCoeruleus';
                        case 'VTA'
                            roitac=TACDATA.VTA;
                            savestr = 'VentralTegmentalArea';
                    end
                    
                    mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2];
                    ASRTM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mroitac(1:k).*dt(1:k));
                    end
                    %LSQ-estimation using lscov
                    [parest se_srtm mse_srtm]   = lscov(ASRTM,roitac);
                    fittac=ASRTM*parest;
                    BP=parest(2)/parest(3)-1;
                    k2p=parest(2)/parest(1);
                    
                    %%%%% DO real SRTM
                    options = optimset('MaxFunEvals',1000);
                    weighs=[0.25*ones(30,1); ones(t_points-30,1)];
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-roitac).*weighs);
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    %                     modfit_esrtm=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-roitac(1:end-11)).*weighs(1:end-11));
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1_bl=parest_srtm(1);
                    k2_bl=parest_srtm(2);
                    BP_bl=parest_srtm(3);
                    modfit_esrtm_bl=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    %                     modfit_esrtm_bl=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    %     weighs=[0.25*ones(30,1); ones(15+46,1); 100*ones(11,1)];
                    %     fobj = @(x) norm((simESRTM_1_0_0(tmidMin([1:76 92:t_points]),reftac([1:76 92:t_points]),t_points-15,x(1),x(2),x(3)*ones(t_points-15,1))-roitac([1:76 92:t_points])).*weighs([1:76 92:t_points]));
                    %     [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    %     R1_task=parest_srtm(1);
                    %     k2_task=parest_srtm(2);
                    %     BP_task=parest_srtm(3);
                    %     modfit_esrtm_task=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    
                    %%% Do lp-ntPET
                    Alpntpet=zeros(t_points,4);
                    Alpntpet(:,1:3)=ASRTM;
                    best_mse=10^20;
                    best_parest=[];
                    best_se=[];
                    
                    % for estimating the optimal peak and alpha parameters
                    for point_rise=break_point:length(tmid)-1
                        t2p_index = find(tmid>tmid(point_rise)); %times(point_rise,1)
                        if length(t2p_index)>1
                            for t_ind=t2p_index(1):t2p_index(end-1)  %t=1:0.1:20 %
                                for alpha=[0.25 1 4]
                                    t_peak=tmid(t_ind)-tmid(point_rise);   %  times(point_rise,1)
                                    p = [1 alpha tmid(point_rise)+t_peak tmid(point_rise)];
                                    actfun = zeros(size(tmid));
                                    actfun(point_rise:t_points) = gamma_variate_madsen(p,tmid(point_rise:t_points));
                                    roitac_gamma = roitac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(roitac))+roitac_gamma(1:length(roitac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,roitac);
                                    
                                    %estimated model TAC
                                    modelfit = Alpntpet*parest;
                                    
                                    %%%%%% magic %%%%%
                                    if (best_mse > mse)
                                        best_mse = mse;
                                        best_parest=parest; % R1 - k2 - k2a - gamma
                                        best_modelfit=modelfit;
                                        best_actfun=actfun;
                                        best_se=se;
                                    end
                                end
                            end                                        %                        breakpoint2=breakpoint;
                        end
                    end
                    
                    % save these and work out receptor occupancy later in
                    % the script
                    
                    k2=best_parest(2);
                    k2a=best_parest(3);
                    BP_lp=k2/k2a-1;
                    BP_srtm=BP__;
                    BP_srtm_bsl=BP_bl;
                    G=best_parest(4);
                    BPND = ((R1__*k2p)/k2a)-1;
                    DBP =  k2./(k2a + G*best_actfun) - 1; % dynamic binding potential
                    OCC = 100*(1-DBP/BP_lp); % receptor occupancy
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.' savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.' savestr '=OCC;' ])
                    eval(['BPND_save{id,d}.' savestr '=BPND;' ])
                    eval(['BP_srtm_save{id,d}.' savestr '=BP_srtm;' ])
                    eval(['BP_srtm_Bsl_save{id,d}.' savestr '=BP_srtm_bsl;' ])
                    
                    %                     BP_lp=best_parest(2)/best_parest(3)-1; % is this the BP from lpntPET?
                    %                     BP_lpntPET=(best_parest(2)/(best_parest(3)+(best_parest(4)).*best_actfun))-1; % isn't this the way?
                    %                     BP_bsl=BP_bl; % baseline from SRTM
                    % save BP for later analysis
                    %                     eval(['BP_lp_orig{id,d}.' savestr '=BP_lp;' ])
                    %                     eval(['BP_lp_save{id,d}.' savestr '=BP_lpntPET;' ])
                    %                     eval(['BP_bsl_save{id,d}.' savestr '=BP_bsl;' ])
                    eval(['gamma_lp_save{id,d}.' savestr '=best_parest(4);' ])
                    %                     BPdata.BP_mrtm(Subj{r})=BP;
                    %                     BPdata.BP_srtm(Subj{r})=BP__;
                    %                     BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                    %                     BPdata.BP_lpnt(Subj{r})=BP_lp;
                    %                     BPdata.BP_logan(Subj{r})=pp(1)-1;
                    
                    %                     fprintf(1,'Here\n')
                    
                    if PlotStrFit
                        %                         legendtext{end+1}=[reg{1} ' raw'];
                        %                         legendtext{end+1}=[reg{1} ' fit SRTM_{Baseline}'];
                        %                         legendtext{end+1}=[reg{1} ' fit SRTM_{All}'];
                        %                         legendtext{end+1}=[reg{1} ' fit lpnt-PET_{All}'];
                        %                         legendtext{end+1}=['Activation function'];
                        SE=abs(BP)*sqrt((se_srtm(2)/parest(2))^2+(se_srtm(3)/parest(3))^2);
                        %                         subplot(3,2,[1 2]);
                        %                         yyaxis left;
                        %                         plot(tmidMin,roitac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
                        %xlim([0 60]);
                        %                         xlabel('Time (min)');
                        %                         ylabel('Radioactivity concentration');
                        %                         legend(legendtext,'Location','south');
                        %                         yyaxis right;
                        %                         plot(tmidMin,best_parest(4)*best_actfun); % plot gamma
                        %                         ylabel('Compensatory function');
                        %ylim([0 4*10^(-4)]);
                        %                         title([Subj{r} ' compartmental fits: BP_{Baseline}=' num2str(BP_bl,'%1.2f') ', BP_{All}=' num2str(BP__,'%1.2f') ', BP_{lp-nt}=' num2str(BP_lp,'%1.2f')]);
                        %                         subplot(3,2,3);
                        %         plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'go',tmidMin,roitac-modfit_esrtm_bl,'co',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        %                         plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'bo',tmidMin,roitac-modfit_esrtm_bl,'ro',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        [h p]=runstest(roitac-fittac);
                        [h1 p1]=runstest(roitac-modfit_esrtm);
                        [h2 p2]=runstest(roitac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
                        %                         title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
                        %                         xlabel('Time (min)');
                        
                        %                         subplot(3,2,4);
                        y=-ASRTM(:,3)./roitac;
                        x=ASRTM(:,2)./roitac;
                        [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                        [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                        %                         plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
                        %                         title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        %                         xlabel(['\int REF/ROI']);
                        %                         ylabel('\int ROI/ROI')
                        %                         subplot(3,2,[5 6]);
                        %                         plot(tmidMin,roitac./reftac,'ko',[0 180],[0 0],'k--');
                        %ylim([0.5 5]);
                        %                         title('Target to reference ratio');
                        %                         xlabel('Time (min)');
                        
                        %                         print('-dpsc2','-append','-bestfit',fullfile('/Users/alex/Documents/RewardTask_ANTscoreg/', [ num2str(IDs(id)) num2str(d) '_TAC_Fit_lpntpet_logan_framesDropped.ps']));
                        %                         close(gcf)
                        BPdata.BP_mrtm(Subj{r})=BP;
                        BPdata.BP_srtm(Subj{r})=BP__;
                        BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                        BPdata.BP_lpnt(Subj{r})=BP_lp;
                        BPdata.BP_logan(Subj{r})=pp(1)-1;
                        
                        % save
                        BPdataSave{id,d}.(reg{1}).BP_mrtm=BP;
                        BPdataSave{id,d}.(reg{1}).BP_srtm=BP__;
                        BPdataSave{id,d}.(reg{1}).BP_srtm_bl=BP_bl;
                        BPdataSave{id,d}.(reg{1}).BP_lpnt=BP_lp;
                        BPdataSave{id,d}.(reg{1}).BP_logan=pp(1)-1;
                        
                        close all
                        continue;
                    end
                    close all
                end
                close all
                
            end
            close all
            keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI BP_srtm_save BPND_save BPdataSave
        end
    end
end


disp('all ROIs done')

%% make maps - task

clc;
diary srtm_task.txt

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            clear SRTMimg baseimage FSlabel existingLabels
            load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')
            
            baseimage = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']));
            SRTMimg = zeros(size(baseimage));
            existingLabels = unique(baseimage);
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    
                    fprintf(['**' FSlabels1{labels,1}{1} '**\n'])
                    
                    % left
                    disp('left')
                    clear indL
                    indL=find(baseimage==FSlabels1{labels,2});
                    SRTMimg(indL)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;
                    
                    % right
                    disp('right')
                    clear indR
                    indR=find(baseimage==FSlabels1{labels,3});
                    SRTMimg(indR)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;
                    
                    disp(['both, ' num2str(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)])
                    
                end
            end
            
            
            %SRTM
            hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
            hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_SRTM_task_onTemplate.nii'];
            hdr.dim = size(SRTMimg);
            hdr = rmfield(hdr,'pinfo');
            hdr.nii = spm_write_vol(hdr,SRTMimg);
            
            fprintf(['\n\n\n'])
        end
    end
end


diary off
type srtm_task.txt

%% make maps - baseline

clc;
diary srtm_bsl.txt

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            clear SRTMimg baseimage FSlabel existingLabels
            load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')
            
            baseimage = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']));
            SRTMimg = zeros(size(baseimage));
            existingLabels = unique(baseimage);
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    
                    fprintf(['**' FSlabels1{labels,1}{1} '**\n'])
                    
                    % left
                    disp('left')
                    clear indL
                    indL=find(baseimage==FSlabels1{labels,2});
                    SRTMimg(indL)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;
                    
                    % right
                    disp('right')
                    clear indR
                    indR=find(baseimage==FSlabels1{labels,3});
                    SRTMimg(indR)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;
                    
                    disp(['both, ' num2str(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)])
                    
                end
            end
            
            
            %SRTM
            hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' num2str(IDs(id)) num2str(d) '/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
            hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_SRTM_baseline_onTemplate.nii'];
            hdr.dim = size(SRTMimg);
            hdr = rmfield(hdr,'pinfo');
            hdr.nii = spm_write_vol(hdr,SRTMimg);
            
            fprintf(['\n\n\n'])
        end
    end
end


diary off
type srtm_bsl.txt

%% mean of baseline SRTM BP

% collapsed imcalc codes
for i1=1
% strtmp=['i1+'];
% for k=2:40
%     strtmp=strcat(strtmp,['i' num2str(k) '+']);
% end
% 
% matlabbatch{1}.spm.util.imcalc.input = {
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40011/40011_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40012/40012_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40021/40021_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40022/40022_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40031/40031_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40041/40041_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40042/40042_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40051/40051_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40052/40052_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40062/40062_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40071/40071_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40081/40081_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40082/40082_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40092/40092_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40101/40101_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40102/40102_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40111/40111_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40112/40112_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40121/40121_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40131/40131_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40132/40132_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40142/40142_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40151/40151_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40152/40152_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40161/40161_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40162/40162_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40171/40171_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40172/40172_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40181/40181_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40182/40182_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40191/40191_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40201/40201_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40202/40202_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40211/40211_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40212/40212_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40222/40222_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40231/40231_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40241/40241_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40262/40262_SRTM_task_onTemplate.nii,1'
%                                         };
% matlabbatch{1}.spm.util.imcalc.output = 'meanBPtask-SRTM_on_Template';
% matlabbatch{1}.spm.util.imcalc.outdir = {'/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi'};
% matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26+i27+i28+i29+i30+i31+i32+i33+i34+i35+i36+i37+i38+i39)/39';
% matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
end



% mean of BPtask SRTM

cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')
            cnt=cnt+1;
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                                        
                    eval(['BPbsl_SRTM_' FSlabels1{labels,1}{1} '(cnt,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;'])
                  
                end
            end
        end
    end
end





for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

        eval(['meanBPbsl_SRTM_' FSlabels1{labels,1}{1} '=nanmean(BPbsl_SRTM_' FSlabels1{labels,1}{1} ');'])

    end
end


clc;
% diary srtm_task.txt

%             fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))

clear SRTMimg baseimage FSlabel existingLabels
load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')

baseimage = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']));
SRTMimg = zeros(size(baseimage));
existingLabels = unique(baseimage);

for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

%         fprintf(['**' FSlabels1{labels,1}{1} '**\n'])

        % left
%         disp('left')
        clear indL
        indL=find(baseimage==FSlabels1{labels,2});
        SRTMimg(indL)=eval(strcat('meanBPbsl_SRTM_',FSlabels1{labels,1}{1}));

        % right
%         disp('right')
        clear indR
        indR=find(baseimage==FSlabels1{labels,3});
        SRTMimg(indR)=eval(strcat('meanBPbsl_SRTM_',FSlabels1{labels,1}{1}));

%         disp(['both, ' num2str(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)])

    end
end


%SRTM
hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/meanBPbsl_SRTM_onTemplate.nii'];
hdr.dim = size(SRTMimg);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,SRTMimg);


%% mean of task SRTM BP

% collapsed imcalc codes
for i1=1
% strtmp=['i1+'];
% for k=2:40
%     strtmp=strcat(strtmp,['i' num2str(k) '+']);
% end
% 
% matlabbatch{1}.spm.util.imcalc.input = {
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40011/40011_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40012/40012_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40021/40021_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40022/40022_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40031/40031_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40041/40041_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40042/40042_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40051/40051_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40052/40052_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40062/40062_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40071/40071_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40081/40081_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40082/40082_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40092/40092_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40101/40101_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40102/40102_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40111/40111_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40112/40112_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40121/40121_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40131/40131_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40132/40132_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40142/40142_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40151/40151_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40152/40152_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40161/40161_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40162/40162_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40171/40171_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40172/40172_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40181/40181_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40182/40182_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40191/40191_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40201/40201_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40202/40202_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40211/40211_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40212/40212_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40222/40222_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40231/40231_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40241/40241_SRTM_task_onTemplate.nii,1'
%                                         '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/40262/40262_SRTM_task_onTemplate.nii,1'
%                                         };
% matlabbatch{1}.spm.util.imcalc.output = 'meanBPtask-SRTM_on_Template';
% matlabbatch{1}.spm.util.imcalc.outdir = {'/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi'};
% matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26+i27+i28+i29+i30+i31+i32+i33+i34+i35+i36+i37+i38+i39)/39';
% matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
end



% mean of BPtask SRTM

cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')
            cnt=cnt+1;
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                                        
                    eval(['BPtask_SRTM_' FSlabels1{labels,1}{1} '(cnt,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;'])
                  
                end
            end
        end
    end
end





for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

        eval(['meanBPtask_SRTM_' FSlabels1{labels,1}{1} '=nanmean(BPtask_SRTM_' FSlabels1{labels,1}{1} ');'])

    end
end


clc;
% diary srtm_task.txt

%             fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))

clear SRTMimg baseimage FSlabel existingLabels
load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')

baseimage = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']));
SRTMimg = zeros(size(baseimage));
existingLabels = unique(baseimage);

for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

%         fprintf(['**' FSlabels1{labels,1}{1} '**\n'])

        % left
%         disp('left')
        clear indL
        indL=find(baseimage==FSlabels1{labels,2});
        SRTMimg(indL)=eval(strcat('meanBPtask_SRTM_',FSlabels1{labels,1}{1}));

        % right
%         disp('right')
        clear indR
        indR=find(baseimage==FSlabels1{labels,3});
        SRTMimg(indR)=eval(strcat('meanBPtask_SRTM_',FSlabels1{labels,1}{1}));

%         disp(['both, ' num2str(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)])

    end
end


%SRTM
hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/meanBPtask_SRTM_onTemplate.nii'];
hdr.dim = size(SRTMimg);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,SRTMimg);


%% make a table, with everything as single sessions

tablestr_task='BPtask_SRTM_bankssts';

for labels=3:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else
        tablestr_task=strcat(tablestr_task,',BPtask_SRTM_',FSlabels1{labels,1});

    end
end

tablestr_bsl='BPbsl_SRTM_bankssts';

for labels=3:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else
        tablestr_bsl=strcat(tablestr_bsl,',BPbsl_SRTM_',FSlabels1{labels,1});

    end
end



% ID=IDs';
eval(strcat('MRPET_BP_SRTM_Table = table(ID,', tablestr_bsl, ',',  tablestr_task, ')'))
writetable(MRPET_BP_SRTM_Table,'MRPET_BPtable_SRTM_single.xls')



%% calculate occupancy per ROI

% occupancies = [];
reg={'Striatum','Caudate','Putamen','Accumbensarea','Thalamus','Hippocampus'};
for id = 1:length(IDs)
    for d=1:2
        if days(id,d) == 0
            disp('no measurement for this date')
            %             occupancies{id,d} = [];
        else
            figure;
            for r=[1 4 5]%:length(reg)
                eval(['plot(Occupancy{id,d}.' reg{r} ',''LineWidth'',3);hold on'])
            end
            title(['ID: ' num2str(IDs(id)) ' Session ' num2str(d)],'FontSize',30); ylim([-5 100]); xlim([70 102]);
            xline(102-15,'--','Linewidth',2) % task onset
            xlabel('Frame','Fontsize',20);
            ylabel('Percentage (%)','Fontsize',20)
            legend({'Striatum','Accumbensarea','Hippocampus','Task Onset'},'Location','northwest','Fontsize',15);
            set(gca,'fontsize',20)
        end
    end
end


%% peak DBP

reg={'Striatum','Caudate','Putamen','Accumbensarea','Thalamus','Hippocampus'};
peakDBP=[];
for id = 1:length(IDs)
    for d=1:2
        if days(id,d) == 0
            disp('no measurement for this date')
        else
            for r=[1 4 5]
                eval(['peakDBP{id,d}.' reg{r} '=max(diff(DBP_save{id,d}.' reg{r} ')/DBP_save{id,d}.' reg{r} '(1))']);
            end
        end
    end
end

%



%% normalise the values

% load('/Users/alex/Dropbox/literatures_IKND/FS_labels_LCSNVTA.mat')

% ROIwise normalisation

ZScored=[];
for id=1:length(IDs)
    
    for d=1:2
        if days(id,d) == 0
            disp('no measurement for this date')
        else
            clear TaskPET WBmask numvoxWB tmp clear refmean
            
            PETtask = [paths.parent num2str(IDs(id)) num2str(days(id,d)) '/coreg_MTonT1.nii'];
            maskpt2  = [paths.parent num2str(IDs(id)) num2str(days(id,d)) '/BrainMask_pt2.nii'];
            
            WBmask = load_nii(maskpt2);
            numvoxWB = sum(WBmask.img(:)>0);
            tmp = load_nii(PETtask);
            WBmaskind=WBmask.img(:)==0;
            
            dummy=[];
            for l1=1:size(tmp.img,4)
                clear tmp2
                tmp2=tmp.img(:,:,:,l1);
                tmp2(WBmaskind)=NaN;
                dummy(:,:,:,l1)=tmp2;
            end
            temp=size(dummy);
            
            TaskPET=reshape(dummy,prod(temp(1:3)),temp(4));
            
            WB_sd = nanstd(TaskPET)';
            refmean = TACs{id,d}.TACDATA_Task.CerC.Bilateral.tac;
            
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    ZScored{id,d}.(FSlabels1{labels,1})=...
                        (TACs{id,d}.TACDATA_Task.(FSlabels1{labels,1}).Bilateral.tac - refmean )./ WB_sd;
                end
            end
            
        end
    end
end

save('/Users/alex/Dropbox/paperwriting/MRPET/data/ZscoredMax.mat','ZScored')

%%

% make maps
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            clear ZScoreImg baseimage FSlabel existingLabels
            load('/Users/alex/Dropbox/literatures_IKND/FS_labels_LCSNVTA.mat')
            
            baseimage = spm_read_vols(spm_vol(['/Volumes/ALEX_DATA1/PET_data/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii']));
            ZScoreImg = zeros(size(baseimage));
            existingLabels = unique(baseimage);
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    
                    fprintf(['**' FSlabels1{labels,1}{1} '**\n'])
                    
                    % left
                    disp('left')
                    clear indL
                    indL=find(baseimage==FSlabels1{labels,2});
                    ZScoreImg(indL)=max(ZScored{id,d}.(FSlabels1{labels,1}));
                    
                    % right
                    disp('right')
                    clear indR
                    indR=find(baseimage==FSlabels1{labels,3});
                    ZScoreImg(indR)=max(ZScored{id,d}.(FSlabels1{labels,1}));
                    
                    
                end
            end
            
            
            % Z Scored, MAX value
            hdr = spm_vol(['/Volumes/ALEX_DATA1/PET_data/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii']); % pick just any header from a file
            hdr.fname = ['/Volumes/ALEX_DATA1/PET_data/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_ZScoreMax.nii'];
            hdr.dim = size(ZScoreImg);
            hdr = rmfield(hdr,'pinfo');
            hdr.nii = spm_write_vol(hdr,ZScoreImg);
            
            fprintf(['\n\n\n'])
        end
    end
end


%% calculate variance

% basal ganglia: caudate, putamen, pallidus, NAcc, SN, VTA, thalamus ,  *LC


load('/Users/alex/Dropbox/paperwriting/MRPET/data/ZscoredMax.mat')
load('/Users/alex/Dropbox/literatures_IKND/FS_labels_LCSNVTA.mat')

% gather the data
Caud_raw=[]; Putamen_raw=[]; Pallidus_raw=[]; NAcc_raw=[];
SN_raw=[]; VTA_raw=[]; Thalamus_raw=[];
Caud_BP=[]; Putamen_BP=[]; Pallidus_BP=[]; NAcc_BP=[];
SN_BP=[]; VTA_BP=[]; Thalamus_BP=[];
Caud_Z=[]; Putamen_Z=[]; Pallidus_Z=[]; NAcc_Z=[];
SN_Z=[]; VTA_Z=[]; Thalamus_Z=[];

Subj_raw=[]; Subj_BP=[]; Subj_Z=[];

for id=1:length(IDs)
    
    for d = 1:2
        if days(id,d) == 0
                        
            Caud_raw(id,d)=NaN; Putamen_raw(id,d)=NaN; Pallidus_raw(id,d)=NaN; NAcc_raw(id,d)=NaN;
            SN_raw(id,d)=NaN; VTA_raw(id,d)=NaN; Thalamus_raw(id,d)=NaN;
            Caud_BP(id,d)=NaN; Putamen_BP(id,d)=NaN; Pallidus_BP(id,d)=NaN; NAcc_BP(id,d)=NaN;
            SN_BP(id,d)=NaN; VTA_BP(id,d)=NaN; Thalamus_BP(id,d)=NaN;
            Caud_Z(id,d)=NaN; Putamen_Z(id,d)=NaN; Pallidus_Z(id,d)=NaN; NAcc_Z(id,d)=NaN;
            SN_Z(id,d)=NaN; VTA_Z(id,d)=NaN; Thalamus_Z(id,d)=NaN;
            
        else
            
            % raw TAC
            Caud_raw(id,d)        = max(TACs{id,d}.TACDATA_Task.Caud.Bilateral.tac);
            Putamen_raw(id,d)     = max(TACs{id,d}.TACDATA_Task.Put.Bilateral.tac);
            Pallidus_raw(id,d)    = max(TACs{id,d}.TACDATA_Task.Pall.Bilateral.tac);
            NAcc_raw(id,d)        = max(TACs{id,d}.TACDATA_Task.Nac.Bilateral.tac);
            SN_raw(id,d)          = max(TACs{id,d}.TACDATA_Task.SN.Bilateral.tac);
            VTA_raw(id,d)         = max(TACs{id,d}.TACDATA_Task.VTA.Bilateral.tac);
            Thalamus_raw(id,d)    = max(TACs{id,d}.TACDATA_Task.ThalProp.Bilateral.tac);
            
            % BP_srtm
            Caud_BP(id,d)        = max(BPdataSave{id,d}.Caud.BP_srtm);
            Putamen_BP(id,d)     = max(BPdataSave{id,d}.Put.BP_srtm);
            Pallidus_BP(id,d)    = max(BPdataSave{id,d}.Pall.BP_srtm);
            NAcc_BP(id,d)        = max(BPdataSave{id,d}.Nac.BP_srtm);
            SN_BP(id,d)          = max(BPdataSave{id,d}.SN.BP_srtm);
            VTA_BP(id,d)         = max(BPdataSave{id,d}.VTA.BP_srtm);
            Thalamus_BP(id,d)    = max(BPdataSave{id,d}.ThalProp.BP_srtm);
            
            % Zscored with cerebellar ctx mean
            Caud_Z(id,d)        = max(ZScored{id,d}.Caud);
            Putamen_Z(id,d)     = max(ZScored{id,d}.Put);
            Pallidus_Z(id,d)    = max(ZScored{id,d}.Pall);
            NAcc_Z(id,d)        = max(ZScored{id,d}.Nac);
            SN_Z(id,d)          = max(ZScored{id,d}.SN);
            VTA_Z(id,d)         = max(ZScored{id,d}.VTA);
            Thalamus_Z(id,d)    = max(ZScored{id,d}.ThalProp);
            
            
        end
    end
    
end

Subj_raw.highDA=[]; Subj_BP.highDA=[]; Subj_Z.highDA=[];
Subj_raw.lowDA=[]; Subj_BP.lowDA=[]; Subj_Z.lowDA=[];

Subj_raw.highDA   =[Caud_raw(:,1) Putamen_raw(:,1) Pallidus_raw(:,1) NAcc_raw(:,1) SN_raw(:,1) VTA_raw(:,1) Thalamus_raw(:,1)]
Subj_raw.lowDA    =[Caud_raw(:,2) Putamen_raw(:,2) Pallidus_raw(:,2) NAcc_raw(:,2) SN_raw(:,2) VTA_raw(:,2) Thalamus_raw(:,2)]

Subj_BP.highDA   =[Caud_BP(:,1) Putamen_BP(:,1) Pallidus_BP(:,1) NAcc_BP(:,1) SN_BP(:,1) VTA_BP(:,1) Thalamus_BP(:,1)]
Subj_BP.lowDA    =[Caud_BP(:,2) Putamen_BP(:,2) Pallidus_BP(:,2) NAcc_BP(:,2) SN_BP(:,2) VTA_BP(:,2) Thalamus_BP(:,2)]

Subj_Z.highDA   =[Caud_Z(:,1) Putamen_Z(:,1) Pallidus_Z(:,1) NAcc_Z(:,1) SN_Z(:,1) VTA_Z(:,1) Thalamus_Z(:,1)]
Subj_Z.lowDA    =[Caud_Z(:,2) Putamen_Z(:,2) Pallidus_Z(:,2) NAcc_Z(:,2) SN_Z(:,2) VTA_Z(:,2) Thalamus_Z(:,2)]


% calculate variance

% within subject
for id=1:11
    
    withinSubjVariance_BP.highDA(id,1) = nanvar([Caud_BP(id,1) Putamen_BP(id,1) Pallidus_BP(id,1) NAcc_BP(id,1) SN_BP(id,1) VTA_BP(id,1) Thalamus_BP(id,1)]);
    withinSubjVariance_BP.lowDA(id,1)  = nanvar([Caud_BP(id,2) Putamen_BP(id,2) Pallidus_BP(id,2) NAcc_BP(id,2) SN_BP(id,2) VTA_BP(id,2) Thalamus_BP(id,2)]);

    withinSubjVariance_Z.highDA(id,1) = nanvar([Caud_Z(id,1) Putamen_Z(id,1) Pallidus_Z(id,1) NAcc_Z(id,1) SN_Z(id,1) VTA_Z(id,1) Thalamus_Z(id,1)]);
    withinSubjVariance_Z.lowDA(id,1)  = nanvar([Caud_Z(id,2) Putamen_Z(id,2) Pallidus_Z(id,2) NAcc_Z(id,2) SN_Z(id,2) VTA_Z(id,2) Thalamus_Z(id,2)]);

    withinSubjSD_BP.highDA(id,1) = nanstd([Caud_BP(id,1) Putamen_BP(id,1) Pallidus_BP(id,1) NAcc_BP(id,1) SN_BP(id,1) VTA_BP(id,1) Thalamus_BP(id,1)]);
    withinSubjSD_BP.lowDA(id,1)  = nanstd([Caud_BP(id,2) Putamen_BP(id,2) Pallidus_BP(id,2) NAcc_BP(id,2) SN_BP(id,2) VTA_BP(id,2) Thalamus_BP(id,2)]);

    withinSubjSD_Z.highDA(id,1) = nanstd([Caud_Z(id,1) Putamen_Z(id,1) Pallidus_Z(id,1) NAcc_Z(id,1) SN_Z(id,1) VTA_Z(id,1) Thalamus_Z(id,1)]);
    withinSubjSD_Z.lowDA(id,1)  = nanstd([Caud_Z(id,2) Putamen_Z(id,2) Pallidus_Z(id,2) NAcc_Z(id,2) SN_Z(id,2) VTA_Z(id,2) Thalamus_Z(id,2)]);
    
    withinSubjVariance_BP.all(id,1) = nanvar([Caud_BP(id,1) Putamen_BP(id,1) Pallidus_BP(id,1) NAcc_BP(id,1) SN_BP(id,1) VTA_BP(id,1) Thalamus_BP(id,1) ...
        [Caud_BP(id,2) Putamen_BP(id,2) Pallidus_BP(id,2) NAcc_BP(id,2) SN_BP(id,2) VTA_BP(id,2) Thalamus_BP(id,2)]]);
    withinSubjSD_BP.all(id,1) = nanstd([Caud_BP(id,1) Putamen_BP(id,1) Pallidus_BP(id,1) NAcc_BP(id,1) SN_BP(id,1) VTA_BP(id,1) Thalamus_BP(id,1) ...
        [Caud_BP(id,2) Putamen_BP(id,2) Pallidus_BP(id,2) NAcc_BP(id,2) SN_BP(id,2) VTA_BP(id,2) Thalamus_BP(id,2)]]);

    withinSubjVariance_Z.all(id,1) = nanvar([Caud_Z(id,1) Putamen_Z(id,1) Pallidus_Z(id,1) NAcc_Z(id,1) SN_Z(id,1) VTA_Z(id,1) Thalamus_Z(id,1) ...
        [Caud_Z(id,2) Putamen_Z(id,2) Pallidus_Z(id,2) NAcc_Z(id,2) SN_Z(id,2) VTA_Z(id,2) Thalamus_Z(id,2)]]);
    withinSubjSD_Z.all(id,1) = nanstd([Caud_Z(id,1) Putamen_Z(id,1) Pallidus_Z(id,1) NAcc_Z(id,1) SN_Z(id,1) VTA_Z(id,1) Thalamus_Z(id,1) ...
        [Caud_Z(id,2) Putamen_Z(id,2) Pallidus_Z(id,2) NAcc_Z(id,2) SN_Z(id,2) VTA_Z(id,2) Thalamus_Z(id,2)]]);
    
end


% across subjects

% variances
acrossSubjVar_BP.highDA.allROI   = nanvar(Subj_BP.highDA(:));
acrossSubjVar_BP.highDA.caudate  = nanvar(Caud_BP(:,1));
acrossSubjVar_BP.highDA.putamen  = nanvar(Putamen_BP(:,1));
acrossSubjVar_BP.highDA.pallidus = nanvar(Pallidus_BP(:,1));
acrossSubjVar_BP.highDA.NAcc     = nanvar(NAcc_BP(:,1));
acrossSubjVar_BP.highDA.SN       = nanvar(SN_BP(:,1));
acrossSubjVar_BP.highDA.VTA      = nanvar(VTA_BP(:,1));
acrossSubjVar_BP.highDA.thalamus = nanvar(Thalamus_BP(:,1));

acrossSubjVar_BP.lowDA.allROI   = nanvar(Subj_BP.lowDA(:));
acrossSubjVar_BP.lowDA.caudate  = nanvar(Caud_BP(:,2));
acrossSubjVar_BP.lowDA.putamen  = nanvar(Putamen_BP(:,2));
acrossSubjVar_BP.lowDA.pallidus = nanvar(Pallidus_BP(:,2));
acrossSubjVar_BP.lowDA.NAcc     = nanvar(NAcc_BP(:,2));
acrossSubjVar_BP.lowDA.SN       = nanvar(SN_BP(:,2));
acrossSubjVar_BP.lowDA.VTA      = nanvar(VTA_BP(:,2));
acrossSubjVar_BP.lowDA.thalamus = nanvar(Thalamus_BP(:,2));

acrossSubjVar_BP.All.allROI   = nanvar([Subj_BP.highDA(:); Subj_BP.lowDA(:)]);
acrossSubjVar_BP.All.caudate  = nanvar([Caud_BP(:,1); Caud_BP(:,2)]);
acrossSubjVar_BP.All.putamen  = nanvar([Putamen_BP(:,1);Putamen_BP(:,2)]);
acrossSubjVar_BP.All.pallidus = nanvar([Pallidus_BP(:,1);Pallidus_BP(:,2)]);
acrossSubjVar_BP.All.NAcc     = nanvar([NAcc_BP(:,1);NAcc_BP(:,2)]);
acrossSubjVar_BP.All.SN       = nanvar([SN_BP(:,1); SN_BP(:,2)]);
acrossSubjVar_BP.All.VTA      = nanvar([VTA_BP(:,1); VTA_BP(:,2)]);
acrossSubjVar_BP.All.thalamus = nanvar([Thalamus_BP(:,1);Thalamus_BP(:,2)]);


acrossSubjVar_Z.highDA.allROI   = nanvar(Subj_Z.highDA(:));
acrossSubjVar_Z.highDA.caudate  = nanvar(Caud_Z(:,1));
acrossSubjVar_Z.highDA.putamen  = nanvar(Putamen_Z(:,1));
acrossSubjVar_Z.highDA.pallidus = nanvar(Pallidus_Z(:,1));
acrossSubjVar_Z.highDA.NAcc     = nanvar(NAcc_Z(:,1));
acrossSubjVar_Z.highDA.SN       = nanvar(SN_Z(:,1));
acrossSubjVar_Z.highDA.VTA      = nanvar(VTA_Z(:,1));
acrossSubjVar_Z.highDA.thalamus = nanvar(Thalamus_Z(:,1));

acrossSubjVar_Z.lowDA.allROI   = nanvar(Subj_Z.lowDA(:));
acrossSubjVar_Z.lowDA.caudate  = nanvar(Caud_Z(:,2));
acrossSubjVar_Z.lowDA.putamen  = nanvar(Putamen_Z(:,2));
acrossSubjVar_Z.lowDA.pallidus = nanvar(Pallidus_Z(:,2));
acrossSubjVar_Z.lowDA.NAcc     = nanvar(NAcc_Z(:,2));
acrossSubjVar_Z.lowDA.SN       = nanvar(SN_Z(:,2));
acrossSubjVar_Z.lowDA.VTA      = nanvar(VTA_Z(:,2));
acrossSubjVar_Z.lowDA.thalamus = nanvar(Thalamus_Z(:,2));

acrossSubjVar_Z.All.allROI   = nanvar([Subj_Z.highDA(:); Subj_Z.lowDA(:)]);
acrossSubjVar_Z.All.caudate  = nanvar([Caud_Z(:,1); Caud_Z(:,2)]);
acrossSubjVar_Z.All.putamen  = nanvar([Putamen_Z(:,1);Putamen_Z(:,2)]);
acrossSubjVar_Z.All.pallidus = nanvar([Pallidus_Z(:,1);Pallidus_Z(:,2)]);
acrossSubjVar_Z.All.NAcc     = nanvar([NAcc_Z(:,1);NAcc_Z(:,2)]);
acrossSubjVar_Z.All.SN       = nanvar([SN_Z(:,1); SN_Z(:,2)]);
acrossSubjVar_Z.All.VTA      = nanvar([VTA_Z(:,1); VTA_Z(:,2)]);
acrossSubjVar_Z.All.thalamus = nanvar([Thalamus_Z(:,1);Thalamus_Z(:,2)]);


% standard deviations
acrossSubjSD_BP.highDA.allROI   = nanstd(Subj_BP.highDA(:));
acrossSubjSD_BP.highDA.caudate  = nanstd(Caud_BP(:,1));
acrossSubjSD_BP.highDA.putamen  = nanstd(Putamen_BP(:,1));
acrossSubjSD_BP.highDA.pallidus = nanstd(Pallidus_BP(:,1));
acrossSubjSD_BP.highDA.NAcc     = nanstd(NAcc_BP(:,1));
acrossSubjSD_BP.highDA.SN       = nanstd(SN_BP(:,1));
acrossSubjSD_BP.highDA.VTA      = nanstd(VTA_BP(:,1));
acrossSubjSD_BP.highDA.thalamus = nanstd(Thalamus_BP(:,1));

acrossSubjSD_BP.lowDA.allROI   = nanstd(Subj_BP.lowDA(:));
acrossSubjSD_BP.lowDA.caudate  = nanstd(Caud_BP(:,2));
acrossSubjSD_BP.lowDA.putamen  = nanstd(Putamen_BP(:,2));
acrossSubjSD_BP.lowDA.pallidus = nanstd(Pallidus_BP(:,2));
acrossSubjSD_BP.lowDA.NAcc     = nanstd(NAcc_BP(:,2));
acrossSubjSD_BP.lowDA.SN       = nanstd(SN_BP(:,2));
acrossSubjSD_BP.lowDA.VTA      = nanstd(VTA_BP(:,2));
acrossSubjSD_BP.lowDA.thalamus = nanstd(Thalamus_BP(:,2));

acrossSubjSD_BP.All.allROI   = nanstd([Subj_BP.highDA(:); Subj_BP.lowDA(:)]);
acrossSubjSD_BP.All.caudate  = nanstd([Caud_BP(:,1); Caud_BP(:,2)]);
acrossSubjSD_BP.All.putamen  = nanstd([Putamen_BP(:,1);Putamen_BP(:,2)]);
acrossSubjSD_BP.All.pallidus = nanstd([Pallidus_BP(:,1);Pallidus_BP(:,2)]);
acrossSubjSD_BP.All.NAcc     = nanstd([NAcc_BP(:,1);NAcc_BP(:,2)]);
acrossSubjSD_BP.All.SN       = nanstd([SN_BP(:,1); SN_BP(:,2)]);
acrossSubjSD_BP.All.VTA      = nanstd([VTA_BP(:,1); VTA_BP(:,2)]);
acrossSubjSD_BP.All.thalamus = nanstd([Thalamus_BP(:,1);Thalamus_BP(:,2)]);


acrossSubjSD_Z.highDA.allROI   = nanstd(Subj_Z.highDA(:));
acrossSubjSD_Z.highDA.caudate  = nanstd(Caud_Z(:,1));
acrossSubjSD_Z.highDA.putamen  = nanstd(Putamen_Z(:,1));
acrossSubjSD_Z.highDA.pallidus = nanstd(Pallidus_Z(:,1));
acrossSubjSD_Z.highDA.NAcc     = nanstd(NAcc_Z(:,1));
acrossSubjSD_Z.highDA.SN       = nanstd(SN_Z(:,1));
acrossSubjSD_Z.highDA.VTA      = nanstd(VTA_Z(:,1));
acrossSubjSD_Z.highDA.thalamus = nanstd(Thalamus_Z(:,1));

acrossSubjSD_Z.lowDA.allROI   = nanstd(Subj_Z.lowDA(:));
acrossSubjSD_Z.lowDA.caudate  = nanstd(Caud_Z(:,2));
acrossSubjSD_Z.lowDA.putamen  = nanstd(Putamen_Z(:,2));
acrossSubjSD_Z.lowDA.pallidus = nanstd(Pallidus_Z(:,2));
acrossSubjSD_Z.lowDA.NAcc     = nanstd(NAcc_Z(:,2));
acrossSubjSD_Z.lowDA.SN       = nanstd(SN_Z(:,2));
acrossSubjSD_Z.lowDA.VTA      = nanstd(VTA_Z(:,2));
acrossSubjSD_Z.lowDA.thalamus = nanstd(Thalamus_Z(:,2));

acrossSubjSD_Z.All.allROI   = nanstd([Subj_Z.highDA(:); Subj_Z.lowDA(:)]);
acrossSubjSD_Z.All.caudate  = nanstd([Caud_Z(:,1); Caud_Z(:,2)]);
acrossSubjSD_Z.All.putamen  = nanstd([Putamen_Z(:,1);Putamen_Z(:,2)]);
acrossSubjSD_Z.All.pallidus = nanstd([Pallidus_Z(:,1);Pallidus_Z(:,2)]);
acrossSubjSD_Z.All.NAcc     = nanstd([NAcc_Z(:,1);NAcc_Z(:,2)]);
acrossSubjSD_Z.All.SN       = nanstd([SN_Z(:,1); SN_Z(:,2)]);
acrossSubjSD_Z.All.VTA      = nanstd([VTA_Z(:,1); VTA_Z(:,2)]);
acrossSubjSD_Z.All.thalamus = nanstd([Thalamus_Z(:,1);Thalamus_Z(:,2)]);


%% gather z-scores of each ROI per subject
%  ROIs: Amygdala, LC, SN, VTA, caudate/putamen, Nacc


%% Collect BP_srtm change in each ROIs

BPchange_SRTM_save=[];

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            clear SRTMimg baseimage FSlabel existingLabels
            load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')
            
            baseimage = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii']));
            SRTMimg = zeros(size(baseimage));
            existingLabels = unique(baseimage);
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    
                    fprintf(['**' FSlabels1{labels,1}{1} '**\n'])
                    
                    % left
                    disp('left')
                    clear indL
                    indL=find(baseimage==FSlabels1{labels,2});
                    SRTMimg(indL)=  (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl-BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)  / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);
                   

                    % right
                    disp('right')
                    clear indR
                    indR=find(baseimage==FSlabels1{labels,3});
                    SRTMimg(indR)= (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl - BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm) / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);
                    
                    BPchange_SRTM_save{id,d}.(FSlabels1{labels,1})=(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl - BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm) / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);
%                     disp(['both, ' num2str(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)])
                    
                end
            end
            
            
            %SRTM
            hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii']); % pick just any header from a file
            hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_SRTM_BPchange.nii'];
            hdr.dim = size(SRTMimg);
            hdr = rmfield(hdr,'pinfo');
            hdr.nii = spm_write_vol(hdr,SRTMimg);
            
            fprintf(['\n\ndone\n\n'])
        end
    end
end


%%
BPchange_SRTM_save=[];
load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0

        else

            clear SRTMimg baseimage FSlabel existingLabels

            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    BPchange_SRTM_save{id,d}.(FSlabels1{labels,1}) ...
                        = (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl - BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm) / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);


                end
            end
        end

    end
end


disp('done')

%% make bpchange map

cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            eval(['BPchange_SRTM_' FSlabels1{labels,1}{1} '_HiLo(id,d)=NaN;'])
        else
            load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')
            cnt=cnt+1;
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                                        
                    eval(['BPchange_SRTM_' FSlabels1{labels,1}{1} '(cnt,1)=BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                    eval(['BPchange_SRTM_HiLo_' FSlabels1{labels,1}{1} '(id,d)=BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                  
                end
            end
        end
    end
end



for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

        eval(['meanBPchange_SRTM_' FSlabels1{labels,1}{1} '=nanmean(BPchange_SRTM_' FSlabels1{labels,1}{1} ');'])

    end
end



clear SRTMimg baseimage FSlabel existingLabels
% load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/FS_labels_LCSNVTA.mat')

baseimage = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']));
SRTMimg = zeros(size(baseimage));
existingLabels = unique(baseimage);

for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

        fprintf(['**' FSlabels1{labels,1}{1} '**\n'])

        % left
%         disp('left')
        clear indL
        indL=find(baseimage==FSlabels1{labels,2});
        SRTMimg(indL)=eval(strcat('meanBPchange_SRTM_',FSlabels1{labels,1}{1}));

        % right
%         disp('right')
        clear indR
        indR=find(baseimage==FSlabels1{labels,3});
        SRTMimg(indR)=eval(strcat('meanBPchange_SRTM_',FSlabels1{labels,1}{1}));

%         disp(['both, ' num2str(BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm)])

    end
end


%SRTM
hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/aparc+aseg_on_template_labelled.nii']); % pick just any header from a file
hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/meanBPchange_SRTM_onTemplate.nii'];
hdr.dim = size(SRTMimg);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,SRTMimg);





%%

hd=0; ld=0; ad=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    eval( strcat('BPchange_all_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                    eval( strcat('BP_SRTM_BSL_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                    eval( strcat('BP_SRTM_TASK_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                else
                    eval( strcat('BPchange_all_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                    eval( strcat('BP_SRTM_BSL_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                    eval( strcat('BP_SRTM_TASK_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
%                     if d==1
%                         hd=hd+1;
%                      eval(  strcat('BPchange_highDA_', (FSlabels1{labels,1}), '(id,1)=BPchange_SRTM_save{id,d}.', (FSlabels1{labels,1}) ) );
%                     else
%                         ld=ld+1;
%                         eval(  strcat('BPchange_lowDA_', (FSlabels1{labels,1}), '(id,1)=BPchange_SRTM_save{id,d}.', (FSlabels1{labels,1}) ) );
%                     end
                end

            end
        else
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    eval( strcat('BPchange_all_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                    eval( strcat('BP_SRTM_BSL_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                    eval( strcat('BP_SRTM_TASK_', (FSlabels1{labels,1}), '(id,d)=NaN;' ) );
                else
                    eval( strcat('BPchange_all_', (FSlabels1{labels,1}), '(id,d)=BPchange_SRTM_save{id,d}.', (FSlabels1{labels,1}),';' ) );
                    eval( strcat('BP_SRTM_BSL_', (FSlabels1{labels,1}), '(id,d)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;' ) );
                    eval( strcat('BP_SRTM_TASK_', (FSlabels1{labels,1}), '(id,d)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;' ) );
%                     if d==1
%                         hd=hd+1;
%                      eval(  strcat('BPchange_highDA_', (FSlabels1{labels,1}), '(id,1)=BPchange_SRTM_save{id,d}.', (FSlabels1{labels,1}) ) );
%                     else
%                         ld=ld+1;
%                         eval(  strcat('BPchange_lowDA_', (FSlabels1{labels,1}), '(id,1)=BPchange_SRTM_save{id,d}.', (FSlabels1{labels,1}) ) );
%                     end
                end

            end
        end
    end
end
disp('done')

% make variables
for labels=2:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else

        eval( strcat('BP_SRTM_BSL_highDA_', FSlabels1{labels,1}, '=BP_SRTM_BSL_', FSlabels1{labels,1}, '(:,1);') )
        eval( strcat('BP_SRTM_BSL_lowDA_', FSlabels1{labels,1}, '=BP_SRTM_BSL_', FSlabels1{labels,1}, '(:,2);') )
        eval( strcat('BP_SRTM_TASK_highDA_', FSlabels1{labels,1}, '=BP_SRTM_TASK_', FSlabels1{labels,1}, '(:,1);') )
        eval( strcat('BP_SRTM_TASK_lowDA_', FSlabels1{labels,1}, '=BP_SRTM_TASK_', FSlabels1{labels,1}, '(:,2);') )

        eval( strcat('BP_SRTM_change_highDA_', FSlabels1{labels,1}, '=BPchange_all_', FSlabels1{labels,1}, '(:,1);') )
        eval( strcat('BP_SRTM_change_lowDA_', FSlabels1{labels,1}, '=BPchange_all_', FSlabels1{labels,1}, '(:,2);') )


    end
end

tablestr_high_bsl='BP_SRTM_BSL_highDA_bankssts';
tablestr_low_bsl='BP_SRTM_BSL_lowDA_bankssts';

for labels=3:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else
        tablestr_high_bsl=strcat(tablestr_high_bsl,',BP_SRTM_BSL_highDA_',FSlabels1{labels,1});
        tablestr_low_bsl=strcat(tablestr_low_bsl,',BP_SRTM_BSL_lowDA_',FSlabels1{labels,1});

    end
end

tablestr_high_TASK='BP_SRTM_TASK_highDA_bankssts';
tablestr_low_TASK='BP_SRTM_TASK_lowDA_bankssts';
for labels=3:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else
        tablestr_high_TASK=strcat(tablestr_high_TASK,',BP_SRTM_TASK_highDA_',FSlabels1{labels,1});
        tablestr_low_TASK=strcat(tablestr_low_TASK,',BP_SRTM_TASK_lowDA_',FSlabels1{labels,1});
    end
end

tablestr_high_change='BP_SRTM_change_highDA_bankssts';
tablestr_low_change='BP_SRTM_change_lowDA_bankssts';

for labels=3:length(FSlabels1)
    if FSlabels1{labels,2}==0
    else
        tablestr_high_change=strcat(tablestr_high_change,',BP_SRTM_change_highDA_',FSlabels1{labels,1});
        tablestr_low_change=strcat(tablestr_low_change,',BP_SRTM_change_lowDA_',FSlabels1{labels,1});

    end
end


% ID=IDs';
% eval(strcat('MRPET_BP_SRTM_Table = table(ID,', tablestr_high_bsl, ',',  tablestr_high_TASK, ',',tablestr_low_bsl, ',', tablestr_low_TASK, ')'))
% writetable(MRPET_BP_SRTM_Table,'MRPET_BPtable_SRTM.xls')
% eval(strcat('MRPET_BPchange_SRTM_Table = table(ID,', tablestr_high_change,',', tablestr_low_change, ')'))
% writetable(MRPET_BPchange_SRTM_Table,'MRPET_BPchange_SRTM.xls')




%% scatter plots

% load data from here: /Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi
        
    close all
    figure,scatter(LC_CRmean_keren(:,2),BPchange_SRTM_HiLo_LC(:,1),120,'MarkerFaceColor',[1,0.2,0.2],'MarkerEdgeColor',[1,0.2,0.2],...
        'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
    x = -0.25 : 0.05;
    m = -0.08;
    b = 0.44;
    y = b*x + m;
    plot(x, y,'LineWidth',4)%,'Color',[1,0.2,0.2],'LineStyle','--')
    set(gca,'FontSize',25)
    ylabel('LC Contrast ratio','FontWeight','bold','FontSize',30)
    xlabel('BP change in LC during high DA session','FontWeight','bold','FontSize',30)

    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    Slope = B(2)
    Intercept = B(1)
    ylim([min(LC_CRmean_keren(:,2))-0.01 max(LC_CRmean_keren(:,2))+0.1])
    %     xlim([min(BPchange_SRTM_HiLo_LC(:,2))-0.002 max(BPchange_SRTM_HiLo_LC(:,2))+0.001])


%%
    close all
    figure,scatter(LC_CRmax_keren(:,2),BPchange_SRTM_HiLo_Amy(:,2),120,'MarkerFaceColor',[1,0.2,0.2],'MarkerEdgeColor',[1,0.2,0.2],...
        'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
    x = -0.25 : 0.05;
    m = -0.08;
    b = 0.44;
    y = b*x + m;
    plot(x, y,'LineWidth',4)%,'Color',[1,0.2,0.2],'LineStyle','--')
    set(gca,'FontSize',25)
    ylabel('LC Contrast ratio','FontWeight','bold','FontSize',30)
    xlabel('BP change in Amygdala during high DA session','FontWeight','bold','FontSize',30)

    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    Slope = B(2)
    Intercept = B(1)
%     ylim([min(LC_CRmean_keren(:,2))-0.01 max(LC_CRmean_keren(:,2))+0.1])
        xlim([0 2])
        ylim([-0.1 0.1])



%%

close all
    figure,scatter(Dprime_rew_highDA_delayed(:,1),BPchange_SRTM_HiLo_VTA(:,1),120,'MarkerFaceColor',[1,0.2,0.2],'MarkerEdgeColor',[1,0.2,0.2],...
        'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
    plot(x, y,'LineWidth',4)%,'Color',[1,0.2,0.2],'LineStyle','--')
    set(gca,'FontSize',25)
    ylabel('Delayed recall D'' in reward trials','FontWeight','bold','FontSize',30)
    xlabel('BP change in VTA during high DA session','FontWeight','bold','FontSize',30)

    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    Slope = B(2);
    Intercept = B(1);
%     ylim([min(Dprime_rew_highDA_delayed(:,1))-0.01 max(Dprime_rew_highDA_delayed(:,1))+0.1])

    %%




scatter(mean_mni(:,2),mean_nat(:,2),120,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);
x = -0.2 : 1;
m = 0.08;
b = 0.44;
y = m*x + b;
plot(x, y,'LineWidth',4,'Color',[0.2,0.2,1],'LineStyle','--')
set(gca,'FontSize',25); hold on
 
x = 3 : 6;
m = 0.12;
b = 4.02;
y = m*x + b;
plot(x, y,'LineWidth',6,'Color',[0 0 0])
set(gca,'FontSize',25)
yticks([0 1 2  3 4 5 6])
xticks([3 4 5 6])
 
legend({'Left LC', 'Left LC fit line', 'Right LC', 'Right LC fit line','Both LC fit line'});




