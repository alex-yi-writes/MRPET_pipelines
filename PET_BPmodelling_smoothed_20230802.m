%% modelling

clear; clc; close all
warning('off','all');

% set paths
paths = [];
paths.parent  = '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/original/';
paths.TACs    = '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/Smoothed/';
paths.TACs_new= '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/Smoothed/';
paths.ROImask = '/Volumes/ALEX3/MRPET/coreg_roi/';
paths.figure  = '/Users/alex/Dropbox/paperwriting/MRPET/figures/modellingFigures/';
paths.rawImg  = '/Volumes/ALEX3/MRPET/img/';
paths.exports = '/Users/alex/Dropbox/paperwriting/MRPET/data/';

addpath(genpath('/Users/alex/Dropbox/paperwriting/MRPET/scripts/modelling'))

% set envs
% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];


cnt=0;
for i1=1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
        else
            cnt=cnt+1;
            ID{cnt}=[num2str(IDs(i1)) num2str(d)];
        end
    end
end

% load TACs
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            TACs{i1,d} = [];
        else
            TACs{i1,d} = load([paths.TACs_new num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr_s3_all.mat']);
        end
    end
end

load('/Users/alex/Dropbox/paperwriting/MRPET/data/FS_labels_LCSNVTA.mat')
load('/Users/alex/Dropbox/paperwriting/MRPET/data/ROIs_LCSNVTA.mat')
ROI=ROI';

disp('prep done')

%% merge-hippocampus, amygdala, SN, VTA, LC

for collapseLoops=1
% for id = 1:length(IDs)
%     
%     for d = 1:2
%         
%         if days(id,d) == 0
%             warning('Skipped')
%         else
%             
%             % pt 1
% %             clear imgtmp indL indR names hdr
% %             names=['ROIs_HPC-AMG-SNVTA-LC.nii'];
% %             
% %             clear imgbase
% %             imgbase=zeros(size(spm_read_vols(spm_vol([ paths.ROImask num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii' ]))));
% %             
% %             imgtmp=spm_read_vols(spm_vol([ paths.ROImask num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii' ]));
% %             indL=( imgtmp==17 | imgtmp==18 | imgtmp==991 | imgtmp==993 | imgtmp==995 );
% %             imgbase(indL)=997;
% %             indR=( imgtmp==53 | imgtmp==54 | imgtmp==992 | imgtmp==994 | imgtmp==996 );
% %             imgbase(indR)=998;
% %             
% %             hdr = spm_vol([ paths.ROImask num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled.nii' ]); % pick just any header from a file
% %             hdr.fname = [paths.ROImask num2str(IDs(id)) num2str(d) '/ROIs_HPC-AMG-SNVTA-LC_pt1.nii'];
% %             hdr.dim = size(imgbase);
% %             hdr = rmfield(hdr,'pinfo');
% %             hdr.nii = spm_write_vol(hdr,imgbase);
%             
%             % pt 2
% %             clear imgtmp indL indR names hdr
% %             names=['ROIs_HPC-AMG-SNVTA-LC.nii'];
% %             
% %             clear imgbase
% %             imgbase=zeros(size(spm_read_vols(spm_vol([ paths.ROImask num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii' ]))));
% %             
% %             imgtmp=spm_read_vols(spm_vol([ paths.ROImask num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii' ]));
% %             indL=( imgtmp==17 | imgtmp==18 | imgtmp==991 | imgtmp==993 | imgtmp==995 );
% %             imgbase(indL)=997;
% %             indR=( imgtmp==53 | imgtmp==54 | imgtmp==992 | imgtmp==994 | imgtmp==996 );
% %             imgbase(indR)=998;
% %             
% %             hdr = spm_vol([ paths.ROImask num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii' ]); % pick just any header from a file
% %             hdr.fname = [paths.ROImask num2str(IDs(id)) num2str(d) '/ROIs_HPC-AMG-SNVTA-LC_pt2.nii'];
% %             hdr.dim = size(imgbase);
% %             hdr = rmfield(hdr,'pinfo');
% %             hdr.nii = spm_write_vol(hdr,imgbase);
%             
%             
%         end
%     end
% end
end

%% calculate TACs for merged mask

% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
TrimPerc=15;

% INFLOW
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask1   =[ paths.ROImask num2str(IDs(id)) num2str(d) '/ROIs_HPC-AMG-SNVTA-LC_pt1.nii'];
            
            % 4-D PET file
            PETflow = [paths.rawImg num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_s3_InFlow' num2str(d) '_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask1);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data

            maskidx.SubcorticalROIs.LongName=strtrim('HPC-AMG-SNVTA-LC');
            for Hemi={'Left','Right','Bilateral'}
                switch Hemi{1}
                    case 'Left'
                        CurIdx=997;
                    case 'Right'
                        CurIdx=998;
                    case 'Bilateral'
                        CurIdx=[997, 997];
                end
                IndList=[];
                for ROIInd=CurIdx
                    if ~isnan(ROIInd)
                        IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                    end
                end
                maskidx.SubcorticalROIs.(Hemi{1})=IndList;
            end
            
            %% Read in 4D-PET data and extract ROI-averages in each frame
            
            clear TACDATA
            
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
            mkdir([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' ])
            save([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_tmp.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end

% BSL & TASK
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask2   =[ paths.ROImask num2str(IDs(id)) num2str(days(id,d)) '/ROIs_HPC-AMG-SNVTA-LC_pt1.nii']; % task and the baseline
            
            % 4-D PET file
            PETtask = [paths.rawImg num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_s3_MT' num2str(d) '_on_T1.nii'];
            PETbsl  = [paths.rawImg num2str(IDs(id)) '_' num2str(days(id,d)) '/coreg_s3_Baseline' num2str(d) '_on_T1.nii'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask2);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            maskidx=[];
            
            ROIDef = {'997998', 'SubcorticalROIs', 'HPC-AMG-SNVTA-LC'};
            maskidx.SubcorticalROIs.LongName=strtrim('HPC-AMG-SNVTA-LC');
            
            for Hemi={'Left','Right','Bilateral'}
                switch Hemi{1}
                    case 'Left'
                        CurIdx=997;
                    case 'Right'
                        CurIdx=998;
                    case 'Bilateral'
                        CurIdx=[997, 997];
                end
                IndList=[];
                for ROIInd=CurIdx
                    if ~isnan(ROIInd)
                        IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                    end
                end
                maskidx.SubcorticalROIs.(Hemi{1})=IndList;
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
            save([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_tmp.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            
            
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
            save([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_tmp.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
        end
    end
end

%% retro-correct the decay

t0_frame_bsl    = 95;
t0_frame_task   = 115;

for CollapseLoops=1
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_tmp.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_tmp.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_tmp.mat']);
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
            
            save([paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_tmp.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
        end
    end
end
end

%% merge this info to the other TACs

for CollapseLoops=1
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            warning('Skipped')
        else
            
            clear TAC_old TAC_new TACDATA_Baseline TACDATA_InFlow TACDATA_Task
            
            TAC_new = load([paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_tmp.mat']);
            TAC_old = load([paths.TACs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_s3.mat']);
            
            TACDATA_Baseline = TAC_old.TACDATA_Baseline;
            TACDATA_Baseline.SubcorticalROIs=TAC_new.TACDATA_Baseline.SubcorticalROIs;
            TACDATA_Task = TAC_old.TACDATA_Task;
            TACDATA_Task.SubcorticalROIs=TAC_new.TACDATA_Task.SubcorticalROIs;
            TACDATA_InFlow = TAC_old.TACDATA_InFlow;
            TACDATA_InFlow.SubcorticalROIs=TAC_new.TACDATA_InFlow.SubcorticalROIs;
            
            save([paths.TACs_new num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_s3_all.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            
        end
    end
end

end

%% run the model: condensed ROIs

set(0, 'DefaultLineLineWidth', 1);

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg={'CerC','Put','Caud','Nac','Hipp','ThalProp','SubcorticalROIs','Amy','LC','SN','VTA'}
                if IDs(id)==4006 && d==2 && strcmp(reg{1},'SubcorticalROIs')
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
                else
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
                end
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
            break_point=find(times(:,1)>=115*60,1,'first'); %% Time of activation start
            
            %%%%%%%%%%%
            PlotStrFit=0;
            %%%%%%%%%%%
            
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg={'Striatum','Caudate','Putamen','Accumbensarea','Thalamus','Hippocampus','Amygdala','SubcorticalROIs'}
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
                            ROItac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                        case 'Caudate'
                            ROItac=TACDATA.Caud;
                            savestr = 'Caudate';
                        case 'Putamen'
                            ROItac=TACDATA.Put;
                            savestr = 'Putamen';
                        case 'Accumbensarea'
                            ROItac=TACDATA.Nac;
                            savestr = 'Accumbensarea';
                        case 'Thalamus'
                            ROItac=TACDATA.ThalProp;
                            savestr = 'Thalamus';
                        case 'Hippocampus'
                            ROItac=TACDATA.Hipp;
                            savestr = 'Hippocampus';
                        case 'SubcorticalROIs'
                            ROItac=TACDATA.SubcorticalROIs;
                            savestr = 'SubcorticalROIs';
                        case 'Amygdala'
                            ROItac=TACDATA.Amy;
                            savestr = 'Amygdala';
                    end
                    
                    mROItac  = [ROItac(1)/2; (ROItac(2:end)+ROItac(1:end-1))/2];
                    ASRTM(:,3)=zeros(t_points,1);
%                     ASTRM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mROItac(1:k).*dt(1:k));
                    end
                    %LSQ-estimation using lscov
                    [parest se_srtm mse_srtm]   = lscov(ASRTM,ROItac);
                    fittac=ASRTM*parest;
                    BP=parest(2)/parest(3)-1;
                    k2p=parest(2)/parest(1);
                    
                    %%%%% DO real SRTM
                    options = optimset('MaxFunEvals',1000);
                    weighs=[0.25*ones(30,1); ones(t_points-30,1)];
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-ROItac).*weighs);
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    fobj = @(x) norm((simESRTM_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-ROItac(1:end-11)).*weighs(1:end-11));
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1_bl=parest_srtm(1);
                    k2_bl=parest_srtm(2);
                    BP_bl=parest_srtm(3);
                    modfit_esrtm_bl=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
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
                                    roitac_gamma = ROItac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(ROItac))+roitac_gamma(1:length(ROItac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,ROItac);
                                    
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
                            end
        %                        breakpoint2=breakpoint;
                        end
                    end
                    
                    y=-ASRTM(:,3)./ROItac;
                    x=ASRTM(:,2)./ROItac;
                    [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                    [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                    
                    % save these and work out receptor occupancy later in
                    % the script

                    k2=best_parest(2);
                    k2a=best_parest(3);
                    BP_lp=k2/k2a-1; % baseline binding potential of lpntPET
                    BP_srtm=BP__; % binding potential until the end of the task
                    BP_srtm_bsl=BP_bl; % binding potential until the end of baseline
                    G=best_parest(4); % gamma, obviously...
                    BPND = ((R1__*k2p)/k2a)-1; % correct? - not sure
%                     DBP =  k2./(k2a + best_actfun) - 1; % dynamic binding potential
                    DBP =  k2./(k2a + G*best_actfun) - 1; % is it this way...? confused
                    OCC = 100*(1-DBP/BP_lp);% receptor occupancy
%                     OCC2 = (1-DBP2/BP_lp);% receptor occupancy
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.' savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.' savestr '=OCC;' ])
%                     eval(['Occupancy2{id,d}.' savestr '=OCC2;' ])
                    eval(['BPND_save{id,d}.' savestr '=BPND;' ])
                    eval(['BP_srtm_save{id,d}.' savestr '=BP_srtm;' ])
                    eval(['BP_srtm_Bsl_save{id,d}.' savestr '=BP_srtm_bsl;' ])
                    eval(['BestActivationFunction{id,d}.' savestr '=best_actfun;'])
                    eval(['BestModelFit{id,d}.' savestr '=best_modelfit;'])
                    eval(['gamma_lp_save{id,d}.' savestr '=best_parest(4);' ])
                    eval(['Residuals{id,d}.' savestr '.BaselineFit=ROItac-fittac;' ])
                    eval(['Residuals{id,d}.' savestr '.AllFit_ESRTM=ROItac-modfit_esrtm;' ])
                    eval(['Residuals{id,d}.' savestr '.lpntPETFit=ROItac-best_modelfit;' ])
                    eval(['CompensatoryFunction{id,d}.' savestr '=best_parest(4)*best_actfun;' ])
                    
                    
                    % Compute residuals
%                     residuals_lpnt = ROItac-best_modelfit;
%                     
%                     % Apply weights
%                     Weight = [ones(1,100) ones(1,11).*2];
%                     weights = 1./sqrt(Weight);
%                     weighted_resid_lpnt = residuals_lpnt .* weights;
%                     
%                     % Square and sum weighted residuals to obtain WRSS
%                     WRSS_lpnt = sum(weighted_resid_lpnt.^2);
                    

%                     BP_lp=best_parest(2)/best_parest(3)-1; % is this the BP from lpntPET?
%                     BP_lpntPET=(best_parest(2)/(best_parest(3)+(best_parest(4)).*best_actfun))-1; % isn't this the way?
%                     BP_bsl=BP_bl; % baseline from SRTM
                    % save BP for later analysis
%                     eval(['BP_lp_orig{id,d}.' savestr '=BP_lp;' ])
%                     eval(['BP_lp_save{id,d}.' savestr '=BP_lpntPET;' ])
%                     eval(['BP_bsl_save{id,d}.' savestr '=BP_bsl;' ])
                    
                    
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
                    BPdataSave{id,d}.(reg{1}).BPND=BPND;
                    BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
%                     BPdataSave{id,d}.(reg{1}).Occupancy2=OCC2;
                    BPdataSave{id,d}.(reg{1}).DBP=DBP;
                    BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);


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
                        plot(tmidMin,ROItac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
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
                        plot(tmidMin,ROItac-fittac,'bo',tmidMin,ROItac-modfit_esrtm,'go',tmidMin,ROItac-modfit_esrtm_bl,'co',tmidMin,ROItac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        plot(tmidMin,ROItac-fittac,'bo',tmidMin,ROItac-modfit_esrtm,'bo',tmidMin,ROItac-modfit_esrtm_bl,'ro',tmidMin,ROItac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        [h p]=runstest(ROItac-fittac);
                        [h1 p1]=runstest(ROItac-modfit_esrtm);
                        [h2 p2]=runstest(ROItac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
                        title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
                        xlabel('Time (min)');

                        subplot(3,2,4);
                        y=-ASRTM(:,3)./ROItac;
                        x=ASRTM(:,2)./ROItac;
                        [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                        [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                        plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
                        title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        xlabel(['\int REF/ROI']);
                        ylabel('\int ROI/ROI')
                        subplot(3,2,[5 6]);
                        plot(tmidMin,ROItac./reftac,'ko',[0 180],[0 0],'k--');
                        ylim([0.5 30]); xlim([0 190])
                        title('Target to reference ratio');
                        xlabel('Time (min)');

                        print('-dpsc2','-append','-bestfit',fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_TAC_Fit_lpntpet_logan_' date '.ps']));
                        close(gcf)
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
                        BPdataSave{id,d}.(reg{1}).BPND=BPND;
                        BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
                        BPdataSave{id,d}.(reg{1}).Occupancy2=OCC2;
                        BPdataSave{id,d}.(reg{1}).DBP=DBP;
                        BPdataSave{id,d}.(reg{1}).DBP2=DBP2;
                        BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);
                        
                        
                        continue;
                    end
                end
                
                
            end
            close all
            keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI BPdataSave BP_lp_save DBP_save Occupancy BPND_save BP_srtm_save BP_srtm_Bsl_save BestActivationFunction...
    BestModelFit gamma_lp_save Residuals CompensatoryFunction
        end
    end
end

disp('done')

save(['/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/MRPET_BPpackage_condensed_smoothed3mm_' date '.mat'],...
    'BPdataSave', 'BP_lp_save', 'DBP_save', 'Occupancy', 'BPND_save', 'BP_srtm_save', 'BP_srtm_Bsl_save', 'BestActivationFunction',...
    'BestModelFit', 'gamma_lp_save', 'Residuals', 'CompensatoryFunction')

                    
%% plot target to reference ratio
set(0, 'DefaultAxesFontSize',10);

regcount=0;
for id = 2:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            figure('Position',[100 100 800 1200],'Visible','on'); hold on;
            spidx=1;
            sgtitle(['ID ' num2str(IDs(id)) num2str(d)],'FontSize',24,'FontWeight','bold')
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg={'CerC','Put','Caud','Nac','Hipp','ThalProp','SubcorticalROIs'}
                if IDs(id)==4006 && d==2 && strcmp(reg{1},'SubcorticalROIs')
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac(10:end); ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
                else
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
                end
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
            break_point=find(times(:,1)>=115*60,1,'first'); %% Time of activation start
            
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
            
            regcount=0;
           

            for r=1
                for reg={'Striatum','Caudate','Putamen','Accumbens-area','Thalamus','Hippocampus','SubcorticalROIs'}
                    regcount=regcount+1;
                    reftac=TACDATA.CerC;
                    mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
                    
                    if PlotStrFit
                        subplot(4,2,1)
                        plot(tmidMin,reftac,'k-'); hold on;
                        title(['TAC in cerebellum'],'FontSize',10);
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
                            ROItac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                        case 'Caudate'
                            ROItac=TACDATA.Caud;
                            savestr = 'Caudate';
                        case 'Putamen'
                            ROItac=TACDATA.Put;
                            savestr = 'Putamen';
                        case 'Accumbens-area'
                            ROItac=TACDATA.Nac;
                            savestr = 'Accumbensarea';
                        case 'Thalamus'
                            ROItac=TACDATA.ThalProp;
                            savestr = 'Thalamus';
                        case 'Hippocampus'
                            ROItac=TACDATA.Hipp;
                            savestr = 'Hippocampus';
                        case 'SubcorticalROIs'
                            ROItac=TACDATA.SubcorticalROIs;
                            savestr = 'SubcorticalROIs';
                    end
                    
                    mROItac  = [ROItac(1)/2; (ROItac(2:end)+ROItac(1:end-1))/2];
                    
                    if PlotStrFit
                        subplot(4,2,regcount+1)
                        plot(tmidMin,ROItac./reftac,'ko',[0 180],[0 0],'k--'); hold on
                        ylim([0.5 30]); xlim([0 190])
                        title(['Target to reference ratio: ' reg{1}],'FontSize',10);
                        xlabel('Time (min)','FontSize',10);
           
                        continue;
                    end
                    
                end
                
                
            end
            print('-dpsc2','-append','-bestfit',fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_T2R_' date '.ps']));
            close all
%             keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI
        end
    end
end


%% all ROIs

load('/Users/alex/Dropbox/paperwriting/MRPET/data/FS_labels_LCSNVTA.mat')
load('/Users/alex/Dropbox/paperwriting/MRPET/data/ROIs_LCSNVTA.mat')


for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            for reg=ROI'%{'CerC','Put','Caud','Nac','Hipp','ThalProp'}
                if IDs(id)==4006 && d==2 && strcmp(reg{1},'SubcorticalROIs')
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
                else
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
                end
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
            PlotStrFit=0;
            %%%%%%%%%%%
            
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg=ROI'
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
                            ROItac=TACDATA.Caud;
                            savestr = 'Caud';
                        case 'Put'
                            ROItac=TACDATA.Put;
                            savestr = 'Put';
                        case 'Striatum' 
                            tempTac=[];
                            tempVol=[];
                            for subReg={'Caud','Put'}
                                tempTac=[tempTac, TACDATA.(subReg{1})];
                                tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol];
                            end
                            ROItac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                        case 'Nac'
                            ROItac=TACDATA.Nac;
                            savestr = 'Nac';
                        case 'ThalProp'
                            ROItac=TACDATA.ThalProp;
                            savestr = 'ThalProp';
                        case 'Hipp'
                            ROItac=TACDATA.Hipp;
                            savestr = 'Hipp';
                        case 'bankssts'
                            ROItac=TACDATA.bankssts;
                            savestr = 'bankssts';
                        case 'caudalanteriorcingulate'
                            ROItac=TACDATA.caudalanteriorcingulate;
                            savestr = 'caudalanteriorcingulate';
                        case 'caudalmiddlefrontal'
                            ROItac=TACDATA.caudalmiddlefrontal;
                            savestr = 'caudalmiddlefrontal';
                        case 'cuneus'
                            ROItac=TACDATA.cuneus;
                            savestr = 'cuneus';
                        case 'entorhinal'
                            ROItac=TACDATA.entorhinal;
                            savestr = 'entorhinal';
                        case 'fusiform'
                            ROItac=TACDATA.fusiform;
                            savestr = 'fusiform';
                        case 'inferiorparietal'
                            ROItac=TACDATA.inferiorparietal;
                            savestr = 'inferiorparietal';
                        case 'inferiortemporal'
                            ROItac=TACDATA.inferiortemporal;
                            savestr = 'inferiortemporal';
                        case 'isthmuscingulate'
                            ROItac=TACDATA.isthmuscingulate;
                            savestr = 'isthmuscingulate';
                        case 'lateraloccipital'
                            ROItac=TACDATA.lateraloccipital;
                            savestr = 'lateraloccipital';
                        case 'lateralorbitofrontal'
                            ROItac=TACDATA.lateralorbitofrontal;
                            savestr = 'lateralorbitofrontal';
                        case 'lingual'
                            ROItac=TACDATA.lingual;
                            savestr = 'lingual';
                        case 'medialorbitofrontal'
                            ROItac=TACDATA.medialorbitofrontal;
                            savestr = 'medialorbitofrontal';
                        case 'parahippocampal'
                            ROItac=TACDATA.parahippocampal;
                            savestr = 'parahippocampal';
                        case 'paracentral'
                            ROItac=TACDATA.paracentral;
                            savestr = 'paracentral';
                        case 'precuneus'
                            ROItac=TACDATA.precuneus;
                            savestr = 'precuneus';
                        case 'parsopercularis'
                            ROItac=TACDATA.parsopercularis;
                            savestr = 'parsopercularis';
                        case 'parsorbitalis'
                            ROItac=TACDATA.parsorbitalis;
                            savestr = 'parsorbitalis';
                        case 'parstriangularis'
                            ROItac=TACDATA.parstriangularis;
                            savestr = 'parstriangularis';
                        case 'pericalcarine'
                            ROItac=TACDATA.pericalcarine;
                            savestr = 'pericalcarine';
                        case 'postcentral'
                            ROItac=TACDATA.postcentral;
                            savestr = 'postcentral';
                        case 'posteriorcingulate'
                            ROItac=TACDATA.posteriorcingulate;
                            savestr = 'posteriorcingulate';
                        case 'precentral'
                            ROItac=TACDATA.precentral;
                            savestr = 'precentral';
                        case 'rostralanteriorcingulate'
                            ROItac=TACDATA.rostralanteriorcingulate;
                            savestr = 'rostralanteriorcingulate';
                        case 'rostralmiddlefrontal'
                            ROItac=TACDATA.rostralmiddlefrontal;
                            savestr = 'rostralmiddlefrontal';
                        case 'superiorfrontal'
                            ROItac=TACDATA.superiorfrontal;
                            savestr = 'superiorfrontal';
                        case 'superiorparietal'
                            ROItac=TACDATA.superiorparietal;
                            savestr = 'superiorparietal';
                        case 'superiortemporal'
                            ROItac=TACDATA.superiortemporal;
                            savestr = 'superiortemporal';
                        case 'supramarginal'
                            ROItac=TACDATA.supramarginal;
                            savestr = 'supramarginal';
                        case 'frontalpole'
                            ROItac=TACDATA.frontalpole;
                            savestr = 'frontalpole';
                        case 'temporalpole'
                            ROItac=TACDATA.temporalpole;
                            savestr = 'temporalpole';
                        case 'transversetemporal'
                            ROItac=TACDATA.transversetemporal;
                            savestr = 'transversetemporal';
                        case 'middletemporal'
                            ROItac=TACDATA.middletemporal;
                            savestr = 'middletemporal';
                        case 'insula'
                            ROItac=TACDATA.insula;
                            savestr = 'insula';
                        case 'CerC'
                            ROItac=TACDATA.CerC;
                            savestr = 'CerC';
                        case 'Pall'
                            ROItac=TACDATA.Pall;
                            savestr = 'Pall';
                        case 'Amy'
                            ROItac=TACDATA.Amy;
                            savestr = 'Amy';
                        case 'SN'
                            ROItac=TACDATA.SN;
                            savestr = 'SN';
                        case 'LC'
                            ROItac=TACDATA.LC;
                            savestr = 'LC';
                        case 'VTA'
                            ROItac=TACDATA.VTA;
                            savestr = 'VTA';
                        case 'SubcorticalROIs'
                            ROItac=TACDATA.SubcorticalROIs;
                            savestr = 'SubcorticalROIs';
                    end
                    
                    mROItac  = [ROItac(1)/2; (ROItac(2:end)+ROItac(1:end-1))/2];
                    ASRTM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mROItac(1:k).*dt(1:k));
                    end
                    %LSQ-estimation using lscov
                    [parest se_srtm mse_srtm]   = lscov(ASRTM,ROItac);
                    fittac=ASRTM*parest;
                    BP=parest(2)/parest(3)-1;
                    k2p=parest(2)/parest(1);
                    
                    %%%%% DO real SRTM
                    options = optimset('MaxFunEvals',1000);
                    weighs=[0.25*ones(30,1); ones(t_points-30,1)];
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-ROItac).*weighs);
                    [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
%                     modfit_esrtm=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
                    
                    %%%%% Do real SRTM up to end of Baseline
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-ROItac(1:end-11)).*weighs(1:end-11));
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
                                    roitac_gamma = ROItac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(ROItac))+roitac_gamma(1:length(ROItac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,ROItac);
                                    
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
                            end%                        breakpoint2=breakpoint;
                        end
                    end
                    
                    % save these and work out receptor occupancy later in
                    % the script
                    
                    y=-ASRTM(:,3)./ROItac;
                    x=ASRTM(:,2)./ROItac;
                    [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                    [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                    
                    k2=best_parest(2);
                    k2a=best_parest(3);
                    BP_lp=k2/k2a-1; % baseline binding potential of lpntPET
                    BP_srtm=BP__; % binding potential until the end of the task
                    BP_srtm_bsl=BP_bl; % binding potential until the end of baseline
                    G=best_parest(4); % gamma, obviously...
                    BPND = ((R1__*k2p)/k2a)-1; % correct? - not sure
%                     DBP =  k2./(k2a + best_actfun) - 1; % dynamic binding potential
                    DBP =  k2./(k2a + G*best_actfun) - 1; % is it this way...? confused
                    OCC = 100*(1-DBP/BP_lp);% receptor occupancy
%                     OCC2 = (1-DBP2/BP_lp);% receptor occupancy
%                     OCC = 100*(1-DBP/BP_lp); % receptor occupancy                    
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.' savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.' savestr '=OCC;' ])
%                     eval(['Occupancy2{id,d}.' savestr '=OCC2;' ])
                    eval(['BPND_save{id,d}.' savestr '=BPND;' ])
                    eval(['BP_srtm_save{id,d}.' savestr '=BP_srtm;' ])
                    eval(['BP_srtm_Bsl_save{id,d}.' savestr '=BP_srtm_bsl;' ])
                    eval(['BestActivationFunction{id,d}.' savestr '=best_actfun;'])
                    eval(['BestModelFit{id,d}.' savestr '=best_modelfit;'])
                    eval(['gamma_lp_save{id,d}.' savestr '=best_parest(4);' ])
                    eval(['Residuals{id,d}.' savestr '.BaselineFit=ROItac-fittac;' ])
                    eval(['Residuals{id,d}.' savestr '.AllFit_ESRTM=ROItac-modfit_esrtm;' ])
                    eval(['Residuals{id,d}.' savestr '.lpntPETFit=ROItac-best_modelfit;' ])
                    eval(['CompensatoryFunction{id,d}.' savestr '=best_parest(4)*best_actfun;' ])
                    
%                     BP_lp=best_parest(2)/best_parest(3)-1; % is this the BP from lpntPET?
%                     BP_lpntPET=(best_parest(2)/(best_parest(3)+(best_parest(4)).*best_actfun))-1; % isn't this the way?
%                     BP_bsl=BP_bl; % baseline from SRTM
                    % save BP for later analysis
%                     eval(['BP_lp_orig{id,d}.' savestr '=BP_lp;' ])
%                     eval(['BP_lp_save{id,d}.' savestr '=BP_lpntPET;' ])
%                     eval(['BP_bsl_save{id,d}.' savestr '=BP_bsl;' ])
                    
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
                    BPdataSave{id,d}.(reg{1}).BPND=BPND;
                    BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
%                     BPdataSave{id,d}.(reg{1}).Occupancy2=OCC2;
                    BPdataSave{id,d}.(reg{1}).DBP=DBP;
                    BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);
                    
                    
%                     BP_lp=best_parest(2)/best_parest(3)-1; % is this the BP from lpntPET?
%                     BP_lpntPET=(best_parest(2)/(best_parest(3)+(best_parest(4)).*best_actfun))-1; % isn't this the way?
%                     BP_bsl=BP_bl; % baseline from SRTM
                    % save BP for later analysis
%                     eval(['BP_lp_orig{id,d}.' savestr '=BP_lp;' ])
%                     eval(['BP_lp_save{id,d}.' savestr '=BP_lpntPET;' ])
%                     eval(['BP_bsl_save{id,d}.' savestr '=BP_bsl;' ])
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
                        [h p]=runstest(ROItac-fittac);
                        [h1 p1]=runstest(ROItac-modfit_esrtm);
                        [h2 p2]=runstest(ROItac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
    %                         title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
    %                         xlabel('Time (min)');
                        
    %                         subplot(3,2,4);
                        y=-ASRTM(:,3)./ROItac;
                        x=ASRTM(:,2)./ROItac;
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
                        BPdataSave{id,d}.(reg{1}).BPND=BPND;
                        BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
%                         BPdataSave{id,d}.(reg{1}).Occupancy2=OCC2;
                        BPdataSave{id,d}.(reg{1}).DBP=DBP;
                        BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);
                        
                        
                        close all
                        continue;
                    end
                    close all
                end
                close all
                
            end
            close all
            keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI BPdataSave BP_lp_save DBP_save Occupancy BPND_save BP_srtm_save BP_srtm_Bsl_save BestActivationFunction...
    BestModelFit gamma_lp_save Residuals CompensatoryFunction
        end
    end
end

disp('done')

% save(['/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/MRPET_BPpackage_All_smoothed3mm_' date '.mat'],...
%     'BPdataSave', 'BP_lp_save', 'DBP_save', 'Occupancy', 'BPND_save', 'BP_srtm_save', 'BP_srtm_Bsl_save', 'BestActivationFunction',...
%     'BestModelFit', 'gamma_lp_save', 'Residuals', 'CompensatoryFunction')
                    
                    
disp('all ROIs done')

%% plot occupancy per ROI

% occupancies = [];
reg={'Striatum','Caudate','Putamen','Accumbensarea','Thalamus','Hippocampus','SubcorticalROIs'};
for id = 1:length(IDs)
    for d=1:2
        if days(id,d) == 0
            disp('no measurement for this date')
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

%% mean of baseline SRTM BP
load('/Users/alex/Dropbox/paperwriting/MRPET/data/FS_labels_LCSNVTA.mat')
cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    
                else
                    % two-session
                    if d==1
                        disp('no data')
                        eval(['BPbsl_SRTM_' FSlabels1{labels,1}{1} '_highDA(id,1)=NaN;'])
                        eval(['BPall_SRTM_' FSlabels1{labels,1}{1} '_highDA(id,1)=NaN;'])
                        eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '_highDA(id,1)=NaN;'])
                        eval(['Gamma_' FSlabels1{labels,1}{1} '_highDA(id,1)=NaN;'])
                    elseif d==2
                        disp('no data')
                        eval(['BPbsl_SRTM_' FSlabels1{labels,1}{1} '_lowDA(id,1)=NaN;'])
                        eval(['BPall_SRTM_' FSlabels1{labels,1}{1} '_lowDA(id,1)=NaN;'])
                        eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '_lowDA(id,1)=NaN;'])
                        eval(['Gamma_' FSlabels1{labels,1}{1} '_lowDA(id,1)=NaN;'])
                    end
                end
            end
            
        else
            cnt=cnt+1;
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    
                else
                    
                    % singlesub
                    eval(['BPbsl_SRTM_' FSlabels1{labels,1}{1} '(cnt,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;'])
                    eval(['BPall_SRTM_' FSlabels1{labels,1}{1} '(cnt,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;'])
                    eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '(cnt,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_lpnt;'])
                    eval(['Gamma_' FSlabels1{labels,1}{1} '(cnt,1)=double(BPdataSave{id,d}.(FSlabels1{labels,1}).gamma);'])
                    
                    % two-session
                    if d==1
                    eval(['BPbsl_SRTM_' FSlabels1{labels,1}{1} '_highDA(id,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;'])
                    eval(['BPall_SRTM_' FSlabels1{labels,1}{1} '_highDA(id,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;'])
                    eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '_highDA(id,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_lpnt;'])
                    eval(['Gamma_' FSlabels1{labels,1}{1} '_highDA(id,1)=double(BPdataSave{id,d}.(FSlabels1{labels,1}).gamma);'])
                    elseif d==2
                    eval(['BPbsl_SRTM_' FSlabels1{labels,1}{1} '_lowDA(id,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl;'])
                    eval(['BPall_SRTM_' FSlabels1{labels,1}{1} '_lowDA(id,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm;'])
                    eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '_lowDA(id,1)=BPdataSave{id,d}.(FSlabels1{labels,1}).BP_lpnt;'])
                    eval(['Gamma_' FSlabels1{labels,1}{1} '_lowDA(id,1)=double(BPdataSave{id,d}.(FSlabels1{labels,1}).gamma);'])
                    end
                    
                end
            end
        end
    end
end

for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    
                else
                    eval(['Gamma_' FSlabels1{labels,1}{1} '_highDA=double(Gamma_' FSlabels1{labels,1}{1} '_highDA);'])
                    eval(['Gamma_' FSlabels1{labels,1}{1} '_lowDA=double(Gamma_' FSlabels1{labels,1}{1} '_lowDA);'])
                    eval(['Gamma_' FSlabels1{labels,1}{1} '=double(Gamma_' FSlabels1{labels,1}{1} ');'])
                    
                    eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '_highDA=double(BP_lpntPET_' FSlabels1{labels,1}{1} '_highDA);'])
                    eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '_lowDA=double(BP_lpntPET_' FSlabels1{labels,1}{1} '_lowDA);'])
                    eval(['BP_lpntPET_' FSlabels1{labels,1}{1} '=double(BP_lpntPET_' FSlabels1{labels,1}{1} ');'])
                    
                end
end


%% calculate occupancy per ROI

% occupancies = [];
reg={'Caud','Put','Nac','ThalProp','Hipp','SubcorticalROIs'};
for id = 1:length(IDs)
    for d=1:2
        if days(id,d) == 0
            disp('no measurement for this date')
            %             occupancies{id,d} = [];
        else
            figure;
            for r=[1 4 5 6]%:length(reg)
                eval(['size(Occupancy{id,d}.' reg{r} ')'])
                disp([num2str(IDs(id)) num2str(d)])
                eval(['plot(Occupancy{id,d}.' reg{r} ',''LineWidth'',3);hold on'])
            end
            title(['ID: ' num2str(IDs(id)) ' Session ' num2str(d)],'FontSize',30); ylim([-0.1 1]); xlim([60 115]);
            xline(100,'--','Linewidth',2) % task onset
            xlabel('Frame','Fontsize',20);
            ylabel('Percentage (%)','Fontsize',20)
            legend({'Striatum','Accumbensarea','Hippocampus','Task Onset'},'Location','northwest','Fontsize',15);
            set(gca,'fontsize',20)

            print('-dpsc2','-append','-bestfit',fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_occupancy_' date '.ps']));

            close all

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


%% output image

load('/Users/alex/Dropbox/literatures_IKND/FS_labels_LCSNVTA.mat')

% make maps
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            clear ZScoreImg baseimage FSlabel existingLabels
            
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

%% Collect BP_srtm change in each ROIs

BPchange_SRTM_save=[];
BPchange_lpnt_save=[];
cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    BPchange_SRTM_save{id,d}.(FSlabels1{labels,1}) = NaN;
                    if d==1
                        eval(['BPchange_SRTM_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                    else
                        eval(['BPchange_SRTM_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                    end

                end
            end

        else
            clear SRTMimg baseimage FSlabel existingLabels
                    cnt=cnt+1;
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    BPchange_SRTM_save{id,d}.(FSlabels1{labels,1}) ...
                        = (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl - BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm) / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);
                    eval(['BPchange_SRTM_' (FSlabels1{labels,1}{1}) '(cnt,1) = BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                    if d==1
                        eval(['BPchange_SRTM_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                    else
                        eval(['BPchange_SRTM_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = BPchange_SRTM_save{id,d}.(FSlabels1{labels,1});'])
                    end

                end
            end
        end

    end
end

cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    BPchange_lpnt_save{id,d}.(FSlabels1{labels,1}) = NaN;
                    if d==1
                        eval(['BPchange_lpnt_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = BPchange_lpnt_save{id,d}.(FSlabels1{labels,1});'])
                    else
                        eval(['BPchange_lpnt_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = BPchange_lpnt_save{id,d}.(FSlabels1{labels,1});'])
                    end
 
                end
            end
 
        else
            clear lpntimg baseimage FSlabel existingLabels
                    cnt=cnt+1;
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                else
                    BPchange_lpnt_save{id,d}.(FSlabels1{labels,1}) ...
                        = (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl - BPdataSave{id,d}.(FSlabels1{labels,1}).BP_lpnt) / (BPdataSave{id,d}.(FSlabels1{labels,1}).BP_srtm_bl);
                    eval(['BPchange_lpnt_' (FSlabels1{labels,1}{1}) '(cnt,1) = BPchange_lpnt_save{id,d}.(FSlabels1{labels,1});'])
                    if d==1
                        eval(['BPchange_lpnt_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = BPchange_lpnt_save{id,d}.(FSlabels1{labels,1});'])
                    else
                        eval(['BPchange_lpnt_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = BPchange_lpnt_save{id,d}.(FSlabels1{labels,1});'])
                    end
 
                end
            end
        end
 
    end
end
 
 


disp('done')


%% calculate area under the curve for PET data

load('/Users/alex/Dropbox/paperwriting/MRPET/data/FS_labels_LCSNVTA.mat')

cnt=0;
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    
                    % Calculate the area under the curve using the trapezoidal rule
                    eval(['Occupancy_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = NaN;']);
                    if d==1
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = NaN;'])
                    else
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = NaN;'])
                    end
                    
                else
                    % Calculate the area under the curve using the trapezoidal rule
                    eval(['Occupancy_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = NaN;']);
                    if d==1
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = NaN;'])
                    else
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = NaN;'])
                    end
                    
                end
            end
        else
            cnt=cnt+1;
            for labels=2:length(FSlabels1)
                if FSlabels1{labels,2}==0
                    
                    % Calculate the area under the curve using the trapezoidal rule
                    eval(['Occupancy_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = NaN;']);
                    eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '(cnt,1) = NaN;'])
                    if d==1
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = NaN;'])
                    else
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = NaN;'])
                    end
                    
                else
                    clear dat
                    dat=eval(['Occupancy{id,d}.' (FSlabels1{labels,1}{1}) ';']);
                    
                    % Calculate the area under the curve using the trapezoidal rule
                    eval(['Occupancy_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = trapz(dat);']);
                    eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '(cnt,1) = trapz(dat);'])
                    if d==1
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = trapz(dat);'])
                    else
                        eval(['Occupancy_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = trapz(dat);'])
                    end
                    
                end
            end
        end
    end
end


%%

cnt=0;
% for id = 1:length(IDs)
%     for d = 1:2
%         if days(id,d) == 0
%             for labels=2:length(FSlabels1)
%                 if FSlabels1{labels,2}==0
%                     
%                     % Calculate the area under the curve using the trapezoidal rule
%                     eval(['Occupancy2_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = NaN;']);
%                     if d==1
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = NaN;'])
%                     else
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = NaN;'])
%                     end
%                     
%                 else
%                     % Calculate the area under the curve using the trapezoidal rule
%                     eval(['Occupancy2_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = NaN;']);
%                     if d==1
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = NaN;'])
%                     else
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = NaN;'])
%                     end
%                     
%                 end
%             end
%         else
%             cnt=cnt+1;
%             for labels=2:length(FSlabels1)
%                 if FSlabels1{labels,2}==0
%                     
%                     % Calculate the area under the curve using the trapezoidal rule
%                     eval(['Occupancy2_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = NaN;']);
%                     eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '(cnt,1) = NaN;'])
%                     if d==1
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = NaN;'])
%                     else
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = NaN;'])
%                     end
%                     
%                 else
%                     clear dat
%                     dat=eval(['Occupancy2{id,d}.' (FSlabels1{labels,1}{1}) ';']);
%                     
%                     % Calculate the area under the curve using the trapezoidal rule
%                     eval(['Occupancy2_AUC{id,d}.' (FSlabels1{labels,1}{1}) ' = trapz(dat);']);
%                     eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '(cnt,1) = trapz(dat);'])
%                     if d==1
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_highDA(id,1) = trapz(dat);'])
%                     else
%                         eval(['Occupancy2_AUC_' (FSlabels1{labels,1}{1}) '_lowDA(id,1) = trapz(dat);'])
%                     end
%                     
%                 end
%             end
%         end
%     end
% end

%% create table - occupancy AUC single datasets

% Define the words to exclude
excludeWords = {'highDA', 'lowDA', 'tablestring', 'Singles', 'save'};

% Get all variable names in the workspace
allVarNames = who('Occupancy_AUC_*');

% Initialize an empty list for variable names to include
OccupancyAUCSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        OccupancyAUCSingles = [OccupancyAUCSingles; allVarNames{i}];
    end
end

% make the vars into a string list
OccupancyAUCSingles_tablestring=[];
for s1=1:length(OccupancyAUCSingles)
    
    if s1==1
        OccupancyAUCSingles_tablestring=OccupancyAUCSingles{s1};
    else
        OccupancyAUCSingles_tablestring=[OccupancyAUCSingles_tablestring ',' OccupancyAUCSingles{s1}];
    end

end

SubjID=ID';

eval(['MRPET_OccupancyAUC_singleSession_Table = table(SubjID,' OccupancyAUCSingles_tablestring ')'])
writetable(MRPET_OccupancyAUC_singleSession_Table,[ paths.exports 'MRPET_OccupancyAUC_SingleSession_smoothed3mm_' date '.xls'])

disp('done')
%%

% Define the words to exclude
excludeWords = {'tablestring', 'Singles', 'save'};

% Get all variable names in the workspace
allVarNames = who('Occupancy_AUC_*_highDA');

% Initialize an empty list for variable names to include
OccupancyAUChighDA = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        OccupancyAUChighDA = [OccupancyAUChighDA; allVarNames{i}];
    end
end

% make the vars into a string list
OccupancyAUChighDA_tablestring=[];
for s1=1:length(OccupancyAUChighDA)
    
    if s1==1
        OccupancyAUChighDA_tablestring=OccupancyAUChighDA{s1};
    else
        OccupancyAUChighDA_tablestring=[OccupancyAUChighDA_tablestring ',' OccupancyAUChighDA{s1}];
    end

end

% Define the words to exclude
excludeWords = {'tablestring', 'Singles', 'save'};

% Get all variable names in the workspace
allVarNames = who('Occupancy_AUC_*_lowDA');

% Initialize an empty list for variable names to include
OccupancyAUClowDA = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        OccupancyAUClowDA = [OccupancyAUClowDA; allVarNames{i}];
    end
end

% make the vars into a string list
OccupancyAUClowDA_tablestring=[];
for s1=1:length(OccupancyAUClowDA)
    
    if s1==1
        OccupancyAUClowDA_tablestring=OccupancyAUClowDA{s1};
    else
        OccupancyAUClowDA_tablestring=[OccupancyAUClowDA_tablestring ',' OccupancyAUClowDA{s1}];
    end

end

SubjID=IDs';

Tablestring_all=strcat(OccupancyAUChighDA_tablestring,',',OccupancyAUClowDA_tablestring);
eval(['MRPET_OccupancyAuC_bothSessions_Table = table(SubjID,' Tablestring_all ')'])
writetable(MRPET_OccupancyAuC_bothSessions_Table,[ paths.exports 'MRPET_MRPET_OccupancyAUC_bothSessions_Table_smoothed3mm_' date '.xls'])

disp('done')

%%


%% create table - single datasets

% Define the words to exclude
excludeWords = {'highDA', 'lowDA', 'tablestring', 'Singles', 'save'};

% Get all variable names in the workspace
allVarNames = who('BPchange_SRTM*');

% Initialize an empty list for variable names to include
BPchangeSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPchangeSingles = [BPchangeSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BPchangeSingles_tablestring=[];
for s1=1:length(BPchangeSingles)
    
    if s1==1
        BPchangeSingles_tablestring=BPchangeSingles{s1};
    else
        BPchangeSingles_tablestring=[BPchangeSingles_tablestring ',' BPchangeSingles{s1}];
    end

end




clear allVarNames

% Get all variable names in the workspace
allVarNames = who('BP_lpntPET*');

% Initialize an empty list for variable names to include
BP_lpntPETSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BP_lpntPETSingles = [BP_lpntPETSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BP_lpntPETSingles_tablestring=[];
for s1=1:length(BP_lpntPETSingles)
    
    if s1==1
        BP_lpntPETSingles_tablestring=BP_lpntPETSingles{s1};
    else
        BP_lpntPETSingles_tablestring=[BP_lpntPETSingles_tablestring ',' BP_lpntPETSingles{s1}];
    end

end




clear allVarNames

% Get all variable names in the workspace
allVarNames = who('BPbsl_SRTM*');

% Initialize an empty list for variable names to include
BPbslSRTMSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPbslSRTMSingles = [BPbslSRTMSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BPbslSingles_tablestring=[];
for s1=1:length(BPbslSRTMSingles)
    
    if s1==1
        BPbslSingles_tablestring=BPbslSRTMSingles{s1};
    else
        BPbslSingles_tablestring=[BPbslSingles_tablestring ',' BPbslSRTMSingles{s1}];
    end

end



clear allVarNames
 
% Get all variable names in the workspace
allVarNames = who('BPall_SRTM*');
 
% Initialize an empty list for variable names to include
BPallSRTMSingles = {};
 
% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPallSRTMSingles = [BPallSRTMSingles; allVarNames{i}];
    end
end
 
% make the vars into a string list
BPallSingles_tablestring=[];
for s1=1:length(BPallSRTMSingles)
    
    if s1==1
        BPallSingles_tablestring=BPallSRTMSingles{s1};
    else
        BPallSingles_tablestring=[BPallSingles_tablestring ',' BPallSRTMSingles{s1}];
    end
 
end




clear allVarNames
 
% Get all variable names in the workspace
allVarNames = who('Gamma_*');
 
% Initialize an empty list for variable names to include
GammaSingles = {};
 
% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        GammaSingles = [GammaSingles; allVarNames{i}];
    end
end
 
% make the vars into a string list
GammaSingles_tablestring=[];
for s1=1:length(GammaSingles)
    
    if s1==1
        GammaSingles_tablestring=GammaSingles{s1};
    else
        GammaSingles_tablestring=[GammaSingles_tablestring ',' GammaSingles{s1}];
    end
 
end

SubjID=ID';

Tablestring_all=strcat(BPchangeSingles_tablestring,',',BPbslSingles_tablestring,',', BPallSingles_tablestring,',',BP_lpntPETSingles_tablestring,',',GammaSingles_tablestring);
eval(['MRPET_SingleSessions_Table = table(SubjID,' Tablestring_all ')'])
writetable(MRPET_SingleSessions_Table,[ paths.exports 'MRPET_BPtable_SingleSession_smoothed3mm_' date '.xls'])

disp('done')

%% creat tables - high-low

% Define the words to exclude
excludeWords = {'tablestring', 'Singles', 'save'};

% Get all variable names in the workspace
allVarNames = who('BPchange_*highDA');

% Initialize an empty list for variable names to include
BPchangeSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPchangeSingles = [BPchangeSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BPchangeSingles_tablestring=[];
for s1=1:length(BPchangeSingles)
    
    if s1==1
        BPchangeSingles_tablestring=BPchangeSingles{s1};
    else
        BPchangeSingles_tablestring=[BPchangeSingles_tablestring ',' BPchangeSingles{s1}];
    end

end




clear allVarNames

% Get all variable names in the workspace
allVarNames = who('BP_lpntPET*highDA');

% Initialize an empty list for variable names to include
BP_lpntPETSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BP_lpntPETSingles = [BP_lpntPETSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BP_lpntPETSingles_tablestring=[];
for s1=1:length(BP_lpntPETSingles)
    
    if s1==1
        BP_lpntPETSingles_tablestring=BP_lpntPETSingles{s1};
    else
        BP_lpntPETSingles_tablestring=[BP_lpntPETSingles_tablestring ',' BP_lpntPETSingles{s1}];
    end

end




clear allVarNames

% Get all variable names in the workspace
allVarNames = who('BPbsl_SRTM*highDA');

% Initialize an empty list for variable names to include
BPbslSRTMSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPbslSRTMSingles = [BPbslSRTMSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BPbslSingles_tablestring=[];
for s1=1:length(BPbslSRTMSingles)
    
    if s1==1
        BPbslSingles_tablestring=BPbslSRTMSingles{s1};
    else
        BPbslSingles_tablestring=[BPbslSingles_tablestring ',' BPbslSRTMSingles{s1}];
    end

end



clear allVarNames
 
% Get all variable names in the workspace
allVarNames = who('BPall_SRTM*highDA');
 
% Initialize an empty list for variable names to include
BPallSRTMSingles = {};
 
% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPallSRTMSingles = [BPallSRTMSingles; allVarNames{i}];
    end
end
 
% make the vars into a string list
BPallSingles_tablestring=[];
for s1=1:length(BPallSRTMSingles)
    
    if s1==1
        BPallSingles_tablestring=BPallSRTMSingles{s1};
    else
        BPallSingles_tablestring=[BPallSingles_tablestring ',' BPallSRTMSingles{s1}];
    end
 
end




clear allVarNames
 
% Get all variable names in the workspace
allVarNames = who('Gamma_*highDA');
 
% Initialize an empty list for variable names to include
GammaSingles = {};
 
% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        GammaSingles = [GammaSingles; allVarNames{i}];
    end
end
 
% make the vars into a string list
GammaSingles_tablestring=[];
for s1=1:length(GammaSingles)
    
    if s1==1
        GammaSingles_tablestring=GammaSingles{s1};
    else
        GammaSingles_tablestring=[GammaSingles_tablestring ',' GammaSingles{s1}];
    end
 
end

SubjID=IDs';

% Tablestring_all=strcat(BPchangeSingles_tablestring,',',BPbslSingles_tablestring,',', BPallSingles_tablestring,',',BP_lpntPETSingles_tablestring,',',GammaSingles_tablestring);
Tablestring_all=strcat(BPchangeSingles_tablestring,',',BPbslSingles_tablestring,',', BPallSingles_tablestring,',',BP_lpntPETSingles_tablestring);
eval(['MRPET_highSession_Table = table(SubjID,' Tablestring_all ')'])
writetable(MRPET_highSession_Table,[ paths.exports 'MRPET_BPtable_highSession_smoothed3mm_' date '.xls'])

disp('done')


% Define the words to exclude
excludeWords = {'tablestring', 'Singles', 'save'};

% Get all variable names in the workspace
allVarNames = who('BPchange_*lowDA');

% Initialize an empty list for variable names to include
BPchangeSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPchangeSingles = [BPchangeSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BPchangeSingles_tablestring=[];
for s1=1:length(BPchangeSingles)
    
    if s1==1
        BPchangeSingles_tablestring=BPchangeSingles{s1};
    else
        BPchangeSingles_tablestring=[BPchangeSingles_tablestring ',' BPchangeSingles{s1}];
    end

end




clear allVarNames

% Get all variable names in the workspace
allVarNames = who('BP_lpntPET*lowDA');

% Initialize an empty list for variable names to include
BP_lpntPETSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BP_lpntPETSingles = [BP_lpntPETSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BP_lpntPETSingles_tablestring=[];
for s1=1:length(BP_lpntPETSingles)
    
    if s1==1
        BP_lpntPETSingles_tablestring=BP_lpntPETSingles{s1};
    else
        BP_lpntPETSingles_tablestring=[BP_lpntPETSingles_tablestring ',' BP_lpntPETSingles{s1}];
    end

end




clear allVarNames

% Get all variable names in the workspace
allVarNames = who('BPbsl_SRTM*lowDA');

% Initialize an empty list for variable names to include
BPbslSRTMSingles = {};

% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPbslSRTMSingles = [BPbslSRTMSingles; allVarNames{i}];
    end
end

% make the vars into a string list
BPbslSingles_tablestring=[];
for s1=1:length(BPbslSRTMSingles)
    
    if s1==1
        BPbslSingles_tablestring=BPbslSRTMSingles{s1};
    else
        BPbslSingles_tablestring=[BPbslSingles_tablestring ',' BPbslSRTMSingles{s1}];
    end

end



clear allVarNames
 
% Get all variable names in the workspace
allVarNames = who('BPall_SRTM*lowDA');
 
% Initialize an empty list for variable names to include
BPallSRTMSingles = {};
 
% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        BPallSRTMSingles = [BPallSRTMSingles; allVarNames{i}];
    end
end
 
% make the vars into a string list
BPallSingles_tablestring=[];
for s1=1:length(BPallSRTMSingles)
    
    if s1==1
        BPallSingles_tablestring=BPallSRTMSingles{s1};
    else
        BPallSingles_tablestring=[BPallSingles_tablestring ',' BPallSRTMSingles{s1}];
    end
 
end




clear allVarNames
 
% Get all variable names in the workspace
allVarNames = who('Gamma_*lowDA');
 
% Initialize an empty list for variable names to include
GammaSingles = {};
 
% Loop through each variable name
for i = 1:length(allVarNames)
    % Check if the variable name contains any of the exclude words
    containsExcludeWord = false;
    for j = 1:length(excludeWords)
        if contains(allVarNames{i}, excludeWords{j})
            containsExcludeWord = true;
            break;
        end
    end
    
    % If the variable name does not contain any of the exclude words, add it to the list of variables to include
    if ~containsExcludeWord
        GammaSingles = [GammaSingles; allVarNames{i}];
    end
end
 
% make the vars into a string list
GammaSingles_tablestring=[];
for s1=1:length(GammaSingles)
    
    if s1==1
        GammaSingles_tablestring=GammaSingles{s1};
    else
        GammaSingles_tablestring=[GammaSingles_tablestring ',' GammaSingles{s1}];
    end
 
end

SubjID=IDs';

% Tablestring_all=strcat(BPchangeSingles_tablestring,',',BPbslSingles_tablestring,',', BPallSingles_tablestring,',',BP_lpntPETSingles_tablestring,',',GammaSingles_tablestring);
Tablestring_all=strcat(BPchangeSingles_tablestring,',',BPbslSingles_tablestring,',', BPallSingles_tablestring,',',BP_lpntPETSingles_tablestring);
eval(['MRPET_lowSession_Table = table(SubjID,' Tablestring_all ')'])
writetable(MRPET_lowSession_Table,[ paths.exports 'MRPET_BPtable_lowSession_smoothed3mm_' date '.xls'])

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


IDnew=IDs';
eval(strcat('MRPET_BP_SRTM_Table = table(IDnew,', tablestr_high_bsl, ',',  tablestr_high_TASK, ',',tablestr_low_bsl, ',', tablestr_low_TASK, ')'))
writetable(MRPET_BP_SRTM_Table,['/Users/yeojin/Desktop/E_data/EB_cleaned/EBD_mrpet/RewardTask/PET/MRPET_BPtable_SRTM_' date '.xls'])
eval(strcat('MRPET_BPchange_SRTM_Table = table(IDnew,', tablestr_high_change,',', tablestr_low_change, ')'))
writetable(MRPET_BPchange_SRTM_Table,['/Users/yeojin/Desktop/E_data/EB_cleaned/EBD_mrpet/RewardTask/PET/MRPET_BPchange_SRTM_' date '.xls'])

disp('done')


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




