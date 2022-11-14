clc;clear


addpath('/Users/alex/Documents/MATLAB/NIfTI_20140122')
opengl hardwarebasic

% IDs
IDs = [4001];
days = [1 2];

path_savetac = ['/Users/alex/Documents/PET/4001_1/searchlight/tacs/'];
path_roitac  = ['/Users/alex/Documents/PET/4001_1/searchlight/ROItacs/'];

% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average
TrimPerc=15;

clear id d

% get cortex mask if you want to only include voxels within the cortex in
% same space as pet data

maskmat = spm_read_vols(spm_vol(['/Users/alex/Documents/PET/4001_1/brainonly.nii'])); % try using different ROI

CortexMask = maskmat;clear maskmat

xl=size(CortexMask,1);
yl=size(CortexMask,2);
zl=size(CortexMask,3);

%% define radius of SearchLight

MaxEuclidianDistanceToSeed=4; % this in voxels in our case, will be a sphere after all

%% first define volume of interest

VoxelsInCortex=find(CortexMask); % get number of voxel in mask / our data have 1mm voxel spacing

%% you need a template to retrieve distances later - make sure it has the same dims as your volume[x, y, z]=ndgrid(1:xl, 1:yl, 1:zl);%create a grid containing entries corresponding to coordinates - this will allow us to find coordinates of voxel indices later
[x, y, z]=ndgrid(1:xl, 1:yl, 1:zl);%create a grid containing entries corresponding to coordinates - this will allow us to find coordinates of voxel indices later

%% now loop through this using every voxel as seed voxel

PETflow = ['/Users/alex/Documents/PET/4001_1/c4001_PET_4D_InFlow1.nii'];
PETtask = ['/Users/alex/Documents/PET/4001_1/c4001_PET_4D_MT1.nii'];
PETbsl  = ['/Users/alex/Documents/PET/4001_1/c4001_PET_4D_Baseline1.nii'];

load([ path_roitac '40011_TACDATA_InFlow.mat']);
TACDATA_InFlow=TACDATA; clear TACDATA
load([ path_roitac '40011_TACDATA_Baseline.mat']);
TACDATA_Baseline=TACDATA; clear TACDATA
load([ path_roitac '40011_TACDATA_Task.mat']);
TACDATA_Task=TACDATA; clear TACDATA

% inFlow
clear DynPET temp ImgData_inflow
DynPET=load_nii(PETflow);
temp=size(DynPET.img);
ImgData_inflow=reshape(DynPET.img,prod(temp(1:3)),temp(4));
% BSL
clear DynPET temp ImgData_bsl
DynPET=load_nii(PETbsl);
temp=size(DynPET.img);
ImgData_bsl=reshape(DynPET.img,prod(temp(1:3)),temp(4));
% task
clear DynPET temp ImgData_task
DynPET=load_nii(PETtask);
temp=size(DynPET.img);
ImgData_task=reshape(DynPET.img,prod(temp(1:3)),temp(4));

%%
for SeedVoxelNo=1:length(VoxelsInCortex) % now we loop through seedvoxels (and have one pet model per voxel
    SeedVoxel=VoxelsInCortex(SeedVoxelNo); %we use only seedvoxels within cortex and retrieve index of SeedVoxel here
    %     try
    [SeedVoxelCoords]=[x(SeedVoxel), y(SeedVoxel), z(SeedVoxel)]; %retrieve coords from template which givex x,y, and z coordinate of this particular voxel
    % calculate index of voxels in searchlight
    xDistanceImage=x-SeedVoxelCoords(1);
    yDistanceImage=y-SeedVoxelCoords(2);
    zDistanceImage=z-SeedVoxelCoords(3); %one distance image for each dimension
    EuclidianDistanceImage=sqrt((xDistanceImage).^2 + (yDistanceImage).^2 + (zDistanceImage).^2); %sqrt of sum of squares gives you euclidean distance
    SearchlightMask=EuclidianDistanceImage<MaxEuclidianDistanceToSeed; %define searchlight by it's radius
    
    % you can check mask position by e.g. figure,imagesc(squeeze(SearchlightMask(:,:,12))),colorbar
    
    %now combine searchlight mask with cortical mask
    clear EffectiveMask tmpmask
    EffectiveMask=SearchlightMask&CortexMask;
    tmpmask = EffectiveMask; % write this one as niftii and use this for TACs
    NoOfEffectiveVoxels=length(find(EffectiveMask));%how many voxels are we effectively using for classification?
    
    %     % write a mask image
    %     hdr = spm_vol('/Users/alex/Documents/PET/4001_1/brainonly.nii'); % pick just any header from a file
    %     hdr.fname = ['/Users/alex/Documents/PET/4001_1/searchlight/sl' num2str(SeedVoxelNo) '.nii'];
    %     hdr.dim = size(tmpmask);
    %     hdr = rmfield(hdr,'pinfo');
    %     hdr.nii = spm_write_vol(hdr,tmpmask);
    
    SearchlightResults.SeedVoxelIndices(SeedVoxelNo)=SeedVoxel ;%please remember these things
    SearchlightResults.SeedVoxelCoords(SeedVoxelNo,:)=SeedVoxelCoords;
    SearchlightResults.NoOfEffectiveVoxels(SeedVoxelNo)=NoOfEffectiveVoxels; % possible to exclude those seed voxel results that are not completely in the brain later / % this is where you decide the number of searchlights
    
    % append more information about modelling results to this
    % SearchlightResults variable and use it for making images (same index!)
    
    % make the name for the searchlight
    clear name
    name = ['sl_' num2str(SeedVoxelNo)];
    %     name = ['sl_' num2str(SearchlightResults.SeedVoxelCoords(1)) '_' ...
    %         num2str(SearchlightResults.SeedVoxelCoords(2)) '_' ...
    %         num2str(SearchlightResults.SeedVoxelCoords(3)) ]
    
    %% calc TACs
    
    %     Mask1   =[ '/Users/alex/Documents/PET/4001_1/searchlight/sl' num2str(SeedVoxelNo) '.nii'];
    
    for ROI={name}
        TACDATA.(ROI{1}).coords=SearchlightResults.SeedVoxelCoords(SeedVoxelNo,:);
        TACDATA.(ROI{1}).NoOfEffectiveVoxels=SearchlightResults.NoOfEffectiveVoxels(SeedVoxelNo);
        TACDATA.(ROI{1}).SeedVoxelIndices=SearchlightResults.SeedVoxelIndices(SeedVoxelNo);
        
        % Read in searchlight image
        clear  ROIMask
        %         searchlightIMG=EffectiveMask;%load_nii(Mask1);
        ROIMask=single(round(EffectiveMask));
        %         ROIMask=round(searchlightIMG.img);
        IndList=[];
        for ROIInd=1
            if ~isnan(ROIInd)
                IndList=unique(union(IndList,find(ROIMask == ROIInd)));
            end
        end
        maskidx.(ROI{1})=IndList;
        
        maskidxcur = []; % input the index of the searchlight?
        maskidxcur = maskidx.(ROI{1});
        TACDATA.(ROI{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
        TACDATA.(ROI{1}).tac_inflow=trimmean(ImgData_inflow(maskidxcur,:),TrimPerc)';
        TACDATA.(ROI{1}).tac_bsl=trimmean(ImgData_bsl(maskidxcur,:),TrimPerc)';
        TACDATA.(ROI{1}).tac_task=trimmean(ImgData_task(maskidxcur,:),TrimPerc)';
        %         TACDATA.(ROI{1}).meantac=nanmean(TACDATA.(ROI{1}).tac);
        
        TACDATA_model.(ROI{1}) = [TACDATA.(ROI{1}).tac_inflow(10:end); ...
            TACDATA.(ROI{1}).tac_bsl;...
            TACDATA.(ROI{1}).tac_task];
        
        disp(ROI{1})
        %% calculate occupancy
        
        Lengths=[10*ones(30,1); 60*ones(55,1)]; % time windows
        tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
        Lengths=60*ones(15,1);
        tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
        Lengths=300*ones(11,1);
        tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
        %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
        times=[tt1(1:length(TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
        
        tmid=mean(times,2);
        tmidMin=tmid/60;
        t_points    = length(tmid);
        dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
        break_point=find(times(:,1)>=115*60,1,'first'); %% TIme of activation start
        
        for r=1
            for reg={'searchlight'}
                reftac=[TACDATA_InFlow.CerC.Bilateral.tac;...
                    TACDATA_Baseline.CerC.Bilateral.tac;...
                    TACDATA_Task.CerC.Bilateral.tac];
                mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
                %% Set the SRTM part of A
                ASRTM = zeros(t_points ,3);
                ASRTM(:,1)  = reftac;
                for k = 1:t_points
                    ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
                end
                refauc=sum(mreftac.*dt);
                roitac=TACDATA_model.(ROI{1});
                savestr = name;
                    
                mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2];
                ASTRM(:,3)=zeros(t_points,1);
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
                DBP =  k2./(k2a + G*best_actfun) - 1;
                OCC = 100*(1-DBP/BP_lp);
                eval(['BP_lp_save.' savestr '=BP_lp;' ])
                eval(['DBP_save.' savestr '=DBP;' ])
                eval(['Occupancy.' savestr '=OCC;' ])
                eval(['gamma_lp_save.' savestr '=best_parest(4);' ])
                
                TACDATA.(ROI{1}).occupancy=OCC;
                TACDATA.(ROI{1}).BP_lp=BP_lp;
                TACDATA.(ROI{1}).DBP=DBP;
                TACDATA.(ROI{1}).gamma_lp=best_parest(4);
                
            end
            
        end
    end
end


% save all
save([ path_savetac '40011_TACDATA_searchlight.mat'],'TACDATA');
save([ path_savetac 'searchlightResults.mat'],'SearchlightResults');

%% EffectiveMask is your mask for the voxels across which you calculate the TAC average.

%     TACDATA=[];
%             for reg={'CerC','Put','Caud','Nac','Hipp','ThalProp'}
%                 TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
%                     TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
%                     TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
%             end
%
%             Lengths=[10*ones(30,1); 60*ones(55,1)]; % time windows
%             tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
%             Lengths=60*ones(15,1);
%             tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
%             Lengths=300*ones(11,1);
%             tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];
%             %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
%             times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
%
%             tmid=mean(times,2);
%             tmidMin=tmid/60;
%             t_points    = length(tmid);
%             dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)]; % size of the time bin per frame
%             break_point=find(times(:,1)>=115*60,1,'first'); %% TIme of activation start
%
%
%             % start plotting the models
%             PlotStrFit=1; % will you plot the models?
%             Tthr=2; % something with threshold?
%             badcases={};
%             badcases2={};
%             BPdata=array2table(NaN*zeros(1,5)); % pre-assign binding potential data variable
%             Subj={[num2str(IDs(id)) num2str(d)]}; % subject id
%             BPdata.Properties.RowNames=Subj; % specify plot properties : subject ids
%             BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'}; % specify plot properties : plot names
%             for r=1 % for this subject
%                 for reg={'Striatum','Caudate','Putamen','Accumbens-area','Thalamus','Hippocampus'} % ROIs as regressors
%                     if PlotStrFit % if you decided to plot the models
%                         figure('Position',[100 100 800 1200],'Visible','on'); hold on;
%                         spidx=1;
%                     end
%                     reftac=TACDATA.CerC; % reference TAC data
%                     mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2]; % average of two hemispheres
%                     %% Set the SRTM part of A
%                     ASRTM = zeros(t_points ,3); % pre-assign SRTM variables
%                     ASRTM(:,1)  = reftac; % first column reference
%                     for k = 1:t_points % for the length of the entire scan frames
%                         ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k)); % the second column is mean of reference * time bin in seconds
%                     end
%                     refauc=sum(mreftac.*dt); % sum of the second column values
%
%                     if PlotStrFit
%                         % draw SRTM model of reference region
%                         subplot(3,2,[1 2]);
%                         plot(tmidMin,reftac,'k-'); hold on;
%                         legendtext={'Cerebellum'};
%                     end
%                     switch(reg{1}) % for each regressor, define TAC values of that ROI
%                         case 'Striatum'
%                             tempTac=[];
%                             tempVol=[];
%                             for subReg={'Caud','Put'}
%                                 tempTac=[tempTac, TACDATA.(subReg{1})];
%                                 tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol];
%                             end
%                             roitac=sum(tempTac.*tempVol,2)./sum(tempVol);
%                         case 'Caudate'
%                             roitac=TACDATA.Caud;
%                         case 'Putamen'
%                             roitac=TACDATA.Put;
%                         case 'Accumbens-area'
%                             roitac=TACDATA.Nac;
%                         case 'Thalamus'
%                             roitac=TACDATA.ThalProp;
%                         case 'Hippocampus'
%                             roitac=TACDATA.Hipp;
%
%                     end
%
%                     mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2]; % mean of ROI TAC between two hemispheres
%                     ASTRM(:,3)=zeros(t_points,1);
%                     for k = 1:t_points
%                         ASRTM(k,3)  = -sum(mroitac(1:k).*dt(1:k));
%                     end
%                     %LSQ-estimation using lscov
%                     [parest se_srtm mse_srtm]   = lscov(ASRTM,roitac); % calculate ordinary least squares
%                     fittac=ASRTM*parest; % calculate the fit curve
%                     % $parest(parameter estimate) is the N-by-1 vector that minimises the sum of squared errors
%                         % parest(1): Ratio K1/K1' how much tracers travel
%                         %       from one compartment to another
%                         % parest(2): Rate constant
%                         % parest(3): BP
%                     BP=parest(2)/parest(3)-1; %
%                     k2p=parest(2)/parest(1);
%
%                     %%%%% DO real SRTM
%                     options = optimset('MaxFunEvals',1000);
%                     weighs=[0.25*ones(30,1); ones(t_points-30,1)];
%                     fobj = @(x) norm((simESRTM_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-roitac).*weighs); % blue and red line, lamm
%                     [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
%                     R1__=parest_srtm(1);
%                     k2__=parest_srtm(2);
%                     BP__=parest_srtm(3);
%                     modfit_esrtm=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
%
%                     %%%%% Do real SRTM up to end of Baseline
%                     fobj = @(x) norm((simESRTM_1_0_0(tmidMin(1:end-11),reftac(1:end-11),t_points-11,x(1),x(2),x(3)*ones(t_points-11,1))-roitac(1:end-11)).*weighs(1:end-11));
%                     [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
%                     R1_bl=parest_srtm(1);
%                     k2_bl=parest_srtm(2);
%                     BP_bl=parest_srtm(3);
%                     modfit_esrtm_bl=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
%
%                     %%%%% Do real SRTM up to end of Baseline
%                     %     weighs=[0.25*ones(30,1); ones(15+46,1); 100*ones(11,1)];
%                     %     fobj = @(x) norm((simESRTM_1_0_0(tmidMin([1:76 92:t_points]),reftac([1:76 92:t_points]),t_points-15,x(1),x(2),x(3)*ones(t_points-15,1))-roitac([1:76 92:t_points])).*weighs([1:76 92:t_points]));
%                     %     [parest_srtm minnorm]=fminsearch(@(x) fobj(x),[1 .3 2],options);
%                     %     R1_task=parest_srtm(1);
%                     %     k2_task=parest_srtm(2);
%                     %     BP_task=parest_srtm(3);
%                     %     modfit_esrtm_task=simESRTM_1_0_0(tmidMin,reftac,t_points,parest_srtm(1),parest_srtm(2),parest_srtm(3)*ones(t_points,1));
%
%
%                     %%% Do lp-ntPET
%                     Alpntpet=zeros(t_points,4);
%                     Alpntpet(:,1:3)=ASRTM;
%                     best_mse=10^20;
%                     best_parest=[];
%                     best_se=[];
%                     % for estimating the optimal peak and alpha parameters
%                     for point_rise=break_point:length(tmid)-1
%                         t2p_index = find(tmid>tmid(point_rise)); %times(point_rise,1)
%                         if length(t2p_index)>1
%                             for t_ind=t2p_index(1):t2p_index(end-1)  %t=1:0.1:20 %
%                                 for alpha=[0.25 1 4]
%                                     t_peak=tmid(t_ind)-tmid(point_rise);   %  times(point_rise,1)
%                                     p = [1 alpha tmid(point_rise)+t_peak tmid(point_rise)];
%                                     actfun = zeros(size(tmid));
%                                     actfun(point_rise:t_points) = gamma_variate_madsen(p,tmid(point_rise:t_points)); % dispersion of bolus as it passes each compartment
%                                     roitac_gamma = roitac.*actfun;
%                                     mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(roitac))+roitac_gamma(1:length(roitac)-1))/2];
%
%                                     Alpntpet(:,4)=0;
%                                     for k = break_point:t_points
%                                         Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
%                                     end
%
%                                     %LSQ-estimation using lscov
%                                     [parest se mse]   = lscov(Alpntpet,roitac);
%
%                                     %estimated model TAC
%                                     modelfit = Alpntpet*parest;
%
%                                     %%%%%% magic %%%%%
%                                     if (best_mse > mse)
%                                         best_mse = mse;
%                                         best_parest=parest;  % R1 - k2 - k2a - gamma
%                                         best_modelfit=modelfit;
%                                         best_actfun=actfun;
%                                         best_se=se;
%                                     end
%                                 end
%                             end                                        %                        breakpoint2=breakpoint;
%                         end
%                     end
%                     BP_lp=best_parest(2)/best_parest(3)-1;
%
%                     fprintf(1,'Here\n')
%
%                     if PlotStrFit
%                         legendtext{end+1}=[reg{1} ' raw'];
%                         legendtext{end+1}=[reg{1} ' fit SRTM_{Baseline}'];
%                         legendtext{end+1}=[reg{1} ' fit SRTM_{All}'];
%                         legendtext{end+1}=[reg{1} ' fit lpnt-PET_{All}'];
%                         legendtext{end+1}=['Activation function'];
%                         SE=abs(BP)*sqrt((se_srtm(2)/parest(2))^2+(se_srtm(3)/parest(3))^2);
%                         subplot(3,2,[1 2]);
%                         yyaxis left;
%                         plot(tmidMin,roitac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
%                         %xlim([0 60]);
%                         xlabel('Time (min)');
%                         ylabel('Radioactivity concentration');
%                         legend(legendtext,'Location','south');
%                         yyaxis right;
%                         plot(tmidMin,best_parest(4)*best_actfun);
%                         ylabel('Compensatory function');
%                         %ylim([0 4*10^(-4)]);
%                         title([Subj{r} ' compartmental fits: BP_{Baseline}=' num2str(BP_bl,'%1.2f') ', BP_{All}=' num2str(BP__,'%1.2f') ', BP_{lp-nt}=' num2str(BP_lp,'%1.2f')]);
%                         subplot(3,2,3);
%                         %         plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'go',tmidMin,roitac-modfit_esrtm_bl,'co',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
%                         plot(tmidMin,roitac-fittac,'bo',tmidMin,roitac-modfit_esrtm,'bo',tmidMin,roitac-modfit_esrtm_bl,'ro',tmidMin,roitac-best_modelfit,'ko',[0 180],[0 0],'k--');
%                         [h p]=runstest(roitac-fittac);
%                         [h1 p1]=runstest(roitac-modfit_esrtm);
%                         [h2 p2]=runstest(roitac-best_modelfit);
%                         if p1<0.05
%                             badcases{end+1}=Subj{r};
%                         end
%                         if p2<0.05
%                             badcases2{end+1}=Subj{r};
%                         end
%                         title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
%                         xlabel('Time (min)');
%
%                         subplot(3,2,4);
%                         y=-ASRTM(:,3)./roitac;
%                         x=ASRTM(:,2)./roitac;
%                         [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
%                         [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
%                         plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
%                         title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
%                         xlabel(['\int REF/ROI']);
%                         ylabel('\int ROI/ROI')
%                         subplot(3,2,[5 6]);
%                         plot(tmidMin,roitac./reftac,'ko',[0 180],[0 0],'k--');
%                         %ylim([0.5 5]);
%                         title('Target to reference ratio');
%                         xlabel('Time (min)');
%
% %                         print('-dpsc2','-append','-bestfit',fullfile(paths.figures, [ num2str(IDs(id)) num2str(d) '_TAC_Fit_lpntpet_logan.ps']));
%                         %         pause
%                         %close(gcf)
%                         BPdata.BP_mrtm(Subj{r})=BP;
%                         BPdata.BP_srtm(Subj{r})=BP__;
%                         BPdata.BP_srtm_bl(Subj{r})=BP_bl;
%                         BPdata.BP_lpnt(Subj{r})=BP_lp;
%                         BPdata.BP_logan(Subj{r})=pp(1)-1;
%                         %return
%
%                         %%% F calc (model fit)
%
%
%                         continue;
%                     end
%                 end
%
%
%             end
% %             close all
%             keep IDs days paths id d TACs

%%%%%%%%%%%%
%     SearchlightResults.occupancy = output from pet model  % save!
%%%%%%%%%%%%%%
% end
