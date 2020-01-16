%% Load simulated data (a struct with many sizes and shapes)
load('NtSimul_RACBI_May2018.mat');

%% Set number of function evaluations
options = optimset('MaxFunEvals',1000);

%% Set common timing information
tmid=mean(NtSimul.plasmatimes,2);
t_points    = length(tmid);
dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
Times=NtSimul.plasmatimes;

%% Set baseline (no activation) and reference region TACs
SA=NtSimul.sa;
ENBasal=NtSimul.ENbasal;
Bmax=NtSimul.Bmax;
Baseline_TAC=NtSimul.ROIdata_Baseline.Tissue(1:end);
Baseline_FreeRAC=SA*NtSimul.ROIdata_Baseline.FreeRAC(1:end-1);
Baseline_BoundRAC=SA*NtSimul.ROIdata_Baseline.BoundRAC(1:end-1);
Baseline_BoundDA=NtSimul.ROIdata_Baseline.BoundDA(1:end-1);
Reference_TAC=NtSimul.REFdata_simul(1:end);
AllFields=fieldnames(NtSimul);
%return
%% Set analysis methods
AnaMethods={'LSRTM','LSRTM2','ESRTM','ESRTM2','lpntPET','lpntPET2'};%,'Logan','Ratio'
OCCDATab=array2table(zeros(t_points,length(AnaMethods)));
OCCDATab.Properties.VariableNames=AnaMethods;
nNoiseReal=1;

%% Filter parameters for a gaussian
Fa = 1;
Fb = [0.25, 0.5, 0.25]; 
for Sd=12%[0:2 4:2:12]
    
    Sd
    %% Loop over various shapes and sizes of activation
    for Activation=AllFields(end)%'Low_gammavar_early_fourth_alpha100','Intermediate_gammavar_early_fourth_alpha100',
        %Activation{1}
        %% Set activation TAC
        Activation_TAC=NtSimul.(Activation{1}).Tissue(1:end);
        Activation_FreeRAC=SA*NtSimul.(Activation{1}).FreeRAC(1:end-1);
        Activation_BoundRAC=SA*NtSimul.(Activation{1}).BoundRAC(1:end-1);
        Activation_BoundDA=NtSimul.(Activation{1}).BoundDA(1:end-1);
        %% Set task specific timing information
        if contains(Activation,'dual')
            n_tasks=2;
            StartTimes=[20,33];
            EndTimes=[26.5,40];
            BaseInd=find(Times==StartTimes(1),1)-1;
        else
            StartTimes=Times(NtSimul.(Activation{1}).StartIndt,1);
            EndTimes=Times(NtSimul.(Activation{1}).EndIndt,2);
            BaseInd=find(Times==StartTimes(1),1)-1;
            n_tasks=1;
        end
        
        %% Tabulate true values for this activation
        TrueParams.R1(Activation{1})=NtSimul.R1;
        TrueParams.BP(Activation{1})=(1-Activation_BoundDA(BaseInd)/NtSimul.Bmax)*NtSimul.Bmax/NtSimul.Kd;
        TrueDArat=(Activation_BoundDA - Activation_BoundDA(BaseInd))./NtSimul.Bmax; %% Subtract before activation ratio
        mDA  = [TrueDArat(1)/2; (TrueDArat(2:end)+TrueDArat(1:end-1))/2];
        TrueParams.AUC(Activation{1})=sum(mDA.*dt);
        [maxDA maxDAIdx]=max(TrueDArat);
        TrueParams.Peak(Activation{1})=maxDA;
        TrueParams.PeakTime(Activation{1})=tmid(maxDAIdx);
        
        %% Generate basis functions
        Options.FunctionName='gamma';
        Options.gamma_alpha=[0.1,.25,0.5,1,2,4];
        %     Options.FunctionName='exp';
        %     Options.exp_tau=[0.1,0.5,1];
        BF=generateBasisFunctions(Times,StartTimes,EndTimes,Options);
        for nn=1:nNoiseReal
            %% Generate noise
            roinoise=Sd*rand(t_points,1)-Sd*rand(t_points,1);
            refnoise=0.25*rand(t_points,1)-0.25*rand(t_points,1);
            if Sd==0
                refnoise=0; 
            end
            %% Set terms for integral (Baseline data without activation, used in fixing k2p)
            roitac=Baseline_TAC+roinoise;
            ntac=roitac;
            roitac(1:end)= filtfilt(Fb,Fa,roitac(1:end));
            mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2];
            reftac=Reference_TAC+refnoise;
            %reftac=filtfilt(Fb,Fa,reftac);
            mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
            
            %%%%%%%%%%%%%% Solving SRTM for baseline (no activation data)
            %% Set the SRTM part of A
            ASRTM = zeros(t_points ,3);
            ASRTM(:,1)  = reftac;
            for k = 1:t_points
                ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
            end
            for k = 1:t_points
                ASRTM(k,3)  = -sum(mroitac(1:k).*dt(1:k));
            end
            
            %LSQ-estimation using lscov
            [parest se mse]   = lscov(ASRTM,roitac);
            k2p=parest(2)/parest(1); %% Set fixed k2 ref region
            
            
            %% Set roitac (ACTIVATION), terms for trapezoidal integrals
            roinoise=Sd*rand(t_points,1)-Sd*rand(t_points,1);
            
            roitac=Activation_TAC+roinoise;
            mroitac  = [roitac(1)/2; (roitac(2:end)+roitac(1:end-1))/2];
            
            %%%%%%%%%%%%%%%%%% Solving with ESRTM
            for Mod={'ESRTM','ESRTM2'}
               switch Mod{1}
                   case 'ESRTM'
                       fobj = @(x) norm(simESRTM_1_0_0(tmid,reftac,t_points,x(1),x(2),[x(3)*ones(BaseInd,1); x(4)*ones(length(tmid)-BaseInd,1)])-roitac);
                       [parest minnorm]=fminsearch(@(x) fobj(x),[1 .3 3 3],options);
                       R1=parest(1);
                       k2=parest(2);
                       BP=parest(3);
                       DBP=[parest(3)*ones(BaseInd,1); parest(4)*ones(length(tmid)-BaseInd,1)];
                   case 'ESRTM2'
                       fobj = @(x) norm(simESRTMfixk2p_1_0_0(tmid,reftac,t_points,x(1),k2p,[x(2)*ones(BaseInd,1); x(3)*ones(length(tmid)-BaseInd,1)])-roitac);
                       [parest minnorm]=fminsearch(@(x) fobj(x),[1 3 3],options);
                       R1=parest(1);
                       BP=parest(2);
                       DBP=[parest(2)*ones(BaseInd,1); parest(3)*ones(length(tmid)-BaseInd,1)];
               end
               
               if BP>0 && BP<10
                   CurResults.BP.(Mod{1})(nn) = BP;
                   CurResults.R1.(Mod{1})(nn)=R1;
                   CurResults.mse.(Mod{1})(nn)=minnorm^2/t_points;
                   
                   % Calculate DA occupancy from dynamic BP
                   OCC = (1-DBP/BP);
                   BoundDA=OCC*(1-Activation_BoundDA(BaseInd)/Bmax);
                   OCCDATab.(Mod{1})=BoundDA;
               end
               
            end
            
            
            %%%%%%%%%%%%%%%%%% Solving with linearized methods
            %% Set up the linear system of equations (SRTM):
            % C(t0)=R1*CR(t0)+k2*int(CR(0-t0))-k2a*int(C(0-t0))-g*int(C(0-t0)*h(0-t0))
            % C(t1)=R1*CR(t1)+k2*int(CR(0-t1))-k2a*int(C(0-t1))-g*int(C(0-t1)*h(0-t1))
            % ...
            % C(tn)=R1*CR(tn)+k2*int(CR(0-tn))-k2a*int(C(0-tn))-g*int(C(0-tn)*h(0-tn))
            % i.e. C=AX, where A(:,1)=CR(t), A(:,2)=int(CR(t)), A(:,3)=int(C(t))
            % A(:,4+k)=int(C(t)*hk(t)), where hk is the k:th activation function
            % c.f. multiple challenges
            
            %% Set the SRTM part of A
            ASRTM = zeros(t_points ,3);
            ASRTM(:,1)  = reftac;
            for k = 1:t_points
                ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
            end
            for k = 1:t_points
                ASRTM(k,3)  = -sum(mroitac(1:k).*dt(1:k));
            end
            
            %LSQ-estimation using lscov (SRTM)
            %---------------------------------
            [parest se mse]   = lscov(ASRTM,roitac);
            CurResults.R1.SRTM(nn)=parest(1);
            CurResults.BP.SRTM(nn) = parest(2)/parest(3)-1;
            CurResults.mse.SRTM(nn)=mse;

            
            %% Set up the linear system of equations (SRTM2):
            % C(t0)=R1*[CR(t0)+k2p*int(CR(0-t0))]-k2a*int(C(0-t0))-g*int(C(0-t0)*h(0-t0))
            % C(t1)=R1*[CR(t1)+k2p*int(CR(0-t1))]-k2a*int(C(0-t1))-g*int(C(0-t1)*h(0-t1))
            % ...
            % C(tn)=R1*[CR(tn)+k2p*int(CR(0-tn))]-k2a*int(C(0-tn))-g*int(C(0-tn)*h(0-tn))
            % i.e. C=AX, where A(:,1)=[CR(t)+k2p*int(CR(0-tn))],  A(:,2)=int(C(t))
            % A(:,3+k)=int(C(t)*hk(t)), where hk is the k:th activation function
            % c.f. multiple challenges
            
            %% Set the SRTM part of A
            ASRTM2 = zeros(t_points ,2);
            for k = 1:t_points
                ASRTM2(k,1)  = reftac(k) +  k2p*sum(mreftac(1:k).*dt(1:k));
            end
            for k = 1:t_points
                ASRTM2(k,2)  = -sum(mroitac(1:k).*dt(1:k));
            end
              
            
            %% Prepare for LSRTM
            ActFun=BF(:,1); %% Set a fixed activation shape
            AAct=[];
            roitac_act=roitac.*ActFun;
            mroitac_act=[roitac_act(1)/2; (roitac_act(2:end)+roitac_act(1:end-1))/2];
            for k=1:t_points
                AAct(k)=-sum(mroitac_act(1:k).*dt(1:k));
            end
            for Mod={'LSRTM','LSRTM2'}
                %LSQ-estimation using lscov (LSRTM or LSRTM2)
                % -----------------------------
                switch Mod{1}
                    case 'LSRTM'
                        [parest se mse]   = lscov([ASRTM AAct'],roitac);
                        R1=parest(1);
                        k2=parest(2);
                        k2a=parest(3);
                        G=parest(4);
                    case 'LSRTM2'
                        [parest se mse]   = lscov([ASRTM2 AAct'],roitac);
                        R1=parest(1);
                        k2=parest(1)*k2p;
                        k2a=parest(2);
                        G=parest(3);
                end
                
                %baseline binding potential
                BP= k2/k2a-1;

                if BP>0 && BP<10
                    % Calculate DA occupancy from dynamic BP
                    DBP =  k2./(k2a + G*ActFun) - 1;
                    OCC = (1-DBP/BP);
                    BoundDA=OCC*(1-Activation_BoundDA(BaseInd)/Bmax);
                    OCCDATab.(Mod{1})=BoundDA;
                end
            end
            
            %% Prepare for lp-ntPET
            % Set the Activation part of A
            % Analyze basis function set
            ActFuns=[];
            AAct=[];
            Aidx=1;
            for taskIdx=1:n_tasks
                for columnIdx=1:size(BF(:,:,taskIdx),2)
                    if ~isempty(find(BF(:,columnIdx,taskIdx)))
                        %% Set current task activation function (null for others)
                        ActFuns(:,taskIdx,Aidx)=BF(:,columnIdx,taskIdx);
                        ActFuns(:,setdiff(1:n_tasks,taskIdx),Aidx)=0;
                        roitac_act=roitac.*BF(:,columnIdx,taskIdx);
                        mroitac_act=[roitac_act(1)/2; (roitac_act(2:end)+roitac_act(1:end-1))/2];
                        AAct(:,:,Aidx)=zeros(t_points,n_tasks);
                        for k=1:t_points
                            AAct(k,taskIdx,Aidx)=-sum(mroitac_act(1:k).*dt(1:k));
                        end
                        Aidx=Aidx+1;
                        
                        %% Set combined activation functions
                        if n_tasks>=taskIdx+1
                            AAct(:,:,Aidx)=zeros(t_points,n_tasks);
                            for taskIdx2=taskIdx:n_tasks
                                ActFuns(:,taskIdx2,Aidx)=BF(:,columnIdx,taskIdx2);
                                roitac_act=roitac.*BF(:,columnIdx,taskIdx2);
                                mroitac_act=[roitac_act(1)/2; (roitac_act(2:end)+roitac_act(1:end-1))/2];
                                for k=1:t_points
                                    AAct(k,taskIdx2,Aidx)=-sum(mroitac_act(1:k).*dt(1:k));
                                end
                                
                            end
                            Aidx=Aidx+1;
                            
                        end
                    end
                end
            end
            
            %% Calculate fits for all basis function combinations
            warning off MATLAB:lscov:RankDefDesignMat
            for Mod={'lpntPET','lpntPET2'}
                parest_all=[];
                se_all=[];
                mse_all=[];
                
                for ai=1:Aidx-1
                    
                    %LSQ-estimation using lscov (SRTM)
                    switch Mod{1}
                        case 'lpntPET'
                            [parest se mse]   = lscov([ASRTM, AAct(:,:,ai)],roitac);
                        case 'lpntPET2'
                            [parest se mse]   = lscov([ASRTM2, AAct(:,:,ai)],roitac);
                    end
                    parest_all(:,ai)=parest;
                    mse_all(ai)=mse;                 
                end
                
                [min_mse, opti]=min(mse_all);
                best_parest=parest_all(:,opti);
                
                %parameter estimates
                switch Mod{1}
                    case 'lpntPET'
                        k2=best_parest(2);
                        k2a=best_parest(3);
                        par_shift=3;
                    case 'lpntPET2'
                        k2=best_parest(1)*k2p;
                        k2a=best_parest(2);
                        par_shift=2;
                end
                
                %baseline binding potential
                BP= k2/k2a-1;
                if BP>0 && BP<10
                    
                    best_actfun=zeros(t_points,1);
                    for idx=1:n_tasks
                        G = best_parest(par_shift+idx);
                        best_actfun=best_actfun+G*ActFuns(:,idx,opti);
                    end
                    % Calculate DA occupancy from dynamic BP
                    DBP =  k2./(k2a + best_actfun) - 1;
                    OCC = (1-DBP/BP);
                    BoundDA=OCC*(1-Activation_BoundDA(BaseInd)/Bmax);
                    OCCDATab.(Mod{1})=BoundDA;
                end
                
                
            end
            

        end
        figure('position',[100 100 1000 800],'paperorientation','landscape');
                                   
        [hAx,hLine1,hLine2]=plotyy(tmid,[Activation_TAC,ntac,roitac,reftac],tmid,[TrueDArat,table2array(OCCDATab)]);
        ylabel(hAx(1),'Radioactivity concentration (kBq mL^{-1})');
        ylabel(hAx(2),'Ratio D2R occupied by DA');
        set(hAx(1),'FontSize',15);
        set(hAx(2),'FontSize',15);
        set(hLine1(1:end),'LineWidth',1.5);
        set(hLine2(2:end),'LineStyle','--');
        set(hLine2(1),'LineWidth',2);
        set(hLine2(2:end),'LineWidth',1.5);
        
        %set(hLine1(4),'LineWidth',0.5);
        %    set(hLine1(2),'LineWidth',2);
        set(hLine2(1),'LineWidth',2);
        %  set(hLine2(2),'LineWidth',2);
        legend('TRUE Tissue','Noisy ROI-tac','Filtered Noisy ROI-tac','Noisy CER-tac','TRUE OCC DA',...
            ['EST OCC DA (' AnaMethods{1} ')'], ...
            ['EST OCC DA (' AnaMethods{2} ')'], ...
            ['EST OCC DA (' AnaMethods{3} ')'], ...
            ['EST OCC DA (' AnaMethods{4} ')'], ...
            ['EST OCC DA (' AnaMethods{5} ')'], ...
            ['EST OCC DA (' AnaMethods{6} ')'], ...
            'Location','NorthWest');
        switch Activation{1}
            case 'Low_exp1'
                ttext=['One phase exponential decay: l=0.22 min^{-1}'];
            case 'Low_box1'
                ttext='Box';
            otherwise
                ttext=['Gamma-variate: t_{peak}=' num2str(NtSimul.(Activation{1}).t_peak) 'min, a=' num2str(NtSimul.(Activation{1}).alpha)];
        end
        xlim(hAx(1:2),[0 61]);
        ylim(hAx(1),[0 max(roitac)]);
        ylim(hAx(2),[0 .6]);
        %        title(ttext);
        title({['Noise-level=' num2str(Sd) ' (std)']; strrep(Activation{1},'_',' ')});
        %print('-dpsc','-append','-bestfit',['Figs/Validity_TACs_GFiltered.ps']);

    end
    %     return
    %save(['Data/DataFitResults_NoiseLevel' num2str(Sd) '.mat'],'FitData');
    %save(['Data/DataFitResults_TrueParams.mat'],'TrueParams');
end
