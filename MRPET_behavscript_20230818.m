%% do ROI extraction

clc;
clear;

% paths
paths = [];
paths.parent  = '/Users/alex/Dropbox/paperwriting/MRPET/data/';
paths.behav   = '/Users/alex/Dropbox/paperwriting/MRPET/data/behav/';
paths.ROI     = ['/Users/alex/Dropbox/Masks/Masks/mni_icbm152/'];

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

TR = 3.6;

% load ROIs
LCmask      = spm_read_vols(spm_vol([paths.ROI 'mni_icbm152_LCmetaMask_final.nii']));
SNVTAmask   = spm_read_vols(spm_vol(['/Users/alex/Dropbox/Masks/Masks/mni_icbm152/SN_VTA/alt/SNVTA.nii']));
BrainstemMask=spm_read_vols(spm_vol([paths.ROI 'brainstemMask_short.nii']));

% load experimental details
expdat = [];
for id = 1:length(IDs)
    
    % set up workspace
    %     mkdir([paths.taskonly num2str(IDs(id))]);
    
    %% gather info
    
    for d = 1:2
        if days(id,d) == 0
            fname_beh{id,d} = {NaN};
            expdat{id,d} = {NaN};
            contingency{id,d} = [];
        else
            fname_beh{id,d}     = [ num2str(IDs(id)) '_' num2str(days(id,d)) '.mat' ];
            expdat{id,d} = load([paths.behav fname_beh{id,d}]);
            
            % define contingency
            eval(['contingency{id,d}(1) = str2num(expdat{id,d}.dat.day' num2str(d) '.RewardCategory(9))']);
            if contingency{id,d}(1) == 1
                contingency{id,d}(2) = 2;
            elseif contingency{id,d}(1) == 3
                contingency{id,d}(2) = 4;
            elseif contingency{id,d}(1) == 2
                contingency{id,d}(2) = 1;
            elseif contingency{id,d}(1) == 4
                contingency{id,d}(2) = 3;
            end
            
            % task info
            clear stim_task stim_del stim_del_rem accuracies
            stim_task = eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.config.stim.fname']);
            if eval(['d' num2str(d) 'm(id,2)==0'])
                disp('no task data')
            else
                stim_del  = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.trl(:,[1 4])']);
                accuracies= eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.accu']);
                stim_del_rem  = stim_del(cell2mat(stim_del(:,2))==1 & accuracies==1,1); % old items and remembered
                stim_del_for  = stim_del(cell2mat(stim_del(:,2))==1 & accuracies==0,1); % old items and forgotten
                
                for q = 1:length(stim_task)
                    %                 idx = find(strcmp([C{:}], 'a'))
                    ind_remembered(q,1) = {find(strcmp(stim_task{q,1},stim_del_rem))};
                end
                tmp = cellfun(@isempty,ind_remembered); indx_remembered = tmp==0;
                
                for q = 1:length(stim_task)
                    %                 idx = find(strcmp([C{:}], 'a'))
                    ind_forgotten(q,1) = {find(strcmp(stim_task{q,1},stim_del_for))};
                end
                tmp = cellfun(@isempty,ind_forgotten); indx_forgotten= tmp==0;
            end
            
            
        end
    end
    
end

% 4008_2 had trials from 51-180
triallist_imm=expdat{8, 2}.dat.day2.memorytest.immediate.results.trl;
triallist_del=expdat{8, 2}.dat.day2.memorytest.delayed.results.trl;
triallist_main=expdat{8, 2}.dat.day2.maintask.results.trl(:,1);
cnt=0;
for trl=1:length(triallist_imm)
    
    if sum(strcmpi(triallist_imm{trl,1},triallist_main))==0 && triallist_imm{trl,4}==1
        trlind_imm(trl,1)=1; % delete
    elseif sum(strcmpi(triallist_imm{trl,1},triallist_main))==1 && triallist_imm{trl,4}==1
        trlind_imm(trl,1)=0; % keep
    elseif sum(strcmpi(triallist_imm{trl,1},triallist_main))==0 && triallist_imm{trl,4}==2
        trlind_imm(trl,1)=0; % keep
    end
    
end
expdat{8, 2}.dat.day2.memorytest.immediate.results.trl(logical(trlind_imm),:)=[];
expdat{8, 2}.dat.day2.memorytest.immediate.results.accu(logical(trlind_imm),:)=[];
expdat{8, 2}.dat.day2.memorytest.immediate.results.all(logical(trlind_imm),:)=[];
expdat{8, 2}.dat.day2.memorytest.immediate.results.confi(logical(trlind_imm),:)=[];
expdat{8, 2}.dat.day2.memorytest.immediate.results.resp(logical(trlind_imm),:)=[];
expdat{8, 2}.dat.day2.memorytest.immediate.results.rt(logical(trlind_imm),:)=[];
expdat{8, 2}.dat.day2.memorytest.immediate.config.stim.stimlist_all(logical(trlind_imm),:)=[];

for trl=1:length(triallist_del)
    
    if sum(strcmpi(triallist_del{trl,1},triallist_main))==0 && triallist_del{trl,4}==1
        trlind_del(trl,1)=1; % delete
    elseif sum(strcmpi(triallist_del{trl,1},triallist_main))==1 && triallist_del{trl,4}==1
        trlind_del(trl,1)=0; % keep
    elseif sum(strcmpi(triallist_del{trl,1},triallist_main))==0 && triallist_del{trl,4}==2
        trlind_del(trl,1)=0; % keep
    end
    
end
expdat{8, 2}.dat.day2.memorytest.delayed.results.trl(logical(trlind_del),:)=[];
expdat{8, 2}.dat.day2.memorytest.delayed.results.accu(logical(trlind_del),:)=[];
expdat{8, 2}.dat.day2.memorytest.delayed.results.all(logical(trlind_del),:)=[];
expdat{8, 2}.dat.day2.memorytest.delayed.results.confi(logical(trlind_del),:)=[];
expdat{8, 2}.dat.day2.memorytest.delayed.results.resp(logical(trlind_del),:)=[];
expdat{8, 2}.dat.day2.memorytest.delayed.results.rt(logical(trlind_del),:)=[];
expdat{8, 2}.dat.day2.memorytest.delayed.config.stim.stimlist_all(logical(trlind_imm),:)=[];


%% calc behavstats

behdat=[]; imms=[]; dels=[];
for id=1:length(IDs)
    
    for d=1:2
        
        if days(id,d) == 0
            behdat{id,d}={NaN};
            imms{id,d}={NaN};
            dels{id,d}={NaN};
        else
            
            behdat{id,d}=expdat{id,d};%load([paths.behav num2str(IDs(id)) '_' num2str(d) '.mat']);
            imms{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate;']);
            if eval(['d' num2str(d) 'm(id,2)==0'])
                disp('no task data')
            else
                dels{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.delayed;']);
            end
        end
    end
    
end


rews=[];neus=[];
for id=1:length(IDs)
    
    for d=1:2
        if days(id,d) == 0
            disp('skipped')
        else
            
            rews(id,d)=str2num(eval(['behdat{id,d}.dat.day' num2str(d) '.RewardCategory(9);']));
            if rews(id,d)==1
                neus(id,d)=2;
            elseif rews(id,d)==2
                neus(id,d)=1;
            elseif rews(id,d)==3
                neus(id,d)=4;
            elseif rews(id,d)==4
                neus(id,d)=3;
            end
            
        end
    end
    
end


%% calculate stats

% immediate recall
overallAccuracy_imm=[]; FAs_imm=[]; hit_FAs_imm=[]; overallAccuracy_imm=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            
            overallAccuracy_imm(id,d)=NaN; 
            FAs_imm(id,d)=NaN; 
%             hit_FAs_imm(id,d)=NaN;
            hits_imm(id,d)= NaN;
            hits_rew_imm(id,d) = NaN;
            hits_neu_imm(id,d) = NaN;
            
            FAs_imm(id,d) = NaN;
            FAs_rew_imm(id,d) = NaN;
            FAs_neu_imm(id,d) = NaN;
%             hit_FAs_rew_imm(id,d) = NaN;
%             hit_FAs_neu_imm(id,d) = NaN;
            
            % number of trials for d-prime calculation
            n_hit_imm(id,d) = NaN;
            n_fa_imm(id,d)  = NaN;
            n_miss_imm(id,d)= NaN;
            n_cr_imm(id,d)  = NaN;
            
            n_hit_rew_imm(id,d) = NaN;
            n_fa_rew_imm(id,d)  = NaN;
            n_miss_rew_imm(id,d)= NaN;
            n_cr_rew_imm(id,d)  = NaN;
            
            n_hit_neu_imm(id,d) = NaN;
            n_fa_neu_imm(id,d)  = NaN;
            n_miss_neu_imm(id,d)= NaN;
            n_cr_neu_imm(id,d)  = NaN;

            items_remembered_imm(id,d) = NaN;
            items_forgotten_imm(id,d) = NaN;
            items_old_total_imm(id,d) = NaN;
            
            confidence_imm(id,d) = NaN;
            confidence_rew_imm(id,d) = NaN;
            confidence_neu_imm(id,d) = NaN;
            conf_corrects_imm(id,d) = NaN;
            conf_rew_corrects_imm(id,d) = NaN;
            conf_neu_corrects_imm(id,d) = NaN;
            conf_wrongs_imm(id,d) = NaN;
            conf_rew_wrongs_imm(id,d)=NaN;
            conf_neu_wrongs_imm(id,d)=NaN;

        else
            
            % tidy the trials
            clear accu trls oldnew missed
            noresp=isnan(imms{id,d}.results.resp);
            accu=imms{id,d}.results.accu;
            trls=cell2mat(imms{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(imms{id,d}.results.trl(:,4))==1;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            
            overallAccuracy_imm(id,d)=nanmean(accu);
            hits_imm(id,d)= ( sum(oldnew==1 & accu==1) )/ sum(oldnew==1);
            hits_rew_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==1) )/ sum(oldnew==1 & trls==1);
            hits_neu_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==0) )/ sum(oldnew==1 & trls==0);
            FAs_imm(id,d)=( sum(oldnew==0 & accu==0) )/ sum(oldnew==0);
            FAs_rew_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==1) )/ sum(oldnew==0 & trls==1);
            FAs_neu_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==0) )/ sum(oldnew==0 & trls==0);
            
            % number of trials for d-prime calculation
            n_hit_imm(id,d) = sum(oldnew==1 & accu==1);
            n_fa_imm(id,d)  = sum(oldnew==0 & accu==0);
            n_miss_imm(id,d)= sum(oldnew==1 & accu==0);
            n_cr_imm(id,d)  = sum(oldnew==0 & accu==1);
            
            n_hit_rew_imm(id,d) = sum(oldnew==1 & accu==1 & trls==1);
            n_fa_rew_imm(id,d)  = sum(oldnew==0 & accu==0 & trls==1);
            n_miss_rew_imm(id,d)= sum(oldnew==1 & accu==0 & trls==1);
            n_cr_rew_imm(id,d)  = sum(oldnew==0 & accu==1 & trls==1);
            
            n_hit_neu_imm(id,d) = sum(oldnew==1 & accu==1 & trls==0);
            n_fa_neu_imm(id,d)  = sum(oldnew==0 & accu==0 & trls==0);
            n_miss_neu_imm(id,d)= sum(oldnew==1 & accu==0 & trls==0);
            n_cr_neu_imm(id,d)  = sum(oldnew==0 & accu==1 & trls==0);  
                
            items_remembered_imm(id,d) = sum(oldnew==1 & accu==1);
            items_forgotten_imm(id,d) = sum(oldnew==1 & accu==0);
            items_old_total_imm(id,d) = sum(oldnew==1)
            
            confidence_imm(id,d)=nanmean(imms{id,d}.results.confi);
            confidence_rew_imm(id,d)=nanmean(imms{id,d}.results.confi(trls==1));
            confidence_neu_imm(id,d)=nanmean(imms{id,d}.results.confi(trls==0));
            conf_corrects_imm(id,d)=nanmean(imms{id,d}.results.confi(accu==1));
            conf_rew_corrects_imm(id,d)=nanmean(imms{id,d}.results.confi(accu==1&trls==1));
            conf_neu_corrects_imm(id,d)=nanmean(imms{id,d}.results.confi(accu==1&trls==0));
            conf_wrongs_imm(id,d)=nanmean(imms{id,d}.results.confi(accu==0));
            conf_rew_wrongs_imm(id,d)=nanmean(imms{id,d}.results.confi(accu==0&trls==1));
            conf_neu_wrongs_imm(id,d)=nanmean(imms{id,d}.results.confi(accu==0&trls==0));
        end
        
    end
end

% delayed recall
hits_del=[]; FAs_del=[]; hit_FAs_del=[]; items_remembered_del=[]; items_forgotten_del=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            overallAccuracy_del(id,d)=NaN;
            hits_del(id,d)=NaN;
            hits_rew_del(id,d) = NaN;
            hits_neu_del(id,d) = NaN;
            FAs_del(id,d)=NaN;
            FAs_rew_del(id,d) = NaN;
            FAs_neu_del(id,d) = NaN;
            
            % number of trials for d-prime calculation
            n_hit_del(id,d) = NaN;
            n_fa_del(id,d)  = NaN;
            n_miss_del(id,d)= NaN;
            n_cr_del(id,d)  = NaN;
            
            n_hit_rew_del(id,d) = NaN;
            n_fa_rew_del(id,d)  = NaN;
            n_miss_rew_del(id,d)= NaN;
            n_cr_rew_del(id,d)  = NaN;
            
            n_hit_neu_del(id,d) = NaN;
            n_fa_neu_del(id,d)  = NaN;
            n_miss_neu_del(id,d)= NaN;
            n_cr_neu_del(id,d)  = NaN;  
            
            items_remembered_del(id,d) = NaN;
            items_forgotten_del(id,d) = NaN;
            items_old_total_del(id,d) = NaN;
            
            confidence_del(id,d)=NaN;
            confidence_rew_del(id,d)=NaN;
            confidence_neu_del(id,d)=NaN;
            conf_corrects_del(id,d)=NaN;
            conf_rew_corrects_del(id,d)=NaN;
            conf_neu_corrects_del(id,d)=NaN;
            conf_wrongs_del(id,d)=NaN;
            conf_rew_wrongs_del(id,d)=NaN;
            conf_neu_wrongs_del(id,d)=NaN;
        else
            % tidy the trials
            clear accu trls oldnew missed
            noresp=isnan(dels{id,d}.results.resp);
            accu=dels{id,d}.results.accu;
            trls=cell2mat(dels{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(dels{id,d}.results.trl(:,4))==1;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            
            overallAccuracy_del(id,d)=nanmean(accu);
            hits_del(id,d)=( sum(oldnew==1 & accu==1) )/ sum(oldnew==1);
            hits_rew_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==1) )/ sum(oldnew==1 & trls==1);
            hits_neu_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==0) )/ sum(oldnew==1 & trls==0);
            FAs_del(id,d)=( sum(oldnew==0 & accu==0) )/ sum(oldnew==0);
            FAs_rew_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==1) )/ sum(oldnew==0 & trls==1);
            FAs_neu_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==0) )/ sum(oldnew==0 & trls==0);

            
            % number of trials for d-prime calculation
            n_hit_del(id,d) = sum(oldnew==1 & accu==1);
            n_fa_del(id,d)  = sum(oldnew==0 & accu==0);
            n_miss_del(id,d)= sum(oldnew==1 & accu==0);
            n_cr_del(id,d)  = sum(oldnew==0 & accu==1);
            
            n_hit_rew_del(id,d) = sum(oldnew==1 & accu==1 & trls==1);
            n_fa_rew_del(id,d)  = sum(oldnew==0 & accu==0 & trls==1);
            n_miss_rew_del(id,d)= sum(oldnew==1 & accu==0 & trls==1);
            n_cr_rew_del(id,d)  = sum(oldnew==0 & accu==1 & trls==1);
            
            n_hit_neu_del(id,d) = sum(oldnew==1 & accu==1 & trls==0);
            n_fa_neu_del(id,d)  = sum(oldnew==0 & accu==0 & trls==0);
            n_miss_neu_del(id,d)= sum(oldnew==1 & accu==0 & trls==0);
            n_cr_neu_del(id,d)  = sum(oldnew==0 & accu==1 & trls==0);

            items_remembered_del(id,d) = sum(oldnew==1 & accu==1);
            items_forgotten_del(id,d) = sum(oldnew==1 & accu==0);
            items_old_total_del(id,d) = sum(oldnew==1);
            
            confidence_del(id,d)=nanmean(dels{id,d}.results.confi);
            confidence_rew_del(id,d)=nanmean(dels{id,d}.results.confi(trls==1));
            confidence_neu_del(id,d)=nanmean(dels{id,d}.results.confi(trls==0));
            conf_corrects_del(id,d)=nanmean(dels{id,d}.results.confi(accu==1));
            conf_rew_corrects_del(id,d)=nanmean(dels{id,d}.results.confi(accu==1&trls==1));
            conf_neu_corrects_del(id,d)=nanmean(dels{id,d}.results.confi(accu==1&trls==0));
            conf_wrongs_del(id,d)=nanmean(dels{id,d}.results.confi(accu==0));
            conf_rew_wrongs_del(id,d)=nanmean(dels{id,d}.results.confi(accu==0&trls==1));
            conf_neu_wrongs_del(id,d)=nanmean(dels{id,d}.results.confi(accu==0&trls==0));
        end
    end
end

%% rts

for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            
            disp('skipped')
            RT_all_immediate{id,d} = NaN;
            RTmean_immediate_rew(id,d) = NaN;
            RTmean_immediate_neu(id,d) = NaN;
            
            RTmean_immediate_rew_FA(id,d)   = NaN;
            RTmean_immediate_rew_hit(id,d)  = NaN;
            RTmean_immediate_rew_CR(id,d)   = NaN;
            RTmean_immediate_rew_miss(id,d) = NaN;
            
            RTmean_immediate_neu_FA(id,d)   = NaN;
            RTmean_immediate_neu_hit(id,d)  = NaN;
            RTmean_immediate_neu_CR(id,d)   = NaN;
            RTmean_immediate_neu_miss(id,d) = NaN;

            RT_all_delayed{id,d} = NaN;
            RTmean_delayed_rew(id,d) = NaN;
            RTmean_delayed_neu(id,d) = NaN;
            
            RTmean_delayed_rew_FA(id,d)   = NaN;
            RTmean_delayed_rew_hit(id,d)  = NaN;
            RTmean_delayed_rew_CR(id,d)   = NaN;
            RTmean_delayed_rew_miss(id,d) = NaN;
            
            RTmean_delayed_neu_FA(id,d)   = NaN;
            RTmean_delayed_neu_hit(id,d)  = NaN;
            RTmean_delayed_neu_CR(id,d)   = NaN;
            RTmean_delayed_neu_miss(id,d) = NaN;
            
        else
            
            % --- immediate --- %
            
            trlinfo{id,d} =  eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.config.stim.stimlist_all(:,2))==rews(id,d);']);
            
            % tidy the trials
            clear accu trls oldnew missed tmptrl noresp
            noresp  = isnan(imms{id,d}.results.resp);
            accu    = imms{id,d}.results.accu;
            trls    = cell2mat(imms{id,d}.results.trl(:,2))==rews(id,d);
            oldnew  = cell2mat(imms{id,d}.results.trl(:,4))==1;
            tmptrl  = trlinfo{id,d};
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];     
            tmptrl(noresp)=[];
            
            RT_all_immediate{id,d} = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.rt;']);
            RT_all_immediate{id,d}(noresp)=[];
            
            RTmean_immediate_rew(id,d)      = nanmean(RT_all_immediate{id,d}(tmptrl==1));
            RTmean_immediate_rew_FA(id,d)   = nanmean(RT_all_immediate{id,d}(tmptrl==1 & oldnew==0 & accu==0));
            RTmean_immediate_rew_hit(id,d)  = nanmean(RT_all_immediate{id,d}(tmptrl==1 & oldnew==1 & accu==1));
            RTmean_immediate_rew_CR(id,d)   = nanmean(RT_all_immediate{id,d}(tmptrl==1 & oldnew==0 & accu==1));
            RTmean_immediate_rew_miss(id,d) = nanmean(RT_all_immediate{id,d}(tmptrl==1 & oldnew==1 & accu==0));
            
            RTmean_immediate_neu(id,d)      = nanmean(RT_all_immediate{id,d}(tmptrl==0));
            RTmean_immediate_neu_FA(id,d)   = nanmean(RT_all_immediate{id,d}(tmptrl==0 & oldnew==0 & accu==0));
            RTmean_immediate_neu_hit(id,d)  = nanmean(RT_all_immediate{id,d}(tmptrl==0 & oldnew==1 & accu==1));
            RTmean_immediate_neu_CR(id,d)   = nanmean(RT_all_immediate{id,d}(tmptrl==0 & oldnew==0 & accu==1));
            RTmean_immediate_neu_miss(id,d) = nanmean(RT_all_immediate{id,d}(tmptrl==0 & oldnew==1 & accu==0));
            
            
            
            % --- delayed --- %
            
            trlinfo{id,d} =  eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.config.stim.stimlist_all(:,2))==rews(id,d);']);
            
            % tidy the trials
            clear accu trls oldnew missed tmptrl noresp
            noresp  = isnan(dels{id,d}.results.resp);
            accu    = dels{id,d}.results.accu;
            trls    = cell2mat(dels{id,d}.results.trl(:,2))==rews(id,d);
            oldnew  = cell2mat(dels{id,d}.results.trl(:,4))==1;
            tmptrl  = trlinfo{id,d};
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];     
            tmptrl(noresp)=[];
            
            RT_all_delayed{id,d} = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.rt;']);
            RT_all_delayed{id,d}(noresp)=[];
                        
            RTmean_delayed_rew(id,d) = nanmean(RT_all_delayed{id,d}(trls==1));
            RTmean_delayed_rew_FA(id,d)   = nanmean(RT_all_delayed{id,d}(trls==1 & oldnew==0 & accu==0));
            RTmean_delayed_rew_hit(id,d)  = nanmean(RT_all_delayed{id,d}(trls==1 & oldnew==1 & accu==1));
            RTmean_delayed_rew_CR(id,d)   = nanmean(RT_all_delayed{id,d}(trls==1 & oldnew==0 & accu==1));
            RTmean_delayed_rew_miss(id,d) = nanmean(RT_all_delayed{id,d}(trls==1 & oldnew==1 & accu==0));
            
            RTmean_delayed_neu(id,d) = nanmean(RT_all_delayed{id,d}(trls==0));
            RTmean_delayed_neu_FA(id,d)   = nanmean(RT_all_delayed{id,d}(trls==0 & oldnew==0 & accu==0));
            RTmean_delayed_neu_hit(id,d)  = nanmean(RT_all_delayed{id,d}(trls==0 & oldnew==1 & accu==1));
            RTmean_delayed_neu_CR(id,d)   = nanmean(RT_all_delayed{id,d}(trls==0 & oldnew==0 & accu==1));
            RTmean_delayed_neu_miss(id,d) = nanmean(RT_all_delayed{id,d}(trls==0 & oldnew==1 & accu==0));
            
            
            % --- main task --- %
            trlinfom{id,d} =  eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,2))==rews(id,d);']);
            RT_all_maintask{id,d} = eval(['expdat{id,d}.dat.day' num2str(d) '.maintask.results.rt;']);
            RTmean_maintask_rew(id,d) = nanmean(RT_all_maintask{id,d}(trlinfom{id,d}==1));
            RTmean_maintask_neu(id,d) = nanmean(RT_all_maintask{id,d}(trlinfom{id,d}==0));
            
            
        end
    end
end

MT_RT_lowDA_Rew     = RTmean_maintask_rew(:,2);
MT_RT_highDA_Rew    = RTmean_maintask_rew(:,1);
MT_RT_highDA_Neu     = RTmean_maintask_neu(:,1);
MT_RT_lowDA_Neu      = RTmean_maintask_neu(:,2);


%% miss rate

for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
        else
            clear missedResponse
            missedResponse = isnan(imms{id,d}.results.resp);
            missedResponse_rate(id,d) = sum(missedResponse)/numel(missedResponse);
            missedResponseRate_reward(id,d) = sum((missedResponse==1 & trlinfo{id,d}==1))/sum(trlinfo{id,d}==1);
            missedResponseRate_neutral(id,d) = sum((missedResponse==1 & trlinfo{id,d}==0))/sum(trlinfo{id,d}==0);
        end
    end
end

%% correct rejection & miss

for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            
            % miss
            Miss_all_imm(id,d) = NaN;
            Miss_rew_imm(id,d) = NaN;
            Miss_neu_imm(id,d) = NaN;
            
            % correct rejection
            CorrectRejection_all_imm(id,d) = NaN;
            CorrectRejection_rew_imm(id,d) = NaN;
            CorrectRejection_neu_imm(id,d) = NaN;
            
        else
            
            clear accu trls oldnew noresp
            noresp=isnan(imms{id,d}.results.resp);
            accu=imms{id,d}.results.accu;
            trls=cell2mat(imms{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(imms{id,d}.results.trl(:,4))==1;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            
            % miss
            Miss_all_imm(id,d) = sum(oldnew==1 & accu==0) / sum(oldnew==1);
            Miss_rew_imm(id,d) = sum(oldnew==1 & accu==0 & trls==1) / sum(oldnew==1 & trls==1);
            Miss_neu_imm(id,d) = sum(oldnew==1 & accu==0 & trls==0) / sum(oldnew==1 & trls==0);
            
            % correct rejection
            CorrectRejection_all_imm(id,d) = sum(oldnew==0 & accu==1) / sum(oldnew==0);
            CorrectRejection_rew_imm(id,d) = sum(oldnew==0 & accu==1 & trls==1) / sum(oldnew==0 & trls==1);
            CorrectRejection_neu_imm(id,d) = sum(oldnew==0 & accu==1 & trls==0) / sum(oldnew==0 & trls==0);
            
        end
    end
end



for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            
            % miss
            Miss_all_del(id,d) = NaN;
            Miss_rew_del(id,d) = NaN;
            Miss_neu_del(id,d) = NaN;
            
            % correct rejection
            CorrectRejection_all_del(id,d) = NaN;
            CorrectRejection_rew_del(id,d) = NaN;
            CorrectRejection_neu_del(id,d) = NaN;
            
        else
            
            clear accu trls oldnew noresp
            noresp=isnan(dels{id,d}.results.resp);
            accu=dels{id,d}.results.accu;
            trls=cell2mat(dels{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(dels{id,d}.results.trl(:,4))==1;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            
            % miss
            Miss_all_del(id,d) = sum(oldnew==1 & accu==0) / sum(oldnew==1);
            Miss_rew_del(id,d) = sum(oldnew==1 & accu==0 & trls==1) / sum(oldnew==1 & trls==1);
            Miss_neu_del(id,d) = sum(oldnew==1 & accu==0 & trls==0) / sum(oldnew==1 & trls==0);
            
            % correct rejection
            CorrectRejection_all_del(id,d) = sum(oldnew==0 & accu==1) / sum(oldnew==0);
            CorrectRejection_rew_del(id,d) = sum(oldnew==0 & accu==1 & trls==1) / sum(oldnew==0 & trls==1);
            CorrectRejection_neu_del(id,d) = sum(oldnew==0 & accu==1 & trls==0) / sum(oldnew==0 & trls==0);
            
        end
    end
end


%% confident trials only

% immediate recall
overallAccuracy_highconf_imm=[]; FAs_highconf_imm=[]; hit_FAs_highconf_imm=[]; overallAccuracy_highconf_imm=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            
            overallAccuracy_highconf_imm(id,d)=NaN; 
            FAs_highconf_imm(id,d)=NaN; 
%             hit_FAs_highconf_imm(id,d)=NaN;
            hits_highconf_imm(id,d)= NaN;
            hits_rew_highconf_imm(id,d) = NaN;
            hits_neu_highconf_imm(id,d) = NaN;
            
            FAs_highconf_imm(id,d) = NaN;
            FAs_rew_highconf_imm(id,d) = NaN;
            FAs_neu_highconf_imm(id,d) = NaN;
            
            % number of trials for d-prime calculation
            n_hit_highconf_imm(id,d) = NaN;
            n_fa_highconf_imm(id,d)  = NaN;
            n_miss_highconf_imm(id,d)= NaN;
            n_cr_highconf_imm(id,d)  = NaN;
            
            n_hit_rew_highconf_imm(id,d) = NaN;
            n_fa_rew_highconf_imm(id,d)  = NaN;
            n_miss_rew_highconf_imm(id,d)= NaN;
            n_cr_rew_highconf_imm(id,d)  = NaN;
            
            n_hit_neu_highconf_imm(id,d) = NaN;
            n_fa_neu_highconf_imm(id,d)  = NaN;
            n_miss_neu_highconf_imm(id,d)= NaN;
            n_cr_neu_highconf_imm(id,d)  = NaN;

            items_remembered_highconf_imm(id,d) = NaN;
            items_forgotten_highconf_imm(id,d) = NaN;
            items_old_total_highconf_imm(id,d) = NaN;
            
        else
            
            % tidy the trials
            clear accu trls oldnew missed confs
            noresp=isnan(imms{id,d}.results.resp);
            accu=imms{id,d}.results.accu;
            trls=cell2mat(imms{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(imms{id,d}.results.trl(:,4))==1;
            confs=imms{id,d}.results.confi;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            confs(noresp)=[];
            
            overallAccuracy_highconf_imm(id,d)=nanmean(accu(confs==1));
            hits_highconf_imm(id,d)= ( sum(oldnew==1 & accu==1 & confs==1) )/ sum(oldnew==1 & confs==1);
            hits_rew_highconf_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==1 & confs==1) )/ sum(oldnew==1 & trls==1 & confs==1);
            hits_neu_highconf_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==0 & confs==1) )/ sum(oldnew==1 & trls==0 & confs==1);
            FAs_highconf_imm(id,d)=( sum(oldnew==0 & accu==0 & confs==1) )/ sum(oldnew==0 & confs==1);
            FAs_rew_highconf_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==1 & confs==1) )/ sum(oldnew==0 & trls==1 & confs==1);
            FAs_neu_highconf_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==0 & confs==1) )/ sum(oldnew==0 & trls==0 & confs==1);
            
            % number of trials for d-prime calculation
            n_hit_highconf_imm(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            n_fa_highconf_imm(id,d)  = sum(oldnew==0 & accu==0 & confs==1);
            n_miss_highconf_imm(id,d)= sum(oldnew==1 & accu==0 & confs==1);
            n_cr_highconf_imm(id,d)  = sum(oldnew==0 & accu==1 & confs==1);
            
            n_hit_rew_highconf_imm(id,d) = sum(oldnew==1 & accu==1 & trls==1 & confs==1);
            n_fa_rew_highconf_imm(id,d)  = sum(oldnew==0 & accu==0 & trls==1 & confs==1);
            n_miss_rew_highconf_imm(id,d)= sum(oldnew==1 & accu==0 & trls==1 & confs==1);
            n_cr_rew_highconf_imm(id,d)  = sum(oldnew==0 & accu==1 & trls==1 & confs==1);
            
            n_hit_neu_highconf_imm(id,d) = sum(oldnew==1 & accu==1 & trls==0 & confs==1);
            n_fa_neu_highconf_imm(id,d)  = sum(oldnew==0 & accu==0 & trls==0 & confs==1);
            n_miss_neu_highconf_imm(id,d)= sum(oldnew==1 & accu==0 & trls==0 & confs==1);
            n_cr_neu_highconf_imm(id,d)  = sum(oldnew==0 & accu==1 & trls==0 & confs==1);  
                
            items_remembered_highconf_imm(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            items_forgotten_highconf_imm(id,d) = sum(oldnew==1 & accu==0 & confs==1);
            items_old_total_highconf_imm(id,d) = sum(oldnew==1 & confs==1);
            
        end
        
    end
end

% delayed recall
hits_highconf_del=[]; FAs_highconf_del=[]; hit_FAs_highconf_del=[]; items_remembered_highconf_del=[]; items_forgotten_highconf_del=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            overallAccuracy_highconf_del(id,d)=NaN;
            hits_highconf_del(id,d)=NaN;
            hits_rew_highconf_del(id,d) = NaN;
            hits_neu_highconf_del(id,d) = NaN;
            FAs_highconf_del(id,d)=NaN;
            FAs_rew_highconf_del(id,d) = NaN;
            FAs_neu_highconf_del(id,d) = NaN;
            
            % number of trials for d-prime calculation
            n_hit_highconf_del(id,d) = NaN;
            n_fa_highconf_del(id,d)  = NaN;
            n_miss_highconf_del(id,d)= NaN;
            n_cr_highconf_del(id,d)  = NaN;
            
            n_hit_rew_highconf_del(id,d) = NaN;
            n_fa_rew_highconf_del(id,d)  = NaN;
            n_miss_rew_highconf_del(id,d)= NaN;
            n_cr_rew_highconf_del(id,d)  = NaN;
            
            n_hit_neu_highconf_del(id,d) = NaN;
            n_fa_neu_highconf_del(id,d)  = NaN;
            n_miss_neu_highconf_del(id,d)= NaN;
            n_cr_neu_highconf_del(id,d)  = NaN;  
            
            items_remembered_highconf_del(id,d) = NaN;
            items_forgotten_highconf_del(id,d) = NaN;
            items_old_total_highconf_del(id,d) = NaN;
            
            confidence_highconf_del(id,d)=NaN;
            confidence_rew_highconf_del(id,d)=NaN;
            confidence_neu_highconf_del(id,d)=NaN;
            conf_corrects_highconf_del(id,d)=NaN;
            conf_rew_corrects_highconf_del(id,d)=NaN;
            conf_neu_corrects_highconf_del(id,d)=NaN;
            conf_wrongs_highconf_del(id,d)=NaN;
            conf_rew_wrongs_highconf_del(id,d)=NaN;
            conf_neu_wrongs_highconf_del(id,d)=NaN;
        else
            % tidy the trials
            clear accu trls oldnew missed confs
            noresp=isnan(dels{id,d}.results.resp);
            accu=dels{id,d}.results.accu;
            trls=cell2mat(dels{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(dels{id,d}.results.trl(:,4))==1;
            confs=dels{id,d}.results.confi;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            confs(noresp)=[];
            
            overallAccuracy_highconf_del(id,d)=nanmean(accu(confs==1));
            hits_highconf_del(id,d)=( sum(oldnew==1 & accu==1 & confs==1) )/ sum(oldnew==1 & confs==1);
            hits_rew_highconf_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==1 & confs==1) )/ sum(oldnew==1 & trls==1 & confs==1);
            hits_neu_highconf_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==0 & confs==1) )/ sum(oldnew==1 & trls==0 & confs==1);
            FAs_highconf_del(id,d)=( sum(oldnew==0 & accu==0 & confs==1) )/ sum(oldnew==0 & confs==1);
            FAs_rew_highconf_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==1 & confs==1) )/ sum(oldnew==0 & trls==1 & confs==1);
            FAs_neu_highconf_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==0 & confs==1) )/ sum(oldnew==0 & trls==0 & confs==1);

            
            % number of trials for d-prime calculation
            n_hit_highconf_del(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            n_fa_highconf_del(id,d)  = sum(oldnew==0 & accu==0 & confs==1);
            n_miss_highconf_del(id,d)= sum(oldnew==1 & accu==0 & confs==1);
            n_cr_highconf_del(id,d)  = sum(oldnew==0 & accu==1 & confs==1);
            
            n_hit_rew_highconf_del(id,d) = sum(oldnew==1 & accu==1 & trls==1 & confs==1);
            n_fa_rew_highconf_del(id,d)  = sum(oldnew==0 & accu==0 & trls==1 & confs==1);
            n_miss_rew_highconf_del(id,d)= sum(oldnew==1 & accu==0 & trls==1 & confs==1);
            n_cr_rew_highconf_del(id,d)  = sum(oldnew==0 & accu==1 & trls==1 & confs==1);
            
            n_hit_neu_highconf_del(id,d) = sum(oldnew==1 & accu==1 & trls==0 & confs==1);
            n_fa_neu_highconf_del(id,d)  = sum(oldnew==0 & accu==0 & trls==0 & confs==1);
            n_miss_neu_highconf_del(id,d)= sum(oldnew==1 & accu==0 & trls==0 & confs==1);
            n_cr_neu_highconf_del(id,d)  = sum(oldnew==0 & accu==1 & trls==0 & confs==1);

            items_remembered_highconf_del(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            items_forgotten_highconf_del(id,d) = sum(oldnew==1 & accu==0 & confs==1);
            items_old_total_highconf_del(id,d) = sum(oldnew==1 & confs==1);
            
        end
    end
end


% immediate recall
overallAccuracy_lowconf_imm=[]; FAs_lowconf_imm=[]; hit_FAs_lowconf_imm=[]; overallAccuracy_lowconf_imm=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            
            overallAccuracy_lowconf_imm(id,d)=NaN; 
            FAs_lowconf_imm(id,d)=NaN; 
%             hit_FAs_lowconf_imm(id,d)=NaN;
            hits_lowconf_imm(id,d)= NaN;
            hits_rew_lowconf_imm(id,d) = NaN;
            hits_neu_lowconf_imm(id,d) = NaN;
            
            FAs_lowconf_imm(id,d) = NaN;
            FAs_rew_lowconf_imm(id,d) = NaN;
            FAs_neu_lowconf_imm(id,d) = NaN;
            
            % number of trials for d-prime calculation
            n_hit_lowconf_imm(id,d) = NaN;
            n_fa_lowconf_imm(id,d)  = NaN;
            n_miss_lowconf_imm(id,d)= NaN;
            n_cr_lowconf_imm(id,d)  = NaN;
            
            n_hit_rew_lowconf_imm(id,d) = NaN;
            n_fa_rew_lowconf_imm(id,d)  = NaN;
            n_miss_rew_lowconf_imm(id,d)= NaN;
            n_cr_rew_lowconf_imm(id,d)  = NaN;
            
            n_hit_neu_lowconf_imm(id,d) = NaN;
            n_fa_neu_lowconf_imm(id,d)  = NaN;
            n_miss_neu_lowconf_imm(id,d)= NaN;
            n_cr_neu_lowconf_imm(id,d)  = NaN;

            items_remembered_lowconf_imm(id,d) = NaN;
            items_forgotten_lowconf_imm(id,d) = NaN;
            items_old_total_lowconf_imm(id,d) = NaN;
            
        else
            
            % tidy the trials
            clear accu trls oldnew missed confs
            noresp=isnan(imms{id,d}.results.resp);
            accu=imms{id,d}.results.accu;
            trls=cell2mat(imms{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(imms{id,d}.results.trl(:,4))==1;
            confs=imms{id,d}.results.confi;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            confs(noresp)=[];
            
            overallAccuracy_lowconf_imm(id,d)=nanmean(accu(confs==1));
            hits_lowconf_imm(id,d)= ( sum(oldnew==1 & accu==1 & confs==1) )/ sum(oldnew==1 & confs==1);
            hits_rew_lowconf_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==1 & confs==1) )/ sum(oldnew==1 & trls==1 & confs==1);
            hits_neu_lowconf_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==0 & confs==1) )/ sum(oldnew==1 & trls==0 & confs==1);
            FAs_lowconf_imm(id,d)=( sum(oldnew==0 & accu==0 & confs==1) )/ sum(oldnew==0 & confs==1);
            FAs_rew_lowconf_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==1 & confs==1) )/ sum(oldnew==0 & trls==1 & confs==1);
            FAs_neu_lowconf_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==0 & confs==1) )/ sum(oldnew==0 & trls==0 & confs==1);
            
            % number of trials for d-prime calculation
            n_hit_lowconf_imm(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            n_fa_lowconf_imm(id,d)  = sum(oldnew==0 & accu==0 & confs==1);
            n_miss_lowconf_imm(id,d)= sum(oldnew==1 & accu==0 & confs==1);
            n_cr_lowconf_imm(id,d)  = sum(oldnew==0 & accu==1 & confs==1);
            
            n_hit_rew_lowconf_imm(id,d) = sum(oldnew==1 & accu==1 & trls==1 & confs==1);
            n_fa_rew_lowconf_imm(id,d)  = sum(oldnew==0 & accu==0 & trls==1 & confs==1);
            n_miss_rew_lowconf_imm(id,d)= sum(oldnew==1 & accu==0 & trls==1 & confs==1);
            n_cr_rew_lowconf_imm(id,d)  = sum(oldnew==0 & accu==1 & trls==1 & confs==1);
            
            n_hit_neu_lowconf_imm(id,d) = sum(oldnew==1 & accu==1 & trls==0 & confs==1);
            n_fa_neu_lowconf_imm(id,d)  = sum(oldnew==0 & accu==0 & trls==0 & confs==1);
            n_miss_neu_lowconf_imm(id,d)= sum(oldnew==1 & accu==0 & trls==0 & confs==1);
            n_cr_neu_lowconf_imm(id,d)  = sum(oldnew==0 & accu==1 & trls==0 & confs==1);  
                
            items_remembered_lowconf_imm(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            items_forgotten_lowconf_imm(id,d) = sum(oldnew==1 & accu==0 & confs==1);
            items_old_total_lowconf_imm(id,d) = sum(oldnew==1 & confs==1);
            
        end
        
    end
end

% delayed recall
hits_lowconf_del=[]; FAs_lowconf_del=[]; hit_FAs_lowconf_del=[]; items_remembered_lowconf_del=[]; items_forgotten_lowconf_del=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            overallAccuracy_lowconf_del(id,d)=NaN;
            hits_lowconf_del(id,d)=NaN;
            hits_rew_lowconf_del(id,d) = NaN;
            hits_neu_lowconf_del(id,d) = NaN;
            FAs_lowconf_del(id,d)=NaN;
            FAs_rew_lowconf_del(id,d) = NaN;
            FAs_neu_lowconf_del(id,d) = NaN;
            
            % number of trials for d-prime calculation
            n_hit_lowconf_del(id,d) = NaN;
            n_fa_lowconf_del(id,d)  = NaN;
            n_miss_lowconf_del(id,d)= NaN;
            n_cr_lowconf_del(id,d)  = NaN;
            
            n_hit_rew_lowconf_del(id,d) = NaN;
            n_fa_rew_lowconf_del(id,d)  = NaN;
            n_miss_rew_lowconf_del(id,d)= NaN;
            n_cr_rew_lowconf_del(id,d)  = NaN;
            
            n_hit_neu_lowconf_del(id,d) = NaN;
            n_fa_neu_lowconf_del(id,d)  = NaN;
            n_miss_neu_lowconf_del(id,d)= NaN;
            n_cr_neu_lowconf_del(id,d)  = NaN;  
            
            items_remembered_lowconf_del(id,d) = NaN;
            items_forgotten_lowconf_del(id,d) = NaN;
            items_old_total_lowconf_del(id,d) = NaN;
            
            confidence_lowconf_del(id,d)=NaN;
            confidence_rew_lowconf_del(id,d)=NaN;
            confidence_neu_lowconf_del(id,d)=NaN;
            conf_corrects_lowconf_del(id,d)=NaN;
            conf_rew_corrects_lowconf_del(id,d)=NaN;
            conf_neu_corrects_lowconf_del(id,d)=NaN;
            conf_wrongs_lowconf_del(id,d)=NaN;
            conf_rew_wrongs_lowconf_del(id,d)=NaN;
            conf_neu_wrongs_lowconf_del(id,d)=NaN;
        else
            % tidy the trials
            clear accu trls oldnew missed confs
            noresp=isnan(dels{id,d}.results.resp);
            accu=dels{id,d}.results.accu;
            trls=cell2mat(dels{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(dels{id,d}.results.trl(:,4))==1;
            confs=dels{id,d}.results.confi;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            confs(noresp)=[];
            
            overallAccuracy_lowconf_del(id,d)=nanmean(accu(confs==1));
            hits_lowconf_del(id,d)=( sum(oldnew==1 & accu==1 & confs==1) )/ sum(oldnew==1 & confs==1);
            hits_rew_lowconf_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==1 & confs==1) )/ sum(oldnew==1 & trls==1 & confs==1);
            hits_neu_lowconf_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==0 & confs==1) )/ sum(oldnew==1 & trls==0 & confs==1);
            FAs_lowconf_del(id,d)=( sum(oldnew==0 & accu==0 & confs==1) )/ sum(oldnew==0 & confs==1);
            FAs_rew_lowconf_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==1 & confs==1) )/ sum(oldnew==0 & trls==1 & confs==1);
            FAs_neu_lowconf_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==0 & confs==1) )/ sum(oldnew==0 & trls==0 & confs==1);

            
            % number of trials for d-prime calculation
            n_hit_lowconf_del(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            n_fa_lowconf_del(id,d)  = sum(oldnew==0 & accu==0 & confs==1);
            n_miss_lowconf_del(id,d)= sum(oldnew==1 & accu==0 & confs==1);
            n_cr_lowconf_del(id,d)  = sum(oldnew==0 & accu==1 & confs==1);
            
            n_hit_rew_lowconf_del(id,d) = sum(oldnew==1 & accu==1 & trls==1 & confs==1);
            n_fa_rew_lowconf_del(id,d)  = sum(oldnew==0 & accu==0 & trls==1 & confs==1);
            n_miss_rew_lowconf_del(id,d)= sum(oldnew==1 & accu==0 & trls==1 & confs==1);
            n_cr_rew_lowconf_del(id,d)  = sum(oldnew==0 & accu==1 & trls==1 & confs==1);
            
            n_hit_neu_lowconf_del(id,d) = sum(oldnew==1 & accu==1 & trls==0 & confs==1);
            n_fa_neu_lowconf_del(id,d)  = sum(oldnew==0 & accu==0 & trls==0 & confs==1);
            n_miss_neu_lowconf_del(id,d)= sum(oldnew==1 & accu==0 & trls==0 & confs==1);
            n_cr_neu_lowconf_del(id,d)  = sum(oldnew==0 & accu==1 & trls==0 & confs==1);

            items_remembered_lowconf_del(id,d) = sum(oldnew==1 & accu==1 & confs==1);
            items_forgotten_lowconf_del(id,d) = sum(oldnew==1 & accu==0 & confs==1);
            items_old_total_lowconf_del(id,d) = sum(oldnew==1 & confs==1);
            
        end
    end
end

%% main task accuracy

for id=1:length(IDs)
    
    for d=1:2
        
        if days(id,d) == 0
            
        else
            clear trials
            trials = eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.maintask.results.trl(:,2))']);
            MT_accuracy_all(id,d)=eval(['nanmean(expdat{id,d}.dat.day' num2str(d) '.maintask.results.accuracy(:))']);
            MT_accuracy_rewards(id,d)=eval(['nanmean(expdat{id,d}.dat.day' num2str(d) '.maintask.results.accuracy(trials==contingency{id,d}(1)))']);
            MT_accuracy_neutrals(id,d)=eval(['nanmean(expdat{id,d}.dat.day' num2str(d) '.maintask.results.accuracy(trials==contingency{id,d}(2)))']);          
        end
    end
end

MT_accuracy_lowDA_Rew  = MT_accuracy_rewards(:,2);
MT_accuracy_highDA_Rew = MT_accuracy_rewards(:,1);
MT_accuracy_highDA_Neu = MT_accuracy_neutrals(:,1);
MT_accuracy_lowDANeu   = MT_accuracy_neutrals(:,2);

%% make an excel table

Hits_highDA_immediate = hits_imm(:,1);
Hits_lowDA_immediate  = hits_imm(:,2);
Hits_highDA_delayed   = hits_del(:,1);
Hits_lowDA_delayed    = hits_del(:,2);

FAs_highDA_immediate = FAs_imm(:,1);
FAs_lowDA_immediate  = FAs_imm(:,2);
FAs_highDA_delayed   = FAs_del(:,1);
FAs_lowDA_delayed    = FAs_del(:,2);

% Dprime_highDA_immediate = dprime_simple(Hits_highDA_immediate, FAs_highDA_immediate);%dprime_loglinear(n_hit_imm(:,1),n_fa_imm(:,1),n_miss_imm(:,1),n_cr_imm(:,1));%dprime_simple(Hits_highDA_immediate,FAs_highDA_immediate);
% Dprime_highDA_delayed   = dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);%dprime_loglinear(n_hit_del(:,1),n_fa_del(:,1),n_miss_del(:,1),n_cr_del(:,1));%dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);
% Dprime_lowDA_immediate  = dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);%dprime_loglinear(n_hit_imm(:,2),n_fa_imm(:,2),n_miss_imm(:,2),n_cr_imm(:,2));%dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);
% Dprime_lowDA_delayed    = dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);%dprime_loglinear(n_hit_del(:,2),n_fa_del(:,2),n_miss_del(:,2),n_cr_del(:,2));%dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);

Dprime_highDA_immediate = dprime_loglinear(n_hit_imm(:,1),n_fa_imm(:,1),n_miss_imm(:,1),n_cr_imm(:,1));%dprime_simple(Hits_highDA_immediate,FAs_highDA_immediate);
Dprime_highDA_delayed   = dprime_loglinear(n_hit_del(:,1),n_fa_del(:,1),n_miss_del(:,1),n_cr_del(:,1));%dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);
Dprime_lowDA_immediate  = dprime_loglinear(n_hit_imm(:,2),n_fa_imm(:,2),n_miss_imm(:,2),n_cr_imm(:,2));%dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);
Dprime_lowDA_delayed    = dprime_loglinear(n_hit_del(:,2),n_fa_del(:,2),n_miss_del(:,2),n_cr_del(:,2));%dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);

Dprime_highDA_highconf_immediate = dprime_loglinear(n_hit_highconf_imm(:,1),n_fa_highconf_imm(:,1),n_miss_highconf_imm(:,1),n_cr_highconf_imm(:,1));%dprime_simple(Hits_highDA_immediate,FAs_highDA_immediate);
Dprime_highDA_highconf_delayed   = dprime_loglinear(n_hit_highconf_del(:,1),n_fa_highconf_del(:,1),n_miss_highconf_del(:,1),n_cr_highconf_del(:,1));%dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);
Dprime_lowDA_highconf_immediate  = dprime_loglinear(n_hit_highconf_imm(:,2),n_fa_highconf_imm(:,2),n_miss_highconf_imm(:,2),n_cr_highconf_imm(:,2));%dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);
Dprime_lowDA_highconf_delayed    = dprime_loglinear(n_hit_highconf_del(:,2),n_fa_highconf_del(:,2),n_miss_highconf_del(:,2),n_cr_highconf_del(:,2));%dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);

Dprime_highDA_lowconf_immediate = dprime_loglinear(n_hit_lowconf_imm(:,1),n_fa_lowconf_imm(:,1),n_miss_lowconf_imm(:,1),n_cr_lowconf_imm(:,1));%dprime_simple(Hits_highDA_immediate,FAs_highDA_immediate);
Dprime_highDA_lowconf_delayed   = dprime_loglinear(n_hit_lowconf_del(:,1),n_fa_lowconf_del(:,1),n_miss_lowconf_del(:,1),n_cr_lowconf_del(:,1));%dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);
Dprime_lowDA_lowconf_immediate  = dprime_loglinear(n_hit_lowconf_imm(:,2),n_fa_lowconf_imm(:,2),n_miss_lowconf_imm(:,2),n_cr_lowconf_imm(:,2));%dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);
Dprime_lowDA_lowconf_delayed    = dprime_loglinear(n_hit_lowconf_del(:,2),n_fa_lowconf_del(:,2),n_miss_lowconf_del(:,2),n_cr_lowconf_del(:,2));%dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);

% Dprime_highDA_immediate = Hits_highDA_immediate-FAs_highDA_immediate;%dprime_loglinear(n_hit_imm(:,1),n_fa_imm(:,1),n_miss_imm(:,1),n_cr_imm(:,1));%dprime_simple(Hits_highDA_immediate,FAs_highDA_immediate);
% Dprime_highDA_delayed   = Hits_highDA_delayed-FAs_highDA_delayed;%dprime_loglinear(n_hit_del(:,1),n_fa_del(:,1),n_miss_del(:,1),n_cr_del(:,1));%dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);
% Dprime_lowDA_immediate  = Hits_lowDA_immediate-FAs_lowDA_immediate;%dprime_loglinear(n_hit_imm(:,2),n_fa_imm(:,2),n_miss_imm(:,2),n_cr_imm(:,2));%dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);
% Dprime_lowDA_delayed    = Hits_lowDA_delayed-FAs_lowDA_delayed;%dprime_loglinear(n_hit_del(:,2),n_fa_del(:,2),n_miss_del(:,2),n_cr_del(:,2));%dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);

Hits_rew_highDA_immediate = hits_rew_imm(:,1);
Hits_rew_lowDA_immediate  = hits_rew_imm(:,2);
Hits_rew_highDA_delayed   = hits_rew_del(:,1);
Hits_rew_lowDA_delayed    = hits_rew_del(:,2);

FAs_rew_highDA_immediate = FAs_rew_imm(:,1);
FAs_rew_lowDA_immediate  = FAs_rew_imm(:,2);
FAs_rew_highDA_delayed   = FAs_rew_del(:,1);
FAs_rew_lowDA_delayed    = FAs_rew_del(:,2);

Hits_neu_highDA_immediate = hits_neu_imm(:,1);
Hits_neu_lowDA_immediate  = hits_neu_imm(:,2);
Hits_neu_highDA_delayed   = hits_neu_del(:,1);
Hits_neu_lowDA_delayed    = hits_neu_del(:,2);

FAs_neu_highDA_immediate = FAs_neu_imm(:,1);
FAs_neu_lowDA_immediate  = FAs_neu_imm(:,2);
FAs_neu_highDA_delayed   = FAs_neu_del(:,1);
FAs_neu_lowDA_delayed    = FAs_neu_del(:,2);

% Dprime_rew_highDA_immediate = dprime_simple(Hits_rew_highDA_immediate,FAs_rew_highDA_immediate);%dprime_loglinear(n_hit_rew_imm(:,1),n_fa_rew_imm(:,1),n_miss_rew_imm(:,1),n_cr_rew_imm(:,1));
% Dprime_rew_highDA_delayed   = dprime_simple(Hits_rew_highDA_delayed,FAs_rew_highDA_delayed);%dprime_loglinear(n_hit_rew_del(:,1),n_fa_rew_del(:,1),n_miss_rew_del(:,1),n_cr_rew_del(:,1));
% Dprime_rew_lowDA_immediate  = dprime_simple(Hits_rew_lowDA_immediate,FAs_rew_lowDA_immediate);%dprime_loglinear(n_hit_rew_imm(:,2),n_fa_rew_imm(:,2),n_miss_rew_imm(:,2),n_cr_rew_imm(:,2));
% Dprime_rew_lowDA_delayed    = dprime_simple(Hits_rew_lowDA_delayed,FAs_rew_lowDA_delayed);%dprime_loglinear(n_hit_rew_del(:,2),n_fa_rew_del(:,2),n_miss_rew_del(:,2),n_cr_rew_del(:,2));
% 
% Dprime_neu_highDA_immediate = dprime_simple(Hits_neu_highDA_immediate,FAs_neu_highDA_immediate);%dprime_loglinear(n_hit_neu_imm(:,1),n_fa_neu_imm(:,1),n_miss_neu_imm(:,1),n_cr_neu_imm(:,1));
% Dprime_neu_highDA_delayed   = dprime_simple(Hits_neu_highDA_delayed,FAs_neu_highDA_delayed);%dprime_loglinear(n_hit_neu_del(:,1),n_fa_neu_del(:,1),n_miss_neu_del(:,1),n_cr_neu_del(:,1));
% Dprime_neu_lowDA_immediate  = dprime_simple(Hits_neu_lowDA_immediate,FAs_neu_lowDA_immediate);%dprime_loglinear(n_hit_neu_imm(:,2),n_fa_neu_imm(:,2),n_miss_neu_imm(:,2),n_cr_neu_imm(:,2));
% Dprime_neu_lowDA_delayed    = dprime_simple(Hits_neu_lowDA_delayed,FAs_neu_lowDA_delayed);%dprime_loglinear(n_hit_neu_del(:,2),n_fa_neu_del(:,2),n_miss_neu_del(:,2),n_cr_neu_del(:,2));

Dprime_rew_highDA_immediate = dprime_loglinear(n_hit_rew_imm(:,1),n_fa_rew_imm(:,1),n_miss_rew_imm(:,1),n_cr_rew_imm(:,1));
Dprime_rew_highDA_delayed   = dprime_loglinear(n_hit_rew_del(:,1),n_fa_rew_del(:,1),n_miss_rew_del(:,1),n_cr_rew_del(:,1));
Dprime_rew_lowDA_immediate  = dprime_loglinear(n_hit_rew_imm(:,2),n_fa_rew_imm(:,2),n_miss_rew_imm(:,2),n_cr_rew_imm(:,2));
Dprime_rew_lowDA_delayed    = dprime_loglinear(n_hit_rew_del(:,2),n_fa_rew_del(:,2),n_miss_rew_del(:,2),n_cr_rew_del(:,2));

Dprime_neu_highDA_immediate = dprime_loglinear(n_hit_neu_imm(:,1),n_fa_neu_imm(:,1),n_miss_neu_imm(:,1),n_cr_neu_imm(:,1));
Dprime_neu_highDA_delayed   = dprime_loglinear(n_hit_neu_del(:,1),n_fa_neu_del(:,1),n_miss_neu_del(:,1),n_cr_neu_del(:,1));
Dprime_neu_lowDA_immediate  = dprime_loglinear(n_hit_neu_imm(:,2),n_fa_neu_imm(:,2),n_miss_neu_imm(:,2),n_cr_neu_imm(:,2));
Dprime_neu_lowDA_delayed    = dprime_loglinear(n_hit_neu_del(:,2),n_fa_neu_del(:,2),n_miss_neu_del(:,2),n_cr_neu_del(:,2));

Dprime_rew_highDA_highconf_immediate = dprime_loglinear(n_hit_rew_highconf_imm(:,1),n_fa_rew_highconf_imm(:,1),n_miss_rew_highconf_imm(:,1),n_cr_rew_highconf_imm(:,1));
Dprime_rew_highDA_highconf_delayed   = dprime_loglinear(n_hit_rew_highconf_del(:,1),n_fa_rew_highconf_del(:,1),n_miss_rew_highconf_del(:,1),n_cr_rew_highconf_del(:,1));
Dprime_rew_lowDA_highconf_immediate  = dprime_loglinear(n_hit_rew_highconf_imm(:,2),n_fa_rew_highconf_imm(:,2),n_miss_rew_highconf_imm(:,2),n_cr_rew_highconf_imm(:,2));
Dprime_rew_lowDA_highconf_delayed    = dprime_loglinear(n_hit_rew_highconf_del(:,2),n_fa_rew_highconf_del(:,2),n_miss_rew_highconf_del(:,2),n_cr_rew_highconf_del(:,2));
 
Dprime_neu_highDA_highconf_immediate = dprime_loglinear(n_hit_neu_highconf_imm(:,1),n_fa_neu_highconf_imm(:,1),n_miss_neu_highconf_imm(:,1),n_cr_neu_highconf_imm(:,1));
Dprime_neu_highDA_highconf_delayed   = dprime_loglinear(n_hit_neu_highconf_del(:,1),n_fa_neu_highconf_del(:,1),n_miss_neu_highconf_del(:,1),n_cr_neu_highconf_del(:,1));
Dprime_neu_lowDA_highconf_immediate  = dprime_loglinear(n_hit_neu_highconf_imm(:,2),n_fa_neu_highconf_imm(:,2),n_miss_neu_highconf_imm(:,2),n_cr_neu_highconf_imm(:,2));
Dprime_neu_lowDA_highconf_delayed    = dprime_loglinear(n_hit_neu_highconf_del(:,2),n_fa_neu_highconf_del(:,2),n_miss_neu_highconf_del(:,2),n_cr_neu_highconf_del(:,2));

Dprime_rew_highDA_lowconf_immediate = dprime_loglinear(n_hit_rew_lowconf_imm(:,1),n_fa_rew_lowconf_imm(:,1),n_miss_rew_lowconf_imm(:,1),n_cr_rew_lowconf_imm(:,1));
Dprime_rew_highDA_lowconf_delayed   = dprime_loglinear(n_hit_rew_lowconf_del(:,1),n_fa_rew_lowconf_del(:,1),n_miss_rew_lowconf_del(:,1),n_cr_rew_lowconf_del(:,1));
Dprime_rew_lowDA_lowconf_immediate  = dprime_loglinear(n_hit_rew_lowconf_imm(:,2),n_fa_rew_lowconf_imm(:,2),n_miss_rew_lowconf_imm(:,2),n_cr_rew_lowconf_imm(:,2));
Dprime_rew_lowDA_lowconf_delayed    = dprime_loglinear(n_hit_rew_lowconf_del(:,2),n_fa_rew_lowconf_del(:,2),n_miss_rew_lowconf_del(:,2),n_cr_rew_lowconf_del(:,2));
 
Dprime_neu_highDA_lowconf_immediate = dprime_loglinear(n_hit_neu_lowconf_imm(:,1),n_fa_neu_lowconf_imm(:,1),n_miss_neu_lowconf_imm(:,1),n_cr_neu_lowconf_imm(:,1));
Dprime_neu_highDA_lowconf_delayed   = dprime_loglinear(n_hit_neu_lowconf_del(:,1),n_fa_neu_lowconf_del(:,1),n_miss_neu_lowconf_del(:,1),n_cr_neu_lowconf_del(:,1));
Dprime_neu_lowDA_lowconf_immediate  = dprime_loglinear(n_hit_neu_lowconf_imm(:,2),n_fa_neu_lowconf_imm(:,2),n_miss_neu_lowconf_imm(:,2),n_cr_neu_lowconf_imm(:,2));
Dprime_neu_lowDA_lowconf_delayed    = dprime_loglinear(n_hit_neu_lowconf_del(:,2),n_fa_neu_lowconf_del(:,2),n_miss_neu_lowconf_del(:,2),n_cr_neu_lowconf_del(:,2));

% Dprime_rew_highDA_immediate = Hits_rew_highDA_immediate-FAs_rew_highDA_immediate;%dprime_loglinear(n_hit_rew_imm(:,1),n_fa_rew_imm(:,1),n_miss_rew_imm(:,1),n_cr_rew_imm(:,1));
% Dprime_rew_highDA_delayed   = Hits_rew_highDA_delayed-FAs_rew_highDA_delayed;%dprime_loglinear(n_hit_rew_del(:,1),n_fa_rew_del(:,1),n_miss_rew_del(:,1),n_cr_rew_del(:,1));
% Dprime_rew_lowDA_immediate  = Hits_rew_lowDA_immediate-FAs_rew_lowDA_immediate;%dprime_loglinear(n_hit_rew_imm(:,2),n_fa_rew_imm(:,2),n_miss_rew_imm(:,2),n_cr_rew_imm(:,2));
% Dprime_rew_lowDA_delayed    = Hits_rew_lowDA_delayed-FAs_rew_lowDA_delayed;%dprime_loglinear(n_hit_rew_del(:,2),n_fa_rew_del(:,2),n_miss_rew_del(:,2),n_cr_rew_del(:,2));
% 
% Dprime_neu_highDA_immediate = Hits_neu_highDA_immediate-FAs_neu_highDA_immediate;%dprime_loglinear(n_hit_neu_imm(:,1),n_fa_neu_imm(:,1),n_miss_neu_imm(:,1),n_cr_neu_imm(:,1));
% Dprime_neu_highDA_delayed   = Hits_neu_highDA_delayed-FAs_neu_highDA_delayed;%dprime_loglinear(n_hit_neu_del(:,1),n_fa_neu_del(:,1),n_miss_neu_del(:,1),n_cr_neu_del(:,1));
% Dprime_neu_lowDA_immediate  = Hits_neu_lowDA_immediate-FAs_neu_lowDA_immediate;%dprime_loglinear(n_hit_neu_imm(:,2),n_fa_neu_imm(:,2),n_miss_neu_imm(:,2),n_cr_neu_imm(:,2));
% Dprime_neu_lowDA_delayed    = Hits_neu_lowDA_delayed-FAs_neu_lowDA_delayed;%dprime_loglinear(n_hit_neu_del(:,2),n_fa_neu_del(:,2),n_miss_neu_del(:,2),n_cr_neu_del(:,2));

confidence_correctTrls_highDA_immdediate = conf_corrects_imm(:,1);
confidence_correctTrls_lowDA_immdediate  = conf_corrects_imm(:,2);
confidence_correctTrls_highDA_delayed    = conf_corrects_del(:,1);
confidence_correctTrls_lowDA_delayed     = conf_corrects_del(:,2);

confidence_correctTrls_rew_highDA_immdediate = conf_rew_corrects_imm(:,1);
confidence_correctTrls_rew_lowDA_immdediate = conf_rew_corrects_imm(:,2);
confidence_correctTrls_rew_highDA_delayed    = conf_rew_corrects_del(:,1);
confidence_correctTrls_rew_lowDA_delayed    = conf_rew_corrects_del(:,2);

confidence_correctTrls_neu_highDA_immdediate = conf_neu_corrects_imm(:,1);
confidence_correctTrls_neu_lowDA_immdediate = conf_neu_corrects_imm(:,2);
confidence_correctTrls_neu_highDA_delayed    = conf_neu_corrects_del(:,1);
confidence_correctTrls_neu_lowDA_delayed    = conf_neu_corrects_del(:,2);

confidence_wrongTrls_highDA_immdediate = conf_wrongs_imm(:,1);
confidence_wrongTrls_lowDA_immdediate  = conf_wrongs_imm(:,2);
confidence_wrongTrls_highDA_delayed    = conf_wrongs_del(:,1);
confidence_wrongTrls_lowDA_delayed     = conf_wrongs_del(:,2);

confidence_wrongTrls_rew_highDA_immdediate = conf_rew_wrongs_imm(:,1);
confidence_wrongTrls_rew_lowDA_immdediate = conf_rew_wrongs_imm(:,2);
confidence_wrongTrls_rew_highDA_delayed    = conf_rew_wrongs_del(:,1);
confidence_wrongTrls_rew_lowDA_delayed    = conf_rew_wrongs_del(:,2);

confidence_wrongTrls_neu_highDA_immdediate = conf_neu_wrongs_imm(:,1);
confidence_wrongTrls_neu_lowDA_immdediate = conf_neu_wrongs_imm(:,2);
confidence_wrongTrls_neu_highDA_delayed    = conf_neu_wrongs_del(:,1);
confidence_wrongTrls_neu_lowDA_delayed    = conf_neu_wrongs_del(:,2);

Misses_highDA_immediate = Miss_all_imm(:,1);
Misses_lowDA_immediate  = Miss_all_imm(:,2);
Misses_highDA_delayed   = Miss_all_del(:,1);
Misses_lowDA_delayed    = Miss_all_del(:,2);

Misses_rew_highDA_immediate = Miss_rew_imm(:,1);
Misses_rew_lowDA_immediate  = Miss_rew_imm(:,2);
Misses_rew_highDA_delayed   = Miss_rew_del(:,1);
Misses_rew_lowDA_delayed    = Miss_rew_del(:,2);

Misses_neu_highDA_immediate = Miss_neu_imm(:,1);
Misses_neu_lowDA_immediate  = Miss_neu_imm(:,2);
Misses_neu_highDA_delayed   = Miss_neu_del(:,1);
Misses_neu_lowDA_delayed    = Miss_neu_del(:,2);

CorrectRejection_highDA_immediate = CorrectRejection_all_imm(:,1);
CorrectRejection_lowDA_immediate  = CorrectRejection_all_imm(:,2);
CorrectRejection_highDA_delayed   = CorrectRejection_all_del(:,1);
CorrectRejection_lowDA_delayed    = CorrectRejection_all_del(:,2);

CorrectRejection_rew_highDA_immediate = CorrectRejection_rew_imm(:,1);
CorrectRejection_rew_lowDA_immediate  = CorrectRejection_rew_imm(:,2);
CorrectRejection_rew_highDA_delayed   = CorrectRejection_rew_del(:,1);
CorrectRejection_rew_lowDA_delayed    = CorrectRejection_rew_del(:,2);

CorrectRejection_neu_highDA_immediate = CorrectRejection_neu_imm(:,1);
CorrectRejection_neu_lowDA_immediate  = CorrectRejection_neu_imm(:,2);
CorrectRejection_neu_highDA_delayed   = CorrectRejection_neu_del(:,1);
CorrectRejection_neu_lowDA_delayed    = CorrectRejection_neu_del(:,2);

RTmean_rew_highDA_immediate = RTmean_immediate_rew(:,1);
RTmean_rew_lowDA_immediate  = RTmean_immediate_rew(:,2);
RTmean_rew_highDA_delayed   = RTmean_delayed_rew(:,1);
RTmean_rew_lowDA_delayed    = RTmean_delayed_rew(:,2);

RTmean_neu_highDA_immediate = RTmean_immediate_neu(:,1);
RTmean_neu_lowDA_immediate  = RTmean_immediate_neu(:,2);
RTmean_neu_highDA_delayed   = RTmean_delayed_neu(:,1);
RTmean_neu_lowDA_delayed    = RTmean_delayed_neu(:,2);

RTmean_FA_rew_highDA_immediate = RTmean_immediate_rew_FA(:,1);
RTmean_HIT_rew_highDA_immediate = RTmean_immediate_rew_hit(:,1);
RTmean_CR_rew_highDA_immediate = RTmean_immediate_rew_CR(:,1);
RTmean_MISS_rew_highDA_immediate = RTmean_immediate_rew_miss(:,1);

RTmean_FA_rew_lowDA_immediate = RTmean_immediate_rew_FA(:,2);
RTmean_HIT_rew_lowDA_immediate = RTmean_immediate_rew_hit(:,2);
RTmean_CR_rew_lowDA_immediate = RTmean_immediate_rew_CR(:,2);
RTmean_MISS_rew_lowDA_immediate = RTmean_immediate_rew_miss(:,2);

RTmean_FA_neu_highDA_immediate = RTmean_immediate_neu_FA(:,1);
RTmean_HIT_neu_highDA_immediate = RTmean_immediate_neu_hit(:,1);
RTmean_CR_neu_highDA_immediate = RTmean_immediate_neu_CR(:,1);
RTmean_MISS_neu_highDA_immediate = RTmean_immediate_neu_miss(:,1);

RTmean_FA_neu_lowDA_immediate = RTmean_immediate_neu_FA(:,2);
RTmean_HIT_neu_lowDA_immediate = RTmean_immediate_neu_hit(:,2);
RTmean_CR_neu_lowDA_immediate = RTmean_immediate_neu_CR(:,2);
RTmean_MISS_neu_lowDA_immediate = RTmean_immediate_neu_miss(:,2);

RTmean_FA_rew_highDA_delayed = RTmean_delayed_rew_FA(:,1);
RTmean_HIT_rew_highDA_delayed = RTmean_delayed_rew_hit(:,1);
RTmean_CR_rew_highDA_delayed = RTmean_delayed_rew_CR(:,1);
RTmean_MISS_rew_highDA_delayed = RTmean_delayed_rew_miss(:,1);
 
RTmean_FA_rew_lowDA_delayed = RTmean_delayed_rew_FA(:,2);
RTmean_HIT_rew_lowDA_delayed = RTmean_delayed_rew_hit(:,2);
RTmean_CR_rew_lowDA_delayed = RTmean_delayed_rew_CR(:,2);
RTmean_MISS_rew_lowDA_delayed = RTmean_delayed_rew_miss(:,2);
 
RTmean_FA_neu_highDA_delayed = RTmean_delayed_neu_FA(:,1);
RTmean_HIT_neu_highDA_delayed = RTmean_delayed_neu_hit(:,1);
RTmean_CR_neu_highDA_delayed = RTmean_delayed_neu_CR(:,1);
RTmean_MISS_neu_highDA_delayed = RTmean_delayed_neu_miss(:,1);
 
RTmean_FA_neu_lowDA_delayed = RTmean_delayed_neu_FA(:,2);
RTmean_HIT_neu_lowDA_delayed = RTmean_delayed_neu_hit(:,2);
RTmean_CR_neu_lowDA_delayed = RTmean_delayed_neu_CR(:,2);
RTmean_MISS_neu_lowDA_delayed = RTmean_delayed_neu_miss(:,2);

ID=IDs';
MRPET_memoryTable = table(...
    ID, ...
    Hits_highDA_immediate, Hits_highDA_delayed, Hits_lowDA_immediate, Hits_lowDA_delayed,...
    FAs_highDA_immediate,  FAs_highDA_delayed,  FAs_lowDA_immediate,  FAs_lowDA_delayed,...
    Hits_rew_highDA_immediate, Hits_rew_highDA_delayed, Hits_rew_lowDA_immediate, Hits_rew_lowDA_delayed,...
    FAs_rew_highDA_immediate,  FAs_rew_highDA_delayed,  FAs_rew_lowDA_immediate,  FAs_rew_lowDA_delayed,...
    Hits_neu_highDA_immediate, Hits_neu_highDA_delayed, Hits_neu_lowDA_immediate, Hits_neu_lowDA_delayed,...
    FAs_neu_highDA_immediate,  FAs_neu_highDA_delayed,  FAs_neu_lowDA_immediate,  FAs_neu_lowDA_delayed,...
    Dprime_highDA_immediate, Dprime_highDA_delayed, Dprime_lowDA_immediate, Dprime_lowDA_delayed,...
    Dprime_rew_highDA_immediate, Dprime_rew_highDA_delayed, Dprime_rew_lowDA_immediate, Dprime_rew_lowDA_delayed,...
    Dprime_neu_highDA_immediate, Dprime_neu_highDA_delayed, Dprime_neu_lowDA_immediate, Dprime_neu_lowDA_delayed,...
    Dprime_rew_highDA_highconf_immediate, Dprime_rew_highDA_highconf_delayed, Dprime_rew_lowDA_highconf_immediate, Dprime_rew_lowDA_highconf_delayed,...
    Dprime_neu_highDA_highconf_immediate, Dprime_neu_highDA_highconf_delayed, Dprime_neu_lowDA_highconf_immediate, Dprime_neu_lowDA_highconf_delayed,...
    Dprime_rew_highDA_lowconf_immediate, Dprime_rew_highDA_lowconf_delayed, Dprime_rew_lowDA_lowconf_immediate, Dprime_rew_lowDA_lowconf_delayed,...
    Dprime_neu_highDA_lowconf_immediate, Dprime_neu_highDA_lowconf_delayed, Dprime_neu_lowDA_lowconf_immediate, Dprime_neu_lowDA_lowconf_delayed,...
    confidence_correctTrls_highDA_immdediate, confidence_correctTrls_highDA_delayed, confidence_correctTrls_lowDA_immdediate, confidence_correctTrls_lowDA_delayed,...
    confidence_correctTrls_rew_highDA_immdediate, confidence_correctTrls_rew_highDA_delayed,...
    confidence_correctTrls_rew_lowDA_immdediate, confidence_correctTrls_rew_lowDA_delayed,...
    confidence_correctTrls_neu_highDA_immdediate, confidence_correctTrls_neu_highDA_delayed,...
    confidence_correctTrls_neu_lowDA_immdediate, confidence_correctTrls_neu_lowDA_delayed,...
    confidence_wrongTrls_highDA_immdediate, confidence_wrongTrls_highDA_delayed, confidence_wrongTrls_lowDA_immdediate, confidence_wrongTrls_lowDA_delayed,...
    confidence_wrongTrls_rew_highDA_immdediate, confidence_wrongTrls_rew_highDA_delayed,...
    confidence_wrongTrls_rew_lowDA_immdediate, confidence_wrongTrls_rew_lowDA_delayed,...
    confidence_wrongTrls_neu_highDA_immdediate, confidence_wrongTrls_neu_highDA_delayed,...
    confidence_wrongTrls_neu_lowDA_immdediate, confidence_wrongTrls_neu_lowDA_delayed,...
    Misses_highDA_immediate, Misses_highDA_delayed, Misses_lowDA_immediate, Misses_lowDA_delayed,...
    Misses_rew_highDA_immediate, Misses_rew_highDA_delayed, Misses_rew_lowDA_immediate, Misses_rew_lowDA_delayed,...
    Misses_neu_highDA_immediate, Misses_neu_highDA_delayed, Misses_neu_lowDA_immediate, Misses_neu_lowDA_delayed,...
    CorrectRejection_highDA_immediate, CorrectRejection_highDA_delayed, CorrectRejection_lowDA_immediate, CorrectRejection_lowDA_delayed,...
    CorrectRejection_rew_highDA_immediate, CorrectRejection_rew_highDA_delayed, CorrectRejection_rew_lowDA_immediate, CorrectRejection_rew_lowDA_delayed,...
    CorrectRejection_neu_highDA_immediate, CorrectRejection_neu_highDA_delayed, CorrectRejection_neu_lowDA_immediate, CorrectRejection_neu_lowDA_delayed,...
    RTmean_rew_highDA_immediate, RTmean_rew_highDA_delayed, RTmean_rew_lowDA_immediate, RTmean_rew_lowDA_delayed,...
    RTmean_neu_highDA_immediate, RTmean_neu_highDA_delayed, RTmean_neu_lowDA_immediate, RTmean_neu_lowDA_delayed,...
    RTmean_FA_rew_highDA_immediate, RTmean_HIT_rew_highDA_immediate,RTmean_CR_rew_highDA_immediate, RTmean_MISS_rew_highDA_immediate,...
    RTmean_FA_rew_lowDA_immediate, RTmean_HIT_rew_lowDA_immediate, RTmean_CR_rew_lowDA_immediate, RTmean_MISS_rew_lowDA_immediate,...
    RTmean_FA_neu_highDA_immediate, RTmean_HIT_neu_highDA_immediate, RTmean_CR_neu_highDA_immediate, RTmean_MISS_neu_highDA_immediate,...
    RTmean_FA_neu_lowDA_immediate, RTmean_HIT_neu_lowDA_immediate, RTmean_CR_neu_lowDA_immediate, RTmean_MISS_neu_lowDA_immediate, ...
    RTmean_FA_rew_highDA_delayed, RTmean_HIT_rew_highDA_delayed, RTmean_CR_rew_highDA_delayed, RTmean_MISS_rew_highDA_delayed, ...
    RTmean_FA_rew_lowDA_delayed, RTmean_HIT_rew_lowDA_delayed, RTmean_CR_rew_lowDA_delayed, RTmean_MISS_rew_lowDA_delayed,...
    RTmean_FA_neu_highDA_delayed, RTmean_HIT_neu_highDA_delayed, RTmean_CR_neu_highDA_delayed, RTmean_MISS_neu_highDA_delayed,...
    RTmean_FA_neu_lowDA_delayed, RTmean_HIT_neu_lowDA_delayed, RTmean_CR_neu_lowDA_delayed, RTmean_MISS_neu_lowDA_delayed)

writetable(MRPET_memoryTable,[ paths.parent 'MRPET_memoryTable_loglinear_' date '.xls'])

%% immediate and delayed, collapsed

nonan_n_hit_rew_imm=n_hit_rew_imm; nonan_n_hit_rew_imm(isnan(n_hit_rew_imm))=0;
nonan_n_hit_rew_del=n_hit_rew_del; nonan_n_hit_rew_del(isnan(n_hit_rew_del))=0;
nonan_n_hit_neu_imm=n_hit_neu_imm; nonan_n_hit_neu_imm(isnan(n_hit_neu_imm))=0;
nonan_n_hit_neu_del=n_hit_neu_del; nonan_n_hit_neu_del(isnan(n_hit_neu_del))=0;

nonan_n_fa_rew_imm=n_fa_rew_imm; nonan_n_fa_rew_imm(isnan(n_fa_rew_imm))=0;
nonan_n_fa_rew_del=n_fa_rew_del; nonan_n_fa_rew_del(isnan(n_fa_rew_del))=0;
nonan_n_fa_neu_imm=n_fa_neu_imm; nonan_n_fa_neu_imm(isnan(n_fa_neu_imm))=0;
nonan_n_fa_neu_del=n_fa_neu_del; nonan_n_fa_neu_del(isnan(n_fa_neu_del))=0;

nonan_n_miss_rew_imm=n_miss_rew_imm; nonan_n_miss_rew_imm(isnan(n_miss_rew_imm))=0;
nonan_n_miss_rew_del=n_miss_rew_del; nonan_n_miss_rew_del(isnan(n_miss_rew_del))=0;
nonan_n_miss_neu_imm=n_miss_neu_imm; nonan_n_miss_neu_imm(isnan(n_miss_neu_imm))=0;
nonan_n_miss_neu_del=n_miss_neu_del; nonan_n_miss_neu_del(isnan(n_miss_neu_del))=0;

nonan_n_cr_rew_imm=n_cr_rew_imm; nonan_n_cr_rew_imm(isnan(n_cr_rew_imm))=0;
nonan_n_cr_rew_del=n_cr_rew_del; nonan_n_cr_rew_del(isnan(n_cr_rew_del))=0;
nonan_n_cr_neu_imm=n_cr_neu_imm; nonan_n_cr_neu_imm(isnan(n_cr_neu_imm))=0;
nonan_n_cr_neu_del=n_cr_neu_del; nonan_n_cr_neu_del(isnan(n_cr_neu_del))=0;

nonan_n_hit_imm=n_hit_imm; nonan_n_hit_imm(isnan(n_hit_imm))=0;
nonan_n_hit_del=n_hit_del; nonan_n_hit_del(isnan(n_hit_del))=0;
nonan_n_hit_imm=n_hit_imm; nonan_n_hit_imm(isnan(n_hit_imm))=0;
nonan_n_hit_del=n_hit_del; nonan_n_hit_del(isnan(n_hit_del))=0;
 
nonan_n_fa_imm=n_fa_imm; nonan_n_fa_imm(isnan(n_fa_imm))=0;
nonan_n_fa_del=n_fa_del; nonan_n_fa_del(isnan(n_fa_del))=0;
nonan_n_fa_imm=n_fa_imm; nonan_n_fa_imm(isnan(n_fa_imm))=0;
nonan_n_fa_del=n_fa_del; nonan_n_fa_del(isnan(n_fa_del))=0;
 
nonan_n_miss_imm=n_miss_imm; nonan_n_miss_imm(isnan(n_miss_imm))=0;
nonan_n_miss_del=n_miss_del; nonan_n_miss_del(isnan(n_miss_del))=0;
nonan_n_miss_imm=n_miss_imm; nonan_n_miss_imm(isnan(n_miss_imm))=0;
nonan_n_miss_del=n_miss_del; nonan_n_miss_del(isnan(n_miss_del))=0;
 
nonan_n_cr_imm=n_cr_imm; nonan_n_cr_imm(isnan(n_cr_imm))=0;
nonan_n_cr_del=n_cr_del; nonan_n_cr_del(isnan(n_cr_del))=0;
nonan_n_cr_imm=n_cr_imm; nonan_n_cr_imm(isnan(n_cr_imm))=0;
nonan_n_cr_del=n_cr_del; nonan_n_cr_del(isnan(n_cr_del))=0;

nonan_n_hit_rew_highconf_imm=n_hit_rew_highconf_imm; nonan_n_hit_rew_highconf_imm(isnan(n_hit_rew_highconf_imm))=0;
nonan_n_hit_rew_highconf_del=n_hit_rew_highconf_del; nonan_n_hit_rew_highconf_del(isnan(n_hit_rew_highconf_del))=0;
nonan_n_hit_neu_highconf_imm=n_hit_neu_highconf_imm; nonan_n_hit_neu_highconf_imm(isnan(n_hit_neu_highconf_imm))=0;
nonan_n_hit_neu_highconf_del=n_hit_neu_highconf_del; nonan_n_hit_neu_highconf_del(isnan(n_hit_neu_highconf_del))=0;
 
nonan_n_fa_rew_highconf_imm=n_fa_rew_highconf_imm; nonan_n_fa_rew_highconf_imm(isnan(n_fa_rew_highconf_imm))=0;
nonan_n_fa_rew_highconf_del=n_fa_rew_highconf_del; nonan_n_fa_rew_highconf_del(isnan(n_fa_rew_highconf_del))=0;
nonan_n_fa_neu_highconf_imm=n_fa_neu_highconf_imm; nonan_n_fa_neu_highconf_imm(isnan(n_fa_neu_highconf_imm))=0;
nonan_n_fa_neu_highconf_del=n_fa_neu_highconf_del; nonan_n_fa_neu_highconf_del(isnan(n_fa_neu_highconf_del))=0;
 
nonan_n_miss_rew_highconf_imm=n_miss_rew_highconf_imm; nonan_n_miss_rew_highconf_imm(isnan(n_miss_rew_highconf_imm))=0;
nonan_n_miss_rew_highconf_del=n_miss_rew_highconf_del; nonan_n_miss_rew_highconf_del(isnan(n_miss_rew_highconf_del))=0;
nonan_n_miss_neu_highconf_imm=n_miss_neu_highconf_imm; nonan_n_miss_neu_highconf_imm(isnan(n_miss_neu_highconf_imm))=0;
nonan_n_miss_neu_highconf_del=n_miss_neu_highconf_del; nonan_n_miss_neu_highconf_del(isnan(n_miss_neu_highconf_del))=0;
 
nonan_n_cr_rew_highconf_imm=n_cr_rew_highconf_imm; nonan_n_cr_rew_highconf_imm(isnan(n_cr_rew_highconf_imm))=0;
nonan_n_cr_rew_highconf_del=n_cr_rew_highconf_del; nonan_n_cr_rew_highconf_del(isnan(n_cr_rew_highconf_del))=0;
nonan_n_cr_neu_highconf_imm=n_cr_neu_highconf_imm; nonan_n_cr_neu_highconf_imm(isnan(n_cr_neu_highconf_imm))=0;
nonan_n_cr_neu_highconf_del=n_cr_neu_highconf_del; nonan_n_cr_neu_highconf_del(isnan(n_cr_neu_highconf_del))=0;

nonan_n_hit_highconf_imm=n_hit_highconf_imm; nonan_n_hit_highconf_imm(isnan(n_hit_highconf_imm))=0;
nonan_n_fa_highconf_imm=n_fa_highconf_imm; nonan_n_fa_highconf_imm(isnan(n_fa_highconf_imm))=0;
nonan_n_miss_highconf_imm=n_miss_highconf_imm; nonan_n_miss_highconf_imm(isnan(n_miss_highconf_imm))=0;
nonan_n_cr_highconf_imm=n_cr_highconf_imm; nonan_n_cr_highconf_imm(isnan(n_cr_highconf_imm))=0;

nonan_n_hit_highconf_del=n_hit_highconf_del; nonan_n_hit_highconf_del(isnan(n_hit_highconf_del))=0;
nonan_n_fa_highconf_del=n_fa_highconf_del; nonan_n_fa_highconf_del(isnan(n_fa_highconf_del))=0;
nonan_n_miss_highconf_del=n_miss_highconf_del; nonan_n_miss_highconf_del(isnan(n_miss_highconf_del))=0;
nonan_n_cr_highconf_del=n_cr_highconf_del; nonan_n_cr_highconf_del(isnan(n_cr_highconf_del))=0;


% calc

Dprime_highDA_immediate = dprime_loglinear(n_hit_imm(:,1),n_fa_imm(:,1),n_miss_imm(:,1),n_cr_imm(:,1));%dprime_simple(Hits_highDA_immediate,FAs_highDA_immediate);
Dprime_highDA_delayed   = dprime_loglinear(n_hit_del(:,1),n_fa_del(:,1),n_miss_del(:,1),n_cr_del(:,1));%dprime_simple(Hits_highDA_delayed,FAs_highDA_delayed);
Dprime_lowDA_immediate  = dprime_loglinear(n_hit_imm(:,2),n_fa_imm(:,2),n_miss_imm(:,2),n_cr_imm(:,2));%dprime_simple(Hits_lowDA_immediate,FAs_lowDA_immediate);
Dprime_lowDA_delayed    = dprime_loglinear(n_hit_del(:,2),n_fa_del(:,2),n_miss_del(:,2),n_cr_del(:,2));%dprime_simple(Hits_lowDA_delayed,FAs_lowDA_delayed);

Dprime_highDA_all       = dprime_loglinear(nonan_n_hit_highconf_imm(:,1)+nonan_n_hit_highconf_del(:,1),...
    nonan_n_fa_highconf_imm(:,1)+nonan_n_fa_highconf_del(:,1),...
    nonan_n_miss_highconf_imm(:,1)+nonan_n_miss_highconf_del(:,1),...
    nonan_n_cr_highconf_imm(:,1)+nonan_n_cr_highconf_del(:,1));
Dprime_lowDA_all       = dprime_loglinear(nonan_n_hit_highconf_imm(:,2)+nonan_n_hit_highconf_del(:,2),...
    nonan_n_fa_highconf_imm(:,2)+nonan_n_fa_highconf_del(:,2),...
    nonan_n_miss_highconf_imm(:,2)+nonan_n_miss_highconf_del(:,2),...
    nonan_n_cr_highconf_imm(:,2)+nonan_n_cr_highconf_del(:,2));

Dprime_highDA_highconf_all       = dprime_loglinear(nonan_n_hit_imm(:,1)+nonan_n_hit_del(:,1),...
    nonan_n_fa_imm(:,1)+nonan_n_fa_del(:,1),...
    nonan_n_miss_imm(:,1)+nonan_n_miss_del(:,1),...
    nonan_n_cr_imm(:,1)+nonan_n_cr_del(:,1));
Dprime_lowDA_highconf_all       = dprime_loglinear(nonan_n_hit_imm(:,2)+nonan_n_hit_del(:,2),...
    nonan_n_fa_imm(:,2)+nonan_n_fa_del(:,2),...
    nonan_n_miss_imm(:,2)+nonan_n_miss_del(:,2),...
    nonan_n_cr_imm(:,2)+nonan_n_cr_del(:,2));

Dprime_rew_highDA_immediate = dprime_loglinear(n_hit_rew_imm(:,1),n_fa_rew_imm(:,1),n_miss_rew_imm(:,1),n_cr_rew_imm(:,1));
Dprime_rew_highDA_delayed   = dprime_loglinear(n_hit_rew_del(:,1),n_fa_rew_del(:,1),n_miss_rew_del(:,1),n_cr_rew_del(:,1));
Dprime_rew_lowDA_immediate  = dprime_loglinear(n_hit_rew_imm(:,2),n_fa_rew_imm(:,2),n_miss_rew_imm(:,2),n_cr_rew_imm(:,2));
Dprime_rew_lowDA_delayed    = dprime_loglinear(n_hit_rew_del(:,2),n_fa_rew_del(:,2),n_miss_rew_del(:,2),n_cr_rew_del(:,2));

Dprime_rew_highDA_all       = dprime_loglinear(nonan_n_hit_rew_imm(:,1)+nonan_n_hit_rew_del(:,1),...
    nonan_n_fa_rew_imm(:,1)+nonan_n_fa_rew_del(:,1),...
    nonan_n_miss_rew_imm(:,1)+nonan_n_miss_rew_del(:,1),...
    nonan_n_cr_rew_imm(:,1)+nonan_n_cr_rew_del(:,1));
Dprime_rew_lowDA_all       = dprime_loglinear(nonan_n_hit_rew_imm(:,2)+nonan_n_hit_rew_del(:,2),...
    nonan_n_fa_rew_imm(:,2)+nonan_n_fa_rew_del(:,2),...
    nonan_n_miss_rew_imm(:,2)+nonan_n_miss_rew_del(:,2),...
    nonan_n_cr_rew_imm(:,2)+nonan_n_cr_rew_del(:,2));

Dprime_neu_highDA_immediate = dprime_loglinear(n_hit_neu_imm(:,1),n_fa_neu_imm(:,1),n_miss_neu_imm(:,1),n_cr_neu_imm(:,1));
Dprime_neu_highDA_delayed   = dprime_loglinear(n_hit_neu_del(:,1),n_fa_neu_del(:,1),n_miss_neu_del(:,1),n_cr_neu_del(:,1));
Dprime_neu_lowDA_immediate  = dprime_loglinear(n_hit_neu_imm(:,2),n_fa_neu_imm(:,2),n_miss_neu_imm(:,2),n_cr_neu_imm(:,2));
Dprime_neu_lowDA_delayed    = dprime_loglinear(n_hit_neu_del(:,2),n_fa_neu_del(:,2),n_miss_neu_del(:,2),n_cr_neu_del(:,2));

Dprime_neu_highDA_all       = dprime_loglinear(nonan_n_hit_neu_imm(:,1)+nonan_n_hit_neu_del(:,1),...
    nonan_n_fa_neu_imm(:,1)+nonan_n_fa_neu_del(:,1),...
    nonan_n_miss_neu_imm(:,1)+nonan_n_miss_neu_del(:,1),...
    nonan_n_cr_neu_imm(:,1)+nonan_n_cr_neu_del(:,1));
Dprime_neu_lowDA_all       = dprime_loglinear(nonan_n_hit_neu_imm(:,2)+nonan_n_hit_neu_del(:,2),...
    nonan_n_fa_neu_imm(:,2)+nonan_n_fa_neu_del(:,2),...
    nonan_n_miss_neu_imm(:,2)+nonan_n_miss_neu_del(:,2),...
    nonan_n_cr_neu_imm(:,2)+nonan_n_cr_neu_del(:,2));

Dprime_rew_highDA_highconf_immediate = dprime_loglinear(n_hit_rew_highconf_imm(:,1),n_fa_rew_highconf_imm(:,1),n_miss_rew_highconf_imm(:,1),n_cr_rew_highconf_imm(:,1));
Dprime_rew_highDA_highconf_delayed   = dprime_loglinear(n_hit_rew_highconf_del(:,1),n_fa_rew_highconf_del(:,1),n_miss_rew_highconf_del(:,1),n_cr_rew_highconf_del(:,1));
Dprime_rew_lowDA_highconf_immediate  = dprime_loglinear(n_hit_rew_highconf_imm(:,2),n_fa_rew_highconf_imm(:,2),n_miss_rew_highconf_imm(:,2),n_cr_rew_highconf_imm(:,2));
Dprime_rew_lowDA_highconf_delayed    = dprime_loglinear(n_hit_rew_highconf_del(:,2),n_fa_rew_highconf_del(:,2),n_miss_rew_highconf_del(:,2),n_cr_rew_highconf_del(:,2));
 
Dprime_rew_highDA_highconf_all       = dprime_loglinear(nonan_n_hit_rew_highconf_imm(:,1)+nonan_n_hit_rew_highconf_del(:,1),...
    nonan_n_fa_rew_highconf_imm(:,1)+nonan_n_fa_rew_highconf_del(:,1),...
    nonan_n_miss_rew_highconf_imm(:,1)+nonan_n_miss_rew_highconf_del(:,1),...
    nonan_n_cr_rew_highconf_imm(:,1)+nonan_n_cr_rew_highconf_del(:,1));
Dprime_rew_lowDA_highconf_all       = dprime_loglinear(nonan_n_hit_rew_highconf_imm(:,2)+nonan_n_hit_rew_highconf_del(:,2),...
    nonan_n_fa_rew_highconf_imm(:,2)+nonan_n_fa_rew_highconf_del(:,2),...
    nonan_n_miss_rew_highconf_imm(:,2)+nonan_n_miss_rew_highconf_del(:,2),...
    nonan_n_cr_rew_highconf_imm(:,2)+nonan_n_cr_rew_highconf_del(:,2));

Dprime_neu_highDA_highconf_immediate = dprime_loglinear(n_hit_neu_highconf_imm(:,1),n_fa_neu_highconf_imm(:,1),n_miss_neu_highconf_imm(:,1),n_cr_neu_highconf_imm(:,1));
Dprime_neu_highDA_highconf_delayed   = dprime_loglinear(n_hit_neu_highconf_del(:,1),n_fa_neu_highconf_del(:,1),n_miss_neu_highconf_del(:,1),n_cr_neu_highconf_del(:,1));
Dprime_neu_lowDA_highconf_immediate  = dprime_loglinear(n_hit_neu_highconf_imm(:,2),n_fa_neu_highconf_imm(:,2),n_miss_neu_highconf_imm(:,2),n_cr_neu_highconf_imm(:,2));
Dprime_neu_lowDA_highconf_delayed    = dprime_loglinear(n_hit_neu_highconf_del(:,2),n_fa_neu_highconf_del(:,2),n_miss_neu_highconf_del(:,2),n_cr_neu_highconf_del(:,2));

Dprime_neu_highDA_highconf_all       = dprime_loglinear(nonan_n_hit_neu_highconf_imm(:,1)+nonan_n_hit_neu_highconf_del(:,1),...
    nonan_n_fa_neu_highconf_imm(:,1)+nonan_n_fa_neu_highconf_del(:,1),...
    nonan_n_miss_neu_highconf_imm(:,1)+nonan_n_miss_neu_highconf_del(:,1),...
    nonan_n_cr_neu_highconf_imm(:,1)+nonan_n_cr_neu_highconf_del(:,1));
Dprime_neu_lowDA_highconf_all       = dprime_loglinear(nonan_n_hit_neu_highconf_imm(:,2)+nonan_n_hit_neu_highconf_del(:,2),...
    nonan_n_fa_neu_highconf_imm(:,2)+nonan_n_fa_neu_highconf_del(:,2),...
    nonan_n_miss_neu_highconf_imm(:,2)+nonan_n_miss_neu_highconf_del(:,2),...
    nonan_n_cr_neu_highconf_imm(:,2)+nonan_n_cr_neu_highconf_del(:,2));

Dprime_rew_highDA_lowconf_immediate = dprime_loglinear(n_hit_rew_lowconf_imm(:,1),n_fa_rew_lowconf_imm(:,1),n_miss_rew_lowconf_imm(:,1),n_cr_rew_lowconf_imm(:,1));
Dprime_rew_highDA_lowconf_delayed   = dprime_loglinear(n_hit_rew_lowconf_del(:,1),n_fa_rew_lowconf_del(:,1),n_miss_rew_lowconf_del(:,1),n_cr_rew_lowconf_del(:,1));
Dprime_rew_lowDA_lowconf_immediate  = dprime_loglinear(n_hit_rew_lowconf_imm(:,2),n_fa_rew_lowconf_imm(:,2),n_miss_rew_lowconf_imm(:,2),n_cr_rew_lowconf_imm(:,2));
Dprime_rew_lowDA_lowconf_delayed    = dprime_loglinear(n_hit_rew_lowconf_del(:,2),n_fa_rew_lowconf_del(:,2),n_miss_rew_lowconf_del(:,2),n_cr_rew_lowconf_del(:,2));
 
Dprime_neu_highDA_lowconf_immediate = dprime_loglinear(n_hit_neu_lowconf_imm(:,1),n_fa_neu_lowconf_imm(:,1),n_miss_neu_lowconf_imm(:,1),n_cr_neu_lowconf_imm(:,1));
Dprime_neu_highDA_lowconf_delayed   = dprime_loglinear(n_hit_neu_lowconf_del(:,1),n_fa_neu_lowconf_del(:,1),n_miss_neu_lowconf_del(:,1),n_cr_neu_lowconf_del(:,1));
Dprime_neu_lowDA_lowconf_immediate  = dprime_loglinear(n_hit_neu_lowconf_imm(:,2),n_fa_neu_lowconf_imm(:,2),n_miss_neu_lowconf_imm(:,2),n_cr_neu_lowconf_imm(:,2));
Dprime_neu_lowDA_lowconf_delayed    = dprime_loglinear(n_hit_neu_lowconf_del(:,2),n_fa_neu_lowconf_del(:,2),n_miss_neu_lowconf_del(:,2),n_cr_neu_lowconf_del(:,2));

MRPET_memoryTable=[];
ID=IDs';
MRPET_memoryTable = table(...
    ID, ...
    Dprime_highDA_all,Dprime_lowDA_all, Dprime_highDA_highconf_all, Dprime_lowDA_highconf_all,...
    Dprime_rew_highDA_all, Dprime_neu_highDA_all, Dprime_rew_lowDA_all, Dprime_neu_lowDA_all,...
    Dprime_rew_highDA_highconf_all, Dprime_neu_highDA_highconf_all, Dprime_rew_lowDA_highconf_all, Dprime_neu_lowDA_highconf_all,...
    Dprime_rew_highDA_immediate, Dprime_rew_highDA_delayed, Dprime_rew_lowDA_immediate, Dprime_rew_lowDA_delayed,...
    Dprime_neu_highDA_immediate, Dprime_neu_highDA_delayed, Dprime_neu_lowDA_immediate, Dprime_neu_lowDA_delayed,...
    Dprime_rew_highDA_highconf_immediate, Dprime_rew_highDA_highconf_delayed, Dprime_rew_lowDA_highconf_immediate, Dprime_rew_lowDA_highconf_delayed,...
    Dprime_neu_highDA_highconf_immediate, Dprime_neu_highDA_highconf_delayed, Dprime_neu_lowDA_highconf_immediate, Dprime_neu_lowDA_highconf_delayed,...
    Dprime_rew_highDA_lowconf_immediate, Dprime_rew_highDA_lowconf_delayed, Dprime_rew_lowDA_lowconf_immediate, Dprime_rew_lowDA_lowconf_delayed,...
    Dprime_neu_highDA_lowconf_immediate, Dprime_neu_highDA_lowconf_delayed, Dprime_neu_lowDA_lowconf_immediate, Dprime_neu_lowDA_lowconf_delayed);

writetable(MRPET_memoryTable,[ paths.parent 'MRPET_memoryTable_collapsed_loglinear_' date '.xls'])

%% check comparability of the measures


% sort stim categories: 1 = nature , 2 = urban, 3 = private, 4 = public
clear stims_imm stims_del missed_imm missed_del accu_imm accu_del
c1=0;c2=0;c3=0;c4=0;
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0'])
            disp('no delay memory test data')
            stims_del{id,d}=NaN;
            missed_del{id,d}=NaN;
            accu_del{id,d}=NaN;
            oldnew_del{id,d}=NaN;
            
            stims_imm{id,d}=eval(['cell2mat(behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl(:,2))']);
            missed_imm{id,d}=isnan(eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.resp']));
            accu_imm{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.accu']);
            accu_imm{id,d}(missed_imm{id,d})=[];
            oldnew_imm{id,d}=eval(['cell2mat(behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl(:,4))==1']);
            oldnew_imm{id,d}(missed_imm{id,d})=[];
            stims_imm{id,d}(missed_imm{id,d})=[];
            
            eval(['c' num2str(contingency{id,d}(1)) '=c' num2str(contingency{id,d}(1)) '+1;'])
            eval(['c' num2str(contingency{id,d}(2)) '=c' num2str(contingency{id,d}(2)) '+1;'])
            
            eval(['hits' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)= ( sum(oldnew_imm{id,d}==1 & accu_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ') )/ sum(oldnew_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ')']);
            eval(['hits' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',1)= ( sum(oldnew_imm{id,d}==1 & accu_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ') )/ sum(oldnew_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ')']);
                        
            eval(['FAs' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)= ( sum(oldnew_imm{id,d}==0 & accu_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ') )/ sum(oldnew_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ')']);
            eval(['FAs' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',1)= ( sum(oldnew_imm{id,d}==0 & accu_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ') )/ sum(oldnew_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ')']);
            
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)= hits' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1) - FAs' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)']);
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',2)= IDs(id)']);
            
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',1)= hits' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(1)) ',1) - FAs' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(1)) ',1)']);
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',2)= IDs(id)']);
               
            
            eval(['hits' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)=NaN;']);
            eval(['FAs' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)=NaN;']);
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)=NaN;']);
            eval(['hits' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',1)=NaN;']);
            eval(['FAs' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',1)=NaN;']);
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',1)=NaN;']);
        
        elseif days(id,d)==0
            
            disp('no delay memory test data')
            stims_del{id,d}=NaN;
            missed_del{id,d}=NaN;
            accu_del{id,d}=NaN;
            oldnew_del{id,d}=NaN;
            
            stims_imm{id,d}=NaN;
            missed_imm{id,d}=NaN;
            accu_imm{id,d}=NaN;
            oldnew_imm{id,d}=NaN;
        
        else
            
            % tidy the trials
            
            stims_imm{id,d}=eval(['cell2mat(behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl(:,2))']);
            stims_del{id,d}=eval(['cell2mat(behdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.trl(:,2))']);
            
            missed_imm{id,d}=isnan(eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.resp']));
            missed_del{id,d}=isnan(eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.resp']));
            
            accu_imm{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.accu']);
            accu_del{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.accu']);
            
            accu_imm{id,d}(missed_imm{id,d})=[];
            accu_del{id,d}(missed_del{id,d})=[];
            
            oldnew_imm{id,d}=eval(['cell2mat(behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl(:,4))==1']);
            oldnew_del{id,d}=eval(['cell2mat(behdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.trl(:,4))==1']);
            
            oldnew_imm{id,d}(missed_imm{id,d})=[];
            oldnew_del{id,d}(missed_del{id,d})=[];
            
            stims_imm{id,d}(missed_imm{id,d})=[];
            stims_del{id,d}(missed_del{id,d})=[];
            
            eval(['c' num2str(contingency{id,d}(1)) '=c' num2str(contingency{id,d}(1)) '+1;'])
            eval(['c' num2str(contingency{id,d}(2)) '=c' num2str(contingency{id,d}(2)) '+1;'])
            
            eval(['hits' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)= ( sum(oldnew_imm{id,d}==1 & accu_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ') )/ sum(oldnew_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ')']);
            eval(['hits' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',1)= ( sum(oldnew_imm{id,d}==1 & accu_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ') )/ sum(oldnew_imm{id,d}==1 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ')']);
            
            eval(['hits' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)= ( sum(oldnew_del{id,d}==1 & accu_del{id,d}==1 & stims_del{id,d}==' num2str(contingency{id,d}(1)) ') )/ sum(oldnew_del{id,d}==1 & stims_del{id,d}==' num2str(contingency{id,d}(1)) ')']);
            eval(['hits' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',1)= ( sum(oldnew_del{id,d}==1 & accu_del{id,d}==1 & stims_del{id,d}==' num2str(contingency{id,d}(2)) ') )/ sum(oldnew_del{id,d}==1 & stims_del{id,d}==' num2str(contingency{id,d}(2)) ')']);
            
            eval(['FAs' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)= ( sum(oldnew_imm{id,d}==0 & accu_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ') )/ sum(oldnew_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(1)) ')']);
            eval(['FAs' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',1)= ( sum(oldnew_imm{id,d}==0 & accu_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ') )/ sum(oldnew_imm{id,d}==0 & stims_imm{id,d}==' num2str(contingency{id,d}(2)) ')']);
            
            eval(['FAs' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)= ( sum(oldnew_del{id,d}==0 & accu_del{id,d}==0 & stims_del{id,d}==' num2str(contingency{id,d}(1)) ') )/ sum(oldnew_del{id,d}==0 & stims_del{id,d}==' num2str(contingency{id,d}(1)) ')']);
            eval(['FAs' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',1)= ( sum(oldnew_del{id,d}==0 & accu_del{id,d}==0 & stims_del{id,d}==' num2str(contingency{id,d}(2)) ') )/ sum(oldnew_del{id,d}==0 & stims_del{id,d}==' num2str(contingency{id,d}(2)) ')']);
            
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)= hits' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1) - FAs' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',1)']);
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_imm(c' num2str(contingency{id,d}(1)) ',2)= IDs(id)']);
            
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',1)= hits' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(1)) ',1) - FAs' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(1)) ',1)']);
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_imm(c' num2str(contingency{id,d}(2)) ',2)= IDs(id)']);
            
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)= hits' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1) - FAs' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',1)']);
            eval(['Hit_FA' num2str(contingency{id,d}(1)) '_del(c' num2str(contingency{id,d}(1)) ',2)= IDs(id)']);
            
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',1)= hits' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(1)) ',1) - FAs' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(1)) ',1)']);
            eval(['Hit_FA' num2str(contingency{id,d}(2)) '_del(c' num2str(contingency{id,d}(2)) ',2)= IDs(id)']);
            
            
        end
        
    end
    
end


% levene's f test: check homogeneity
indoorDs_imm=[Hit_FA1_imm(:,1) ones(length(Hit_FA1_imm(:,1)),1); Hit_FA2_imm(:,1) ones(length(Hit_FA2_imm(:,1)),1)];
outdoorDs_imm=[Hit_FA3_imm(:,1) 2*ones(length(Hit_FA3_imm(:,1)),1); Hit_FA4_imm(:,1) 2*ones(length(Hit_FA4_imm(:,1)),1)];

indoorDs_del=[Hit_FA1_del(:,1) ones(length(Hit_FA1_del(:,1)),1); Hit_FA2_del(:,1) ones(length(Hit_FA2_del(:,1)),1)];
outdoorDs_del=[Hit_FA3_del(:,1) 2*ones(length(Hit_FA3_del(:,1)),1); Hit_FA4_del(:,1) 2*ones(length(Hit_FA4_del(:,1)),1)];

Levenetest([indoorDs_imm;outdoorDs_imm])

%% extract T values

for id=1:length(IDs)
    
    
    stim_rew_neu = spm_read_vols(spm_vol([paths.taskonly num2str(IDs(id)) '/con_0002_mni.nii']));
    fb_rew_neu   = spm_read_vols(spm_vol([paths.taskonly num2str(IDs(id)) '/con_0008_mni.nii']));
    
    
    Tval_LC_stim_rew_neu(id,1)    = nanmean(stim_rew_neu(LCmask==1));
    Tval_SNVTA_stim_rew_neu(id,1) = nanmean(stim_rew_neu(SNVTAmask==1));
    
    Tval_LC_fb_rew_neu(id,1)    = nanmean(fb_rew_neu(LCmask==1));
    Tval_SNVTA_fb_rew_neu(id,1) = nanmean(fb_rew_neu(SNVTAmask==1));
    
    clear stim_rew_neu fb_rew_neu
    
    
end


%% figure

clear data dat
data = readtable('/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_memoryTable_DprimeLogLinear_20-Mar-2023.xls');

%%
% colors = cbrewer2('qual', 'Set1', 10);
colors = [1 0.1 0; 1 0.4 0; 0 0.1 1; 0 0.4 1];
dat{1,1}(:, 1) = data{:,30}; % high rew imm
dat{2,1}(:, 1) = data{:,34}; % high neu imm
dat{3,1}(:, 1) = data{:,32}; % low rew imm
dat{4,1}(:, 1) = data{:,36}; % low neu imm

close all

subplot(1,2,1); % rather than a square plot, make it thinner
hold on;
% subplot(1,2,1)
% if we want each bar to have a different color, loop
for b = 1:size(dat, 1),
    bar(b, nanmean(dat{b,1}(:)), 'FaceColor',  colors(b, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
end
 
stderror = @ (x) nanstd(x)/sqrt(length(x));

% show standard deviation on top
h = ploterr(1:4, cell2mat(cellfun(@nanmean,dat,'UniformOutput',false)), [], cell2mat(cellfun(@ (x) nanstd(x)/sqrt(length(x)),dat,'UniformOutput',false)), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none'); % remove marker
 
% label what we're seeing
% if labels are too long to fit, use the xticklabelrotation with about -30
% to rotate them so they're readable
set(gca, 'xtick', [1 2 3 4], 'xticklabel', {'highDA-Rew-Imm','highDA-Neu-Imm','lowDA-Rew-Imm','lowDA-Neu-Imm'}, ...
    'xlim', [0 5], 'ylim', [0 1.5]);
ylabel('D-prime'); xlabel('Categories');
 
% if these data are paired, show the differences
% plot(dat', '.k-', 'linewidth', 0.2, 'markersize', 2);
 
% significance star for the difference
[~, pval] = ttest2(dat{1,1}, dat{2,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [1 2], 1.4, pval);
end

[~, pval] = ttest2(dat{3,1}, dat{4,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [3 4], 1.4, pval);
end

[~, pval] = ttest2(dat{1,1}, dat{3,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [1 3], 1.5, pval);
end

[~, pval] = ttest2(dat{2,1}, dat{4,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [2 4], 1.6, pval);
end

[~, pval] = ttest2(dat{1,1}, dat{4,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [1 4], 1.7, pval);
end

[~, pval] = ttest2(dat{2,1}, dat{3,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [2 3], 1.8, pval);
end


% add significance stars for each bar
% for b = 1:3,
%     [~, pval] = ttest(dat(:, b));
%     yval = mean(dat(:, b)) * 0.5; % plot this on top of the bar
%     mysigstar(gca, b, yval, pval);
%     % if mysigstar gets just 1 xpos input, it will only plot stars
% end


% manage grids

% grid on
% set(gca,'color',[230/255 230/255 230/255])
% box off
% grid on
% set(gca,'LineWidth',4) 
% ax = gca;
% ax.GridColor = [1 1 1];
% ax.GridAlpha = 1;


set(gca,'LineWidth',3)



%%

colors = [1 0.1 0; 1 0.4 0; 0 0.1 1; 0 0.4 1];
dat{1,1}(:, 1) = data{:,31}; % high rew del
dat{2,1}(:, 1) = data{:,35}; % high neu del
dat{3,1}(:, 1) = data{:,33}; % low rew del
dat{4,1}(:, 1) = data{:,37}; % low neu del

% subplot(4,7,6); % rather than a square plot, make it thinner
% subplot(1,2,2)
% if we want each bar to have a different color, loop
for b = 1:size(dat, 1),
    bar(b, nanmean(dat{b,1}(:)), 'FaceColor',  colors(b, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
end
 
stderror = @ (x) nanstd(x)/sqrt(length(x));

% show standard deviation on top
h = ploterr(1:4, cell2mat(cellfun(@nanmean,dat,'UniformOutput',false)), [], cell2mat(cellfun(@ (x) nanstd(x)/sqrt(length(x)),dat,'UniformOutput',false)), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none'); % remove marker
 
% label what we're seeing
% if labels are too long to fit, use the xticklabelrotation with about -30
% to rotate them so they're readable
set(gca, 'xtick', [1 2 3 4], 'xticklabel', {'highDA-Rew-Del','highDA-Neu-Del','lowDA-Rew-Del','lowDA-Neu-Del'}, ...
    'xlim', [0 5], 'ylim', [0 1.5]);
ylabel('D-prime'); xlabel('Categories');
 
% if these data are paired, show the differences
% plot(dat', '.k-', 'linewidth', 0.2, 'markersize', 2);
 
% significance star for the difference
[~, pval] = ttest2(dat{1,1}, dat{2,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [1 2], 1.4, pval);
end

[~, pval] = ttest2(dat{3,1}, dat{4,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [3 4], 1.4, pval);
end

[~, pval] = ttest2(dat{1,1}, dat{3,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [1 3], 1.5, pval);
end

[~, pval] = ttest2(dat{2,1}, dat{4,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [2 4], 1.6, pval);
end

[~, pval] = ttest2(dat{1,1}, dat{4,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [1 4], 1.7, pval);
end

[~, pval] = ttest2(dat{2,1}, dat{3,1});
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
if pval<0.05
mysigstar(gca, [2 3], 1.8, pval);
end


% add significance stars for each bar
% for b = 1:3,
%     [~, pval] = ttest(dat(:, b));
%     yval = mean(dat(:, b)) * 0.5; % plot this on top of the bar
%     mysigstar(gca, b, yval, pval);
%     % if mysigstar gets just 1 xpos input, it will only plot stars
% end


% manage grids

% grid on
% set(gca,'color',[230/255 230/255 230/255])
% box off
% grid on
% set(gca,'LineWidth',4) 
% ax = gca;
% ax.GridColor = [1 1 1];
% ax.GridAlpha = 1;


set(gca,'LineWidth',3)





%%

close all
figure,scatter(Tval_LC_fb_rew_neu(:,1),mean(hit_FAs_del,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_LC_fb_rew_neu(:,1),mean(hit_FAs_del,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_LC_fb_rew_neu(:,1));
plot(Tval_LC_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('LC activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Hit - FA rate','FontWeight','bold','FontSize',30)


% close all
figure,scatter(Tval_LC_stim_rew_neu(:,1),mean(hit_FAs_del,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_LC_fb_rew_neu(:,1),mean(hit_FAs_del,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_LC_fb_rew_neu(:,1));
plot(Tval_LC_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('LC activation (Stim: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Hit - FA rate','FontWeight','bold','FontSize',30)


%%

close all
figure,scatter(Tval_SNVTA_fb_rew_neu(:,1),mean(hit_FAs_del,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_SNVTA_fb_rew_neu(:,1),mean(hit_FAs_del,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_SNVTA_fb_rew_neu(:,1));
plot(Tval_SNVTA_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('SNVTA activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Hit - FA rate','FontWeight','bold','FontSize',30)


% close all
figure,scatter(Tval_SNVTA_stim_rew_neu(:,1),mean(hit_FAs_del,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_SNVTA_fb_rew_neu(:,1),mean(hit_FAs_del,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_SNVTA_fb_rew_neu(:,1));
plot(Tval_SNVTA_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('SNVTA activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Hit - FA rate','FontWeight','bold','FontSize',30)


%%


close all
figure,scatter(Tval_SNVTA_fb_rew_neu(:,1),mean(conf_corrects,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_SNVTA_fb_rew_neu(:,1),mean(conf_corrects,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_SNVTA_fb_rew_neu(:,1));
plot(Tval_SNVTA_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('SNVTA activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Confidence rating, correct trials','FontWeight','bold','FontSize',30)


% close all
figure,scatter(Tval_SNVTA_stim_rew_neu(:,1),mean(conf_corrects,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_SNVTA_fb_rew_neu(:,1),mean(conf_corrects,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_SNVTA_fb_rew_neu(:,1));
plot(Tval_SNVTA_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('SNVTA activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Confidence rating, correct trials','FontWeight','bold','FontSize',30)



%%
close all
figure,scatter(Tval_LC_fb_rew_neu(:,1),mean(conf_corrects,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_LC_fb_rew_neu(:,1),mean(conf_corrects,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_LC_fb_rew_neu(:,1));
plot(Tval_LC_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('LC activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Confidence rating, correct trials','FontWeight','bold','FontSize',30)


% close all
figure,scatter(Tval_LC_stim_rew_neu(:,1),mean(conf_corrects,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_LC_fb_rew_neu(:,1),mean(conf_corrects,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_LC_fb_rew_neu(:,1));
plot(Tval_LC_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('LC activation (stim: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Confidence rating, correct trials','FontWeight','bold','FontSize',30)


%%

close all
figure,scatter(Tval_LC_fb_rew_neu(:,1),mean(hit_FAs_del,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_LC_fb_rew_neu(:,1),mean(hit_FAs_del,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_LC_fb_rew_neu(:,1));
plot(Tval_LC_fb_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('LC activation (FB: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Hit - FA rate','FontWeight','bold','FontSize',30)


% close all
figure,scatter(Tval_LC_stim_rew_neu(:,1),mean(hit_FAs_del,2),200,'MarkerFaceColor',[0.2,0.2,1],'MarkerEdgeColor',[0.2,0.2,1],...
    'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6); grid on; hold on
% figure;
Fit = polyfit(Tval_LC_stim_rew_neu(:,1),mean(hit_FAs_del,2),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
f = polyval(Fit,Tval_LC_stim_rew_neu(:,1));
plot(Tval_LC_stim_rew_neu(:,1),f,'-','LineWidth',6)
% x = 3 : 6;
% m = 0.1;
% b = 4.29;
% y = m*x + b;
% plot(x, y,'LineWidth',4,'Color',[1,0.2,0.2],'LineStyle','--')
set(gca,'FontSize',25)
xlabel('LC activation (Stim: rew vs neu)','FontWeight','bold','FontSize',30)
ylabel('Hit - FA rate','FontWeight','bold','FontSize',30)
