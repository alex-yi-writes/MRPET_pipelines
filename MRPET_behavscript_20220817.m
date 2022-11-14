%% do ROI extraction

clc;
clear;

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';
paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';
paths.ROI     = '/Users/yeojin/Desktop/E_data/EE_atlases_templates/';

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4026];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0; 0 2];
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 0 0; 0 0; 1 2; 0 0]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 0 0; 1 2; 1 2; 1 2; 0 0; 0 0; 1 2];

TR = 3.6;

% load ROIs
LCmask      = spm_read_vols(spm_vol([paths.ROI 'mni_icbm152_LCmetaMask_final.nii']));
SNVTAmask   = spm_read_vols(spm_vol([paths.ROI 'mni_icbm152_SN_VTA.nii']));
BrainstemMask=spm_read_vols(spm_vol([paths.ROI 'mni_icbm152_brainstemMask_short.nii']));

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


%% calc behavstats


behdat=[]; imms=[]; dels=[];
for id=1:length(IDs)
    
    for d=1:2
        
        if days(id,d) == 0
            behdat{id,d}={NaN};
            imms{id,d}={NaN};
            dels{id,d}={NaN};
        else
            
            behdat{id,d}=load([paths.behav num2str(IDs(id)) '_' num2str(d) '.mat']);
            imms{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.immediate;'])
            if eval(['d' num2str(d) 'm(id,2)==0'])
                disp('no task data')
            else
                dels{id,d}=eval(['behdat{id,d}.dat.day' num2str(d) '.memorytest.delayed;'])
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


% calculate stats

% immediate recall
overallAccuracy_imm=[]; FAs_imm=[]; hit_FAs_imm=[]; overallAccuracy_imm=[];
for id=1:length(IDs)
    
    for d=1:2
        if eval(['d' num2str(d) 'm(id,2)==0 | days(id,d) == 0'])
            disp('no task data')
            
            overallAccuracy_imm(id,d)=NaN; 
            FAs_imm(id,d)=NaN; 
            hit_FAs_imm(id,d)=NaN;
            hits_imm(id,d)= NaN;
            hits_rew_imm(id,d) = NaN;
            hits_neu_imm(id,d) = NaN;
            
            FAs_imm(id,d) = NaN;
            FAs_rew_imm(id,d) = NaN;
            FAs_neu_imm(id,d) = NaN;
            hit_FAs_rew_imm(id,d) = NaN;
            hit_FAs_neu_imm(id,d) = NaN;
            
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

            highConfTrials_OverallAccuracy_imm(id,d)=NaN;
            lowConfTrials_OverallAccuracy_imm(id,d)=NaN;

            highConfTrials_hits_imm(id,d)=NaN;
            lowConfTrials_hits_imm(id,d)=NaN;
            highConfTrials_FA_imm(id,d)=NaN;
            lowConfTrials_FA_imm(id,d)=NaN;

            highConfTrials_Dprime_imm(id,d)=NaN;
            lowConfTrials_Dprime_imm(id,d)=NaN;

            numHighConfTrials_rew_imm(id,d)=NaN;
            numHighConfTrials_neu_imm(id,d)=NaN;

        else
            
            % tidy the trials
            clear accu trls oldnew confi missed
            noresp=isnan(imms{id,d}.results.resp);
            accu=imms{id,d}.results.accu;
            trls=cell2mat(imms{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(imms{id,d}.results.trl(:,4))==1;
            confi=imms{id,d}.results.confi;

            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            confi(noresp)=[];
            
            overallAccuracy_imm(id,d)=nanmean(accu);
            hits_imm(id,d)= ( sum(oldnew==1 & accu==1) )/ sum(oldnew==1);
            hits_rew_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==1) )/ sum(oldnew==1 & trls==1);
            hits_neu_imm(id,d) = ( sum(oldnew==1 & accu==1 & trls==0) )/ sum(oldnew==1 & trls==0);
            FAs_imm(id,d)=( sum(oldnew==0 & accu==0) )/ sum(oldnew==0);
            FAs_rew_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==1) )/ sum(oldnew==0 & trls==1);
            FAs_neu_imm(id,d) = ( sum(oldnew==0 & accu==0 & trls==0) )/ sum(oldnew==0 & trls==0);
            hit_FAs_imm(id,d)=hits_imm(id,d)-FAs_imm(id,d);
            hit_FAs_rew_imm(id,d) = hits_rew_imm(id,d)-FAs_rew_imm(id,d);
            hit_FAs_neu_imm(id,d) = hits_neu_imm(id,d)-FAs_neu_imm(id,d);
            
            items_remembered_imm(id,d) = sum(oldnew==1 & accu==1);
            items_forgotten_imm(id,d) = sum(oldnew==1 & accu==0);
            items_old_total_imm(id,d) = sum(oldnew==1);
            
            confidence_imm(id,d)=nanmean(confi);
            confidence_rew_imm(id,d)=nanmean(confi(trls==1));
            confidence_neu_imm(id,d)=nanmean(confi(trls==0));
            conf_corrects_imm(id,d)=nanmean(confi(accu==1));
            conf_rew_corrects_imm(id,d)=nanmean(confi(accu==1&trls==1));
            conf_neu_corrects_imm(id,d)=nanmean(confi(accu==1&trls==0));
            conf_wrongs_imm(id,d)=nanmean(confi(accu==0));
            conf_rew_wrongs_imm(id,d)=nanmean(confi(accu==0&trls==1));
            conf_neu_wrongs_imm(id,d)=nanmean(confi(accu==0&trls==0));

            % confidence-related parameters
            
            highConfTrials_OverallAccuracy_imm(id,d)=nanmean(accu(confi==1));
            lowConfTrials_OverallAccuracy_imm(id,d)=nanmean(accu(confi==0));

            highConfTrials_hits_imm(id,d)=( sum(oldnew==1 & accu==1 & confi==1) ) / sum(oldnew==1 & confi==1);
            lowConfTrials_hits_imm(id,d)=( sum(oldnew==1 & accu==1 & confi==0) ) / sum(oldnew==1 & confi==0);
            highConfTrials_FA_imm(id,d)=( sum(oldnew==0 & accu==0 & confi==1) ) / sum(oldnew==0 & confi==1);
            lowConfTrials_FA_imm(id,d)=( sum(oldnew==0 & accu==0 & confi==0) ) / sum(oldnew==0 & confi==0);

            highConfTrials_Dprime_imm(id,d)=highConfTrials_hits_imm(id,d)-highConfTrials_FA_imm(id,d);
            lowConfTrials_Dprime_imm(id,d)=lowConfTrials_hits_imm(id,d)-lowConfTrials_FA_imm(id,d);

            numHighConfTrials_rew_imm(id,d) = sum(confi==1 & trls==1)/sum(trls==1);
            numHighConfTrials_neu_imm(id,d) = sum(confi==1 & trls==0)/sum(trls==0);


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
            hit_FAs_del(id,d) =NaN;
            hit_FAs_rew_del(id,d) =NaN;
            hit_FAs_neu_del(id,d) =NaN;
            
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

            highConfTrials_OverallAccuracy_del(id,d)=NaN;
            lowConfTrials_OverallAccuracy_del(id,d)=NaN;

            highConfTrials_hits_del(id,d)=NaN;
            lowConfTrials_hits_del(id,d)=NaN;
            highConfTrials_FA_del(id,d)=NaN;
            lowConfTrials_FA_del(id,d)=NaN;

            highConfTrials_Dprime_del(id,d)=NaN;
            lowConfTrials_Dprime_del(id,d)=NaN;

            numHighConfTrials_rew_del(id,d)=NaN;
            numHighConfTrials_neu_del(id,d)=NaN;

        else
            % tidy the trials
            clear accu trls oldnew missed
            noresp=isnan(dels{id,d}.results.resp);
            accu=dels{id,d}.results.accu;
            trls=cell2mat(dels{id,d}.results.trl(:,2))==rews(id,d);
            oldnew=cell2mat(dels{id,d}.results.trl(:,4))==1;
            confi=dels{id,d}.results.confi;
            
            accu(noresp)=[];
            trls(noresp)=[];
            oldnew(noresp)=[];
            confi(noresp)=[];
            
            overallAccuracy_del(id,d)=nanmean(accu);
            hits_del(id,d)=( sum(oldnew==1 & accu==1) )/ sum(oldnew==1);
            hits_rew_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==1) )/ sum(oldnew==1 & trls==1);
            hits_neu_del(id,d) = ( sum(oldnew==1 & accu==1 & trls==0) )/ sum(oldnew==1 & trls==0);
            FAs_del(id,d)=( sum(oldnew==0 & accu==0) )/ sum(oldnew==0);
            FAs_rew_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==1) )/ sum(oldnew==0 & trls==1);
            FAs_neu_del(id,d) = ( sum(oldnew==0 & accu==0 & trls==0) )/ sum(oldnew==0 & trls==0);
            hit_FAs_del(id,d) = hits_del(id,d)-FAs_del(id,d);
            hit_FAs_rew_del(id,d) = hits_rew_del(id,d)-FAs_rew_del(id,d);
            hit_FAs_neu_del(id,d) = hits_neu_del(id,d)-FAs_neu_del(id,d);
            
            items_remembered_del(id,d) = sum(oldnew==1 & accu==1);
            items_forgotten_del(id,d) = sum(oldnew==1 & accu==0);
            items_old_total_del(id,d) = sum(oldnew==1)
            
            confidence_del(id,d)=nanmean(dels{id,d}.results.confi);
            confidence_rew_del(id,d)=nanmean(dels{id,d}.results.confi(trls==1));
            confidence_neu_del(id,d)=nanmean(dels{id,d}.results.confi(trls==0));
            conf_corrects_del(id,d)=nanmean(dels{id,d}.results.confi(accu==1));
            conf_rew_corrects_del(id,d)=nanmean(dels{id,d}.results.confi(accu==1&trls==1));
            conf_neu_corrects_del(id,d)=nanmean(dels{id,d}.results.confi(accu==1&trls==0));
            conf_wrongs_del(id,d)=nanmean(dels{id,d}.results.confi(accu==0));
            conf_rew_wrongs_del(id,d)=nanmean(dels{id,d}.results.confi(accu==0&trls==1));
            conf_neu_wrongs_del(id,d)=nanmean(dels{id,d}.results.confi(accu==0&trls==0));


            highConfTrials_OverallAccuracy_del(id,d)=nanmean(accu(confi==1));
            lowConfTrials_OverallAccuracy_del(id,d)=nanmean(accu(confi==0));

            highConfTrials_hits_del(id,d)=( sum(oldnew==1 & accu==1 & confi==1) ) / sum(oldnew==1 & confi==1);
            lowConfTrials_hits_del(id,d)=( sum(oldnew==1 & accu==1 & confi==0) ) / sum(oldnew==1 & confi==0);
            highConfTrials_FA_del(id,d)=( sum(oldnew==0 & accu==0 & confi==1) ) / sum(oldnew==0 & confi==1);
            lowConfTrials_FA_del(id,d)=( sum(oldnew==0 & accu==0 & confi==0) ) / sum(oldnew==0 & confi==0);

            highConfTrials_Dprime_del(id,d)=highConfTrials_hits_del(id,d)-highConfTrials_FA_del(id,d);
            lowConfTrials_Dprime_del(id,d)=lowConfTrials_hits_del(id,d)-lowConfTrials_FA_del(id,d);

            numHighConfTrials_rew_del(id,d) = sum(confi==1 & trls==1)/sum(trls==1);
            numHighConfTrials_neu_del(id,d) = sum(confi==1 & trls==0)/sum(trls==0);

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

            RT_all_delayed{id,d} = NaN;
            RTmean_delayed_rew(id,d) = NaN;
            RTmean_delayed_neu(id,d) = NaN;
            
        else
            
            trlinfo{id,d} =  eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.config.stim.stimlist_all(:,2))==rews(id,d);']);
            RT_all_immediate{id,d} = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.immediate.results.rt;']);
            RTmean_immediate_rew(id,d) = nanmean(RT_all_immediate{id,d}(trlinfo{id,d}==1));
            RTmean_immediate_neu(id,d) = nanmean(RT_all_immediate{id,d}(trlinfo{id,d}==0));
            
            trlinfo{id,d} =  eval(['cell2mat(expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.config.stim.stimlist_all(:,2))==rews(id,d);']);
            RT_all_delayed{id,d} = eval(['expdat{id,d}.dat.day' num2str(d) '.memorytest.delayed.results.rt;']);
            RTmean_delayed_rew(id,d) = nanmean(RT_all_delayed{id,d}(trlinfo{id,d}==1));
            RTmean_delayed_neu(id,d) = nanmean(RT_all_delayed{id,d}(trlinfo{id,d}==0));
            
        end
    end
end



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

%% make an excel table

Hits_highDA_immediate = hits_imm(:,1);
Hits_lowDA_immediate  = hits_imm(:,2);
Hits_highDA_delayed   = hits_del(:,1);
Hits_lowDA_delayed    = hits_del(:,2);

FAs_highDA_immediate = FAs_imm(:,1);
FAs_lowDA_immediate  = FAs_imm(:,2);
FAs_highDA_delayed   = FAs_del(:,1);
FAs_lowDA_delayed    = FAs_del(:,2);

Dprime_highDA_immediate = Hits_highDA_immediate-FAs_highDA_immediate;
Dprime_highDA_delayed   = Hits_highDA_delayed-FAs_highDA_delayed;
Dprime_lowDA_immediate  = Hits_lowDA_immediate-FAs_lowDA_immediate;
Dprime_lowDA_delayed    = Hits_lowDA_delayed-FAs_lowDA_delayed;

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

Dprime_rew_highDA_immediate = Hits_rew_highDA_immediate-FAs_rew_highDA_immediate;
Dprime_rew_highDA_delayed   = Hits_rew_highDA_delayed-FAs_rew_highDA_delayed;
Dprime_rew_lowDA_immediate  = Hits_rew_lowDA_immediate-FAs_rew_lowDA_immediate;
Dprime_rew_lowDA_delayed    = Hits_rew_lowDA_delayed-FAs_rew_lowDA_delayed;

Dprime_neu_highDA_immediate = Hits_neu_highDA_immediate-FAs_neu_highDA_immediate;
Dprime_neu_highDA_delayed   = Hits_neu_highDA_delayed-FAs_neu_highDA_delayed;
Dprime_neu_lowDA_immediate  = Hits_neu_lowDA_immediate-FAs_neu_lowDA_immediate;
Dprime_neu_lowDA_delayed    = Hits_neu_lowDA_delayed-FAs_neu_lowDA_delayed;

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


highConfTrials_OverallAccuracy_highDA_delayed=highConfTrials_OverallAccuracy_del(:,1);
highConfTrials_OverallAccuracy_lowDA_delayed=highConfTrials_OverallAccuracy_del(:,2);
lowConfTrials_OverallAccuracy_highDA_delayed=lowConfTrials_OverallAccuracy_del(:,1);
lowConfTrials_OverallAccuracy_lowDA_delayed=lowConfTrials_OverallAccuracy_del(:,2);

highConfTrials_Hits_highDA_delayed=highConfTrials_hits_del(:,1);
highConfTrials_Hits_lowDA_delayed=highConfTrials_hits_del(:,2);
lowConfTrials_Hits_highDA_delayed=lowConfTrials_hits_del(:,1);
lowConfTrials_Hits_lowDA_delayed=lowConfTrials_hits_del(:,2);

highConfTrials_FA_highDA_delayed=highConfTrials_FA_del(:,1);
highConfTrials_FA_lowDA_delayed=highConfTrials_FA_del(:,2);
lowConfTrials_FA_highDA_delayed=lowConfTrials_FA_del(:,1);
lowConfTrials_FA_lowDA_delayed=lowConfTrials_FA_del(:,2);


highConfTrials_Dprime_highDA_delayed=highConfTrials_Dprime_del(:,1);
highConfTrials_Dprime_lowDA_delayed=highConfTrials_Dprime_del(:,2);
lowConfTrials_Dprime_highDA_delayed=lowConfTrials_Dprime_del(:,1);
lowConfTrials_Dprime_lowDA_delayed=lowConfTrials_Dprime_del(:,2);

numHighConfTrials_highDA_rew_del=numHighConfTrials_rew_del(:,1);
numHighConfTrials_lowDA_rew_del=numHighConfTrials_rew_del(:,2);
numHighConfTrials_highDA_neu_del=numHighConfTrials_neu_del(:,1);
numHighConfTrials_lowDA_neu_del=numHighConfTrials_neu_del(:,2);


ID=IDs';
MRPET_memoryTable = table(...
    ID, ...
    Hits_highDA_immediate, Hits_highDA_delayed, Hits_lowDA_immediate, Hits_lowDA_delayed,...
    FAs_highDA_immediate,  FAs_highDA_delayed,  FAs_lowDA_immediate,  FAs_lowDA_delayed,...
    Dprime_highDA_immediate, Dprime_highDA_delayed, Dprime_lowDA_immediate, Dprime_lowDA_delayed,...
    Hits_rew_highDA_immediate, Hits_rew_highDA_delayed, Hits_rew_lowDA_immediate, Hits_rew_lowDA_delayed,...
    FAs_rew_highDA_immediate,  FAs_rew_highDA_delayed,  FAs_rew_lowDA_immediate,  FAs_rew_lowDA_delayed,...
    Hits_neu_highDA_immediate, Hits_neu_highDA_delayed, Hits_neu_lowDA_immediate, Hits_neu_lowDA_delayed,...
    FAs_neu_highDA_immediate,  FAs_neu_highDA_delayed,  FAs_neu_lowDA_immediate,  FAs_neu_lowDA_delayed,...
    Dprime_rew_highDA_immediate, Dprime_rew_highDA_delayed, Dprime_rew_lowDA_immediate, Dprime_rew_lowDA_delayed,...
    Dprime_neu_highDA_immediate, Dprime_neu_highDA_delayed, Dprime_neu_lowDA_immediate, Dprime_neu_lowDA_delayed,...
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
    highConfTrials_OverallAccuracy_highDA_delayed, highConfTrials_OverallAccuracy_lowDA_delayed, ...
    lowConfTrials_OverallAccuracy_highDA_delayed, lowConfTrials_OverallAccuracy_lowDA_delayed, ...
    highConfTrials_Hits_highDA_delayed, highConfTrials_Hits_lowDA_delayed, ...
    lowConfTrials_Hits_highDA_delayed, lowConfTrials_Hits_lowDA_delayed, ...
    highConfTrials_FA_highDA_delayed, highConfTrials_FA_lowDA_delayed, lowConfTrials_FA_highDA_delayed, lowConfTrials_FA_lowDA_delayed, ...
    highConfTrials_Dprime_highDA_delayed, highConfTrials_Dprime_lowDA_delayed, lowConfTrials_Dprime_highDA_delayed, lowConfTrials_Dprime_lowDA_delayed,...    
    numHighConfTrials_highDA_rew_del, numHighConfTrials_lowDA_rew_del, numHighConfTrials_highDA_neu_del, numHighConfTrials_lowDA_neu_del)

writetable(MRPET_memoryTable,'MRPET_memoryTable.xls')

%% single-dataset table


% Dprime_rew_highDA_immediate, Dprime_rew_highDA_delayed, Dprime_rew_lowDA_immediate, Dprime_rew_lowDA_delayed,...
%     Dprime_neu_highDA_immediate, Dprime_neu_highDA_delayed, Dprime_neu_lowDA_immediate, Dprime_neu_lowDA_delayed,...

singlesession_DPrime_highDA_immediate=[]; hc=0;
singlesession_DPrime_lowDA_immediate=[];  lc=0;
singlesession_DPrime_highDA_delayed=[];  
singlesession_DPrime_lowDA_delayed=[];   

singlesession_Dprime_rew_highDA_immediate=[];, singlesession_Dprime_rew_highDA_delayed=[];, singlesession_Dprime_rew_lowDA_immediate=[];, singlesession_Dprime_rew_lowDA_delayed=[];,...
    singlesession_Dprime_neu_highDA_immediate=[];, singlesession_Dprime_neu_highDA_delayed=[];, singlesession_Dprime_neu_lowDA_immediate=[];, singlesession_Dprime_neu_lowDA_delayed=[];

for id=1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
        else
            if d==1
                hc=hc+1;
                singlesession_DPrime_highDA_immediate(hc,1)=Dprime_highDA_immediate(id,1);
                singlesession_DPrime_highDA_delayed(hc,1)=Dprime_highDA_delayed(id,1);

                singlesession_Dprime_rew_highDA_immediate(hc,1)=Dprime_rew_highDA_immediate(id,1);
                singlesession_Dprime_rew_highDA_delayed(hc,1)=Dprime_rew_highDA_delayed(id,1);
                singlesession_Dprime_neu_highDA_immediate(hc,1)=Dprime_neu_highDA_immediate(id,1);
                singlesession_Dprime_neu_highDA_delayed(hc,1)=Dprime_neu_highDA_delayed(id,1);
            else
                lc=lc+1;
                singlesession_DPrime_lowDA_immediate(lc,1)=Dprime_lowDA_immediate(id,1);
                singlesession_DPrime_lowDA_delayed(lc,1)=Dprime_lowDA_delayed(id,1);

                singlesession_Dprime_rew_lowDA_immediate(hc,1)=Dprime_rew_lowDA_immediate(id,1);
                singlesession_Dprime_rew_lowDA_delayed(hc,1)=Dprime_rew_lowDA_delayed(id,1);
                singlesession_Dprime_neu_lowDA_immediate(hc,1)=Dprime_neu_lowDA_immediate(id,1);
                singlesession_Dprime_neu_lowDA_delayed(hc,1)=Dprime_neu_lowDA_delayed(id,1);

            end

        end

    end
end
singlesession_DPrime_lowDA_immediate=[singlesession_DPrime_lowDA_immediate];
singlesession_DPrime_lowDA_delayed=[singlesession_DPrime_lowDA_delayed];
singlesession_DPrime_highDA_immediate=[singlesession_DPrime_lowDA_immediate];
singlesession_DPrime_highDA_delayed=[singlesession_DPrime_lowDA_delayed];
MRPET_memoryTable_singlesession = table(singlesession_DPrime_highDA_immediate,singlesession_DPrime_lowDA_immediate,...
    singlesession_DPrime_highDA_delayed,singlesession_DPrime_lowDA_delayed,...
    singlesession_Dprime_rew_highDA_immediate,singlesession_Dprime_rew_highDA_delayed,singlesession_Dprime_neu_highDA_immediate,singlesession_Dprime_neu_highDA_delayed,...
    singlesession_Dprime_rew_lowDA_immediate,singlesession_Dprime_rew_lowDA_delayed,singlesession_Dprime_neu_lowDA_immediate,singlesession_Dprime_neu_lowDA_delayed)
writetable(MRPET_memoryTable_singlesession,'MRPET_memoryTable_singlesession.xls')


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
