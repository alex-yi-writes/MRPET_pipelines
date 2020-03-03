%% assess movement parameters per subject
clc
clear
close all
addpath('/Users/dorothea/Dropbox/matlab/toolboxes/spm12');
cwd = '/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/';cd(cwd)

% IDs
ids = [4001 4002 4003 4004 4005 4006 4007];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0];
fname = [40011 40012 40021 40022 40031 40032 40041 40042 40051 40052 40061 40062 40071 40072];

mr = nan(length(fname),6,6);

% x larger: more right, y larger, more anterior, z larger higher up
% on average: x around zero, y positive, so moved forward, z largest
% variance, if people sink lower, what happens?
cc = 0;
for id = 1:length(ids)
    for d = 1:2
        if days(id,d) == 0
            cc=cc+1;
            mr(cc,1:85,1:6) = nan(85,6);
        else
            cc=cc+1;
            tmp = spm_select('list',[cwd num2str(ids(id)) '_' num2str(d) '/'], ['rp.*\InFlow' num2str(d) '.txt$']);% 4:6 is rotation, 1:3 is translation, assumable x,y,z
            dum = load([cwd num2str(ids(id)) '_' num2str(d) '/' tmp]);
            if size(dum,1)<85
                mr(cc,1:85,1:size(dum,2)) = [nan(9,size(dum,2));dum];clear dum
            else
                mr(cc,1:size(dum,1),1:size(dum,2)) = dum;clear dum
            end
            
        end
    end
end

%% plot
d1 = [1 3 5 7 9 11 13];
d2 = [2 4 6 8 10 12 14];

for i =1:length(fname)
    mdmove(i,:) = nanmedian(squeeze(mr(i,:,:))); % along scan repetitions
    mdmoveX(i,:) = nanmedian(squeeze(mr(i,:,1))); % along scan repetitions
    mdmoveY(i,:) = nanmedian(squeeze(mr(i,:,2))); % along scan repetitions
    mdmoveZ(i,:) = nanmedian(squeeze(mr(i,:,3))); % along scan repetitions
    maxmoveZ(i,:) = nanmax(squeeze(mr(i,:,3))); % along scan repetitions
    minmoveZ(i,:) = nanmin(squeeze(mr(i,:,3))); % along scan repetitions
    mdabsmove(i,:) = nanmedian(abs(squeeze(mr(i,:,:)))); % along scan repetitions
    mvarmove(i,:) = nanvar(squeeze(mr(i,:,:))); % along scan repetitions
end

%%
% close all
% figure;
% 
% barweb([nanmax(mdabsmove(d1,:)) nanmax(mdabsmove(d2,:))],...
%     [nanstd(mdabsmove(d1,:))./sqrt(length(d1)) nanstd(mdabsmove(d2,:))./sqrt(length(d2))]);
% 
% figure;
% barweb([nanmax(mdmove(d1,:)) nanmax(mdmove(d2,:))],...
%     [nanstd(mdmove(d1,:))./sqrt(length(d1)) nanstd(mdmove(d2,:))./sqrt(length(d2))]);
% 
% figure;
% barweb([nanmax(mvarmove(d1,1:3)) nanmax(mvarmove(d2,1:3))],...
%     [nanstd(mvarmove(d1,1:3))./sqrt(length(d1)) nanstd(mvarmove(d2,1:3))./sqrt(length(d2))]);
% 
% figure;
% barweb([nanmax(mvarmove(d1,1:3)) nanmax(mvarmove(d2,1:3))],...
%     [nanstd(mvarmove(d1,1:3))./sqrt(length(d1)) nanstd(mvarmove(d2,1:3))./sqrt(length(d2))]);
% 
% figure;
% barweb([nanmean(mdmove(d1,1:3)) nanmean(mdmove(d2,1:3))],...
%     [nanstd(mdmove(d1,1:3))./sqrt(length(d1)) nanstd(mdmove(d2,1:3))./sqrt(length(d2))]);

%%
close all
cd('/Users/yeojin/Desktop/C_writings/CB_figures/MRPET/MainTask/etc')
cc2=0;
for id = 1:length(ids)
    for d = 1:2
        if days(id,d) == 0
           cc2=cc2+1;
        else
            cc2=cc2+1;
            figure;
            plot(mr(cc2,:,1),'MarkerSize',10,'LineWidth',5);hold on;set(gca,'FontSize',20, 'FontWeight', 'bold')
            plot(mr(cc2,:,2),'MarkerSize',10,'LineWidth',5);hold on;set(gca,'FontSize',20, 'FontWeight', 'bold')
            plot(mr(cc2,:,3),'MarkerSize',10,'LineWidth',5);hold on;set(gca,'FontSize',20, 'FontWeight', 'bold')
            legend({'X','Y','Z'},'FontSize',10); xlim([0,85])
            title(num2str(ids(id)));xlabel('Frames','FontSize',20, 'FontWeight', 'bold');ylabel('millimetre','FontSize',20, 'FontWeight', 'bold')
            
            fig = gcf; % save
            fig.PaperPositionMode = 'auto';
            print(['movement_' num2str(fname(cc2)) ],'-dpng','-r0')
        end
    end
end

%% see what is related to median value (X, Y or Z)
movedim = minmoveZ; % Z very strongly related to mdmove
figure;scatter(nanmean(mdmove'),movedim)
[r p] = corr(nanmean(mdmove')',movedim)

%% stats
[h c p stats] = ttest2(nanmean(mvarmove(d1,:)'),nanmean(mvarmove(d2,:)'))
[h c p stats] = ttest2(nanmean(mdmove(d1,:)'),nanmean(mdmove(d2,:)'))
[h c p stats] = ttest2(nanmean(mdabsmove(d1,:)'),nanmean(mdabsmove(d2,:)'))

group = nan(1,length(ids));
group(d1) = 1;group(d2) = 2;
[p anovatab stats] = kruskalwallis(nanmean(mvarmove')',group)


%% export
export = [str2num(cell2mat(ids')) nanmean(mdmove')' nanmean(mdabsmove')' nanmean(mvarmove')']