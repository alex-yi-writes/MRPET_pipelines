%% PET data preprocessing pipeline
%  
%% work log

%   18-02-2020      separated this part of the script from the previous
%               pipeline that was both fMRI and PET

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.spm     = '/Users/yeojin/Documents/MATLAB/spm12';
paths.funx_MRI= [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
paths.funx_PET= [paths.parent 'B_scripts/BA_preprocessing/BAC_PET/preproc_functions/'];
paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
paths.preproc = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAB_MRI/EABX_history/MainTask/'];
paths.behav   = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPET/';
paths.seg     = [paths.parent 'E_data/EA_raw/EAD_PET/EADD_segmented/'];
paths.TACs    = [paths.parent 'E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask/'];
paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/']

% add toolboxes and functions
% addpath(genpath('/Users/yeojin/Documents/MATLAB/spm12'))
addpath(paths.funx_MRI)
addpath(paths.funx_PET)

% IDs

% IDs = [4017 4018 4019 4020 4021 4022 4023 4024 4026];
% days = [0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0; 0 2]; 
% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 0];
d1m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 0 0; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 0 0; 1 2; 1 2; 1 2; 1 2; 1 2]; 


% load experimental details
expdat = [];
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            fname_beh{i1,d} = {NaN};
            expdat{i1,d} = {NaN};
        else
            fname_beh{i1,d}     = [ num2str(IDs(i1)) '_' num2str(days(i1,d)) '.mat' ];
            expdat{i1,d} = load([paths.behav fname_beh{i1,d}]);
        end
    end
end

% record history
% flags: 0 not done / 1 done
fg_1_resliced      = 0;
fg_2_realigned     = 0;
fg_3_coregistered  = 0;
fg_4_segmented     = 0;
fg_5_normalised    = 0;
fg_6_smoothed      = 0;

config = [];

fprintf('\n preparation done \n')


%% MRI part

%% run

clear id d

for id = 27:28%1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            fprintf('\n*** ID %d being preprocessed: %2.0f out of %2.0f ***\n', IDs(id), id, length(IDs))
            
            %% manage configuration
            
            config.ID   = IDs(id);
            config.exp  = expdat{id,d}.dat;
            
            %% reslice
            
            if fg_1_resliced == 0
                
                clear nvols
%                 cd([paths.dat_3D num2str(IDs(id)) '_' num2str(days(id,d))])
%                 tmp = dir('*MoCoSeries'); cd(tmp(1).name)
                nvols = length(spm_vol([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']))
                
                for cc = 1:nvols
                    funclist{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end
%                 flist    = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_Main_4D_MT' num2str(days(id,d)) '.nii']);
                nslices  = 51;
                TR       = 3.6;
                refslice = 25;
                
                [config]        = preproc_reslice(IDs(id),funclist,nslices,TR,refslice,[1:2:51,2:2:50]);
                
                fg_1_resliced = 1;
                config.flags.resliced = 1;
                fprintf('\n *** resliced ***\n')
            else
                fprintf('\n Already resliced batch \n')
            end
            
            %% realign
            
            %% calculate VDM, then realign and unwarp
            
            if fg_2_realigned == 0
                
                shortTE         = 10;
                longTE          = 11.02;
                Total_TReadouot = 0.0452196;
                
                clear flist
                for cc = 1:nvols
                    funclist{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/a' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end
                
                
                % ------------------- make mean functional image ------------------- %
                clear matlabbatch
                
                spm_jobman('initcfg')
                matlabbatch{1}.spm.spatial.realign.write.data = funclist;
                matlabbatch{1}.spm.spatial.realign.write.roptions.which = [0 1];
                matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
                matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'mean';
                spm_jobman('run',matlabbatch)
                
                
                
                % ----------------------- calculate VDM map ------------------------ %
                
                meanfMRI        = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']);
                Phase_FMap      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_fieldmap' num2str(d) '_2.nii,1']);
                Magnitude_FMap  = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_fieldmap' num2str(d) '_1.nii,1']);
                
                clear matlabbatch
                
                spm_jobman('initcfg')
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase = Phase_FMap;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude = Magnitude_FMap;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsfile = {'/Users/yeojin/Desktop/B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/pm_defaults_MRPET.m'};
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session.epi = meanfMRI;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 0;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;           
                spm_jobman('run', matlabbatch)
                
                
                
                % ---------------------- realign and unwarp ----------------------- %
                clear flist
                for cc = 1:nvols
                    funclist{cc,1}    = [paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/a' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii,' num2str(cc)];
                end
                VDMmap                = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/vdm5_sc' num2str(IDs(id)) '_MRI_4D_fieldmap' num2str(d) '_2.nii,1']);
                
                clear matlabbatch
                
                spm_jobman('initcfg')
                matlabbatch{1}.spm.spatial.realignunwarp.data.scans = funclist;
                matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = VDMmap;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
                matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
                matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
                spm_jobman('run', matlabbatch)
                
                fg_2_realigned  = 1;
                config.flags.realigned = 1;
                fprintf('\n *** realigned ***\n')
            else
                fprintf('\n Already realigned batch \n')
            end
            
            %% coregister
%             
%             if fg_3_coregistered == 0
%                 
%                 fMRI_mean= cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/meana' num2str(IDs(id)) '_Main_4D_MT' num2str(days(id,d)) '.nii']);
%                 try
%                     sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_Main_4D_MPRAGE' num2str(days(id,d)) '.nii']);
%                     fMRI     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_Main_4D_MT' num2str(days(id,d)) '.nii']);
%                 
%                     [config] = preproc_coregister(IDs(id),fMRI_mean,sMRI,fMRI);
%                 catch
%                     sMRI = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/' num2str(IDs(id)) '_Main_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
%                     fMRI     = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ra' num2str(IDs(id)) '_Main_4D_MT' num2str(days(id,d)) '.nii']);
%                 
%                     [config] = preproc_coregister(IDs(id),fMRI_mean,sMRI,fMRI);
%                 end
%                 fg_3_coregistered = 1;
%                 config.flags.coregistered = 1;
%                 fprintf('\n *** coregistered ***\n')
%             else
%                 fprintf('\n Already realigned batch \n')
%             end
            
            
            %% segment
%             
%             if fg_4_segmented == 0
%                 
%                 [config] = preproc_segment(IDs(id),sMRI); % header information has changed from the previous step
%                 
%                 fg_4_segmented = 1;
%                 config.flags.segmented = 1;
%                 fprintf('\n *** segmented ***\n')
%             else
%                 fprintf('\n Already segmented batch \n')
%             end
            
            %% normalise
            
%             if fg_5_normalised == 0
%                 
%                 try
%                     Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_Main_4D_MPRAGE' num2str(days(id,d)) '.nii']);
%                     [config] = preproc_normalise(IDs(id),Def_field,fMRI);
%                 catch
%                     Def_field = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/y_' num2str(IDs(id)) '_Main_4D_MPRAGE' num2str(days(id,d)) '_1.nii']);
%                     [config] = preproc_normalise(IDs(id),Def_field,fMRI);
%                 end
%                 
%                 fg_5_normalised = 1;
%                 config.flags.normalised = 1;
%                 fprintf('\n *** normalised ***\n')
%             else
%                 fprintf('\n Already normalised batch \n')
%             end
            
            %% smoothe
            
            if fg_6_smoothed == 0
                
                fMRI      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/ua' num2str(IDs(id)) '_MRI_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
%                 fMRI      = cellstr([paths.preproc num2str(IDs(id)) '_' num2str(days(id,d)) '/wra' num2str(IDs(id)) '_Main_4D_MT' num2str(days(id,d)) '.nii']); % normalised and realined
                
                [config] = preproc_smoothe(IDs(id),fMRI,3);
                
                fg_6_smoothed = 1;
                config.flags.smoothed = 1;
                fprintf('\n *** smoothed ***\n')
            else
                fprintf('\n Already smoothed batch \n')
            end
            
            
            %% wrap up
            
            % save configuration
%             save([paths.history num2str(IDs(id)) '_' num2str(days(id,d)) '_config.mat'],'config')
            
            % flags hoisted again
            fg_1_resliced      = 0;
            fg_2_realigned     = 0;
            fg_3_coregistered  = 0;
            fg_4_segmented     = 0;
            fg_5_normalised    = 0;
            fg_6_smoothed      = 0;
            
        end
    end
    fprintf('\n***********\nID %d preprocessing done\n***********\n',IDs(id))
    
end

%% merge series - highDA+lowDA


for id = 1:length(IDs) 
    
    mkdir([paths.preproc num2str(IDs(id))])   
    
    clear matlabbatch fMRI1 fMRI2 volnum1 volnum2 list_scan

    if days(id,1) == 1

        copyfile([paths.preproc num2str(IDs(id)) '_1/ua' num2str(IDs(id)) '_MRI_4D_MT1.nii'],...
            [paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func1.nii'])

        clear matlabbatch
        spm_jobman('initcfg')
        matlabbatch{1}.spm.util.split.vol = {[paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func1.nii']};
        matlabbatch{1}.spm.util.split.outdir = {[paths.preproc num2str(IDs(id))]};
        spm_jobman('run',matlabbatch)

        volnum1 = length(spm_vol([paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func1.nii']));
        list_scan1=[];
        for v1 = 1:volnum1
            list_scan1{v1,1} = [paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func1.nii,' num2str(v1)];
        end
    end

    if days(id,2) == 2

        copyfile([paths.preproc num2str(IDs(id)) '_2/ua' num2str(IDs(id)) '_MRI_4D_MT2.nii'],...
            [paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func2.nii'])

        clear matlabbatch
        spm_jobman('initcfg')
        matlabbatch{1}.spm.util.split.vol = {[paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func2.nii']};
        matlabbatch{1}.spm.util.split.outdir = {[paths.preproc num2str(IDs(id))]};
        spm_jobman('run',matlabbatch)

        volnum2 = length(spm_vol([paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func2.nii']));
        list_scan2=[];
        for v2 = 1:volnum2
            list_scan2{v2,1} = [paths.preproc num2str(IDs(id)) '/ua' num2str(IDs(id)) 'func2.nii,' num2str(v2)];
        end
    end

    if sum(days(id,:)) == 3

        clear files3D
        files3D = [list_scan1;list_scan2];
        cd([paths.preproc num2str(IDs(id))])
        clear matlabbatch
        spm_jobman('initcfg')
        matlabbatch{1}.spm.util.cat.vols = files3D;
        matlabbatch{1}.spm.util.cat.name = [ 'all_ua' num2str(IDs(id)) '.nii'];
        matlabbatch{1}.spm.util.cat.dtype = 4;
        matlabbatch{1}.spm.util.cat.RT = 3.6;
        spm_jobman('run',matlabbatch)
        eval(['!rm ' paths.preproc num2str(IDs(id)) '/ua4*.nii']);

        volnum = length(spm_vol([paths.preproc num2str(IDs(id)) '/all_ua' num2str(IDs(id)) '.nii']));
        list_scan=[];
        for va = 1:volnum
            list_scan{va,1} = [paths.preproc num2str(IDs(id)) '/all_ua' num2str(IDs(id)) '.nii,' num2str(va)];
        end

        clear matlabbatch

        spm_jobman('initcfg')
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_scan}';
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        spm_jobman('run',matlabbatch)

        list_realigned=[];
        for vall=1:length(list_scan)
            list_realigned{vall,1} = [paths.preproc num2str(IDs(id)) '/rall_ua' num2str(IDs(id)) '.nii,' num2str(vall)];
        end

        [config] = preproc_smoothe(IDs(id),list_realigned,3);

    elseif sum(days(id,:)) == 1
        list_scan=list_scan1;
        [config] = preproc_smoothe(IDs(id),list_scan,3);
        eval(['!rm ' paths.preproc num2str(IDs(id)) '/ua4*.nii']);

    elseif sum(days(id,:)) == 2
        list_scan=list_scan2;
        [config] = preproc_smoothe(IDs(id),list_scan,3);
        eval(['!rm ' paths.preproc num2str(IDs(id)) '/ua4*.nii']);
    end

end

