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
paths.parent   = '/Users/yeojin/Desktop/E_data/transmat/kalina/';
paths.spm      = '/please paste your own SPM12 path here/'; % https://www.fil.ion.ucl.ac.uk/spm/software/download/
paths.funx_MRI = [paths.parent 'scripts/preproc_functions/'];
paths.raw      = [paths.parent 'data/raw/mri/'];
paths.preproc  = [paths.parent 'data/preproc/'];
paths.behav    = [paths.parent 'data/raw/behav/'];

% add toolboxes and functions
addpath(paths.funx_MRI)
load([paths.preproc 'SliceOrder.mat'])
load([paths.preproc 'configuration.mat'])


% IDs
IDs  = [4001 4002 4004 4005 4008 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4028 4030 4031 4032 4033];
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];  % 1 or 2 means that the data are there, 0 means no data

% config = [];
% for id=1:length(IDs)
%
%     for d=1:2
%         if days(id,d) == 0
%             warning('Skipped')
%         else
%             config{id,d}.ID = IDs(id);
%             config{id,d}.preproc = [];
%             config{id,d}.preproc.flag_1_slicetimeCorrected = 0;
%             config{id,d}.preproc.flag_2_unwarped = 0;
%             config{id,d}.preproc.flag_3_sessionMerged = 0;
%             config{id,d}.preproc.flag_4_smoothed = 0;
%             config{id,d}.preproc.flag_5_coregistered = 0;
%
%         end
%     end
%
% end

fprintf('\n preparation done \n')

%% run

clear id d

for id = 1:length(IDs)

    for d = 1:2

        if days(id,d) == 0
            warning('Skipped')
        else

            fprintf('\n*** ID %d being preprocessed: %2.0f out of %2.0f ***\n', IDs(id), id, length(IDs))

            %% reslice

            if config{id,d}.preproc.flag_1_slicetimeCorrected == 0

                clear nvols
                nvols = length(spm_vol([paths.raw num2str(IDs(id)) '/sess' num2str(d) '/fMRI.nii']));
                for cc = 1:nvols
                    funclist{cc,1}    = [paths.raw num2str(IDs(id)) '/sess' num2str(d) '/fMRI.nii,' num2str(cc)];
                end

                nslices  = 51;
                TR       = 3.6;
                refslice = 0;

                [configurations]        = preproc_reslice(IDs(id),funclist,nslices,TR,refslice,sliceOrder);

                config{id,d}.preproc.flag_1_slicetimeCorrected = 1;
                save([paths.preproc 'configuration.mat'],'config')

                fprintf('\n *** resliced ***\n')

            else
                fprintf('\n Already slice-time corrected batch \n')
            end

            %% unwarp and realign
            %  calculate VDM, then realign and unwarp

            if config{id,d}.preproc.flag_2_unwarped == 0

                shortTE         = 10;
                longTE          = 11.02;
                Total_TReadouot = 0.0452196;

                clear flist
                for cc = 1:nvols
                    funclist{cc,1}    = [paths.raw num2str(IDs(id)) '/sess' num2str(d) '/afMRI.nii,' num2str(cc)];
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

                meanfMRI        = cellstr([paths.raw num2str(IDs(id)) '/sess' num2str(d) '/meanafMRI.nii']);
                Phase_FMap      = cellstr([paths.raw num2str(IDs(id)) '/sess' num2str(d) '/phase.nii,1']);
                Magnitude_FMap  = cellstr([paths.raw num2str(IDs(id)) '/sess' num2str(d) '/magnitude.nii,1']);

                clear matlabbatch

                spm_jobman('initcfg')
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase = Phase_FMap;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude = Magnitude_FMap;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsfile = {[paths.funx_MRI '/pm_defaults_MRPET.m']};
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
                    funclist{cc,1}    = [paths.raw num2str(IDs(id)) '/sess' num2str(d) '/afMRI.nii,' num2str(cc)];
                end
                VDMmap                = cellstr([paths.raw num2str(IDs(id)) '/sess' num2str(d) '/vdm5_scphase.nii']);

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

                config{id,d}.preproc.flag_2_unwarped  = 1;
                save([paths.preproc 'configuration.mat'],'config')

                fprintf('\n *** unwarped ***\n')
            else
                fprintf('\n Already unwarped batch \n')
            end

        end
    end

    fprintf('\n***********\nID %d preprocessing done\n***********\n',IDs(id))

end

%% merge series - highDA+lowDA


for id = 1:length(IDs)

    mkdir([paths.preproc num2str(IDs(id))])

    clear matlabbatch fMRI1 fMRI2 volnum1 volnum2 list_scan

    if config{id,d}.preproc.flag_3_sessionMerged == 0

        if days(id,1) == 1

            copyfile([paths.raw num2str(IDs(id)) '/sess1/uafMRI.nii'],...
                [paths.preproc num2str(IDs(id)) '/uafunc1.nii'])

            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[paths.preproc num2str(IDs(id)) '/uafunc1.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[paths.preproc num2str(IDs(id))]};
            spm_jobman('run',matlabbatch)

            volnum1 = length(spm_vol([paths.preproc num2str(IDs(id)) '/uafunc1.nii']));
            list_scan1=[];
            for v1 = 1:volnum1
                list_scan1{v1,1} = [paths.preproc num2str(IDs(id)) '/uafunc1.nii,' num2str(v1)];
            end
        end

        if days(id,2) == 2

            copyfile([paths.preproc num2str(IDs(id)) '/sess2/uafMRI.nii'],...
                [paths.preproc num2str(IDs(id)) '/uafunc2.nii'])

            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[paths.preproc num2str(IDs(id)) '/uafunc2.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[paths.preproc num2str(IDs(id))]};
            spm_jobman('run',matlabbatch)

            volnum2 = length(spm_vol([paths.preproc num2str(IDs(id)) '/uafunc2.nii']));
            list_scan2=[];
            for v2 = 1:volnum2
                list_scan2{v2,1} = [paths.preproc num2str(IDs(id)) '/uafunc2.nii,' num2str(v2)];
            end
        end

        if sum(days(id,:)) == 3

            clear files3D
            files3D = [list_scan1;list_scan2];
            cd([paths.preproc num2str(IDs(id))])
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = files3D;
            matlabbatch{1}.spm.util.cat.name = [ 'all_uafMRI.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            matlabbatch{1}.spm.util.cat.RT = 3.6;
            spm_jobman('run',matlabbatch)
            eval(['!rm ' paths.preproc num2str(IDs(id)) '/uafunc*.nii']);

            volnum = length(spm_vol([paths.preproc num2str(IDs(id)) '/all_uafMRI.nii']));
            list_scan=[];
            for va = 1:volnum
                list_scan{va,1} = [paths.preproc num2str(IDs(id)) '/all_uafMRI.nii,' num2str(va)];
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
                list_realigned{vall,1} = [paths.preproc num2str(IDs(id)) '/rall_uafMRI.nii,' num2str(vall)];
            end

            [config] = preproc_smoothe(IDs(id),list_realigned,3);

        elseif sum(days(id,:)) == 1
            list_scan=list_scan1;
            [config] = preproc_smoothe(IDs(id),list_scan,3);
            eval(['!rm ' paths.preproc num2str(IDs(id)) '/uafunc*.nii']);

        elseif sum(days(id,:)) == 2
            list_scan=list_scan2;
            [config] = preproc_smoothe(IDs(id),list_scan,3);
            eval(['!rm ' paths.preproc num2str(IDs(id)) '/uafunc*.nii']);
        end

    else
        fprintf('\n Already session-merged batch \n')
    end

end

