%% Generate physio regressors for MRPET data!
%

%% work log

%   08-04-2020      created the script

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
paths.seg     = [paths.parent 'E_data/EA_raw/EAB_MRI/EABD_segmented/'];
paths.TACs    = [paths.parent 'E_data/EA_raw/EAD_PET/EADC_TACs/RewardTask/'];
paths.figures = [paths.parent 'C_writings/CB_figures/MRPET/MainTask/TACs/'];
paths.physraw = [paths.parent 'E_data/EA_raw/EAE_physio/MRPET/'];

% add toolboxes and functions
% addpath(genpath('/Users/yeojin/Documents/MATLAB/spm12'))
addpath(paths.funx_MRI)
addpath(paths.funx_PET)

% IDs
% IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4026];
% days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0; 0 2];
IDs=[4007];
days=[1 0];


% load experimental details
expdat = [];
for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            fname_beh{i1,d} = {NaN};
            expdat{i1,d} = {NaN};
        else
            mkdir([paths.physraw num2str(IDs(i1)) num2str(d)])
            fname_beh{i1,d}     = [ num2str(IDs(i1)) '_' num2str(days(i1,d)) '.mat' ];
            expdat{i1,d} = load([paths.behav fname_beh{i1,d}]);
        end
    end
end

%% set up pre-defined parameters

Nslices = 51;
TR      = 3.6;
Nscans  = 915;

% pulse file names and paths
for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            puls{id,d} = [];
            resp{id,d} = [];
            scan{id,d} = [];
        else
            if (IDs(id) == 4007 && d==1) ||  (IDs(id) == 4009 && d==2) || ...
                    (IDs(id) == 4012 && d==1) || (IDs(id) == 4013 && d==2)% this id has only resp and info!
            tmp2 = dir([paths.physraw num2str(IDs(id)) num2str(d) '/*/*RESP.log']);
            resp{id,d} = fullfile(tmp2.folder,tmp2.name);
            puls{id,d} = [];
            resp{id,d} = [];
            tmp3 = dir([paths.physraw num2str(IDs(id)) num2str(d) '/*/*Info.log']);
            scan{id,d} = fullfile(tmp3.folder,tmp3.name);

            else
            clear tmp1 tmp2 tmp3
            tmp1 = dir([paths.physraw num2str(IDs(id)) num2str(d) '/*/*PULS.log']);
            puls{id,d} = fullfile(tmp1.folder,tmp1.name);
            tmp2 = dir([paths.physraw num2str(IDs(id)) num2str(d) '/*/*RESP.log']);
            resp{id,d} = fullfile(tmp2.folder,tmp2.name);
            tmp3 = dir([paths.physraw num2str(IDs(id)) num2str(d) '/*/*Info.log']);
            scan{id,d} = fullfile(tmp3.folder,tmp3.name);
            end
        end
    end
end

%% run batch

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            close all
            clear matlabbatch
            
            mkdir(['/Users/yeojin/Desktop/E_data/EB_cleaned/EBD_mrpet/RewardTask/physio/' num2str(IDs(id)) num2str(days(id,d)) '/'])
            
            spm_jobman('initcfg');
            
            matlabbatch{1}.spm.tools.physio.save_dir = {['/Users/yeojin/Desktop/E_data/EB_cleaned/EBD_mrpet/RewardTask/physio/' num2str(IDs(id)) num2str(days(id,d)) '/']};
            matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Siemens_Tics';
            if (IDs(id) == 4007 && d==1) ||  (IDs(id) == 4009 && d==2) || ...
                    (IDs(id) == 4012 && d==1) || (IDs(id) == 4013 && d==2)% this id has only resp and info!
            matlabbatch{1}.spm.tools.physio.log_files.cardiac = '<UNDEFINED>';
            else
            matlabbatch{1}.spm.tools.physio.log_files.cardiac = puls(id,d);
            end
            matlabbatch{1}.spm.tools.physio.log_files.respiration = resp(id,d);
            matlabbatch{1}.spm.tools.physio.log_files.scan_timing = scan(id,d);
            matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = [];
            matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
            matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'last';
            
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = Nslices;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = Nscans;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = 1;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
            
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'ECG';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.filter.no = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
            
            matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'multiple_regressors.txt';
            matlabbatch{1}.spm.tools.physio.model.output_physio = ['physio.mat'];
            matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
            matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
            matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
            
            matlabbatch{1}.spm.tools.physio.verbose.level = 2;
            matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = '';
            matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;
            
            spm_jobman('run', matlabbatch) % run batch
            
        end
    end
end


%% append to the realignment parameters

path_physio = '/Users/yeojin/Desktop/E_data/EB_cleaned/EBD_mrpet/RewardTask/physio/';
path_movemt = paths.preproc;

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0
            
        else
            
            clear filepath_physio filepath_move reg_move reg_all reg_physio
            
            filepath_physio = fullfile(path_physio,[num2str(IDs(id)) num2str(d)],'multiple_regressors.txt');
            reg_physio      = importdata(filepath_physio);
            
            filepath_move   = fullfile(path_movemt,[num2str(IDs(id)) '_' num2str(d)],['rp_a' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.txt']);
            reg_move        = importdata(filepath_move);
            
            reg_all         = [reg_move reg_physio];
            
            cd([path_physio num2str(IDs(id)) num2str(d)])
            dlmwrite('reg_all.txt',reg_all, 'delimiter','	','newline','pc');
            
        end
    end
end