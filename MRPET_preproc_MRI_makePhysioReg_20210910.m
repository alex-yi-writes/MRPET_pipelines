%% make physio data for MRPET

%% prepare

clc;clear

% IDs
IDs         = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015];
days        = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2];

% paths
path_par    ='/Users/yeojin/Desktop/E_data/EA_raw/EAE_physio/physio/MRPET/';
path_save   ='/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/';

% predefined scan parameters
Nslices=51; % how many slices per volume?
TR=3.6; % repetition time
OnsetSlice=23; % slice to which the regressors are temporally aligned. this supposed to be where the most important activation is expected, e.g. where LC lies

%% run

for id=1:length(IDs) % for the length of IDs
    for d=1:2 % for the length of sessions
        if days(id,d) == 0
            disp('no data')
            
        else
            clear matlabbatch name cardiac resp timing
            
            name    = [num2str(IDs(id)) num2str(d)];
            
            % how many files per parameter are there? collect all the names!
            if (IDs(id) == 4007 && d==1) || (IDs(id) == 4009 && d==2) || (IDs(id) == 4012 && d==1) ...
                    || (IDs(id) == 4013&& d==1)
            disp('no cardiac data')
            else
            cardiac =dir([path_par name '/Physio*PULS*']); % change the names here accordingly
            end
            resp    =dir([path_par name '/Physio*RESP*']); % change the names here accordingly
            timing  =dir([path_par name '/Physio*Info*']); % change the names here accordingly
            Nscans = length(spm_vol([path_save num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii'])); % how many volumes are there?
            
            
            spm_jobman('initcfg')
            
            matlabbatch{1}.spm.tools.physio.save_dir = {[path_save num2str(IDs(id)) '_' num2str(d)]}; % where is it saved?
            matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Siemens_Tics'; % which machine is collecting the physio and the scans?
            if (IDs(id) == 4007 && d==1) || (IDs(id) == 4009 && d==2) || (IDs(id) == 4012 && d==1) ...
                    || (IDs(id) == 4013&& d==1) % some IDs were missing the cardio data
            matlabbatch{1}.spm.tools.physio.log_files.cardiac = {''};
            else
            matlabbatch{1}.spm.tools.physio.log_files.cardiac = {[path_par name '/' cardiac.name]}; % each file goes here
            end
            matlabbatch{1}.spm.tools.physio.log_files.respiration = {[path_par name '/' resp.name]}; % each file goes here
            matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {[path_par name '/' timing.name]}; % each file goes here
            matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = [];
            matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
            matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'last'; % align the physio data to the last scan
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = Nslices;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = Nscans;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = OnsetSlice;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.filter.no = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
            matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'multiple_regressors.txt'; % name of the output regressors
            matlabbatch{1}.spm.tools.physio.model.output_physio = 'physio.mat'; % name of the output matrix
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
             
            spm_jobman('run',matlabbatch)
            disp([num2str(IDs(id)) ' done! tada~'])
        end % close session conditionals 
        
    end % close day loop
end % close ID loop