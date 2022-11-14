%% Convert raw dicom files to 4D niftii

%% work log

%   02-12-2019    created the script
%   03-12-2019    conversion via SPM currently doesn't work for some
%                   reason, so i decided to use MRIcron toolbox, dcm2niiX.
%                   for this you have to call in the unix command line
%                   indirectly via matlab... 
%                   *** dependencies: MRIcron or dcm2niiX (curl -fLO https://github.com/rordenlab/dcm2niix/releases/latest/download/dcm2niix_mac.zip)

%% set environmental variables

clear; clc
warning('off','all');

% paths
paths = [];
paths.parent  = '/Users/yeojin/Desktop/';
paths.spm     = ['/Users/yeojin/Documents/MATLAB/spm12/'];
paths.funx    = [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
paths.MT      = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAD_PET/EABX_history/'];
paths.dcm2nii = '/Applications/MRIcron.app/Contents/Resources/dcm2niix'; % specify the path to your mricron toolbox 

% add toolboxes and functions
% addpath(paths.spm)  
addpath(paths.funx)  % those two are useless for now

% IDs & sessions
IDs  = [4027 4028];
days = [1 0; 1 0]; 
load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/ID_KZ.mat')


% file-specicfic variables: choose whatever you'd like to convert
extensions    = '*.IMA';
% filetypes_MRI = {'fieldmap'};
% MRI_filenames = {'fieldmap'};
% filetypes_MRI = {'MoCoSeries';'Hippocampus';'GRE3D';'t2_tse'};
filetypes_MRI = {'fMRI';'Hippocampus';'GRE3D';'t2_tse';'fieldmap'};
MRI_filenames = {'MT';'Hippocampus';'GRE3D';'T2slab';'fieldmap'};
filetypes_PET = {'PET-InFlow*_AC_Images';'PET-Baseline*_AC_'};
PET_filenames = {'InFlow';'Baseline'};

%% MRI part, except MPRAGE

%% run

for id = [50:51]%1:length(ID_KZ)

    for imc = 1:length(filetypes_MRI)

        fprintf('\n %s, %s \n', filetypes_MRI{imc}, MRI_filenames{imc})

        fprintf('\n ID: %s\n', ID_KZ{id,1})


        spm_jobman('initcfg')
        clear matlabbatch

        %% dcm to 3D

        % setup directories and working directories
        cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
        cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2}])
        tmpfname = dir('study*'); cd(tmpfname.name)
        numiter   = length(dir(['*' filetypes_MRI{imc} '*' ]));

        if numiter >= 2

            clear fnames
            cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            templist1   = dir(['*' filetypes_MRI{imc} '*' ]);

            for scans = 1:numiter

                clear fnames
                path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

                cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
                tempdir1    = templist1(scans).name; cd(tempdir1);
                fnames = dir('MR*');

                cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(scans).name);
                savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(scans).name '/'];

                clear flist_dicom flist_dicom_tmp
                for i2 = 1:length(fnames) % assemble into a list
                    flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(scans).name '/' fnames(i2,1).name ];
                end
                flist_dicom = flist_dicom_tmp;

                % run 3D conversion
                clear matlabbatch
                matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
                matlabbatch{1}.spm.util.import.dicom.root = 'flat';
                matlabbatch{1}.spm.util.import.dicom.outdir = {savepath};
                matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
                matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
                matlabbatch{1}.spm.util.import.dicom.convopts.meta = 1;
                matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
                spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

                clear flist_dicom
                %                         [s,w] = unix([paths.dcm2nii... % run unix command
                %                             ' -f "' fname_string '"'... % specify filename
                %                             ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
                %                             paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
                %                             '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
                %                         clear flist_dicom

                fprintf('\n\n 3D CONVERT DONE \n\n')

                %% 3D to 4D

                % setup directories and working directories
                cd(savepath);
                fnames_3D = dir('*.nii');

                clear flist_3D_tmp flist_3D flist_4D flist_4D_alt
                for i3 = 1:length(fnames_3D)
                    flist_3D_tmp{i3,:}  = [pwd '/'  fnames_3D(i3,1).name ', 1'];
                    flist_3D{i3,:}      = [pwd '/'  fnames_3D(i3,1).name];
                    flist_4D{i3,:}      = flist_3D{i3,1}(1,:);
                end

                clear matlabbatch

                % run 4D conversion
                matlabbatch{1}.spm.util.cat.vols = flist_4D;
                fname_string = [ID_KZ{id,1}(1:4) '_MRI_4D_' MRI_filenames{imc} ID_KZ{id,1}(5) '_' num2str(scans) '.nii'];
                matlabbatch{1}.spm.util.cat.name = fname_string;
                matlabbatch{1}.spm.util.cat.dtype = 4;
                matlabbatch{1}.spm.util.cat.RT = NaN;

                spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

                fprintf('\n\n 4D CONVERT DONE \n\n')

                %% move file for preprocessing
                movefile([savepath '/' fname_string], ...
                    [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
                cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
                copyfile([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string], ...
                    [paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
                clear seqname matlabbatch fname_string

            end
        elseif numiter < 2

            clear fnames
            cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            templist1   = dir(['*' filetypes_MRI{imc} '*' ]);

            path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            tempdir1    = templist1(1).name; cd(tempdir1);
            fnames = dir('MR*');

            cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(1).name);
            savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(1).name '/'];

            clear flist_dicom flist_dicom_tmp
            for i2 = 1:length(fnames) % assemble into a list
                flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
            end
            flist_dicom = flist_dicom_tmp;

            % run 3D conversion
            cd(savepath)
            clear matlabbatch
            matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
            matlabbatch{1}.spm.util.import.dicom.root = 'flat';
            matlabbatch{1}.spm.util.import.dicom.outdir = cellstr(savepath);
            matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
            matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
            matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
            matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
            spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

            clear flist_dicom
            %                     [s,w] = unix([paths.dcm2nii... % run unix command
            %                         ' -f "' fname_string '"'... % specify filename
            %                         ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
            %                         paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
            %                         '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
            %                     clear flist_dicom

            fprintf('\n\n 3D CONVERT DONE \n\n')

            %% 3D to 4D

            % setup directories and working directories
            cd(savepath);
            fnames_3D = dir('*.nii');

            clear flist_3D_tmp flist_3D flist_4D flist_4D_alt
            for i3 = 1:length(fnames_3D)
                flist_3D_tmp{i3,:}  = [pwd '/'  fnames_3D(i3,1).name ', 1'];
                flist_3D{i3,:}      = [pwd '/'  fnames_3D(i3,1).name];
                flist_4D{i3,:}      = flist_3D{i3,1}(1,:);
            end

            clear matlabbatch

            % run 4D conversion
            matlabbatch{1}.spm.util.cat.vols = flist_4D;
            fname_string = [ID_KZ{id,1}(1:4) '_MRI_4D_' MRI_filenames{imc} ID_KZ{id,1}(5) '.nii'];
            matlabbatch{1}.spm.util.cat.name = fname_string;
            matlabbatch{1}.spm.util.cat.dtype = 4;
            matlabbatch{1}.spm.util.cat.RT = NaN;

            spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

            fprintf('\n\n 4D CONVERT DONE \n\n')

            %% move file for preprocessing
            movefile([savepath '/' fname_string], ...
                [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
            cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
            copyfile([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string], ...
                [paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
            clear seqname matlabbatch fname_string

        end
    end

end


%% Rest of the MRI part, only MPRAGE
%% run

clear filetypes_MRI MRI_filenames
filetypes_MRI = {'MPRAGE'};
MRI_filenames = {'MPRAGE'};

for id = [50:51]%1:length(ID_KZ)

    for imc = 1:length(filetypes_MRI)
        fprintf('\n %s, %s \n', 'MPRAGE', MRI_filenames{imc})


        fprintf('\n ID: %s\n', ID_KZ{id,1})


        spm_jobman('initcfg')
        clear matlabbatch

        %% dcm to 3D

        % setup directories and working directories
        cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
        cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2}])
        tmpfname = dir('study*'); cd(tmpfname.name)
        numiter   = length(dir(['*' filetypes_MRI{imc} '*' ]));

        if numiter >= 2

            clear fnames
            cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            templist1   = dir(['*' filetypes_MRI{imc} '*' ]);
            [~, reindex] = sort(str2double(regexp({templist1.name}, '\d+', 'match', 'once')))
            templist1 = templist1(reindex)

            for scans = [2 4]

                clear fnames
                path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

                cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
                tempdir1    = templist1(scans).name; cd(tempdir1);
                fnames = dir('MR*');

                cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(scans).name);
                savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(scans).name '/'];

                clear flist_dicom flist_dicom_tmp
                for i2 = 1:length(fnames) % assemble into a list
                    flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(scans).name '/' fnames(i2,1).name ];
                end
                flist_dicom = flist_dicom_tmp;

                % run 3D conversion
                templist1(scans).name
                cd(savepath)
                clear matlabbatch
                matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
                matlabbatch{1}.spm.util.import.dicom.root = 'flat';
                matlabbatch{1}.spm.util.import.dicom.outdir = {savepath};
                matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
                matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
                matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
                matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
                spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

                clear flist_dicom
                %                         [s,w] = unix([paths.dcm2nii... % run unix command
                %                             ' -f "' fname_string '"'... % specify filename
                %                             ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
                %                             paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
                %                             '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
                %                         clear flist_dicom

                fprintf('\n\n 3D CONVERT DONE \n\n')

                %% 3D to 4D
                
                % setup directories and working directories
                cd(savepath);
                fnames_3D = dir('*.nii');

                clear flist_3D_tmp flist_3D flist_4D flist_4D_alt
                for i3 = 1:length(fnames_3D)
                    flist_3D_tmp{i3,:}  = [pwd '/'  fnames_3D(i3,1).name ', 1'];
                    flist_3D{i3,:}      = [pwd '/'  fnames_3D(i3,1).name];
                    flist_4D{i3,:}      = flist_3D{i3,1}(1,:);
                end

                clear matlabbatch

                % run 4D conversion
                matlabbatch{1}.spm.util.cat.vols = flist_4D;
                if scans == 2
                    fname_string = [ID_KZ{id,1}(1:4) '_MRI_4D_' MRI_filenames{imc} ID_KZ{id,1}(5) '_pt1.nii'];
                elseif scans == 4
                    fname_string = [ID_KZ{id,1}(1:4) '_MRI_4D_' MRI_filenames{imc} ID_KZ{id,1}(5) '_pt2.nii'];
                end
                matlabbatch{1}.spm.util.cat.name = fname_string;
                matlabbatch{1}.spm.util.cat.dtype = 4;
                matlabbatch{1}.spm.util.cat.RT = NaN;

                spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

                fprintf('\n\n 4D CONVERT DONE \n\n')

                %% move file for preprocessing
                movefile([savepath '/' fname_string], ...
                    [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
                cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
                copyfile([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string], ...
                    [paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
                clear seqname matlabbatch fname_string

            end
        elseif numiter < 2

            clear fnames
            cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            templist1   = dir(['*' filetypes_MRI{imc} '*' ]);

            path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            tempdir1    = templist1(1).name; cd(tempdir1);
            fnames = dir('MR*');

            cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(1).name);
            savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(1).name '/'];

            clear flist_dicom flist_dicom_tmp
            for i2 = 1:length(fnames) % assemble into a list
                flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
            end
            flist_dicom = flist_dicom_tmp;

            % run 3D conversion
            cd(savepath)
            clear matlabbatch
            matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
            matlabbatch{1}.spm.util.import.dicom.root = 'flat';
            matlabbatch{1}.spm.util.import.dicom.outdir = cellstr(savepath);
            matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
            matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
            matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
            matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
            spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

            clear flist_dicom
            %                     [s,w] = unix([paths.dcm2nii... % run unix command
            %                         ' -f "' fname_string '"'... % specify filename
            %                         ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
            %                         paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
            %                         '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
            %                     clear flist_dicom

            fprintf('\n\n 3D CONVERT DONE \n\n')

            %% 3D to 4D

            % setup directories and working directories
            cd(savepath);
            fnames_3D = dir('*.nii');

            clear flist_3D_tmp flist_3D flist_4D flist_4D_alt
            for i3 = 1:length(fnames_3D)
                flist_3D_tmp{i3,:}  = [pwd '/'  fnames_3D(i3,1).name ', 1'];
                flist_3D{i3,:}      = [pwd '/'  fnames_3D(i3,1).name];
                flist_4D{i3,:}      = flist_3D{i3,1}(1,:);
            end

            clear matlabbatch

            % run 4D conversion
            matlabbatch{1}.spm.util.cat.vols = flist_4D;
            fname_string = [ID_KZ{id,1}(1:4) '_MRI_4D_' MRI_filenames{imc} ID_KZ{id,1}(5) '.nii'];
            matlabbatch{1}.spm.util.cat.name = fname_string;
            matlabbatch{1}.spm.util.cat.dtype = 4;
            matlabbatch{1}.spm.util.cat.RT = NaN;

            spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

            fprintf('\n\n 4D CONVERT DONE \n\n')

            %% move file for preprocessing
            movefile([savepath '/' fname_string], ...
                [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
            cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
            copyfile([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string], ...
                [paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
            clear seqname matlabbatch fname_string

        end
    end

end


%% PET part: other
%% run

for id = [51]%1:length(ID_KZ)

    for imc = 1%1:length(filetypes_PET)
        fprintf('\n %s, %s \n', filetypes_PET{imc}, PET_filenames{imc})


        fprintf('\n ID: %s\n', ID_KZ{id,1})


        spm_jobman('initcfg')
        clear matlabbatch

        %% dcm to 3D

        % setup directories and working directories
        cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
        cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2}])
        tmpfname = dir('study*'); cd(tmpfname.name)
        numiter   = length(dir(['*' filetypes_PET{imc} '*' ]));

        if numiter >= 2

            clear fnames
            cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            templist1   = dir(['*' filetypes_PET{imc} '*' ]);

            for scans = 1:numiter

                clear fnames
                path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

                cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
                tempdir1    = templist1(scans).name; cd(tempdir1);
                fnames = dir('PI*');

                cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(scans).name);
                savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(scans).name '/'];

                clear flist_dicom flist_dicom_tmp
                for i2 = 1:length(fnames) % assemble into a list
                    flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(scans).name '/' fnames(i2,1).name ];
                end
                flist_dicom = flist_dicom_tmp;

                % run 3D conversion
                cd(savepath)
                clear matlabbatch
                matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
                matlabbatch{1}.spm.util.import.dicom.root = 'flat';
                matlabbatch{1}.spm.util.import.dicom.outdir = {savepath};
                matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
                matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
                matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
                matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
                spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

                clear flist_dicom
                %                         [s,w] = unix([paths.dcm2nii... % run unix command
                %                             ' -f "' fname_string '"'... % specify filename
                %                             ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
                %                             paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
                %                             '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
                %                         clear flist_dicom

                fprintf('\n\n 3D CONVERT DONE \n\n')

                %% 3D to 4D

                % setup directories and working directories
                cd(savepath);
                fnames_3D = dir('s*.nii');

                clear flist_3D_tmp flist_3D flist_4D flist_4D_alt
                for i3 = 1:length(fnames_3D)
                    flist_3D_tmp{i3,:}  = [pwd '/'  fnames_3D(i3,1).name ', 1'];
                    flist_3D{i3,:}      = [pwd '/'  fnames_3D(i3,1).name];
                    flist_4D{i3,:}      = flist_3D{i3,1}(1,:);
                end

                clear matlabbatch

                % run 4D conversion
                matlabbatch{1}.spm.util.cat.vols = flist_4D;
                fname_string = [ID_KZ{id,1}(1:4) '_PET_4D_' PET_filenames{imc} ID_KZ{id,1}(5) '_' num2str(scans) '.nii'];
                matlabbatch{1}.spm.util.cat.name = fname_string;
                matlabbatch{1}.spm.util.cat.dtype = 4;
                matlabbatch{1}.spm.util.cat.RT = NaN;

                spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

                fprintf('\n\n 4D CONVERT DONE \n\n')

                %% move file for preprocessing
                movefile([savepath '/' fname_string], ...
                    [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
                cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
                copyfile([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string], ...
                    [paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
                clear seqname matlabbatch fname_string

            end
        elseif numiter < 2

            clear fnames
            cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2}  '/' tmpfname.name])
            templist1   = dir(['*' filetypes_PET{imc} '*' ]);

            path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

            cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
            tempdir1    = templist1(1).name; cd(tempdir1);
            fnames = dir('PI*');

            cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(1).name);
            savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(1).name '/'];

            clear flist_dicom flist_dicom_tmp
            for i2 = 1:length(fnames) % assemble into a list
                flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
            end
            flist_dicom = flist_dicom_tmp;

            % run 3D conversion
            cd(savepath)
            clear matlabbatch
            matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
            matlabbatch{1}.spm.util.import.dicom.root = 'flat';
            matlabbatch{1}.spm.util.import.dicom.outdir = {savepath};
            matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
            matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
            matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
            matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
            spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

            clear flist_dicom
            %                     [s,w] = unix([paths.dcm2nii... % run unix command
            %                         ' -f "' fname_string '"'... % specify filename
            %                         ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
            %                         paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
            %                         '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
            %                     clear flist_dicom

            fprintf('\n\n 3D CONVERT DONE \n\n')

            %% 3D to 4D

            % setup directories and working directories
            cd(savepath);
            fnames_3D = dir('s*.nii');

            clear flist_3D_tmp flist_3D flist_4D flist_4D_alt
            for i3 = 1:length(fnames_3D)
                flist_3D_tmp{i3,:}  = [pwd '/'  fnames_3D(i3,1).name ', 1'];
                flist_3D{i3,:}      = [pwd '/'  fnames_3D(i3,1).name];
                flist_4D{i3,:}      = flist_3D{i3,1}(1,:);
            end

            clear matlabbatch

            % run 4D conversion
            matlabbatch{1}.spm.util.cat.vols = flist_4D;
            fname_string = [ID_KZ{id,1}(1:4) '_PET_4D_' PET_filenames{imc} ID_KZ{id,1}(5) '.nii'];
            matlabbatch{1}.spm.util.cat.name = fname_string;
            matlabbatch{1}.spm.util.cat.dtype = 4;
            matlabbatch{1}.spm.util.cat.RT = NaN;

            spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

            fprintf('\n\n 4D CONVERT DONE \n\n')

            %% move file for preprocessing
            movefile([savepath '/' fname_string], ...
                [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
            cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
            copyfile([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string], ...
                [paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
            clear seqname matlabbatch fname_string

        end
    end

end



%% PET part: task
% run

filetypes_pTask = {'PET-Task-MoCo*T0*_AC_';'PET-Task-MoCo*T300_*_AC_';'PET-Task-MoCo*T600*_AC_';'PET-Task-MoCo*T900*_AC_';...
    'PET-Task-MoCo*T1200*_AC_';'PET-Task-MoCo*T1500*_AC_';'PET-Task-MoCo*T1800*_AC_';'PET-Task-MoCo*T2100*_AC_';...
    'PET-Task-MoCo*T2400*_AC_';'PET-Task-MoCo*T2700*_AC_';'PET-Task-MoCo*T3000_*_AC_'};
binname       = {'T0';'T300';'T600';'T900';'T1200';'T1500';'T1800';'T2100';'T2400';'T2700';'T3000'};
seriesnum     = numel(~isnan(cell2mat(strfind(filetypes_pTask,'Task')))); % see how many binned PET data there are
binsize       = 300; % how long are the time bins for now?

for id = [40]%1:length(ID_KZ)

    for imc = 1:length(filetypes_pTask)
        fprintf('\n %s \n', [filetypes_pTask{imc} num2str((imc-1)*binsize)])

        fprintf('\n ID: %d\n', ID_KZ{id,1})

        spm_jobman('initcfg')
        clear matlabbatch

        %% dcm to 3D

        % setup directories and working directories
        cd(paths.dat_3D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])
        cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2}])
        tmpfname = dir('study*'); cd(tmpfname.name)
        numiter   = length(dir(['*' filetypes_pTask{imc} '*' ]));

        clear fnames
        path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)])

        cd([paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name])
        templist1   = dir(['*' filetypes_pTask{imc} '*' ]);
        tempdir1    = templist1.name; cd(tempdir1);
        fnames = dir('PI*');

        cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]); mkdir(templist1(1).name);
        savepath    = [paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' templist1(1).name '/'];

        clear flist_dicom flist_dicom_tmp
        for i2 = 1:length(fnames) % assemble into a list
            flist_dicom_tmp{i2,:} = [ paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,2} '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
        end
        flist_dicom = flist_dicom_tmp;
        clear fnames

        fname_string = [ID_KZ{id,1}(1:4) '_PET_3D_T' num2str((imc-1)*binsize) '_MT' ID_KZ{id,1}(5) '.nii'];

        % run 3D conversion
        cd(savepath)
        clear matlabbatch
        matlabbatch{1}.spm.util.import.dicom.data = flist_dicom;
        matlabbatch{1}.spm.util.import.dicom.root = 'flat';
        matlabbatch{1}.spm.util.import.dicom.outdir = {savepath};
        matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
        spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

        clear flist_dicom
        %                 [s,w] = unix([paths.dcm2nii... % run unix command
        %                     ' -f "' fname_string '"'... % specify filename
        %                     ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
        %                     paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '" ' ... % then, specify the output path
        %                     '"' paths.raw ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
        %                 clear flist_dicom
        %
        fprintf('\n\n 3D CONVERT DONE \n\n')


        %% move file for preprocessing
        cd(paths.dat_4D); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
        cd([paths.dat_3D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' tempdir1])
        tmp=dir('*.nii');tmp2=tmp.name;
        copyfile(tmp2,[paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' fname_string])
        clear seqname matlabbatch fname_string tmp tmp2

    end
end

% now make 4D volume out of PET-task data

for id = [40]%1:length(ID_KZ)

    %             try
    flist_PET = [];
    for frames = 1:length(binname)
        flist_PET{frames,1} = [paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ...
            ID_KZ{id,1}(1:4) '_PET_3D_T' num2str((frames-1)*binsize) '_MT' ID_KZ{id,1}(5) '.nii,1'];
    end

    % run
    clear matlabbatch
    spm_jobman('initcfg');      % initiate job manager
    matlabbatch{1}.spm.util.cat.vols = flist_PET;
    matlabbatch{1}.spm.util.cat.name = [ ID_KZ{id,1}(1:4) '_PET_4D_MT' ID_KZ{id,1}(5) '.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 4;
    matlabbatch{1}.spm.util.cat.RT = binsize;
    spm_jobman('run', matlabbatch) % 'run'  laeuft sofort los, 'interactive' laedt alles nochmal in den Batcheditor

    fprintf('\n\n 4D CONVERT DONE \n\n')

    %% move file for preprocessing

    fname_string = [ID_KZ{id,1}(1:4) '_PET_4D_MT' ID_KZ{id,1}(5)];

    cd(paths.MT); mkdir([ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5)]);
    cd([paths.dat_4D ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/'])
    copyfile([fname_string '.nii'],[paths.MT ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/'])
    clear seqname matlabbatch fname_string
    %             catch
    %                 fprintf('\n error: ID: %d\n', ID_KZ{id,1})
    %
    %             end

end

%%

c1=0;
for id=[50:51]
    c1=c1+1;
    IDID{c1}=[ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '1'];
    c1=c1+1;
    IDID{c1}=[ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '2'];

    mkdir(['/Users/yeojin/Desktop/E_data/ED_coreg/mrpetseg/' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '1/'])
    copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,1}(1:4) '_MRI_4D_MPRAGE'  ID_KZ{id,1}(5) '_pt1.nii'],...
        ['/Users/yeojin/Desktop/E_data/ED_coreg/mrpetseg/' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '1/T1pt1.nii'])

    mkdir(['/Users/yeojin/Desktop/E_data/ED_coreg/mrpetseg/' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '2/'])
    copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' ID_KZ{id,1}(1:4) '_' ID_KZ{id,1}(5) '/' ID_KZ{id,1}(1:4) '_MRI_4D_MPRAGE'  ID_KZ{id,1}(5) '_pt2.nii'],...
        ['/Users/yeojin/Desktop/E_data/ED_coreg/mrpetseg/' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '2/T1pt2.nii'])


    cd('/Users/yeojin/Desktop/E_data/ED_coreg/scripts')

    insert_here = 9;

    fid = fopen( 'mrpetseg_pt1base.sh' );
    cac = textscan( fid,'%s', 'Delimiter','\n','CollectOutput',true );
    fclose( fid )
    fid = fopen( ['mrpetseg' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '1.sh'], 'w');
    for jj = 1 : insert_here
        fprintf( fid, '%s\n', cac{1}{jj} );
    end
    fprintf( fid, '%s\n', ['ID=' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '1'] );
    for jj = insert_here+1 : length(cac{1})
        fprintf( fid,'%s\n', cac{1}{jj} );
    end
    fclose( fid );

    fid = fopen( 'mrpetseg_pt2base.sh' );
    cac = textscan( fid,'%s', 'Delimiter','\n','CollectOutput',true );
    fclose( fid )
    fid = fopen( ['mrpetseg' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '2.sh'], 'w');
    for jj = 1 : insert_here
        fprintf( fid, '%s\n', cac{1}{jj} );
    end
    fprintf( fid, '%s\n', ['ID=' ID_KZ{id,1}(1:4) ID_KZ{id,1}(5) '2'] );
    for jj = insert_here+1 : length(cac{1})
        fprintf( fid,'%s\n', cac{1}{jj} );
    end
    fclose( fid );

end