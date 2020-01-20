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
paths.spm     = [paths.parent 'B_scripts/BE_toolboxes/spm12/'];
paths.funx    = [paths.parent 'B_scripts/BA_preprocessing/BAB_MRI/preproc_functions/'];
paths.raw     = [paths.parent 'E_data/EA_raw/EAD_PET/EADY_originals/DOPE/'];
paths.dat_3D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/A_3D/'];
paths.dat_4D  = [paths.parent 'E_data/EA_raw/EAD_PET/EADA_converted/RewardTask/B_4D/'];
paths.MT      = [paths.parent 'E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/'];
paths.history = [paths.parent 'E_data/EA_raw/EAD_PET/EABX_history/'];
paths.dcm2nii = '/Volumes/MRIcron/MRIcron.app/Contents/Resources/dcm2niix'; % specify the path to your mricron toolbox 

% add toolboxes and functions
addpath(paths.spm)
addpath(paths.funx)  % those two are useless for now

% IDs & sessions
IDs = [4001];
days = [0 2]; 

% file-specicfic variables: choose whatever you'd like to convert
extensions    = '*.IMA';
filetypes_MRI = {'MoCoSeries';'MPRAGE';'Hippocampus';'GRE3D';'t2_tse'};
filetypes_PET = {'PET-InFlow_HD*_AC_';'PET-Baseline_HD*_AC_'};

%% MRI part

%% run

for id = 1:length(IDs)
    
    for imc = 1:length(filetypes_MRI)
        fprintf('\n %s \n', filetypes_MRI{imc})
        
        for d = 1:2
            
            if days(id,d) == 0
                warning('day %d doesn''t exist for this participant',d)
            else
                
                fprintf('\n ID: %d\n', IDs(id))
                
                
                spm_jobman('initcfg')
                clear matlabbatch
                
                %% dcm to 3D
                
                % setup directories and working directories
                cd(paths.dat_3D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                cd([paths.raw num2str(IDs(id)) '_' num2str(d)])
                tmpfname = dir('study*'); cd(tmpfname.name)
                numiter   = length(dir(['*' filetypes_MRI{imc} '*' ]));
                
                if numiter >= 2
                    
                    clear fnames
                    cd(paths.dat_3D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                    cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                    templist1   = dir(['*' filetypes_MRI{imc} '*' ]);
                    
                    for scans = 1:numiter
                        
                        clear fnames
                        path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                        
                        cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                        tempdir1    = templist1(scans).name; cd(tempdir1);
                        fnames = dir('MR*');
                        
                        cd([paths.dat_3D num2str(IDs(id)) '_' num2str(d)]); mkdir(templist1(scans).name);
                        savepath    = [paths.dat_3D num2str(IDs(id)) '_' num2str(d) '/' templist1(scans).name '/'];
                        
                        clear flist_dicom flist_dicom_tmp
                        for i2 = 1:length(fnames) % assemble into a list
                            flist_dicom_tmp{i2,:} = [ paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(scans).name '/' fnames(i2,1).name ];
                        end
                        flist_dicom = flist_dicom_tmp;
                        
                        % generate headers
%                         hdr = spm_dicom_headers(flist_dicom); % what's happening?!?!?!
                        if imc == 1 % maintask
                            fname_string = [num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '_' num2str(scans)];
                        else
                            fname_string = [num2str(IDs(id)) '_MRI_4D_' filetypes_MRI{imc} num2str(d) '_' num2str(scans)];
                        end
                        
                        % run 3D conversion
                        [s,w] = unix([paths.dcm2nii... % run unix command
                            ' -f "' fname_string '"'... % specify filename
                            ' -p y -z n -ba n -o "' ... % specify options (for help run: [s,w] = unix([paths.dcm2nii ' -h']))
                            paths.dat_4D num2str(IDs(id)) '_' num2str(d) '" ' ... % then, specify the output path
                            '"' paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
                        clear flist_dicom
                        
                        fprintf('\n\n 3D CONVERT DONE \n\n')
                        
                        %% move file for preprocessing
                        cd(paths.MT); mkdir([num2str(IDs(id))  '_' num2str(d)]);
                        copyfile([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/' fname_string '.nii'], ...
                            [paths.MT num2str(IDs(id)) '_' num2str(d) '/' fname_string '.nii'])
                        clear seqname matlabbatch
                        
                    end
                elseif numiter < 2
                    
                    clear fnames
                    path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                    
                    cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                    templist1   = dir(['*' filetypes_MRI{imc} '*' ]);
                    tempdir1    = templist1.name; cd(tempdir1);
                    fnames = dir('MR*');
                    
                    cd([paths.dat_3D num2str(IDs(id)) '_' num2str(d)]); mkdir(templist1(1).name);
                    savepath    = [paths.dat_3D num2str(IDs(id)) '_' num2str(d) '/' templist1(1).name '/'];
                    
                    clear flist_dicom flist_dicom_tmp
                    for i2 = 1:length(fnames) % assemble into a list
                        flist_dicom_tmp{i2,:} = [ paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
                    end
                    flist_dicom = flist_dicom_tmp;
                    clear fnames
                    
                    if imc == 1 % maintask
                        fname_string = [num2str(IDs(id)) '_MRI_4D_MT' num2str(d)];
                    else
                        fname_string = [num2str(IDs(id)) '_MRI_4D_' filetypes_MRI{imc} num2str(d)];
                    end
                    
                    % run 4D conversion
                    [s,w] = unix([paths.dcm2nii... % run unix command
                        ' -f "' fname_string '"'... % specify filename
                        ' -p y -z n -ba n -o "' ... 
                        paths.dat_4D num2str(IDs(id)) '_' num2str(d) '" ' ... % first, specify the output path
                         '"' paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '"']) % now specify the dicom path
                    clear flist_dicom
                    
                    fprintf('\n\n 4D CONVERT DONE \n\n')
                    
                    
                    %% move file for preprocessing
                    cd(paths.MT); mkdir([num2str(IDs(id)) '_' num2str(d)]);
                    copyfile([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/' fname_string '.nii'], ...
                        [paths.MT num2str(IDs(id)) '_' num2str(d) '/' fname_string '.nii'])
                    clear seqname matlabbatch fname_string
                    
                end
            end
        end
        
    end
    
end

%% PET part: other
%% run

% spm fmri
for id = 1:length(IDs)
    
    for imc = 1:length(filetypes_PET)
        fprintf('\n %s \n', filetypes_PET{imc})
        
        for d = 1:2
            
            if days(id,d) == 0
                warning('day %d doesn''t exist for this participant',d)
            else
                
                fprintf('\n ID: %d\n', IDs(id))
                
                
                spm_jobman('initcfg')
                clear matlabbatch
                
                %% dcm to 3D
                
                % setup directories and working directories
                cd(paths.dat_3D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                cd([paths.raw num2str(IDs(id)) '_' num2str(d)])
                tmpfname = dir('study*'); cd(tmpfname.name)
                numiter   = length(dir(['*' filetypes_PET{imc} '*' ]));
                
                if numiter >= 2
                    
                    clear fnames
                    cd(paths.dat_3D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                    cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                    templist1   = dir(['*' filetypes_PET{imc} '*' ]);
                    
                    for scans = 1:numiter
                        
                        clear fnames
                        path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                        
                        cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                        tempdir1    = templist1(scans).name; cd(tempdir1);
                        fnames = dir('PI*');
                        
                        cd([paths.dat_3D num2str(IDs(id)) '_' num2str(d)]); mkdir(templist1(scans).name);
                        savepath    = [paths.dat_3D num2str(IDs(id)) '_' num2str(d) '/' templist1(scans).name '/'];
                        
                        clear flist_dicom flist_dicom_tmp
                        for i2 = 1:length(fnames) % assemble into a list
                            flist_dicom_tmp{i2,:} = [ paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(scans).name '/' fnames(i2,1).name ];
                        end
                        flist_dicom = flist_dicom_tmp;
                        
                        if imc <= seriesnum % maintask
                            fname_string = [num2str(IDs(id)) '_PET_3D_T' num2str((imc-1)*binsize) '_MT' num2str(d) '_' num2str(scans)];
                        else
                            fname_string = [num2str(IDs(id)) '_PET_4D_' filetypes_PET{imc} num2str(d) '_' num2str(scans)];
                        end
                        
                        % run 3D conversion
                        [s,w] = unix([paths.dcm2nii... 
                            ' -f "' fname_string '"'... 
                            ' -p y -z n -ba n -o "' ...
                            paths.dat_4D num2str(IDs(id)) '_' num2str(d) '" ' ... 
                            '"' paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '"']) 
                        clear flist_dicom
                        
                        fprintf('\n\n 3D CONVERT DONE \n\n')
                                               
                        
                        %% move file for preprocessing
                        cd(paths.MT); mkdir([num2str(IDs(id)) '_' num2str(d)]);
                        cd([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/'])
                        copyfile([fname_string '.nii'],[paths.MT num2str(IDs(id)) '_' num2str(d) '/'])
                        clear seqname matlabbatch fname_string
                        
                    end
                elseif numiter < 2
                    
                    clear fnames
                    path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                    
                    cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                    templist1   = dir(['*' filetypes_PET{imc} '*' ]);
                    tempdir1    = templist1.name; cd(tempdir1);
                    fnames = dir('PI*');
                    
                    cd([paths.dat_3D num2str(IDs(id)) '_' num2str(d)]); mkdir(templist1(1).name);
                    savepath    = [paths.dat_3D num2str(IDs(id)) '_' num2str(d) '/' templist1(1).name '/'];
                    
                    clear flist_dicom flist_dicom_tmp
                    for i2 = 1:length(fnames) % assemble into a list
                        flist_dicom_tmp{i2,:} = [ paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
                    end
                    flist_dicom = flist_dicom_tmp;
                    clear fnames
                    
                    fname_string = [num2str(IDs(id)) '_PET_4D_' filetypes_PET{imc} num2str(d)];
                    
                    % run 3D conversion
                    [s,w] = unix([paths.dcm2nii... 
                        ' -f "' fname_string '"'... 
                        ' -p y -z n -ba n -o "' ... 
                        paths.dat_4D num2str(IDs(id)) '_' num2str(d) '" ' ... 
                         '"' paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '"'])
                    clear flist_dicom
                    
                    fprintf('\n\n 3D CONVERT DONE \n\n')
                    
                    
                    %% move file for preprocessing
                    cd(paths.MT); mkdir([num2str(IDs(id)) '_' num2str(d)]);
                    cd([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/'])
                    copyfile([fname_string '.nii'],[paths.MT num2str(IDs(id)) '_' num2str(d) '/'])
                    clear seqname matlabbatch fname_string
                    
                end
            end
        end
        
    end
    
end

%% PET part: task
%% run

filetypes_pTask = {'PET-Task-MoCo*T0*_AC_';'PET-Task-MoCo*T300_*_AC_';'PET-Task-MoCo*T600*_AC_';'PET-Task-MoCo*T900*_AC_';...
    'PET-Task-MoCo*T1200*_AC_';'PET-Task-MoCo*T1500*_AC_';'PET-Task-MoCo*T1800*_AC_';'PET-Task-MoCo*T2100*_AC_';...
    'PET-Task-MoCo*T2400*_AC_';'PET-Task-MoCo*T2700*_AC_';'PET-Task-MoCo*T3000_*_AC_'};
seriesnum     = numel(~isnan(cell2mat(strfind(filetypes_pTask,'Task')))); % see how many binned PET data there are
binsize       = 300; % how long are the time bins for now?

for id = 1:length(IDs)
    
    for imc = 1:length(filetypes_pTask)
        fprintf('\n %s \n', [filetypes_pTask{imc} num2str((imc-1)*binsize)])
        
        for d = 1:2
            
            if days(id,d) == 0
                warning('day %d doesn''t exist for this participant',d)
            else
                
                fprintf('\n ID: %d\n', IDs(id))
                
                
                spm_jobman('initcfg')
                clear matlabbatch
                
                %% dcm to 3D
                
                % setup directories and working directories
                cd(paths.dat_3D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                cd([paths.raw num2str(IDs(id)) '_' num2str(d)])
                tmpfname = dir('study*'); cd(tmpfname.name)
                numiter   = length(dir(['*' filetypes_pTask{imc} '*' ]));
                
                if numiter >= 2
                    
                    clear fnames
                    cd(paths.dat_3D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                    cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                    templist1   = dir(['*' filetypes_pTask{imc} '*' ]);
                    
                    for scans = 1:numiter
                        
                        clear fnames
                        path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                        
                        cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                        tempdir1    = templist1(scans).name; cd(tempdir1);
                        fnames = dir('PI*');
                        
                        cd([paths.dat_3D num2str(IDs(id)) '_' num2str(d)]); mkdir(templist1(scans).name);
                        savepath    = [paths.dat_3D num2str(IDs(id)) '_' num2str(d) '/' templist1(scans).name '/'];
                        
                        clear flist_dicom flist_dicom_tmp
                        for i2 = 1:length(fnames) % assemble into a list
                            flist_dicom_tmp{i2,:} = [ paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(scans).name '/' fnames(i2,1).name ];
                        end
                        flist_dicom = flist_dicom_tmp;
                        
                        if imc <= seriesnum % maintask
                            fname_string = [num2str(IDs(id)) '_PET_3D_T' num2str((imc-1)*binsize) '_MT' num2str(d) '_' num2str(scans)];
                        else
                            fname_string = [num2str(IDs(id)) '_PET_4D_' filetypes_pTask{imc} num2str(d) '_' num2str(scans)];
                        end
                        
                        % run 3D conversion
                        [s,w] = unix([paths.dcm2nii... 
                            ' -f "' fname_string '"'... 
                            ' -p y -z n -ba n -o "' ...
                            paths.dat_4D num2str(IDs(id)) '_' num2str(d) '" ' ... 
                            '"' paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '"']) 
                        clear flist_dicom
                        
                        fprintf('\n\n 3D CONVERT DONE \n\n')
                                               
                        
                        %% move file for preprocessing
                        cd(paths.MT); mkdir([num2str(IDs(id)) '_' num2str(d)]);
                        cd([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/'])
                        copyfile([fname_string '.nii'],[paths.MT num2str(IDs(id)) '_' num2str(d) '/'])
                        clear seqname matlabbatch fname_string
                        
                    end
                elseif numiter < 2
                    
                    clear fnames
                    path_origin = paths.dat_4D; cd(paths.dat_4D); mkdir([num2str(IDs(id)) '_' num2str(d)])
                    
                    cd([paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name])
                    templist1   = dir(['*' filetypes_pTask{imc} '*' ]);
                    tempdir1    = templist1.name; cd(tempdir1);
                    fnames = dir('PI*');
                    
                    cd([paths.dat_3D num2str(IDs(id)) '_' num2str(d)]); mkdir(templist1(1).name);
                    savepath    = [paths.dat_3D num2str(IDs(id)) '_' num2str(d) '/' templist1(1).name '/'];
                    
                    clear flist_dicom flist_dicom_tmp
                    for i2 = 1:length(fnames) % assemble into a list
                        flist_dicom_tmp{i2,:} = [ paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '/' fnames(i2,1).name ];
                    end
                    flist_dicom = flist_dicom_tmp;
                    clear fnames
                    
                    if imc <= seriesnum % maintask
                        fname_string = [num2str(IDs(id)) '_PET_3D_T' num2str((imc-1)*binsize) '_MT' num2str(d)];
                    else
                        fname_string = [num2str(IDs(id)) '_PET_4D_' filetypes_pTask{imc} num2str(d)];
                    end
                    
                    % run 3D conversion
                    [s,w] = unix([paths.dcm2nii... 
                        ' -f "' fname_string '"'... 
                        ' -p y -z n -ba n -o "' ... 
                        paths.dat_4D num2str(IDs(id)) '_' num2str(d) '" ' ... 
                         '"' paths.raw num2str(IDs(id)) '_' num2str(d) '/' tmpfname.name '/' templist1(1).name '"'])
                    clear flist_dicom
                    
                    fprintf('\n\n 3D CONVERT DONE \n\n')
                    
                    
                    %% move file for preprocessing
                    cd(paths.MT); mkdir([num2str(IDs(id)) '_' num2str(d)]);
                    cd([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/'])
                    copyfile([fname_string '.nii'],[paths.MT num2str(IDs(id)) '_' num2str(d) '/'])
                    clear seqname matlabbatch fname_string
                    
                end
            end
        end
        
    end
    
end

%% now make 4D volume out of PET-task data

for id = 1:length(IDs)
    for d = 1:2
        
        if days(id,d) == 0
            warning('day %d doesn''t exist for this participant',d)
            
        else
            flist_PET = [];
            for frames = 1:seriesnum
                flist_PET{frames,1} = [paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/' ...
                    num2str(IDs(id)) '_PET_3D_T' num2str((frames-1)*binsize) '_MT' num2str(d) '.nii,1'];
            end
            
            % run
            clear matlabbatch
            spm_jobman('initcfg');      % initiate job manager
            matlabbatch{1}.spm.util.cat.vols = flist_PET;
            matlabbatch{1}.spm.util.cat.name = [ num2str(IDs(id)) '_PET_4D_MT' num2str(d) '.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            matlabbatch{1}.spm.util.cat.RT = binsize;
            spm_jobman('run', matlabbatch) % run batch
            
            fprintf('\n\n 4D CONVERT DONE \n\n')
            
            %% move file for preprocessing
            
            fname_string = [num2str(IDs(id)) '_PET_4D_MT' num2str(d)];
            
            cd(paths.MT); mkdir([num2str(IDs(id)) '_' num2str(d)]);
            cd([paths.dat_4D num2str(IDs(id)) '_' num2str(d) '/'])
            copyfile([fname_string '.nii'],[paths.MT num2str(IDs(id)) '_' num2str(d) '/'])
            clear seqname matlabbatch fname_string
            
        end
    end
end