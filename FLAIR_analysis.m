%% FLAIR - LST analysis

clc; clear;

% prepare IDs
% load('/Volumes/ALEX3/MRPET/IDs_singleSession.mat')
load('/Volumes/ALEX3/MRPET/ID_KZ.mat')

% prepare paths
path_base       = '/Volumes/ALEX3/MRPET/';
path_original   = [path_base 'original/'];
path_raw        = [path_base 'FLAIR/raw/'];
path_analysis   = [path_base 'FLAIR/analysis/'];
path_app        = '/Applications/MRIcron.app/Contents/Resources/dcm2niix';

% make workspace
for id=1:length(ID_KZ)
    
    mkdir([path_analysis ID_KZ{id,1}])
    mkdir([path_raw ID_KZ{id,1}])
    
end

%% convert images from original zip and tidy up for the pipeline

for id=30:31%1:length(ID_KZ)
    
    % unzip
    clear tmp
    tmp=dir([path_original ID_KZ{id,1} '/' ID_KZ{id,2} '*.zip']);
    
    eval(['!unzip -qq ' path_original ID_KZ{id,1} '/' tmp(1).name...
        ' ''*/*/*dzne_FLAIR*iso/*'' -d ''' path_raw ID_KZ{id,1} ''''])
    eval(['!unzip -qq ' path_original ID_KZ{id,1} '/' tmp(1).name...
        ' ''*/*/*dzne_MPRAGE*iso/*'' -d ''' path_raw ID_KZ{id,1} ''''])
    
    % convert
    clear tmp
    tmp=dir([path_raw ID_KZ{id,1} '/' ID_KZ{id,2} '/study*/*FLAIR*']);
    eval(['!' path_app ' -f "FLAIR" -p y -z n -o "' path_analysis ID_KZ{id,1}...
        '" "' tmp(1).folder '/' tmp(1).name '"'])
    
    clear tmp
    tmp=dir([path_raw ID_KZ{id,1} '/' ID_KZ{id,2} '/study*/*MPRAGE*']);
    eval(['!' path_app ' -f "T1" -p y -z n -o "' path_analysis ID_KZ{id,1}...
        '" "' tmp(1).folder '/' tmp(1).name '"'])
    
    eval(['!rm -r ' path_raw ID_KZ{id,1} '/' ID_KZ{id,2}])
    
    disp([ID_KZ{id,1} ' conversion done'])
        
end

%% first LPA, because we don't have prior / probability map

for id=30%:length(ID_KZ)

    clear FLAIR T1 LPAbatch
    FLAIR=[path_analysis ID_KZ{id,1} '/FLAIR.nii'];
    T1=[path_analysis ID_KZ{id,1} '/T1.nii'];
    LPAbatch{1}.spm.tools.LST.lpa.data_F2 = cellstr(FLAIR);
    LPAbatch{1}.spm.tools.LST.lpa.data_coreg = cellstr(T1);
    LPAbatch{1}.spm.tools.LST.lpa.html_report = 1;
    spm_jobman('run', LPAbatch);

end


