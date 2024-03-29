%% average LC images

clc;clear;

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin']);
setenv('ANTSPATH','/usr/local/bin')

paths=[];
paths.original='/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADY_originals/DOPE/';
paths.converted='/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/LCaverage/';
% 
% IDs = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016];
% days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2];

IDs  = [4017 4018 4019 4020 4021 4022 4023 4024 4026];
days = [1 2; 1 2; 1 0; 1 2; 0 2; 0 2; 1 0; 1 0; 0 2]; 

%% convert images, use dcm2nii


for id=1%:length(IDs)
    
    series=0; % reset the series number
    eval(['!mkdir ' paths.converted '/' num2str(IDs(id))])
    
    for d=1:2
        if days(id,d)==0
            disp('skipped')
            
        else
            
            
            clear tmppath tmp tmp2
            
            try
                tmp=dir([paths.original num2str(IDs(id)) '_' num2str(d) '/study*']);
                tmp=tmp(~cellfun(@(x) x==0, {tmp.isdir}));
                tmp=tmp(~ismember({tmp.name} ,{'.','..','.DS_Store'}));
                tmppath=[tmp(1).folder '/' tmp(1).name];
                
            catch
                tmp=dir([paths.original num2str(IDs(id)) '_' num2str(d)]);
                tmp=tmp(~cellfun(@(x) x==0, {tmp.isdir}));
                tmp=tmp(~ismember({tmp.name} ,{'.','..','.DS_Store'}));
                tmp2=dir([tmp(1).folder '/' tmp(1).name '/study*']);
                tmppath=[tmp2(1).folder '/' tmp2(1).name];
            end
            
            clear tmp3
            tmp3=dir([tmppath '/*GRE3D*BW130']);
            for s1=1:length(tmp3)
                clear GREpath
                series=series+1;
                GREpath=[tmp3(s1).folder '/' tmp3(s1).name];
                
                % convert
                eval(['!/Applications/MRIcron.app/Contents/Resources/dcm2niix -f "GRE' num2str(series) '" -p y -z y -ba n "' GREpath '"'])
                
                % move
                eval(['!mv -v ' GREpath '/GRE' num2str(series) '.nii.gz ' paths.converted num2str(IDs(id)) '/GRE' num2str(series) '.nii.gz'])
            end
        end
    end
end

%% coregister and average

for id=1%:length(IDs)
    
    if sum(days(id,:))==3
        series=2:4;
    elseif sum(days(id,:))==2 || sum(days(id,:))==1
        series=2;
    end
    
    referenceImage=[paths.converted num2str(IDs(id)) '/GRE1.nii.gz'];
    clear regstring
    regstring=[paths.converted num2str(IDs(id)) '/GRE1.nii.gz '];
    for s1=series
        eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -f ' referenceImage ' -m ' paths.converted num2str(IDs(id)) '/GRE' num2str(s1) '.nii.gz -o ' paths.converted num2str(IDs(id)) '/GRE'  num2str(s1) '_' ])
        regstring=[regstring paths.converted num2str(IDs(id)) '/GRE'  num2str(s1) '_Warped.nii.gz '];
    end
    eval(['!AverageImages 3 ' paths.converted num2str(IDs(id)) '/LCslab_averaged.nii 1 ' regstring])
    eval(['!rm ' paths.converted num2str(IDs(id)) '/GRE*.nii.gz'])
    eval(['!rm ' paths.converted num2str(IDs(id)) '/GRE*.mat'])

end


%% copy the files for coregistraiton

for id=1:length(IDs)
            copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/LCaverage/' num2str(IDs(id)) '/LCslab_averaged.nii'],...
               ['/Users/yeojin/Desktop/E_data/ED_coreg/MRPET_LCtemplate/' num2str(IDs(id))  '_LCslab.nii'])
end


%buildtemplateparallel.sh -d 3 -c 1 -r 1 -n 0 -i 8 -t GR -s CC -m 30x90x30 -o zz_ *_LCslab.nii