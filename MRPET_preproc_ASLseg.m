% IDs
IDs         = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017];
days        = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2];



for id=length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('day %d doesn''t exist for this participant',d)
        else
            
            clear matlabbatch
            
            spm_jobman('initcfg')
            
            matlabbatch{1}.spm.spatial.preproc.channel.vols = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d)...
                '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Applications/spm12/tpm/TPM.nii,1'};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Applications/spm12/tpm/TPM.nii,2'};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Applications/spm12/tpm/TPM.nii,3'};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'//Applications/spm12/tpm/TPM.nii,4'};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Applications/spm12/tpm/TPM.nii,5'};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Applications/spm12/tpm/TPM.nii,6'};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
            
            spm_jobman('run',matlabbatch)
            
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.imcalc.input = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d)...
                '/c1' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii,1']};
            matlabbatch{1}.spm.util.imcalc.output = 'grey';
            matlabbatch{1}.spm.util.imcalc.outdir = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/' num2str(IDs(id)) num2str(d) '1']};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.8';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.imcalc.input = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d)...
                '/c2' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii,1']};
            matlabbatch{1}.spm.util.imcalc.output = 'white';
            matlabbatch{1}.spm.util.imcalc.outdir = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/' num2str(IDs(id)) num2str(d) '1']};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.8';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            
        end
    end
end


%% segmentations

load('/Users/yeojin/Desktop/E_data/FSlabels.mat')

for id=1:length(IDs)
    
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('day %d doesn''t exist for this participant',d)
        else
            
            clear aparc
            aparc=spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d)...
                '/aparc+aseg_pt1_nat.nii']));
            
            for l1=2:length(FSlabels)
                
                if isempty(find(aparc==FSlabels{l1,1}))
                    disp(['no ' FSlabels{l1,2} ' seg for ID ' num2str(IDs(id)) '_' num2str(d)])
                else
                    
                    mkdir(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/transmat/katja_segs/' num2str(IDs(id)) num2str(d) '/'])
                    copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/' num2str(IDs(id)) num2str(d) '1/' FSlabels{l1,2} '.nii'],...
                        ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/transmat/katja_segs/' num2str(IDs(id)) num2str(d) '/' FSlabels{l1,2} '.nii'])
%                 clear matlabbatch
%                 spm_jobman('initcfg')
%                 matlabbatch{1}.spm.util.imcalc.input = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' num2str(IDs(id)) '_' num2str(d)...
%                     '/aparc+aseg_pt1_nat.nii,1']};
%                 matlabbatch{1}.spm.util.imcalc.output = FSlabels{l1,2};
%                 matlabbatch{1}.spm.util.imcalc.outdir = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/' num2str(IDs(id)) num2str(d) '1']};
%                 matlabbatch{1}.spm.util.imcalc.expression = ['i1==' num2str(FSlabels{l1,1})];
%                 matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%                 matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
%                 matlabbatch{1}.spm.util.imcalc.options.mask = 0;
%                 matlabbatch{1}.spm.util.imcalc.options.interp = 1;
%                 matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%                 spm_jobman('run',matlabbatch)



                end
            end
            
%             eval(['!gzip /Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/transmat/katja_segs/' num2str(IDs(id)) num2str(d) '/*.nii'])
        end
        
        
        
        
    end
end