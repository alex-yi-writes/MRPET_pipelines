function [config] = preproc_coregister_wrt_PET(ID,sMRI,PET)

spm_jobman('initcfg')
clear matlabbatch

matlabbatch{1}.spm.spatial.coreg.write.ref = sMRI;
matlabbatch{1}.spm.spatial.coreg.write.source = PET;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'k';

spm_jobman('run', matlabbatch);
clear matlabbatch

config.MRI.preproc.coregister = 'estimate and write';

end