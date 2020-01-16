function [config] = preproc_coregister_estwrt_PET(ID,PET_mean,sMRI,PET)

spm_jobman('initcfg')
clear matlabbatch

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = sMRI; 
matlabbatch{1}.spm.spatial.coreg.estwrite.source = PET_mean;
matlabbatch{1}.spm.spatial.coreg.estwrite.other = PET; % functional images that follows
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';

spm_jobman('run', matlabbatch);
clear matlabbatch

config.MRI.preproc.coregister = 'estimate';

end