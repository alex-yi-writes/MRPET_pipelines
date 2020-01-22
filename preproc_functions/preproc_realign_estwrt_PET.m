function [config] = preproc_realign_estwrt_PET(ID,flist)

spm_jobman('initcfg')

%% estimate only
% matlabbatch{1}.spm.spatial.realign.estimate.data = {flist};
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 4;
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
% matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';

%% estimate & write
matlabbatch{1}.spm.spatial.realign.estwrite.data = {flist};
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

spm_jobman('run', matlabbatch);

clear matlabbatch

config.MRI.preproc.realign = 'reslice, 4th-degree B spline';

end