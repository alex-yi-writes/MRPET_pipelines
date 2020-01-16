function [config] = preproc_normalise(ID,Def_Field,scans)

clear matlabbatch
spm_jobman('initcfg')

matlabbatch{1}.spm.spatial.normalise.write.subj.def = Def_Field;
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = scans;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm_jobman('run', matlabbatch);
clear matlabbatch

config.MRI.preproc.normalise = 'normalised, 2*2*2';

end

