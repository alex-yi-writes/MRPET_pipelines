function [config] = preproc_smoothe(ID,scans,varargin)

if ~isempty(varargin)
    n = varargin{1}(1);
    fprintf('\n*** Smoothing with specified kernel size %d X %d X %d ***\n',n,n,n)
else
    n = 8;
    fprintf('\n***Smoothing with default kernel size %d X %d X %d ***\n',n,n,n)
end


spm_jobman('initcfg')
clear matlabbatch

matlabbatch{1}.spm.spatial.smooth.data = scans;
matlabbatch{1}.spm.spatial.smooth.fwhm = [n n n];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';


spm_jobman('run', matlabbatch);
clear matlabbatch

config.MRI.preproc.normalise = ['normalised:' num2str(n)];

end
