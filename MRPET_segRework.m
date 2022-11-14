%% add SNVTALC to the aparc+aseg in native space

clc;clear

load('/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/ID_MRPET.mat')

%%

for id=1:length(ID)
    
    clear SN VTA LC parc
    
    SN.pt1  = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt1.nii']));
    SN.pt2  = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt2.nii']));
    
    VTA.pt1 = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt1.nii']));
    VTA.pt2 = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt2.nii']));
    
    LC.pt1  = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt1.nii']));
    LC.pt2  = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt2.nii']));
    
    parc.pt1  = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt1_nat.nii']));
    parc.pt2  = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt2_nat.nii']));
    
    
    %% SN
    
    
    % ----------- pt 1 ----------- %
    
    clear x y z midpt coord indR indL hdr
    [x,y,z] = ind2sub(size(SN.pt1),find(SN.pt1~=0));
    midpt   = median(unique(z));
    
    coord = [x y z];
    indR  = coord(:,3) < midpt;
    indL  = coord(:,3) > midpt;
    
    coordR = coord(indR,:);
    coordL = coord(indL,:);
    
    for i1=1:length(coordR)
        SN.pt1(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 992;
        parc.pt1(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 992;
    end
    for i1=1:length(coordL)
        SN.pt1(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 991;
        parc.pt1(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 991;
    end
    
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt1.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt1_labelled.nii'];
    hdr.dim = size(SN.pt1);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,SN.pt1);

    
    % ----------- pt 2 ----------- %
    clear x y z midpt coord indR indL hdr
    [x,y,z] = ind2sub(size(SN.pt2),find(SN.pt2~=0));
    midpt   = median(unique(z));
    
    coord = [x y z];
    indR  = coord(:,3) < midpt;
    indL  = coord(:,3) > midpt;
    
    coordR = coord(indR,:);
    coordL = coord(indL,:);
    
    for i1=1:length(coordR)
        SN.pt2(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 992;
        parc.pt2(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 992;
    end
    for i1=1:length(coordL)
        SN.pt2(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 991;
        parc.pt2(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 991;
    end
    
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt2.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt2_labelled.nii'];
    hdr.dim = size(SN.pt2);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,SN.pt2);
    
    disp('SN done')
    
    %% VTA
    
    % ----------- pt 1 ----------- %
    
    clear x y z midpt coord indR indL hdr
    [x,y,z] = ind2sub(size(VTA.pt1),find(VTA.pt1~=0));
    midpt   = median(unique(z));
    
    coord = [x y z];
    indR  = coord(:,3) < midpt;
    indL  = coord(:,3) > midpt;
    
    coordR = coord(indR,:);
    coordL = coord(indL,:);
    
    for i1=1:length(coordR)
        VTA.pt1(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 996;
        parc.pt1(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 996;
    end
    for i1=1:length(coordL)
        VTA.pt1(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 995;
        parc.pt1(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 995;
    end
    
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt1.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt1_labelled.nii'];
    hdr.dim = size(VTA.pt1);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,VTA.pt1);
    
    
    % ----------- pt 2 ----------- %
    
    clear x y z midpt coord indR indL hdr
    [x,y,z] = ind2sub(size(VTA.pt2),find(VTA.pt2~=0));
    midpt   = median(unique(z));
    
    coord = [x y z];
    indR  = coord(:,3) < midpt;
    indL  = coord(:,3) > midpt;
    
    coordR = coord(indR,:);
    coordL = coord(indL,:);
    
    for i1=1:length(coordR)
        VTA.pt2(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 996;
        parc.pt2(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 996;
    end
    for i1=1:length(coordL)
        VTA.pt2(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 995;
        parc.pt2(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 995;
    end
    
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt2.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt2_labelled.nii'];
    hdr.dim = size(VTA.pt2);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,VTA.pt2);
    
    disp('VTA done')
    
    %% LC
    
    % ----------- pt 1 ----------- %
    
    clear x y z midpt coord indR indL hdr
    [x,y,z] = ind2sub(size(LC.pt1),find(LC.pt1~=0));
    midpt   = median(unique(z));
    
    coord = [x y z];
    indR  = coord(:,3) < midpt;
    indL  = coord(:,3) > midpt;
    
    coordR = coord(indR,:);
    coordL = coord(indL,:);
    
    for i1=1:length(coordR)
        LC.pt1(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 994;
        parc.pt1(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 994;
    end
    for i1=1:length(coordL)
        LC.pt1(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 993;
        parc.pt1(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 993;
    end
    
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt1.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt1_labelled.nii'];
    hdr.dim = size(LC.pt1);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,LC.pt1);
    
    
    % ----------- pt 2 ----------- %
    
    clear x y z midpt coord indR indL hdr
    [x,y,z] = ind2sub(size(LC.pt2),find(LC.pt2~=0));
    midpt   = median(unique(z));
    
    coord = [x y z];
    indR  = coord(:,3) < midpt;
    indL  = coord(:,3) > midpt;
    
    coordR = coord(indR,:);
    coordL = coord(indL,:);
    
    for i1=1:length(coordR)
        LC.pt2(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 994;
        parc.pt2(coordR(i1,1),coordR(i1,2),coordR(i1,3)) = 994;
    end
    for i1=1:length(coordL)
        LC.pt2(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 993;
        parc.pt2(coordL(i1,1),coordL(i1,2),coordL(i1,3)) = 993;
    end
    
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt2.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt2_labelled.nii'];
    hdr.dim = size(LC.pt2);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,LC.pt2);
    
    disp('LC done')
    
    
    %% write out the parc images
    
    
    clear hdr
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt1_nat.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt1_nat_labelled.nii'];
    hdr.dim = size(parc.pt1);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,parc.pt1);
    clear hdr
    hdr = spm_vol(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt2_nat.nii']); % pick just any header from a file
    hdr.fname = ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt2_nat_labelled.nii'];
    hdr.dim = size(parc.pt2);
    hdr = rmfield(hdr,'pinfo');
    hdr.nii = spm_write_vol(hdr,parc.pt2);
    
    disp([ID{id} ' done'])
    
    %% now assemble
    %  this ain't working?!
    
%     clear matlabbatch
%     spm_jobman('initcfg')
%     matlabbatch{1}.spm.util.imcalc.input = {
%         ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt1_nat.nii,1']
%         ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/LC_onT1pt1_labelled.nii,1']
%         ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/SN_onT1pt1_labelled.nii,1']
%         ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/data/VTA_onT1pt1_labelled.nii,1']
%         };
%     matlabbatch{1}.spm.util.imcalc.output = 'aparc+aseg_pt1_nat_LCSNVTA';
%     matlabbatch{1}.spm.util.imcalc.outdir = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id}]};
%     matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3+i4';
%     matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%     matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
%     matlabbatch{1}.spm.util.imcalc.options.mask = 0;
%     matlabbatch{1}.spm.util.imcalc.options.interp = 1;
%     matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%     spm_jobman('run',matlabbatch)
        
end

%% copy files

% for id=34%35:length(ID)
% 
%     copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' ID{id}(1:4) '_' ID{id}(5) '/aparc+aseg_pt1_nat.nii' ],...
%         ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt1_nat.nii'])
%     copyfile(['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADB_preprocessed/RewardTask/' ID{id}(1:4) '_' ID{id}(5) '/aparc+aseg_pt2_nat.nii' ],...
%         ['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt2_nat.nii'])
% 
% end


%% brain masks

for id=18%1:length(ID)
    
    clear matlabbatch
    spm_jobman('initcfg')
    matlabbatch{1}.spm.util.imcalc.input = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt1_nat.nii,1']};
    matlabbatch{1}.spm.util.imcalc.output = ['BrainMask_pt1'];
    matlabbatch{1}.spm.util.imcalc.outdir = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/'  ID{id}]};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch)
    
    clear matlabbatch
    spm_jobman('initcfg')
    matlabbatch{1}.spm.util.imcalc.input = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/' ID{id} '/aparc+aseg_pt2_nat.nii,1']};
    matlabbatch{1}.spm.util.imcalc.output = ['BrainMask_pt2'];
    matlabbatch{1}.spm.util.imcalc.outdir = {['/Users/yeojin/Desktop/E_data/EA_raw/EAD_PET/EADD_segmented/coreg_roi/'  ID{id}]};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch)
    
end