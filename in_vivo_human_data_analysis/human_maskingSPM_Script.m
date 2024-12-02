%-----------------------------------------------------------------------
% Job saved on 25-May-2018 12:55:11 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
addpath('freesurfer/6.0.0beta/ubuntu-xenial-amd64/matlab')

spm('defaults','fmri');
spm_jobman('initcfg');
clear jobs
matlabbatch={};

subjects_dir = 'derivatives/MPMs_wOLS-true_rfcorr-percon/';
    
subj = {'sub-01', 'sub-02','sub-03', 'sub-04', 'sub-05', 'sub-06','sub-07', 'sub-08', 'sub-09', 'sub-10','sub-11', 'sub-12','sub-13', 'sub-14', 'sub-15'}

for j = 1:length(subj)    
mri = MRIread([subjects_dir subj{j} '-1mm/mri/T1.mgz']);
err = MRIwrite(mri,[subjects_dir subj{j} '-1mm/mri/T1.nii']);

matlabbatch{1}.spm.spatial.preproc.channel.vols = {[subjects_dir subj{j} '-1mm/mri/T1.nii,1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'Software/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'Software/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'Software/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'Software/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'Software/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'Software/spm12/tpm/TPM.nii,6'};
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
matlabbatch{2}.spm.util.imcalc.input = {
                                        [subjects_dir subj{j} '-1mm/mri/c1T1.nii,1']
                                        [subjects_dir subj{j} '-1mm/mri/c2T1.nii,1']
                                        [subjects_dir subj{j} '-1mm/mri/c3T1.nii,1']
                                        };
                                    
matlabbatch{2}.spm.util.imcalc.output = 'mask_spm';
matlabbatch{2}.spm.util.imcalc.outdir = {[subjects_dir subj{j} '-1mm/mri']};
matlabbatch{2}.spm.util.imcalc.expression = 'i1>0|i2>0|i3>0';
% matlabbatch{2}.spm.util.imcalc.expression = 'i1>0|i2>0';

matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
matlabbatch{3}.spm.util.imcalc.input = {
                                        [subjects_dir subj{j} '-1mm/mri/T1.nii,1']
                                        [subjects_dir subj{j} '-1mm/mri/mask_spm.nii,1']
                                        };
matlabbatch{3}.spm.util.imcalc.output = 'brainmask';
matlabbatch{3}.spm.util.imcalc.outdir = {[subjects_dir subj{j} '-1mm/mri']};
matlabbatch{3}.spm.util.imcalc.expression = 'i1.*i2';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
end
