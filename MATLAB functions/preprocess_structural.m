function preprocess_structural(subject_list,file_id,Data_dir,SPM_dir)
    % preprocess_structural Performs preprocessing of the structural images.
    %   preprocess_structural(subject_list,file_id,Data_dir,SPM_dir)
    %   
    %   subject_list  should be a R-by-C character array, containing R subject IDs of length C.
    %   file_id       should be a 1-by-C character array, containing the file name pattern of the structural images.
    %   Data_dir      should be a 1-by-C character array, containing the MRI data directory.
    %   SPM_dir       should be a 1-by-C character array, containing the SPM directory.
    %   
    %   See also remove_first_five_volumes, preprocess_functional, extract_regressors, extract_time_series.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(file_id) || ~isrow(file_id); error('file_id should be a 1-by-C character array.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    if ~ischar(SPM_dir) || ~isrow(SPM_dir); error('SPM_dir should be a 1-by-C character array.'); end
    if ~exist(SPM_dir,'dir'); error('The given SPM directory does not exist.'); end
    
    tic
    for i = 1:size(subject_list,1) % loop over subjects
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        
        anat = dir(fullfile(subject_folder,'**',file_id)); % find structural image
        anat_file = anat(1).name;
        anat_folder = anat(1).folder;
        disp(['<strong>Preprocessing structural image: ' anat_file '</strong>' ...
            ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
        
        [nifti_image,~]=spm_select('FPList',anat_folder,['^' anat_file '$']); % select .nii file
        clear matlabbatch
        %% segmentation
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {nifti_image};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[SPM_dir '/tpm/TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[SPM_dir '/tpm/TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[SPM_dir '/tpm/TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[SPM_dir '/tpm/TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[SPM_dir '/tpm/TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[SPM_dir '/tpm/TPM.nii,6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
        %% make skull-stripped bias-corrected image
        matlabbatch{2}.spm.util.imcalc.input(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{2}.spm.util.imcalc.input(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{2}.spm.util.imcalc.input(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{2}.spm.util.imcalc.input(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
        matlabbatch{2}.spm.util.imcalc.output = ['Skullstr_biascor_' anat_file(1:end-4)];
        matlabbatch{2}.spm.util.imcalc.outdir = {anat_folder};
        matlabbatch{2}.spm.util.imcalc.expression = 'i1.*(i2+i3+i4)';
        matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = 1;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
        %% write structural to MNI space
        matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep(['Image Calculator: ImCalc Computed Image: ' anat_folder '/Skullstr_biascor_' anat_file(1:end-4)], substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
        matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 
                                                                  78 76 85];
        matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % 1 mm MNI space
        matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %%
        spm_jobman('run',matlabbatch)
        
        % create preproc folder
        if ~exist([anat_folder '/preproc'],'dir'); mkdir([anat_folder '/preproc']); end
        
        % move files to preproc folder
        movefile([anat_folder '/' anat_file(1:end-4) '_seg8.mat'],[anat_folder '/preproc'])
        movefile([anat_folder '/BiasField_' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/m' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/c1' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/c2' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/c3' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/c4' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/c5' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/y_' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/Skullstr_biascor_' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/wSkullstr_biascor_' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/wc2' anat_file],[anat_folder '/preproc'])
        movefile([anat_folder '/wc3' anat_file],[anat_folder '/preproc'])
        
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Preprocessing structural images: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
