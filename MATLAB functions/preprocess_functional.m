function preprocess_functional(subject_list,file_ids,params,multiple_runs,Data_dir,SPM_dir)
    % preprocess_functional Performs preprocessing of the (trimmed) functional images.
    %   preprocess_functional(subject_list,file_ids,params,multiple_runs,Data_dir,SPM_dir)
    %   
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_ids       should be a 1-by-1 structure, containing the file name pattern of the (untrimmed) functional images and 
    %                  the file name patterns of the magnitude and phase difference images.
    %   params         should be a 1-by-1 structure as specified in main_script.m, with
    %                   params.func  should be a R-by-1 map container, containing at least:
    %                     the number of slices ('nslices'), the repetition time in seconds ('TR(s)'), the slice order
    %                     indexes or milliseconds ('sliceorder(idx/ms)'), the reference slice index or milliseconds ('refslice(idx/ms)'),
    %                     the blip direction ('blipdir(1/-1)'), the number of echoes ('nechoes'), the echo spacing
    %                     in milliseconds ('ESP(ms)') and the phase acceleration factor ('phsaccfactor').
    %                   params.fmap  should be a R-by-1 map container, containing at least:
    %                     the short echo time in milliseconds('shortTE(ms)') and the long echo time in milliseconds ('longTE(ms)').
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to use all available runs.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   SPM_dir        should be a 1-by-C character array, containing the SPM directory.
    %   
    %   See also preprocess_structural, remove_first_five_volumes, extract_regressors, extract_time_series.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~isstruct(file_ids) || ~isscalar(file_ids); error('file_ids should be a 1-by-1 structure.'); end
    if ~isstruct(params) || ~isscalar(params); error('params should be a 1-by-1 structure.'); end
    if ~isa(params.func,'containers.Map'); error('params.func should be a R-by-1 map container.'); end
    if ~all(isKey(params.func,{'nslices', 'TR(s)', 'sliceorder(idx/ms)', 'refslice(idx/ms)', 'blipdir(1/-1)', 'nechoes', 'ESP(ms)', 'phsaccfactor'})) 
        error('params.func should contain the keys: ''nslices'', ''TR(s)'', ''sliceorder(idx/ms)'', ''refslice(idx/ms)'', ''blipdir(1/-1)'', ''nechoes'', ''ESP(ms)'' and ''phsaccfactor''.'); end
    if ~all(cellfun(@isrow,values(params.func))); error('All values in params.func should be scalars or row vectors.'); end
    if ~isa(params.fmap,'containers.Map'); error('params.fmap should be a R-by-1 map container.'); end
    if ~all(isKey(params.fmap,{'shortTE(ms)', 'longTE(ms)'})) 
        error('params.fmap should contain the keys: ''shortTE(ms)'', ''longTE(ms)''.'); end
    if ~all(cellfun(@isrow,values(params.fmap))); error('All values in params.fmap should be scalars or row vectors.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    if ~ischar(SPM_dir) || ~isrow(SPM_dir); error('SPM_dir should be a 1-by-C character array.'); end
    if ~exist(SPM_dir,'dir'); error('The given SPM directory does not exist.'); end
    
    tic
    blipdir = params.func('blipdir(1/-1)');
    for i = 1:size(subject_list,1)
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        
        func = dir(fullfile(subject_folder,'**',[file_ids.func(1:end-4) '_trimmed.nii'])); % find trimmed images
        fmap_magn = dir(fullfile(subject_folder,'**',file_ids.fmap.magn)); % find magnitude and phase diff images
        fmap_phsd = dir(fullfile(subject_folder,'**',file_ids.fmap.phsd));
        fmap_folder = fmap_magn(1).folder;
        anat = dir(fullfile(subject_folder,'**','y_*.nii')); % find forward deformation field
        anat_preproc_folder = anat(1).folder;
        
        if multiple_runs; nruns = length(func); else; nruns = 1; end
        for j = 1:nruns
            func_file = func(j).name;
            func_folder = func(j).folder;
            disp(['<strong>Preprocessing functional images: ' func_file '</strong>' ...
                ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
            disp(['   Phase difference image: ' fmap_phsd(1).name])
            disp(['   Magnitude image: ' fmap_magn(1).name])
            
            [nifti_images,~]=spm_select('ExtFPList',func_folder,['^' func_file '$'],inf); % select .nii file (all frames)
            clear matlabbatch
            %% slice time correction
            matlabbatch{1}.spm.temporal.st.scans{1} = cellstr(nifti_images);
            matlabbatch{1}.spm.temporal.st.nslices = params.func('nslices');
            matlabbatch{1}.spm.temporal.st.tr = params.func('TR(s)');
            matlabbatch{1}.spm.temporal.st.ta = params.func('TR(s)')-params.func('TR(s)')/params.func('nslices');
            matlabbatch{1}.spm.temporal.st.so = params.func('sliceorder(idx/ms)');
            matlabbatch{1}.spm.temporal.st.refslice = params.func('refslice(idx/ms)');
            matlabbatch{1}.spm.temporal.st.prefix = 'a';
            %% calculate VDM
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {[fmap_folder '/' fmap_phsd(1).name]};
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {[fmap_folder '/' fmap_magn(1).name]};
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [params.fmap('shortTE(ms)') params.fmap('longTE(ms)')];
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdir(j);
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = params.func('nechoes')*params.func('ESP(ms)')/params.func('phsaccfactor');
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {[SPM_dir '/toolbox/FieldMap/T1.nii']};
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.session.epi(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
            matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
            %% spatial realignment and unwarping
            matlabbatch{3}.spm.spatial.realignunwarp.data.scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
            matlabbatch{3}.spm.spatial.realignunwarp.data.pmscan(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.sep = 4;
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.rtm = 0;
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.einterp = 2;
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
            matlabbatch{3}.spm.spatial.realignunwarp.eoptions.weight = '';
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.jm = 0;
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.sot = [];
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.rem = 1;
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.noi = 5;
            matlabbatch{3}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
            matlabbatch{3}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
            matlabbatch{3}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
            matlabbatch{3}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
            matlabbatch{3}.spm.spatial.realignunwarp.uwroptions.mask = 1;
            matlabbatch{3}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
            %% co-register to structural image
            matlabbatch{4}.spm.spatial.coreg.estimate.ref = {[anat_preproc_folder '/Skullstr_biascor_' anat(1).name(3:end)]};
            matlabbatch{4}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
            matlabbatch{4}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            %% write functional scans to MNI space
            matlabbatch{5}.spm.spatial.normalise.write.subj.def = {[anat_preproc_folder '/' anat(1).name]};
            matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
            matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                      78 76 85];
            matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2]; % 2 mm MNI space
            matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
            %% smoothing
            matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
            matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];
            matlabbatch{6}.spm.spatial.smooth.dtype = 0;
            matlabbatch{6}.spm.spatial.smooth.im = 0;
            matlabbatch{6}.spm.spatial.smooth.prefix = 's';
            %%
            spm_jobman('run',matlabbatch)
            
            % create preproc folders
            if ~exist([func_folder '/preproc'],'dir'); mkdir([func_folder '/preproc']); end
            if ~exist([fmap_folder '/preproc'],'dir'); mkdir([fmap_folder '/preproc']); end
            
            % move files to preproc folders
            movefile([func_folder '/a' func_file],[func_folder '/preproc'])
            movefile([func_folder '/a' func_file(1:end-3) 'mat'],[func_folder '/preproc'])
            movefile([func_folder '/a' func_file(1:end-4) '_uw.mat'],[func_folder '/preproc'])
            movefile([func_folder '/rp_a' func_file(1:end-3) 'txt'],[func_folder '/preproc'])
            movefile([func_folder '/wfmag_a' func_file],[func_folder '/preproc'])
            movefile([func_folder '/meanua' func_file],[func_folder '/preproc'])
            movefile([func_folder '/ua' func_file],[func_folder '/preproc'])
            movefile([func_folder '/ua' func_file(1:end-3) 'mat'],[func_folder '/preproc'])
            movefile([func_folder '/wmeanua' func_file],[func_folder '/preproc'])
            movefile([func_folder '/wua' func_file],[func_folder '/preproc'])
            movefile([func_folder '/swmeanua' func_file],[func_folder '/preproc'])
            movefile([func_folder '/swua' func_file],[func_folder '/preproc'])
            
            movefile([fmap_folder '/sc' fmap_phsd(1).name],[fmap_folder '/preproc'])
            movefile([fmap_folder '/m' fmap_magn(1).name],[fmap_folder '/preproc'])
            movefile([fmap_folder '/bmask' fmap_magn(1).name],[fmap_folder '/preproc'])
            movefile([fmap_folder '/fpm_sc' fmap_phsd(1).name],[fmap_folder '/preproc'])
            movefile([fmap_folder '/vdm5_sc' fmap_phsd(1).name],[fmap_folder '/preproc'])
            
        end
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Preprocessing functional images: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
