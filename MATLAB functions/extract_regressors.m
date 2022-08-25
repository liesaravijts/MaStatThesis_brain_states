function extract_regressors(subject_list,file_id,params,multiple_runs,Data_dir)
    % extract_regressors Extracts the head movement and extra-cerebral component nuisance regressors.
    %   extract_regressors(subject_list,file_id,params,multiple_runs,Data_dir)
    %   
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_id        should be a 1-by-C character array, containing the file name pattern of the (untrimmed) functional images.
    %   params         should be a R-by-1 map container, containing at least:
    %                    the number of slices ('nslices'), the repetition time in seconds ('TR(s)'), the slice times 
    %                    in milliseconds ('sliceorder(idx/ms)'), the reference slice in milliseconds ('refslice(ms)'),
    %                    and the multiband acceleration factor ('mbaccfactor').
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to use all available runs.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   
    %   See also preprocess_structural, remove_first_five_volumes, preprocess_functional, extract_time_series.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(file_id) || ~isrow(file_id); error('file_id should be a 1-by-C character array.'); end
    if ~isa(params,'containers.Map'); error('params should be a R-by-1 map container.'); end
    if ~all(isKey(params,{'nslices','TR(s)','sliceorder(idx/ms)','refslice(idx/ms)','mbaccfactor'}))
        error('params should contain the keys: ''nslices'', ''TR(s)'', ''sliceorder(idx/ms)'', ''refslice(idx/ms)'' and ''mbaccfactor''.'); end
    if ~all(cellfun(@isrow,values(params))); error('All values in params should be scalars or row vectors.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    
    tic
    for i = 1:size(subject_list,1)
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        
        func = dir(fullfile(subject_folder,'**',['swua' file_id(1:end-4) '_trimmed.nii']));
        anat = dir(fullfile(subject_folder,'**','wc*.nii'));
        anat_preproc_folder = anat(1).folder;
        
        if multiple_runs; nruns = length(func); else; nruns = 1; end
        for j = 1:nruns
            func_file = func(j).name;
            func_preproc_folder = func(j).folder; func_folder = fileparts(func_preproc_folder);
            disp(['<strong>Extracting nuisance regressors: ' func_file '</strong>' ...
                ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
            
            % create the folder func/timeseries/run/regressors
            run_name = char(regexp(func_file,'task-.*bold','match'));
            regressors_folder = [func_folder '/timeseries/' run_name(6:end-5) '/regressors'];
            if ~exist(regressors_folder,'dir'); mkdir(regressors_folder); end
            
            %% DCT set
            nvols = params('nvols')-5;
            full_dct_set = spm_dctmtx(nvols,nvols);
            
            LL = 1/128; % lower limit Hz
            UL = 0.1; % upper limit Hz
            
            LL_component = fix(2*nvols*params('TR(s)')*LL+1);
            UL_component = fix(2*nvols*params('TR(s)')*UL+1);
            
            R = full_dct_set(:,(LL_component+1):(UL_component-1)); 
            save([regressors_folder '/dct.mat'],'R')
            
            %% head movement parameters (Friston-24)
            rp24 = load([func_preproc_folder '/rp_' func_file(4:end-4) '.txt']);
            rp24(:,7:12) = [zeros(1,6); rp24(1:end-1,:)]; % historical effects (1 time point)
            rp24(:,13:24) = rp24(:,1:12).^2; % squares
            
            % write to text file
            writematrix(rp24,[regressors_folder '/rp24.txt'],'delimiter','\t')
            
            %% WM and CSF signal
            [nifti_images,~]=spm_select('ExtFPList',func_preproc_folder,['^' func_file '$'],inf); % select .nii file (all frames)
            clear matlabbatch
            %% model specification
            matlabbatch{1}.spm.stats.fmri_spec.dir = {regressors_folder};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = params('TR(s)');
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = params('nslices')/params('mbaccfactor');
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = find(sort(params('sliceorder(idx/ms)'))==params('refslice(idx/ms)'),1,'last')/params('mbaccfactor');
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(nifti_images);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.6; % more lenient threshold
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            %% model estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            %% extract WM signal
            matlabbatch{3}.spm.util.voi.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{3}.spm.util.voi.adjust = NaN;
            matlabbatch{3}.spm.util.voi.session = 1;
            matlabbatch{3}.spm.util.voi.name = 'WM_signal';
            matlabbatch{3}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
            matlabbatch{3}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{3}.spm.util.voi.roi{2}.mask.image = {[anat_preproc_folder '/' anat(1).name]}; % WM mask
            matlabbatch{3}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            matlabbatch{3}.spm.util.voi.roi{3}.sphere.centre = [0 -27 -33];
            matlabbatch{3}.spm.util.voi.roi{3}.sphere.radius = 7;
            matlabbatch{3}.spm.util.voi.roi{3}.sphere.move.fixed = 1;
            matlabbatch{3}.spm.util.voi.expression = 'i1&i2&i3';
            %% extract CSF signal
            matlabbatch{4}.spm.util.voi.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{4}.spm.util.voi.adjust = NaN;
            matlabbatch{4}.spm.util.voi.session = 1;
            matlabbatch{4}.spm.util.voi.name = 'CSF_signal';
            matlabbatch{4}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
            matlabbatch{4}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{4}.spm.util.voi.roi{2}.mask.image = {[anat_preproc_folder '/' anat(2).name]}; % CSF mask
            matlabbatch{4}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            matlabbatch{4}.spm.util.voi.roi{3}.sphere.centre = [0 -40 -5];
            matlabbatch{4}.spm.util.voi.roi{3}.sphere.radius = 5;
            matlabbatch{4}.spm.util.voi.roi{3}.sphere.move.fixed = 1;
            matlabbatch{4}.spm.util.voi.expression = 'i1&i2&i3';
            %%
            spm_jobman('run',matlabbatch)
            
        end
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Extracting nuisance regressors: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
