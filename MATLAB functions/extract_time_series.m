function extract_time_series(subject_list,file_id,params,multiple_runs,ROI_coords,Data_dir)
    % extract_time_series Extracts the BOLD time series from the given ROIs.
    %   extract_time_series(subject_list,file_id,params,multiple_runs,ROI_coords,Data_dir)
    %   
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_id        should be a 1-by-C character array, containing the file name pattern of the (untrimmed) functional images.
    %   params         should be a R-by-1 map container, containing at least:
    %                    the number of slices ('nslices'), the repetition time in seconds ('TR(s)'), the slice times 
    %                    in milliseconds ('sliceorder(idx/ms)'), the reference slice in milliseconds ('refslice(idx/ms)'),
    %                    and the multiband acceleration factor ('mbaccfactor').
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to use all available runs.
    %   ROI_coords     should be a R-by-2 cell array, containing the ROI names and centroid coordinates.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   
    %   See also preprocess_structural, remove_first_five_volumes, preprocess_functional, extract_regressors.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(file_id) || ~isrow(file_id); error('file_id should be a 1-by-C character array.'); end
    if ~isa(params,'containers.Map'); error('params should be a R-by-1 map container.'); end
    if ~all(isKey(params,{'nslices','TR(s)','sliceorder(idx/ms)','refslice(idx/ms)','mbaccfactor','nvols'}))
        error('params should contain the keys: ''nslices'', ''TR(s)'', ''sliceorder(idx/ms)'', ''refslice(idx/ms)'', ''mbaccfactor'' and ''nvols''.'); end
    if ~all(cellfun(@isrow,values(params))); error('All values in params should be scalars or row vectors.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~iscell(ROI_coords) || ~ismatrix(ROI_coords); error('ROI_coords should be a R-by-2 cell matrix.'); end
    if ~all(cellfun(@(x) ischar(x) & isrow(x),ROI_coords(:,1))); error('The 1st column of ROI_coords should contain 1-by-C character arrays.'); end
    if ~all(cellfun(@(x) isnumeric(x) & isrow(x),ROI_coords(:,2))); error('The 2nd column of ROI_coords should contain 1-by-3 numerical arrays.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    
    tic
    for i = 1:size(subject_list,1)
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        
        func = dir(fullfile(subject_folder,'**',['swua' file_id(1:end-4) '_trimmed.nii']));
        if multiple_runs; nruns = length(func); else; nruns = 1; end
        for j = 1:nruns
            func_file = func(j).name;
            func_preproc_folder = func(j).folder; func_folder = fileparts(func_preproc_folder);
            disp(['<strong>Extracting regional time series: ' func_file '</strong>' ...
                ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
            
            run_name = char(regexp(func_file,'task-.*bold','match'));
            timeseries_folder = [func_folder '/timeseries/' run_name(6:end-5)];
            load([timeseries_folder '/regressors/dct.mat'],'R')
             
            [nifti_images,~]=spm_select('ExtFPList',func_preproc_folder,['^' func_file '$'],inf); % select .nii file
            clear matlabbatch
            %% model specification
            matlabbatch{1}.spm.stats.fmri_spec.dir = {timeseries_folder};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = params('TR(s)');
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = params('nslices')/params('mbaccfactor');
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = find(sort(params('sliceorder(idx/ms)'))==params('refslice(idx/ms)'),1,'last')/params('mbaccfactor');
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(nifti_images);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
                                                                 [timeseries_folder '/regressors/dct.mat']
                                                                 [timeseries_folder '/regressors/rp24.txt']
                                                                 [timeseries_folder '/regressors/VOI_WM_signal_1.mat']
                                                                 [timeseries_folder '/regressors/VOI_CSF_signal_1.mat']
                                                                 };
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            %% model estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            %% define contrasts
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'EOI'; % effects of interest
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(size(R,2));
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.delete = 0;
            %% extract regional time series
            for k = 1:length(ROI_coords)
                matlabbatch{3+k}.spm.util.voi.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                matlabbatch{3+k}.spm.util.voi.adjust = 1;
                matlabbatch{3+k}.spm.util.voi.session = 1;
                matlabbatch{3+k}.spm.util.voi.name = char(ROI_coords(k,1));
                matlabbatch{3+k}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
                matlabbatch{3+k}.spm.util.voi.roi{1}.mask.threshold = 0.5;
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.contrast = 1;
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.conjunction = 1;
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.threshdesc = 'none';
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.thresh = 0.001;
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.extent = 0;
                matlabbatch{3+k}.spm.util.voi.roi{2}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
                matlabbatch{3+k}.spm.util.voi.roi{3}.sphere.centre = cell2mat(ROI_coords(k,2));
                matlabbatch{3+k}.spm.util.voi.roi{3}.sphere.radius = 8;
                matlabbatch{3+k}.spm.util.voi.roi{3}.sphere.move.fixed = 1;
                matlabbatch{3+k}.spm.util.voi.expression = 'i1&i2&i3';
            end
            %%
            spm_jobman('run',matlabbatch)
            
            %% save ROI data
            ROI_data = zeros(params('nvols')-5,length(ROI_coords));
            for k = 1:length(ROI_coords)
                load([timeseries_folder '/VOI_' char(ROI_coords(k,1)) '_1.mat'],'Y')
                ROI_data(:,k) = Y;
            end
            save([timeseries_folder '/ROI_data.mat'],'ROI_data')
            
        end
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Extracting regional time series: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
