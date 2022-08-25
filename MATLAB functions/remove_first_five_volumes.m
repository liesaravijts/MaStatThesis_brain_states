function remove_first_five_volumes(subject_list,file_id,multiple_runs,Data_dir)
    % remove_first_five_volumes Removes the first five volumes of the functional scans.
    %   remove_first_five_volumes(subject_list,file_id,multiple_runs,Data_dir)
    %   
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_id        should be a 1-by-C character array, containing the file name pattern of the functional images.
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to use all available runs.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   
    %   See also preprocess_structural, preprocess_functional, extract_regressors, extract_time_series.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(file_id) || ~isrow(file_id); error('file_id should be a 1-by-C character array.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    
    tic
    for i = 1:size(subject_list,1)
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        
        func = dir(fullfile(subject_folder,'**',file_id)); % find functional images
        if multiple_runs; nruns = length(func); else; nruns = 1; end
        for j = 1:nruns
            func_file = func(j).name;
            func_folder = func(j).folder;
            disp(['<strong>Removing first five volumes: ' func_file '</strong>' ...
                ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
            
            [nifti_images,~]=spm_select('ExtFPList',func_folder,['^' func_file '$'],inf); % select .nii file (all frames)
            %%
            clear matlabbatch
            matlabbatch{1}.spm.util.cat.vols = cellstr(nifti_images(6:end,:));
            matlabbatch{1}.spm.util.cat.name = [func_folder '/' func_file(1:end-4) '_trimmed.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 0;
            spm_jobman('run',matlabbatch)
        end
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Removing first five volumes: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
