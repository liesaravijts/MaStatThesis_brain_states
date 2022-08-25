function unzip_files(subject_list,file_ids,multiple_runs,Data_dir)
    % unzip_files Unpacks the .nii.gz files found by file_ids for all subjects in subject_list.
    %   unzip_files(subject_list,file_ids,multiple_runs,Data_dir)
    %   
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_ids       should be a 1-by-1 structure, containing the file name patterns as specified in main_script.m.
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to unpack all available runs.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   
    %   See also read_ids, sort_participants, set_parameters, move_unselected, unzip_all.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~isstruct(file_ids) || ~isscalar(file_ids); error('file_ids should be a 1-by-1 structure as specified in main_script.m.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    
    tic
    for i = 1:size(subject_list,1) % loop over subjects
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        disp(['<strong>Unzipping files: ' subject_id '</strong>' ...
            ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
        
        % find .nii.gz files
        anat = dir(fullfile(subject_folder,'**',[file_ids.anat '.gz']));
        func = dir(fullfile(subject_folder,'**',[file_ids.func '.gz']));
        fmap_magn = dir(fullfile(subject_folder,'**',[file_ids.fmap.magn '.gz']));
        fmap_phsd = dir(fullfile(subject_folder,'**',[file_ids.fmap.phsd '.gz']));
        
        % create cell array containing the file paths
        if multiple_runs == false
            files = {[anat(1).folder '/' anat(1).name]...
                     [func(1).folder '/' func(1).name]...
                     [fmap_magn(1).folder '/' fmap_magn(1).name]...
                     [fmap_phsd(1).folder '/' fmap_phsd(1).name]};
        else
            files = [{[anat(1).folder '/' anat(1).name]}...
                     strcat({func.folder},'/',{func.name})...
                     strcat({fmap_magn.folder},'/',{fmap_magn.name})...
                     strcat({fmap_phsd.folder},'/',{fmap_phsd.name})];
        end
        %%
        clear matlabbatch
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files = files'; % input column vector
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir = {''};
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep = true;
        spm_jobman('run',matlabbatch)
        
        % display file names when unzipped
        disp(char(cellfun(@(x) x(find(x=='/',1,'last')+1:end-3),files,'UniformOutput',false)))
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Unzipping files: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
