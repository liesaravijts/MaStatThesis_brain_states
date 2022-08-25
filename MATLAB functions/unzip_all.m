function unzip_all(subject_list,Data_dir)
    % unzip_all Unpacks all .nii.gz files found in Data_dir for all subjects in subject_list.
    %   unzip_all(subject_list,Data_dir)
    %   
    %   subject_list  should be a R-by-C character array, containing R subject IDs of length C.
    %   Data_dir      should be a 1-by-C character array, containing the MRI data directory.
    %   
    %   See also read_ids, sort_participants, set_parameters, move_unselected, unzip_files.
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    
    tic
    for i = 1:size(subject_list,1) % loop over subjects
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        disp(['<strong>Unzipping files: ' subject_id '</strong>' ...
            ' (' num2str(i) '/' num2str(size(subject_list,1)) '...)'])
        
        files = dir(fullfile(subject_folder,'**','*.nii.gz')); % search for all .nii.gz files
        for j = 1:length(files) % loop over .nii.gz files
            file_name = files(j).name;
            file_folder = files(j).folder; % subfolder
            
            gunzip([file_folder '/' file_name])
            disp(file_name(1:end-3)) % display file name when unzipped
        end
    end
    
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Unzipping files: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
