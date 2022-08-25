function move_unselected(subject_list,Data_dir,quiet)
    % move_unselected Moves unselected participants to a separate folder Data_dir/not-selected.
    %   move_unselected(subject_list,Data_dir,quiet=true)
    %   
    %   subject_list  should be a R-by-C character array, containing R subject IDs of length C.
    %   Data_dir      should be a 1-by-C character array, containing the MRI data directory.
    %   quiet         should be a 1-by-1 logical value, indicating whether to silence the function (default: true).
    %   
    %   See also read_ids, sort_participants, set_parameters, unzip_all, unzip_files.
    
    % set default value
    if ~exist('quiet','var'); quiet = true; end
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    if ~islogical(quiet) || ~isscalar(quiet); error('quiet should be a 1-by-1 logical value (default: true).'); end
    
    % create folder
    if ~exist([Data_dir '/not-selected'],'dir'); mkdir([Data_dir '/not-selected']); end
    
    full_list = dir(fullfile(Data_dir,'sub-*'));
    if all(ismember(char(full_list.name),subject_list,'rows'))
            if ~quiet; disp('No folders were moved'); end
    else
        count = 0;
        for i = 1:length(full_list)
            subject_id = full_list(i).name;
            if ~ismember(subject_id,subject_list,'rows')
                subject_folder = [Data_dir '/' subject_id];
                movefile(subject_folder,[Data_dir '/not-selected'])
                if ~quiet; disp(['Moved ' subject_folder ' to ' Data_dir '/not-selected/' subject_id]); end
                count = count + 1;
            end
        end
        if ~quiet; if count ~= 1; disp([num2str(count) ' folders were moved']);
            else; disp('1 folder was moved'); end; end
    end
    
end
