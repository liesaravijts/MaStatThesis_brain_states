function subject_list=read_ids(file_path,delimiter)
    % read_ids Reads subject IDs from a text file to a character array.
    %   subject_list=read_ids(file_path,delimiter='\r\n')
    %   
    %   file_path     should be a 1-by-C character array, referring to a text file containing subject IDs.
    %   delimiter     should be a 1-by-C character array, specifying how the subject IDs are separated (default: '\r\n').
    %   subject_list  is a R-by-C character array, containing R subject IDs of length C.
    %   
    %   See also sort_participants, set_parameters, move_unselected, unzip_all, unzip_files.
    
    % set default value
    if ~exist('delimiter','var'); delimiter = '\r\n'; end
    
    % check input arguments
    if ~ischar(file_path) || ~isrow(file_path); error('file_path should be a 1-by-C character array.'); end
    if ~exist(file_path,'file'); error('The given file does not exist.'); end
    [~,~,fExt] = fileparts(file_path); if ~strcmpi(fExt,'.txt'); error('Unexpected file extension.'); end
    if ~ischar(delimiter) || ~isrow(delimiter); error('delimiter should be a 1-by-C character array (default: ''\r\n'').'); end
    
    subject_list = strsplit(fileread(file_path),delimiter); % read and split IDs to cell array
    subject_list = char(subject_list(~cellfun(@isempty,subject_list))); % remove empty cells and convert to character array
    
end
