function [LEMON,NC,MPILMBB]=sort_participants(n,Data_dir,Work_dir,write)
    % sort_participants Sorts the subject IDs of the MPI-Leipzig MBB dataset into three character arrays
    % according to the different protocols.
    %   [LEMON,NC,MPILMBB]=sort_participants(n,Data_dir,Work_dir,write=false)
    %   
    %   n         should be a 1-by-1 numerical value, containing the largest subject number.
    %   Data_dir  should be a 1-by-C character array, containing the MRI data directory.
    %   Work_dir  should be a 1-by-C character array, containing the working directory.
    %   write     should be a 1-by-1 logical value, indicating whether to write LEMON and NC to text files (default: false).
    %   LEMON, NC and MPILMBB are R-by-10 character arrays, containing R subject IDs of length 10 (format: 'sub-010xxx').
    %   
    %   LEMON   - 228 subject IDs included in the LEMON protocol (ses-01).
    %   NC      - 199 subject IDs included in the N&C protocol (ses-02).
    %   MPILMBB - all 318 subject IDs in the dataset (109 took part in both protocols).
    %   
    %   If write = true, the arrays are written to Work_dir/sub_ids_LEMON.txt, Work_dir/sub_ids_N&C.txt and Work_dir/sub_ids_MPILMBB.txt.
    %   
    %   See also read_ids, set_parameters, move_unselected, unzip_all, unzip_files.
    
    % set default value
    if ~exist('write','var'); write = false; end
    
    % check input arguments
    if ~isnumeric(n) || ~isscalar(n) || mod(n,1)~=0 || n<=0; error('n should be a whole number > 0.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    if ~ischar(Work_dir) || ~isrow(Work_dir); error('Work_dir should be a 1-by-C character array.'); end
    if ~exist(Work_dir,'dir'); error('The given working directory does not exist.'); end
    if ~islogical(write) || ~isscalar(write); error('write should be a 1-by-1 logical value (default: false).'); end
    
    if write % open files
        fLEMON = fopen([Work_dir '/sub_ids_LEMON.txt'],'w');
        fNC = fopen([Work_dir '/sub_ids_N&C.txt'],'w');
        fMPILMBB = fopen([Work_dir '/sub_ids_MPILMBB.txt'],'w');
    end
    
    % initialise arrays
    LEMON = char(zeros(n,10));
    NC = char(zeros(n,10));
    MPILMBB = char(zeros(n,10));
    
    for subject = 1:n % loop over all possible IDs
        id = ['sub-010' sprintf('%03d',subject)];
        folder = [Data_dir '/' id];
        
        if exist(folder,'dir') % sub-010013, sub-010095 and sub-010248 don't exist
            if exist([folder '/ses-01'],'dir'); LEMON(subject,:) = id;
                if write; fprintf(fLEMON,[id '\r\n']); end; end
            if exist([folder '/ses-02'],'dir'); NC(subject,:) = id;
                if write; fprintf(fNC,[id '\r\n']); end; end
            MPILMBB(subject,:) = id;
            if write; fprintf(fMPILMBB,[id '\r\n']); end
        end
    end
    
    if write; fclose('all'); end % close files
    
    % remove empty rows
    LEMON = LEMON(~cellfun(@isempty,cellstr(LEMON)),:);
    NC = NC(~cellfun(@isempty,cellstr(NC)),:);
    MPILMBB = MPILMBB(~cellfun(@isempty,cellstr(MPILMBB)),:);
    
end
