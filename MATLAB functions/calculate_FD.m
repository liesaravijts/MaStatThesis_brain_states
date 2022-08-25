function FD_summary=calculate_FD(subject_list,file_id,multiple_runs,Data_dir,Work_dir,write_output)
    % calculate_FD Computes the framewise displacement (FD). (Power et al., 2012).
    %   FD_summary=calculate_FD(subject_list,file_id,multiple_runs,Data_dir,Work_dir,write_output)
    %
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_id        should be a 1-by-C character array, containing the file name pattern of the functional images.
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to use all available runs.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   Work_dir       should be a 1-by-C character array, containing the working directory.
    %   write_output   should be a 1-by-1 logical value, indicating whether to write the results to Work_dir/Results.
    %
    % Based on https://github.com/CPernet/spmup/blob/master/QA/spmup_FD.m
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(file_id) || ~isrow(file_id); error('file_id should be a 1-by-C character array.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    if ~ischar(Work_dir) || ~isrow(Work_dir); error('Work_dir should be a 1-by-C character array.'); end
    if ~exist(Work_dir,'dir'); error('The given working directory does not exist.'); end
    if ~islogical(write_output) || ~isscalar(write_output); error('write should be a 1-by-1 logical value.'); end
    
    if write_output; fFD = fopen([Work_dir '/Results/FD_summary.csv'],'w'); % open csv file
        fprintf(fFD,'participant_id,meanFD,maxFD\r\n'); end % set header
    
    FD_summary = array2table(zeros(0,3),'VariableNames',{'participant_id','meanFD','maxFD'});
    for i = 1:size(subject_list,1)
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        
        func = dir(fullfile(subject_folder,'**',['rp_a' file_id(1:end-4) '_trimmed.txt']));
        if multiple_runs; nruns = length(func); else; nruns = 1; end
        for j = 1:nruns
            rp_file = func(j).name;
            func_preproc_folder = func(j).folder;
            disp(rp_file)
            
            rp = load([func_preproc_folder '/' rp_file]);
            rp = spm_detrend(rp,1);
            rp(:,4:6) = rp(:,4:6).*50; % rotational displacement in mm
            
            FD = [0; sum(abs(diff(rp)),2)];
            meanFD = mean(FD); maxFD = max(FD);
            if write_output; save([func_preproc_folder '/FD' rp_file(3:end-4) '.mat'],'FD','meanFD','maxFD'); end
            
            FD_summary = [FD_summary; {subject_id,meanFD,maxFD}]; %#ok<AGROW>
            if write_output; fprintf(fFD,[subject_id ',' num2str(meanFD) ',' num2str(maxFD) '\r\n']); end
        end
    end
    
    if write_output; fclose('all'); end
    if write_output; save([Work_dir '/Results/FD_summary.mat'],'FD_summary'); end % save as mat file
    
end
