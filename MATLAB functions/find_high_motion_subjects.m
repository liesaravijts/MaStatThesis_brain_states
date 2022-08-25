function high_motion_subjects=find_high_motion_subjects(subject_list,file_id,Data_dir,FD_threshold,cutoff)
    % find_high_motion_subjects Finds subjects with high motion as defined by FD_threshold and cutoff.
    %   high_motion_subjects=find_high_motion_subjects(subject_list,file_id,Data_dir,FD_threshold,cutoff)
    
    % set default values
    if ~exist('FD_threshold','var'); FD_threshold = 0.2; end
    if ~exist('cutoff','var'); cutoff = 0.25; end
    
    high_motion = boolean(zeros(size(subject_list,1),1));
    for i = 1:size(subject_list,1)
        subject_id = subject_list(i,:);
        subject_folder = [Data_dir '/' subject_id];
        func = dir(fullfile(subject_folder,'**',['FD_a' file_id(1:end-4) '_trimmed.mat']));
        
        load([func(1).folder '/' func(1).name],'FD','meanFD','maxFD')
        for j = 1:length(FD_threshold)
            if high_motion(i) ~= 1
                if isnumeric(cutoff{j}); high_motion(i) = sum(FD > FD_threshold(j))/length(FD) > cutoff{j};
                elseif ischar(cutoff{j}) && strcmp(cutoff{j},'mean'); high_motion(i) = meanFD > FD_threshold(j);
                elseif ischar(cutoff{j}) && strcmp(cutoff{j},'max'); high_motion(i) = maxFD > FD_threshold(j); 
                end
            end
        end
    end
    
    high_motion_subjects = subject_list(high_motion == 1,:);
    
end
