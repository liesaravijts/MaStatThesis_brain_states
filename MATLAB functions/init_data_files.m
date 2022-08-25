% Initialise data files as is required for the HMM-MAR toolbox

N = size(subject_list,1);

% f{i} contains a list of file names (data file for each subject)
% T{i} contains the number of time points (per run)
f = cell(N,1); T = cell(N,1);
for i = 1:N
    subject_id = subject_list(i,:);
    func = dir(fullfile([Data_dir '/' subject_id],'**',file_id));
    if multiple_runs % temporally concatenated runs (from concat_ROI_data.m)
        f{i} = [Data_dir '/' subject_id '/ROI_data_concat.mat'];
        T{i} = repelem(params('nvols')-5,length(func));
    else
        run_name = char(regexp(func(1).name,'task.*bold','match'));
        f{i} = [func(1).folder '/timeseries/' run_name(6:end-5) '/ROI_data.mat'];
        T{i} = params('nvols')-5;
    end
end
