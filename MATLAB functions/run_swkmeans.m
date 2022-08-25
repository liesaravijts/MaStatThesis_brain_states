function [spath,FC,state_dynamics,sumD,D,run_mintsumD]=run_swkmeans(subject_list,file_id,params,multiple_runs,K,Data_dir,Work_dir,write_output,nreps,window_length,step_size,data)
    % run_swkmeans Runs the sliding window with k-means clustering (SWK).
    %   [spath,FC,state_dynamics,sumD,D,run_mintsumD]=run_swkmeans(subject_list,file_id,params,multiple_runs,K,Data_dir,Work_dir,write_output,nreps=3,window_length=1/fmin,step_size=1(,data))
    %   
    %   subject_list   should be a R-by-C character array, containing R subject IDs of length C.
    %   file_id        should be a 1-by-C character array, containing the file name pattern of the (untrimmed) functional images.
    %   params         should be a R-by-1 map container, containing at least:
    %                    the repetition time in seconds ('TR(s)') and the number of volumes ('nvols').
    %   multiple_runs  should be a 1-by-1 logical value, indicating whether to use all available runs.
    %   K              should be a 1-by-1 numerical value, containing the maximum number of states to be inferred.
    %   Data_dir       should be a 1-by-C character array, containing the MRI data directory.
    %   Work_dir       should be a 1-by-C character array, containing the working directory.
    %   write_output   should be a 1-by-1 logical value, indicating whether to write all output to Work_dir/Results.
    %   nreps          should be a 1-by-1 numerical value (default: 3).
    %   window_length  should be a 1-by-1 numerical value, containing the window size in TRs (default: 1/(1/128)/params('TR(s)')).
    %   step_size      should be a 1-by-1 numerical value, containing the window offset in TRs (default: 1). 
    %   data           should be a R-by-C data matrix (optional).
    %   
    %   Makes use of the <a href="matlab:web('https://github.com/OHBA-analysis/HMM-MAR')">HMM-MAR Toolbox</a> (Vidaurre et al., 2016, 2017).
    %   
    %   See also run_hmm, run_gmm, get_state_dynamics.
    
    % set default values
    if ~exist('nreps','var'); nreps = 3; end
    if ~exist('window_length','var'); window_length = ceil(128/params('TR(s)')); end % smallest possible w following Leonardi et al. (2015)
    if ~mod(window_length,2); window_length = window_length + 1;
        if nargin > 9; warning(['Window length is set to ' num2str(window_length) '.']); end
    end
    if ~exist('step_size','var'); step_size = 1; end % default: 1 TR
    
    % check input arguments
    if ~ischar(subject_list) || ~ismatrix(subject_list); error('subject_list should be a R-by-C character array.'); end
    if ~ischar(file_id) || ~isrow(file_id); error('file_id should be a 1-by-C character array.'); end
    if ~isa(params,'containers.Map'); error('params should be a R-by-1 map container.'); end
    if ~all(isKey(params,{'TR(s)','nvols'}))
        error('params should contain the keys: ''TR(s)'' and ''nvols''.'); end
    if ~all(cellfun(@isrow,values(params))); error('All values in params should be scalars or row vectors.'); end
    if ~islogical(multiple_runs) || ~isscalar(multiple_runs); error('multiple_runs should be a 1-by-1 logical value.'); end
    if ~isnumeric(K) || ~isscalar(K); error('nstates should be a 1-by-1 numerical value'); end
    if ~ischar(Data_dir) || ~isrow(Data_dir); error('Data_dir should be a 1-by-C character array.'); end
    if ~exist(Data_dir,'dir'); error('The given data directory does not exist.'); end
    if ~ischar(Work_dir) || ~isrow(Work_dir); error('Work_dir should be a 1-by-C character array.'); end
    if ~exist(Work_dir,'dir'); error('The given working directory does not exist.'); end
    if ~islogical(write_output) || ~isscalar(write_output); error('write_output should be a 1-by-1 logical value.'); end
    if ~isnumeric(nreps) || ~isscalar(nreps) || mod(nreps,1)~=0 || nreps<=0; error('nreps should be a whole number > 0.'); end
    if ~isnumeric(window_length) || ~isscalar(window_length) || mod(window_length,1)~=0 || window_length<=0; error('window_size should be a whole number > 0.'); end
    if ~isnumeric(step_size) || ~isscalar(step_size) || mod(step_size,1)~=0 || step_size<=0; error('window_offset should be a whole number > 0.'); end
    if exist('data','var') && (~isnumeric(data) || ~ismatrix(data)); error('data should be a R-by-C numerical array.'); end
    
    tic
    N = num2str(size(subject_list,1));
    nTR_scan = (params('nvols')-5); dur_scan = nTR_scan * params('TR(s)');
    nTR_total = nTR_scan * size(subject_list,1); dur_total = nTR_total * params('TR(s)');
    
    disp(['<strong>Running sliding window with K-means clustering</strong> (subjects = ' num2str(N) ', repetitions = ' num2str(nreps) ')'])
    disp(['Number of samples/TRs per subject: ' num2str(nTR_scan) ' TRs (' num2str(dur_scan) ' s)'])
    disp(['Total length of the time series: ' num2str(nTR_total) ' TRs (' num2str(dur_total) ' s)'])
    disp(['Window size: ' num2str(window_length) ' TRs (' num2str(window_length * params('TR(s)')) ' s)'])
    disp(['Step size: ' num2str(step_size) ' TRs (' num2str(step_size * params('TR(s)')) ' s)'])
    
    %% initialise data files and load data
    init_data_files
    
    % load data (part adapted from checkdatacell.m)
    if ~exist('data','var')
        t = 0;
        for i = 1:length(f) %#ok<*USENS>
            load(f{i},'ROI_data')
            if i==1; data = zeros(sum(cell2mat(T)),size(ROI_data,2)); end
            data((1:size(ROI_data,1)) + t,:) = ROI_data; t = t + size(ROI_data,1);
        end
    end
    
    %% standardise ROI data for each subject
    % using standardisedata.m function from HMM-MAR toolbox
    sddata = standardisedata(data,cell2mat(T),1); 
    
    %% run sliding window + K-means    
    W = nTR_total/step_size - window_length/step_size + 1; % total number of windows
    R = size(sddata,2); Rsq = R^2; % number of elements in one corr matrix (R^2 with R = number of ROIs)
    
    [spath,sumD,D,FC,state_dynamics] = deal(cell(nreps,1));
    for rep = 1:nreps
        disp(['<strong>SWK run ' num2str(rep) '/' num2str(nreps) '</strong>'])
        
        corrmat_all = zeros(W,Rsq); % initialise matrix
        t = 0;
        for w = 1:W % sliding window
            % break loop when end of data is reached
            if t + window_length > size(sddata,1); continue; end 
            
            window = sddata((1:window_length) + t,:); % window
            corrmat_window = corr(window); % correlation matrix
            corrmat_all(w,:) = reshape(corrmat_window,1,Rsq); % reshape to vector and add
            t = t + step_size;
        end
        
        % k-means clustering
        [spath{rep},C,sumD{rep},D{rep}] = kmeans(corrmat_all,K,'Distance','cityblock','MaxIter',250,'Display','Iter');
        
        % get FC matrix for each state (correlation matrices)
        FC{rep} = zeros(R,R,K);
        for s = 1:K % loop over states
            %covmat_state = reshape(C(s,:),size(sddata,2),size(sddata,2));
            %FC{rep}(:,:,s) = corrcov(covmat_state);
            FC{rep}(:,:,s) = reshape(C(s,:),size(sddata,2),size(sddata,2));
        end
        
        % state course
        spath{rep} = repelem(spath{rep},step_size); % account for step size
        lost = (window_length-1)/2;
        
        % get state dynamics
        state_dynamics{rep} = get_state_dynamics(T(2:end-1),spath{rep}(nTR_scan-lost+1:end-(nTR_scan-lost))); %#ok<FNCOLND>
        
        % save results
        if write_output
            SWK_results = struct('rep',rep,'files',char(f),'T',cell2mat(T),'window_length',window_length,...
                                 'step_size',step_size,'K',K,'FC',FC{rep},'spath',spath{rep},...
                                 'totsumD',sum(sumD{rep}),'D',D{rep},'state_dynamics',state_dynamics{rep});
            save([Work_dir '/Results/SWK_K' num2str(K) '_rep' num2str(rep,'%02.f') '.mat'],'SWK_results');
        end
    end
    
    %% get run(s) with minimum cost (total sum of distances)
    for rep = 1:nreps
        if rep == 1; minD = sum(sumD{rep}); run_mintsumD = 1;
        elseif sum(sumD{rep}) < minD; minD = sum(sumD{rep}); run_mintsumD = rep;
        elseif sum(sumD{rep}) == minD; run_mintsumD(length(run_mintsumD) + 1) = rep; %#ok<AGROW>
        end
    end
    
    %%
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Running sliding window with K-means clustering: Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
