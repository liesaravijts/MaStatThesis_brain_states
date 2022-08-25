function [gmm,Gamma,spath,state_dynamics,GammaInit,fehist,run_minfe]=run_gmm(subject_list,file_id,params,multiple_runs,K,Data_dir,Work_dir,write_output,nreps,data)
    % run_gmm Runs the Gaussian mixture model (GMM).
    %   [gmm,Gamma,spath,state_dynamics,GammaInit,fehist,run_minfe]=run_gmm(subject_list,file_id,params,multiple_runs,K,Data_dir,Work_dir,write_output,nreps=3(,data))
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
    %   data           should be a R-by-C data matrix (optional).
    %   
    %   Makes use of the <a href="matlab:web('https://github.com/OHBA-analysis/HMM-MAR')">HMM-MAR Toolbox</a> (Vidaurre et al., 2016, 2017).
    %   
    %   See also run_hmm, run_swkmeans, get_state_dynamics.
    %
    % Reference: Vidaurre et al. (2017) Brain network dynamics are hierarchically organized in time, PNAS 114(48):12827–12832.
    
    % set default value
    if ~exist('nreps','var'); nreps = 3; end
    
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
    if ~isnumeric(nreps) || ~isscalar(nreps) || mod(nreps,1)~=0 || nreps<=0; error('n should be a whole number > 0.'); end
    if exist('data','var') && (~isnumeric(data) || ~ismatrix(data)); error('data should be a R-by-C numerical array.'); end
    
    tic
    disp(['<strong>Running Gaussian mixture model (GMM)</strong> (subjects = ' num2str(size(subject_list,1)) ', repetitions = ' num2str(nreps) ')'])
    
    %% initialise data files
    init_data_files
    
    if exist('data','var'); f = data; T = cell2mat(T); end %#ok<NODEF>
    
    %% set options
    options = struct();
    % model specification
    options.K = K; % maximum number of states to infer
    options.order = 0; % no autoregressive components
    options.zeromean = 1; % do not model the mean
    options.covtype = 'full'; % full covariance matrix
    options.id_mixture = 1; % estimate GMM (ignore temporal structure)
    % data options
    options.Fs = 1/params('TR(s)'); % sampling frequency
    options.standardise = 1;
    % training options
    options.useParallel = 1;
    options.inittype = 'hmmmar'; % default
    options.initrep = 5; % default
    options.initcyc = 25; % default
    options.cyc = 500; % default
    options.dropstates = 1;
    % display options
    options.verbose = 1; % show algorithm progress
    options.plotGamma = 0; % show graphical progress
    
    %% run GMM
    [gmm,Gamma,spath,GammaInit,fehist,state_dynamics] = deal(cell(nreps,1));
    for rep = 1:nreps
        disp(['<strong>GMM run ' num2str(rep) '/' num2str(nreps) '</strong>'])
        % run the model
        [gmm{rep},Gamma{rep},~,~,GammaInit{rep},~,fehist{rep}] = hmmmar(f,T,options);
        
        % FC matrix for each state (correlation matrices)
        for s = 1:gmm{rep}.K % loop over states
            [~,gmm{rep}.state(s).FC] = getFuncConn(gmm{rep},s);
        end
        
        % state course
        spath{rep} = zeros(length(Gamma{rep}),1);
        for t = 1:length(Gamma{rep})
            [~,spath{rep}(t)] = max(Gamma{rep}(t,:));
        end
        
        % state dynamics
        state_dynamics{rep} = get_state_dynamics(T,spath{rep},Gamma{rep});
        
        % save results
        if write_output
            GMM_results = struct('rep',rep,'files',char(f),'T',cell2mat(T),'gmm',gmm{rep},'Gamma',Gamma{rep},'spath',spath{rep},...
                                 'GammaInit',GammaInit{rep},'fehist',fehist{rep},'state_dynamics',state_dynamics{rep});
            save([Work_dir '/Results/GMM_K' num2str(K) '_rep' num2str(rep,'%02.f') '.mat'],'GMM_results');
        end
    end
    
    %% get run(s) with the least free energy
    for rep = 1:nreps
        if rep == 1; minfe = fehist{rep}(end); run_minfe = 1;
        elseif fehist{rep}(end) < minfe; minfe = fehist{rep}(end); run_minfe = rep;
        elseif fehist{rep}(end) == minfe; run_minfe(length(run_minfe) + 1) = rep; %#ok<AGROW>
        end
    end
    
    %%
    endtime = seconds(toc); endtime.Format = 'hh:mm:ss';
    disp(['<strong>Running Gaussian mixture model (GMM): Done</strong>  (Time elapsed: ' char(endtime) ')'])
    
end
