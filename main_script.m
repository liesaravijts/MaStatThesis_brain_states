% Main script for "Modelling dynamical changes in brain states with resting-state fMRI: A comparison of methods"
% Author: Liesa Ravijts
%
%
%%% MRI data: MPI-Leipzig Mind-Brain-Body dataset
% https://openneuro.org/datasets/ds000221/versions/1.0.0
% The dataset consists of two protocols:
% LEMON: 1 resting-state run,  ses-01 subfolder, n=228
%        Babayan et al., 2019: https://www.nature.com/articles/sdata2018308
% N&C:   4 resting-state runs, ses-02 subfolder, n=199
%        Mendes et al., 2019: https://www.nature.com/articles/sdata2018307
% 109 participants took part in both protocols (total n=318)
%
%%% MRI data structure (for the files used in this script)
% Data_dir/
%   sub-010xxx/
%       ses-01/
%           anat/
%               sub-010xxx_ses-01_acq-mp2rage_T1w.nii.gz
%           dwi/
%           fmap/
%               sub-010xxx_ses-01_acq-GEfmap_run-01_magnitude1.nii.gz
%               sub-010xxx_ses-01_acq-GEfmap_run-01_phasediff.nii.gz
%           func/
%               sub-010xxx_ses-01_task-rest_acq-AP_run-01_bold.nii.gz
%       ses-02/
%           anat/ (only participants not included in the LEMON protocol)
%               sub-010xxx_ses-02_acq-mp2rage_T1w.nii.gz
%           fmap/
%               sub-010xxx_ses-02_acq-GEfmap_run-01_magnitude1.nii.gz  (for runs 1 and 2)
%               sub-010xxx_ses-02_acq-GEfmap_run-01_phasediff.nii.gz
%               sub-010xxx_ses-02_acq-GEfmap_run-02_magnitude1.nii.gz  (for runs 3 and 4)
%               sub-010xxx_ses-02_acq-GEfmap_run-02_phasediff.nii.gz
%           func/
%               sub-010xxx_ses-02_task-rest_acq-AP_run-01_bold.nii.gz  (run 1)
%               sub-010xxx_ses-02_task-rest_acq-AP_run-02_bold.nii.gz  (run 3)
%               sub-010xxx_ses-02_task-rest_acq-PA_run-01_bold.nii.gz  (run 2)
%               sub-010xxx_ses-02_task-rest_acq-PA_run-02_bold.nii.gz  (run 4)
%
% ! Due to a bug in multiband sequence used in the N&C protocol,
% the echo time for N&C resting-state is longer than in LEMON — 39.4 ms vs 30 ms. (< README file)
%
%%% Behavioural data acquired in the MPILMBB dataset:
% https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VMJ6NV
% Mind wandering questionnaire data are stored at Work_dir/Data
%   i.e. Work_dir/Data/SNYCQ_cleaned.csv
%
%%% Software specifications
% OS: Windows 10 21H2 (build 19044) Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz
% MATLAB: version 9.6.0.1072779 (R2019a)
% SPM: SPM12 version 7771
%
%%% Dependencies
% The script makes use of the following toolboxes (located under Work_dir/Toolboxes):
% HMM-MAR, PermCCA and BrewerMap
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialisation
clear;clc;close all
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

% set working directory
Work_dir = 'D:/Users/Liesa/OneDrive/Universiteit Gent/Statistical Data Analysis/THESIS/Code';
cd(Work_dir)

%help(cd) % uncomment to get an overview of all matlab scripts and functions
% or use 'help function' (e.g. help unzip_files) to get more info on a specific function

% specify the data directory (containing the MRI data) and SPM directory
% my data are stored on an external hard drive E:/
Data_dir = 'E:/ds000221-download';
SPM_dir = 'D:/Users/Liesa/Documents/MATLAB/spm12';

% add SPM12 and the other toolboxes to the search path
addpath(SPM_dir); addpath(genpath('Toolboxes'))

% sort the MPILMBB subject ids
n = 321; % 318 + 3 non-existing ids
[LEMON,NC,MPILMBB]=sort_participants(n,Data_dir,Work_dir,false); % set to true to write arrays to text files

% read the list of selected subject ids
subject_list = read_ids([Work_dir '/sub_ids_selected.txt']);
% to use the code with a new dataset, load a new text file with the subject ids
% and change the data specifications accordingly
% the data structure should be similar as specified above 
% (i.e. in BIDS format; https://bids.neuroimaging.io)

%% specifications
% specify how the nifti files can be recognised
% search starts from subject folder and selects alphabetically
% be as specific as possible (e.g. make sure the fieldmap corresponds to the EPI images)
% * is a wildcard character, do not use regex characters here
files.anat = 'sub*T1w.nii';
files.func = 'sub*ses-02*bold.nii';
files.fmap.magn = 'sub*ses-02*magnitude1.nii'; 
files.fmap.phsd = 'sub*ses-02*phasediff.nii';

multiple_runs = false;
% set to true to loop over all available functional runs
% if set to false, the first run found alphabetically is used

% specify the imaging parameters

% load slice times (convert to ms)
slice_timing = load('slice_timing.txt')'*1000;
% central slice as reference slice
slice_timing_sorted = sort(slice_timing);
refslice = slice_timing_sorted(length(slice_timing)/2);

parameters.anat = set_parameters({'voxsize(mm³)'},...
                                 {[1 1 1]});
parameters.func = set_parameters({'nvols' 'voxsize(mm³)' 'blipdir(1/-1)' 'nechoes' 'nslices' 'TR(s)' 'sliceorder(idx/ms)'...
                                  'refslice(idx/ms)' 'ESP(ms)' 'phsaccfactor' 'mbaccfactor'},...
                                 {657, [2.3 2.3 2.3], -1, 88, 64, 1.4, slice_timing, refslice, 0.67, 8/7, 4});
parameters.fmap = set_parameters({'shortTE(ms)' 'longTE(ms)'},...
                                 {5.19, 7.65});
                     % blip direction: defined positive from posterior to anterior (PA)
                     % if multiple_runs specify blipdir as follows e.g. 2 runs: [-1 1]

% remove high motion subjects?
motion.remove_high_motion_subjects = false; % subjects removed in questionnaire_data.R
motion.FD_threshold = [0.3 0.3]; % FD threshold defining high motion time points
motion.cutoff = {'mean' 0.5}; % proportion of time points, 'mean' or 'max'

% define ROIs (Schaefer et al., 2018; Yeo et al., 2011)
Schaefer100 = readtable('Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
Schaefer100 = [Schaefer100.ROIName num2cell([Schaefer100.R Schaefer100.A Schaefer100.S],2)]; % reformat
ROIs = Schaefer100(load('Schaefer2018_selected_idx.txt'),:); % selected ROIs

% load mind wandering data
% the data were first inspected and cleaned in R (see questionnaire_data.R)
MW_scores = readtable('Data/SNYCQ_cleaned.csv','TreatAsEmpty','NA');
MW_scores = MW_scores(ismember(MW_scores.participant_id,subject_list),[1 14:25]);
MW_scores = table2array(MW_scores(:,2:13));

% model specifications
K = 6; % maximum number of states to infer for each model
nreps = 10; % number of repetitions to run the models


% set to true to write results and figures to Work_dir/Results and Work_dir/Figures resp.
write_output = false;

% from here on function input arguments should not be changed
%%---------------------------------------------------------------------------------------
%% (optional) move unselected participants to Data_dir/not-selected
move_unselected(subject_list,Data_dir,false)

%% unpack .nii.gz files
%unzip_all(subject_list,Data_dir) % unpack all .nii.gz files
unzip_files(subject_list,files,multiple_runs,Data_dir) % unpack only the files used


%% MRI preprocessing
% SPM preprocessing pipelines were adapted from
% https://github.com/halmgren/Pipeline_effect_GSR_effective_connectivity_rsfMRI
% and further based on the SPM12 Manual https://www.fil.ion.ucl.ac.uk/spm/software/spm12

% structural images
preprocess_structural(subject_list,files.anat,Data_dir,SPM_dir)
% functional images
remove_first_five_volumes(subject_list,files.func,multiple_runs,Data_dir)
preprocess_functional(subject_list,files,parameters,multiple_runs,Data_dir,SPM_dir)

%% compute framewise displacement (FD)
FD_summary = calculate_FD(subject_list,files.func,multiple_runs,Data_dir,Work_dir,write_output);

%%
load([Work_dir '/Results/FD_summary.mat'])
%FD_summary(124,:) = []; % with high motion subjects
FD_summary = FD_summary(ismember(FD_summary.participant_id,subject_list),:); % without
FD_summary = table2array(FD_summary(:,2:3));

% histogram
figure
subplot(1,2,1); histogram(FD_summary(:,1)); title('Mean FD')
subplot(1,2,2); histogram(FD_summary(:,2)); title('Maximum FD')
% calculate sample quantiles (built-in matlab function 'quantile' uses interpolation)
round(sample_quantile([FD_summary(:,1) FD_summary(:,2)],[0 0.25 0.5 0.75 0.98 1]),3)

%% subjects with high FD values
high_motion_subjects = find_high_motion_subjects(subject_list,files.func,Data_dir,motion.FD_threshold,motion.cutoff); %#ok<*UNRCH>

if motion.remove_high_motion_subjects % removed in questionnaire_data.R
    subject_list(ismember(subject_list,high_motion_subjects,'rows'),:) = [];
end

%% extract BOLD time series
extract_regressors(subject_list,files.func,parameters.func,multiple_runs,Data_dir)
extract_time_series(subject_list,files.func,parameters.func,multiple_runs,ROIs,Data_dir)
cd(Work_dir)

% if working with multiple runs, temporally concatenate them
if multiple_runs; concat_ROI_data(subject_list,files.func,Data_dir); end

%%-----------------------------------------------------
%% Statistical modelling
% The HMM and GMM are estimated using the HMM-MAR toolbox
% https://github.com/OHBA-analysis/HMM-MAR

nvols = parameters.func('nvols')-5;

%%------------------------------
%% Hidden Markov model (HMM)
[~,~,~,~,~,~,run_minfe]=run_hmm(subject_list,files.func,parameters.func,multiple_runs,K,Data_dir,Work_dir,write_output,nreps);

%%
load([Work_dir '/Results/HMM_K6_rep06.mat'],'HMM_results') % run with least fe

% FC matrices (correlation)
figure
for s = 1:HMM_results.hmm.K
    if ~write_output; subplot(ceil(HMM_results.hmm.K/3),3,s); end
    plot_FCmatrix(HMM_results.hmm,s,Work_dir,write_output); title([' State ' sprintf('%d',s)],'fontsize',20)
end

% max fractional occupancy
figure
subplot(1,2,1); histogram(HMM_results.state_dynamics.maxFO,'FaceColor',[0.6 0.6 0.6]); 
xlabel('Max. FO'); ylabel('Number of subjects')
set(gca,'FontSize',12,'box','off')
% switching rate
subplot(1,2,2); histogram(HMM_results.state_dynamics.SwitchingRate,'FaceColor',[0.6 0.6 0.6]);
xlabel('Switching rate'); ylabel('Number of subjects')
set(gca,'FontSize',12,'box','off')
if write_output; print(gcf,[Work_dir '/Figures/HMM_mFO_SR.png'],'-dpng','-r300'); end

round(sample_quantile(HMM_results.state_dynamics.maxFO,[0 0.25 0.5 0.75 1]),2)
round(sample_quantile(HMM_results.state_dynamics.SwitchingRate,[0 0.25 0.5 0.75 1]),2)

% fractional occupancies
figure
for s = 1:HMM_results.hmm.K
    subplot(ceil(HMM_results.hmm.K/3),3,s)
    histogram(HMM_results.state_dynamics.FO(:,s)); 
    title(['State ' sprintf('%d',s)]); xlabel('Fractional occupancy (FO)'); ylabel('Number of subjects')
end
round(sample_quantile(HMM_results.state_dynamics.FO,[0 0.25 0.5 0.75 1]),2)
% dwell times
figure
for s = 1:HMM_results.hmm.K
    subplot(ceil(HMM_results.hmm.K/3),3,s)
    histogram(HMM_results.state_dynamics.StateLifeTimes{s}); 
    title([' State ' sprintf('%d',s)]); xlabel('Dwell time (number of TRs)'); ylabel('State visits')
    sample_quantile(HMM_results.state_dynamics.StateLifeTimes{s}',[0 0.25 0.5 0.75 1])
end

% time course for one subject
subj = 10; period = (1+nvols*(subj-1)):(nvols*subj);
figure
subplot(2,1,1); plot_prob_time_course(HMM_results.Gamma,period,Work_dir,write_output);
title('State activation probablity time course')
subplot(2,1,2); plot_state_path(HMM_results.vpath,HMM_results.hmm.K,period,Work_dir,write_output);
title('Most likely state path (Viterbi path)')

% transition probabilities between states (without considering persistence probabilities)
figure; plot_TPmatrix(HMM_results.hmm,Work_dir,write_output)
HMM_TP = getTransProbs(HMM_results.hmm);

%%------------------------------
%% Gaussian mixture model (GMM)
[~,~,~,~,~,~,run_minfe]=run_gmm(subject_list,files.func,parameters.func,multiple_runs,K,Data_dir,Work_dir,write_output,nreps);

%%
load([Work_dir '/Results/GMM_K6_rep01.mat'],'GMM_results') % run with least fe

% FC matrices
figure
for s = 1:GMM_results.gmm.K
    if ~write_output; subplot(ceil(GMM_results.gmm.K/3),3,s); end
    plot_FCmatrix(GMM_results.gmm,s,Work_dir,write_output); title([' State ' sprintf('%d',s)],'fontsize',20)
end

% max fractional occupancy
figure
subplot(1,2,1); histogram(GMM_results.state_dynamics.maxFO,'FaceColor',[0.6 0.6 0.6]); 
xlabel('Max. FO'); ylabel('Number of subjects')
set(gca,'FontSize',12,'box','off')
% switching rate
subplot(1,2,2); histogram(GMM_results.state_dynamics.SwitchingRate,'FaceColor',[0.6 0.6 0.6]);
xlabel('Switching rate'); ylabel('Number of subjects')
set(gca,'FontSize',12,'box','off')
if write_output; print(gcf,[Work_dir '/Figures/GMM_mFO_SR.png'],'-dpng','-r300'); end

round(sample_quantile(GMM_results.state_dynamics.maxFO,[0 0.25 0.5 0.75 1]),2)
round(sample_quantile(GMM_results.state_dynamics.SwitchingRate,[0 0.25 0.5 0.75 1]),2)

% fractional occupancies
figure
for s = 1:GMM_results.gmm.K
    subplot(ceil(GMM_results.gmm.K/3),3,s)
    histogram(GMM_results.state_dynamics.FO(:,s)); 
    title(['State ' sprintf('%d',s)]); xlabel('Fractional occupancy (FO)'); ylabel('Number of subjects')
end
round(sample_quantile(GMM_results.state_dynamics.FO,[0 0.25 0.5 0.75 1]),2)
% dwell times
figure
for s = 1:GMM_results.gmm.K
    subplot(ceil(GMM_results.gmm.K/3),3,s)
    histogram(GMM_results.state_dynamics.StateLifeTimes{s}); 
    title([' State ' sprintf('%d',s)]); xlabel('Dwell time (number of TRs)'); ylabel('State visits')
    sample_quantile(GMM_results.state_dynamics.StateLifeTimes{s}',[0 0.25 0.5 0.75 1])
end

% time course for one subject
subj = 10; period = (1+nvols*(subj-1)):(nvols*subj);
figure
subplot(2,1,1); plot_prob_time_course(GMM_results.Gamma,period,Work_dir,write_output);
title('State activation probablity time course')
subplot(2,1,2); plot_state_path(GMM_results.spath,GMM_results.gmm.K,period,Work_dir,write_output);
title('Most likely state path')

% transition probabilities between states (without considering the persistence probabilities)
figure; plot_TPmatrix(GMM_results.gmm,Work_dir,write_output)

%%------------------------------
%% Sliding window with k-means clustering (SWK)
rng(111)
[~,~,~,~,~,run_mintsumD]=run_swkmeans(subject_list,files.func,parameters.func,multiple_runs,K,Data_dir,Work_dir,write_output,nreps);

%%
load([Work_dir '/Results/SWK_K6_rep01.mat'],'SWK_results') % run with smallest total sumD

% FC matrices
figure
for s = 1:SWK_results.K
    if ~write_output; subplot(ceil(SWK_results.K/3),3,s); end
    plot_FCmatrix(SWK_results.FC,s,Work_dir,write_output); title([' State ' sprintf('%d',s)],'fontsize',20)
end

% max fractional occupancy
figure
subplot(1,2,1); histogram(SWK_results.state_dynamics.maxFO,'FaceColor',[0.6 0.6 0.6]); 
xlabel('Max. FO'); ylabel('Number of subjects')
set(gca,'FontSize',12,'box','off')
% switching rate
subplot(1,2,2); histogram(SWK_results.state_dynamics.SwitchingRate,'FaceColor',[0.6 0.6 0.6]);
xlabel('Switching rate'); ylabel('Number of subjects')
set(gca,'FontSize',12,'box','off')
if write_output; print(gcf,[Work_dir '/Figures/HMM_mFO_SR.png'],'-dpng','-r300'); end

round(sample_quantile(SWK_results.state_dynamics.maxFO,[0 0.25 0.5 0.75 1]),2)
round(sample_quantile(SWK_results.state_dynamics.SwitchingRate,[0 0.25 0.5 0.75 1]),2)

% fractional occupancies
figure
for s = 1:SWK_results.K
    subplot(ceil(SWK_results.K/3),3,s)
    histogram(SWK_results.state_dynamics.FO(:,s)); 
    title(['State ' sprintf('%d',s)]); xlabel('Fractional occupancy (FO)'); ylabel('Number of subjects')
end
round(sample_quantile(SWK_results.state_dynamics.FO,[0 0.25 0.5 0.75 1]),2)
% dwell times
figure
for s = 1:SWK_results.K
    subplot(ceil(SWK_results.K/3),3,s)
    histogram(SWK_results.state_dynamics.StateLifeTimes{s}); 
    title([' State ' sprintf('%d',s)]); xlabel('Dwell time (number of TRs)'); ylabel('State visits')
    sample_quantile(SWK_results.state_dynamics.StateLifeTimes{s}',[0 0.25 0.5 0.75 1])
end

% time course for one subject
lost = (SWK_results.window_length-1)/2; % time points lost at both ends
subj = 10; period = ((nvols-lost+1)+nvols*(subj-2)):(nvols*subj-lost); % do not set to 1 or N
figure; plot_state_path(SWK_results.spath,SWK_results.K,period,Work_dir,write_output);
title('Most likely state path')

%%-----------------------------------------------------
%% Correlation with mind wandering scores

%% correlations between the 12 MW items
% full correlations
[MW_corr.r,MW_corr.r_p] = corr(MW_scores,'type','Spearman');
% partial correlations
[MW_corr.pr,MW_corr.pr_p] = partialcorr(MW_scores,'type','Spearman');

figure;
subplot(1,2,1); plot_corrmatrix(MW_corr.r,MW_corr.r_p,'Spearman',false,Work_dir,write_output);
subplot(1,2,2); plot_corrmatrix(MW_corr.pr,MW_corr.pr_p,'Spearman',true,Work_dir,write_output);

%% canonical correlation analysis (CCA) 
% with permutation testing
% https://github.com/andersonwinkler/PermCCA

rng(111)
nperms = 1000;
HMM_CCA.X = MW_scores; GMM_CCA.X = MW_scores; SWK_CCA.X = MW_scores(2:end-1,:);

% set labels
labels.MW = {'Positive','Negative','Future','Past','Myself','Others','Surrounding',...
             'Vigilance','Images','Words','Specific','Intrusive'};
labels.S = {'State 1','State 2','State 3','State 4','State 5','State 6'};

%%------------------------------
%% HMM
HMM_CCA.Y = HMM_results.state_dynamics.FO;
[HMM_CCA.p,HMM_CCA.r,HMM_CCA.A,HMM_CCA.B,HMM_CCA.U,HMM_CCA.V] = permcca(HMM_CCA.X,HMM_CCA.Y,nperms); %#ok<*ASGLU>

figure; bar(HMM_CCA.r); xlabel('Dimension'); ylabel('Canonical correlation')
figure; scatter(HMM_CCA.U(:,1),HMM_CCA.V(:,1),'filled'); xlabel('Canonical variate U'); ylabel('Canonical variate V')

% loadings and cross-loadings
[HMM_CCA.XU,HMM_CCA.XV] = deal(zeros(size(HMM_CCA.X,2),length(HMM_CCA.r)));
[HMM_CCA.YV,HMM_CCA.YU] = deal(zeros(size(HMM_CCA.Y,2),length(HMM_CCA.r)));
for i = 1:length(HMM_CCA.r)
    for j = 1:size(HMM_CCA.X,2)
        HMM_CCA.XU(j,i) = corr(HMM_CCA.X(:,j),HMM_CCA.U(:,i));
        HMM_CCA.XV(j,i) = corr(HMM_CCA.X(:,j),HMM_CCA.V(:,i));
    end
    for j = 1:size(HMM_CCA.Y,2)
        HMM_CCA.YV(j,i) = corr(HMM_CCA.Y(:,j),HMM_CCA.V(:,i));
        HMM_CCA.YU(j,i) = corr(HMM_CCA.Y(:,j),HMM_CCA.U(:,i));
    end
end

figure; 
subplot(1,2,1); barh(HMM_CCA.XU(:,1),'FaceColor',[0.6 0.6 0.6]); 
yticklabels(labels.MW); xlabel('Loadings of canonical variate U')
subplot(1,2,2); barh(HMM_CCA.YV(:,1),'FaceColor',[0.6 0.6 0.6]); 
yticklabels(labels.S); xlabel('Loadings of canonical variate V')

figure; textscatter(HMM_CCA.XU(:,1),HMM_CCA.XU(:,2),labels.MW);
hold on; textscatter(HMM_CCA.YV(:,1),HMM_CCA.YV(:,2),labels.S);
xline(0,'--'); yline(0,'--'); xlabel('Dimension 1'); ylabel('Dimension 2')

if write_output; save([Work_dir '/Results/HMM_CCA.mat'],'HMM_CCA'); end

%%------------------------------
%% GMM
GMM_CCA.Y = GMM_results.state_dynamics.FO;
[GMM_CCA.p,GMM_CCA.r,GMM_CCA.A,GMM_CCA.B,GMM_CCA.U,GMM_CCA.V] = permcca(GMM_CCA.X,GMM_CCA.Y,nperms); %#ok<*ASGLU>

figure; bar(GMM_CCA.r); xlabel('Dimension'); ylabel('Canonical correlation')
figure; scatter(GMM_CCA.U(:,1),GMM_CCA.V(:,1),'filled'); xlabel('Canonical variate U'); ylabel('Canonical variate V')

% loadings and cross-loadings
[GMM_CCA.XU,GMM_CCA.XV] = deal(zeros(size(GMM_CCA.X,2),length(GMM_CCA.r)));
[GMM_CCA.YV,GMM_CCA.YU] = deal(zeros(size(GMM_CCA.Y,2),length(GMM_CCA.r)));
for i = 1:length(GMM_CCA.r)
    for j = 1:size(GMM_CCA.X,2)
        GMM_CCA.XU(j,i) = corr(GMM_CCA.X(:,j),GMM_CCA.U(:,i));
        GMM_CCA.XV(j,i) = corr(GMM_CCA.X(:,j),GMM_CCA.V(:,i));
    end
    for j = 1:size(GMM_CCA.Y,2)
        GMM_CCA.YV(j,i) = corr(GMM_CCA.Y(:,j),GMM_CCA.V(:,i));
        GMM_CCA.YU(j,i) = corr(GMM_CCA.Y(:,j),GMM_CCA.U(:,i));
    end
end

figure; barh(GMM_CCA.XU(:,1)); yticklabels(labels.MW)
figure; barh(GMM_CCA.YV(:,1)); yticklabels(labels.S)

figure; textscatter(GMM_CCA.XU(:,1),GMM_CCA.XU(:,2),labels.MW);
hold on; textscatter(GMM_CCA.YV(:,1),GMM_CCA.YV(:,2),labels.S);
xline(0,'--'); yline(0,'--'); xlabel('Dimension 1'); ylabel('Dimension 2')

if write_output; save([Work_dir '/Results/GMM_CCA.mat'],'GMM_CCA'); end

%%------------------------------
%% SWK
SWK_CCA.Y = SWK_results.state_dynamics.FO;
[SWK_CCA.p,SWK_CCA.r,SWK_CCA.A,SWK_CCA.B,SWK_CCA.U,SWK_CCA.V] = permcca(SWK_CCA.X,SWK_CCA.Y,nperms); %#ok<*ASGLU>

figure; bar(SWK_CCA.r); xlabel('Dimension'); ylabel('Canonical correlation')
figure; scatter(SWK_CCA.U(:,1),SWK_CCA.V(:,1),'filled'); xlabel('Canonical variate U'); ylabel('Canonical variate V')

% loadings and cross-loadings
[SWK_CCA.XU,SWK_CCA.XV] = deal(zeros(size(SWK_CCA.X,2),length(SWK_CCA.r)));
[SWK_CCA.YV,SWK_CCA.YU] = deal(zeros(size(SWK_CCA.Y,2),length(SWK_CCA.r)));
for i = 1:length(SWK_CCA.r)
    for j = 1:size(SWK_CCA.X,2)
        SWK_CCA.XU(j,i) = corr(SWK_CCA.X(:,j),SWK_CCA.U(:,i));
        SWK_CCA.XV(j,i) = corr(SWK_CCA.X(:,j),SWK_CCA.V(:,i));
    end
    for j = 1:size(SWK_CCA.Y,2)
        SWK_CCA.YV(j,i) = corr(SWK_CCA.Y(:,j),SWK_CCA.V(:,i));
        SWK_CCA.YU(j,i) = corr(SWK_CCA.Y(:,j),SWK_CCA.U(:,i));
    end
end

figure; barh(SWK_CCA.XU(:,1)); yticklabels(labels.MW)
figure; barh(SWK_CCA.YV(:,1)); yticklabels(labels.S)

figure; textscatter(SWK_CCA.XU(:,1),SWK_CCA.XU(:,2),labels.MW);
hold on; textscatter(SWK_CCA.YV(:,1),SWK_CCA.YV(:,2),labels.S);
xline(0,'--'); yline(0,'--'); xlabel('Dimension 1'); ylabel('Dimension 2')

if write_output; save([Work_dir '/Results/SWK_CCA.mat'],'SWK_CCA'); end

