%ep_comparePipelines

%For each subject, for each of 3 ROIs common to the two pipelines (A1, angular
%gyrus, precuneus), import (V) x T x cond x rep. 
%Ultimate goal: for each subject, correlation
%matrix of AFNI[A1, cond1, rep1], AFNI[A1, cond1, rep2]... etc. vs.
%Python[A1, cond1, rep1]...etc, so a 12 x 12 matrix of correlations

clear; 
subject = 123; 
n_cropped_TRs = 10; %nTRs to crop from beginning and end
AFNI_data_type = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2'; %v1_original_regressors, v2_jamals_regressors, v3_jamals_regressors_smoothing=1, v4_jamals_regressors_smoothing=1_defaultGMmask, v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3, v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2, v7_15_regressors_no_smoothing_defaultGMmask_polort=2  
Python_data_type = 'HPF=.01Hz'; %HPF=.01Hz, HPF=.03Hz, HPF=.06Hz
all_subjects = [103 105 108 115 117 120 121 122 123]; whichSubject = find(all_subjects==subject);
groups = {'AM', 'M', 'M', 'AM', 'M', 'AM', 'M', 'M', 'AM'}; group = groups{whichSubject};

%Define the 12 cond+rep combinations
scramble_indices = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3; 4 1; 4 2; 4 3]; matrix_length = size(scramble_indices,1);
index_labels = {'1B_1', '1B_2', '1B_3', '2B_1', '2B_2', '2B_3', '8B_1', '8B_2', '8B_3', 'I_1', 'I_2', 'I_3'};
nTRs = 148 - n_cropped_TRs - n_cropped_TRs;

AFNI_ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', ...
    'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC'};
Python_ROIs = {'A1', 'AngularGyrus', 'Erez-DMN', 'Precuneus', 'vmPFC'};

%Load the ROI x TR x cond x rep data
load(['../../common_space_AFNI/reshaped_by_conditions/' AFNI_data_type '/sub-' num2str(subject) '.mat']);
data_AFNI = data_ROIavg_scramble;

load(['../../common_space_Python/reshaped_by_conditions/' Python_data_type '/sub-' num2str(subject) '.mat']);
data_python = data_ROIavg_scramble;

A1_data_AFNI = zeros(nTRs,matrix_length); A1_data_python = zeros(nTRs,matrix_length);
AngG_data_AFNI = zeros(nTRs,matrix_length); AngG_data_python = zeros(nTRs,matrix_length);
Precun_data_AFNI = zeros(nTRs,matrix_length); Precun_data_python = zeros(nTRs,matrix_length);

%Load in the TR x 12 data
for i = 1:matrix_length
    A1_data_AFNI(:,i) = zscore(data_AFNI(9,n_cropped_TRs+1:end-n_cropped_TRs,scramble_indices(i,1),scramble_indices(i,2)));
    A1_data_python(:,i) = zscore(data_python(1,n_cropped_TRs+1:end-n_cropped_TRs,scramble_indices(i,1),scramble_indices(i,2)));
    
    AngG_data_AFNI(:,i) = zscore(data_AFNI(1,n_cropped_TRs+1:end-n_cropped_TRs,scramble_indices(i,1),scramble_indices(i,2)));
    AngG_data_python(:,i) = zscore(data_python(2,n_cropped_TRs+1:end-n_cropped_TRs,scramble_indices(i,1),scramble_indices(i,2)));
    
    Precun_data_AFNI(:,i) = zscore(data_AFNI(8,n_cropped_TRs+1:end-n_cropped_TRs,scramble_indices(i,1),scramble_indices(i,2)));
    Precun_data_python(:,i) = zscore(data_python(4,n_cropped_TRs+1:end-n_cropped_TRs,scramble_indices(i,1),scramble_indices(i,2)));
end

%Compute correlation matrices
A1_corr_matrix = corr(A1_data_AFNI, A1_data_python);
AngG_corr_matrix = corr(AngG_data_AFNI, AngG_data_python);
Precun_corr_matrix = corr(Precun_data_AFNI, Precun_data_python);

figsize = [100 100 2000 350]; figure('Units', 'pixels', 'Position', figsize); 

subplot(1,3,1); imagesc(A1_corr_matrix); colorbar; caxis([-.2 .8]); title('A1'); xlabel('Python pipeline'); ylabel('AFNI pipeline'); 
set(gca, 'XTick', 1:12, 'YTick', 1:12, 'XTickLabel', index_labels, 'YTickLabel', index_labels, 'FontSize', 16, 'FontName', 'Helvetica'); 

subplot(1,3,2); imagesc(AngG_corr_matrix); colorbar; caxis([-.2 .8]); title('Angular Gyrus'); xlabel('Python pipeline'); ylabel('AFNI pipeline'); 
set(gca, 'XTick', 1:12, 'YTick', 1:12, 'XTickLabel', index_labels, 'YTickLabel', index_labels, 'FontSize', 16, 'FontName', 'Helvetica'); 

subplot(1,3,3); imagesc(Precun_corr_matrix); colorbar; caxis([-.2 .8]); title('Precuneus'); xlabel('Python pipeline'); ylabel('AFNI pipeline'); 
set(gca, 'XTick', 1:12, 'YTick', 1:12, 'XTickLabel', index_labels, 'YTickLabel', index_labels, 'FontSize', 16, 'FontName', 'Helvetica'); 


print(gcf, '-dtiff', ['../figures/s' num2str(subject) '_' group '_compare pipelines_' AFNI_data_type '_' Python_data_type '_croppedTRs= ' num2str(n_cropped_TRs) '.tif']);


%Inspect 3 "well-matched" time series, 3 "badly-matched" time series for each subject

%Based on the corr matrices, define which are good and bad
A1_good = [1 3 1 9 8 7 10 10 3];   A1_good_cond = {'1B1', '1B3', '1B1', '8B3', '8B2', '8B1', 'I1', 'I1', '1B3'};
A1_bad =  [12 5 12 5 9 9 3 3 11];  A1_bad_cond = {'I3', '2B2', 'I3', '2B2', '8B3', '8B3', '8B3', '1B3', 'I2'};

AngG_good = [2 8 8 9 8 10 6 7 10];  AngG_good_cond = {'1B2', '8B2', '8B2', '8B3', '8B2', 'I1', '2B3', '8B1', 'I1'};
AngG_bad =  [9 11 12 6 7 12 1 3 2]; AngG_bad_cond = {'8B3', 'I2', 'I3', '2B3', '8B1', 'I3', '1B1', '1B3', '1B2'};

Precun_good = [2 12 4 9 4 2 6 4 10]; Precun_good_cond = {'1B2', 'I3', '2B1', '8B3', '2B1', '1B2', '2B3', '2B1', 'I1'};
Precun_bad =  [9 1 12 6 9 12 1 2 1];  Precun_bad_cond = {'8B3', '1B1', 'I3', '2B3', '8B3', 'I3', '1B1', '1B2', '1B1'};
    

%Plot
figure('Units', 'pixels', 'Position', [0 0 1500 800]);

%Plot a high corr in A1
subplot(3,2,1); plot(A1_data_AFNI(:,A1_good(whichSubject)),'LineWidth',2); hold on; plot(A1_data_python(:,A1_good(whichSubject)),'LineWidth',2);
xlabel('TR'); ylabel('BOLD'); title(['s' num2str(subject) ', group=' group ', cond=' A1_good_cond{whichSubject} ', A1']); set(gca,'FontSize',12); legend({'python', 'AFNI'})

%Plot a low corr in A1
subplot(3,2,2); plot(A1_data_AFNI(:,A1_bad(whichSubject)),'LineWidth',2); hold on; plot(A1_data_python(:,A1_bad(whichSubject)),'LineWidth',2);
xlabel('TR'); ylabel('BOLD'); title(['s' num2str(subject) ', group=' group ', cond=' A1_bad_cond{whichSubject} ', A1']); set(gca,'FontSize',12); 

%Plot a high corr in AngG
subplot(3,2,3); plot(AngG_data_AFNI(:,AngG_good(whichSubject)),'LineWidth',2); hold on; plot(AngG_data_python(:,AngG_good(whichSubject)),'LineWidth',2);
xlabel('TR'); ylabel('BOLD'); title(['s' num2str(subject) ', group=' group ', cond=' AngG_good_cond{whichSubject} ', AngG']); set(gca,'FontSize',12); 
% 
%Plot a low corr in AngG
subplot(3,2,4); plot(AngG_data_AFNI(:,AngG_bad(whichSubject)),'LineWidth',2); hold on; plot(AngG_data_python(:,AngG_bad(whichSubject)),'LineWidth',2);
xlabel('TR'); ylabel('BOLD'); title(['s' num2str(subject) ', group=' group ', cond=' AngG_bad_cond{whichSubject} ', AngG']); set(gca,'FontSize',12); 

%Plot a high corr in PreC
subplot(3,2,5); plot(Precun_data_AFNI(:,Precun_good(whichSubject)),'LineWidth',2); hold on; plot(Precun_data_python(:,Precun_good(whichSubject)),'LineWidth',2);
xlabel('TR'); ylabel('BOLD'); title(['s' num2str(subject) ', group=' group ', cond=' Precun_good_cond{whichSubject} ', Precun']); set(gca,'FontSize',12); 

%Plot a low corr in PreC
subplot(3,2,6); plot(Precun_data_AFNI(:,Precun_bad(whichSubject)),'LineWidth',2); hold on; plot(Precun_data_python(:,Precun_bad(whichSubject)),'LineWidth',2);
xlabel('TR'); ylabel('BOLD'); title(['s' num2str(subject) ', group=' group ', cond=' Precun_bad_cond{whichSubject} ', Precun']); set(gca,'FontSize',12); 

print(gcf, '-dtiff', ['../figures/s' num2str(subject) '_' group '_compare pipelines_timeseries_' AFNI_data_type '_' Python_data_type '_croppedTRs= ' num2str(n_cropped_TRs) '.tif']);

% %Plot difference between python and AFNI time series
% diff_btwn_data_types = A1_data_python - A1_data_AFNI;
% figure('Units', 'pixels', 'Position', [0 0 1500 800]);
% for i = 1:6
%     subplot(6,2,(i*2)-1); plot(diff_btwn_data_types(:,(i*2)-1),'r','LineWidth',2); xlabel('TR'); ylabel('Python - AFNI'); set(gca, 'FontSize', 12); line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k'); 
%     subplot(6,2,i*2); plot(diff_btwn_data_types(:,i*2),'r','LineWidth',2); xlabel('TR'); ylabel('Python - AFNI'); set(gca, 'FontSize', 12); line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k'); 
% end
% print(gcf, '-dtiff', ['../figures/s' num2str(subject) '_' group '_subtract_timeseries_' AFNI_data_type '_' Python_data_type '_croppedTRs= ' num2str(n_cropped_TRs) '.tif']);
% 

