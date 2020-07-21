%ep_computeWindowedPatternCorr
%For each subject, for each ROI, for each of the 4 scramble conditions (and then each of the 3 control conditions), 
%compute a TRxTR pattern corr matrix (compute the correlation between the rep-averaged voxels in each TR and every other TR)

clear;

ROI_to_plot = 2;

group = 'AM';
subjects = [103 115 120 123]; %AM

% group = 'M';
% subjects = [105 108 117 121 122]; %M;
 
n_cropped_TRs = 0;
TRs = 148 - 2*n_cropped_TRs;

nSubs = length(subjects);

all_subjects = [103 105 108 115 117 120 121 122 123]; 

preproc_type = 'AFNI';
preproc_param = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', 'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC'}; nROIs = length(ROIs);

TR_matrix_scramble = zeros(nROIs,TRs,TRs,4,nSubs);
TR_matrix_control = zeros(nROIs,TRs,TRs,3,nSubs);

for s = 1:nSubs
    subject = subjects(s);
    
    load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
    n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
    n_control_cond = size(data_ROIavg_control,3); n_control_reps = size(data_ROIavg_control,4);
    
    for ROI = 1:nROIs
        
        %Extract the full 4D matrices for this ROI (V x T x cond x rep)
        data_scramble_thisROI = data_scramble{ROI};
        data_control_thisROI = data_control{ROI};
        
        %Crop N TRs from beginning and end and average across reps
        data_repAvg_scramble = mean(data_scramble_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:),4);
        data_repAvg_control = mean(data_control_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:),4);
                
        %For each scramble condition, extract the V x T data for that
        %condition, compute the TR x TR correlation matrix, and save that matrix for this condition and subject 
        for scramble_cond = 1:n_scramble_cond           
           scramble_cond_data = data_repAvg_scramble(:,:,scramble_cond);
           TR_matrix_scramble(ROI,:,:,scramble_cond,s) = corr(scramble_cond_data,scramble_cond_data);           
        end  
        
        %Same for each control condition 
        for control_cond = 1:n_control_cond           
           control_cond_data = data_repAvg_control(:,:,control_cond);
           TR_matrix_control(ROI,:,:,control_cond,s) = corr(control_cond_data,control_cond_data);           
        end  
        
    end
    
end

%For each scramble condition, plot the group-averaged pattern correlation matrix
figsize = [100 100 1400 250]; 
figure('Units', 'pixels', 'Position', figsize);

subplot(1,4,1); imagesc(squeeze(mean(TR_matrix_scramble(ROI,:,:,1,:),5))); title('1B'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
subplot(1,4,2); imagesc(squeeze(mean(TR_matrix_scramble(ROI,:,:,2,:),5))); title('2B'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
subplot(1,4,3); imagesc(squeeze(mean(TR_matrix_scramble(ROI,:,:,3,:),5))); title('8B'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
subplot(1,4,4); imagesc(squeeze(mean(TR_matrix_scramble(ROI,:,:,4,:),5))); title('I'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
print(gcf, '-dtiff', ['../figures/Pattern correlation/Pattern correlation (' ROIs{ROI_to_plot} ', scramble, ' group ' group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);


%For each control condition, plot the group-averaged pattern correlation matrix 
figsize = [100 100 1000 250]; 
figure('Units', 'pixels', 'Position', figsize);

subplot(1,3,1); imagesc(squeeze(mean(TR_matrix_control(ROI,:,:,1,:),5))); title('I_N'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
subplot(1,3,2); imagesc(squeeze(mean(TR_matrix_control(ROI,:,:,2,:),5))); title('I_A'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
subplot(1,3,3); imagesc(squeeze(mean(TR_matrix_control(ROI,:,:,3,:),5))); title('I_I'); xlabel('TRs'); ylabel('TRs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
print(gcf, '-dtiff', ['../figures/Pattern correlation/Pattern correlation (' ROIs{ROI_to_plot} ', control, ' group ' group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
