%ep_computeFC
%For each subject, for each of the 4 scramble conditions (and then each of the 3 control conditions), 
%compute the functional connectivity (correlation) between the rep-averaged time course in that ROI and every other ROI. 

%Note: I was going to do a similar thing but load all voxels for each ROI,
%but that's WAY too large (FC matrix has 3e9 elements)

clear;

% group = 'AM';
% subjects = [103 115 120 123]; %AM

group = 'M';
subjects = [105 108 117 121 122]; %M;

n_cropped_TRs = 0;

nSubs = length(subjects);

all_subjects = [103 105 108 115 117 120 121 122 123]; 

preproc_type = 'AFNI';
preproc_param = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', 'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC'}; nROIs = length(ROIs);
ROI_order = [9 3 4 2 5 6 1 7 8 10];

ROI_matrix_scramble = zeros(nROIs,nROIs,4,nSubs);
ROI_matrix_control = zeros(nROIs,nROIs,3,nSubs);

for s = 1:nSubs
    subject = subjects(s);
    
    load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
    n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
    n_control_cond = size(data_ROIavg_control,3); n_control_reps = size(data_ROIavg_control,4);
        
    %Crop N TRs from beginning and end, extract ROIs in more reasonable order, and average across reps
    data_ROIavg_scramble = mean(data_ROIavg_scramble(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,:,:),4);
    data_ROIavg_control = mean(data_ROIavg_control(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,:,:),4);
        
        %For each scramble condition, extract the ROI x TR data for that
        %condition, compute the correlation matrix, and save that matrix for this condition and subject 
        for scramble_cond = 1:n_scramble_cond           
           scramble_cond_data = data_ROIavg_scramble(:,:,scramble_cond);
           ROI_matrix_scramble(:,:,scramble_cond,s) = corr(scramble_cond_data',scramble_cond_data');           
        end  
        
        %Same for each control condition 
        for control_cond = 1:n_control_cond           
           control_cond_data = data_ROIavg_control(:,:,control_cond);
           ROI_matrix_control(:,:,control_cond,s) = corr(control_cond_data',control_cond_data');           
        end  

end


%For each scramble condition, plot the group-averaged FC matrix
figsize = [100 100 1400 250]; 
figure('Units', 'pixels', 'Position', figsize);

subplot(1,4,1); imagesc(mean(ROI_matrix_scramble(:,:,1,:),4)); title('1B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
subplot(1,4,2); imagesc(mean(ROI_matrix_scramble(:,:,2,:),4)); title('2B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
subplot(1,4,3); imagesc(mean(ROI_matrix_scramble(:,:,3,:),4)); title('8B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
subplot(1,4,4); imagesc(mean(ROI_matrix_scramble(:,:,4,:),4)); title('I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
print(gcf, '-dtiff', ['../figures/Functional connectivity/Functional Connectivity (scramble, ' group ' group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);


%For each control condition, plot the group-averaged FC matrix 
figsize = [100 100 1000 250]; 
figure('Units', 'pixels', 'Position', figsize);

subplot(1,3,1); imagesc(mean(ROI_matrix_control(:,:,1,:),4)); title('I_N'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
subplot(1,3,2); imagesc(mean(ROI_matrix_control(:,:,2,:),4)); title('I_A'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
subplot(1,3,3); imagesc(mean(ROI_matrix_control(:,:,3,:),4)); title('I_I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
print(gcf, '-dtiff', ['../figures/Functional connectivity/Functional Connectivity (control, ' group ' group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);

