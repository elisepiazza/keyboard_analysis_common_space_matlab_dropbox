%ep_computeFC
%For each subject, for each of the 4 scramble conditions, compute the functional connectivity (correlation)
%between the rep-averaged time course in that ROI and every other ROI

clear;
% subjects = [105 108 117 121 122]; %M;
subjects = [103 115 120 123]; %AM

nSubs = length(subjects);

all_subjects = [103 105 108 115 117 120 121 122 123]; 

preproc_type = 'AFNI';
preproc_param = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', 'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC'}; nROIs = length(ROIs);
ROI_order = [9 3 4 2 5 6 1 7 8 10];

n_cropped_TRs = 10;

ROI_matrix = zeros(nROIs,nROIs,4,nSubs);

for s = 1:nSubs
    subject = subjects(s);
    
    load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
    n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
        
    %Crop N TRs from beginning and end, extract ROIs in more reasonable order, and average across reps
    data_ROIavg_scramble = mean(data_ROIavg_scramble(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,:,:),4);
        
        %For each scramble condition, extract the ROI x TR data for that
        %condition, compute the correlation matrix, and save that matrix for this condition and subject 
        for cond = 1:n_scramble_cond           
           cond_data = data_ROIavg_scramble(:,:,cond);
           ROI_matrix(:,:,cond,s) = corr(cond_data',cond_data');           
        end  
end


%Plot the correlation matrices for the 4 scramble conditions
figsize = [100 100 1400 300]; 
figure('Units', 'pixels', 'Position', figsize);

subplot(1,4,1); imagesc(mean(ROI_matrix(:,:,1,:),4)); title('1B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
subplot(1,4,2); imagesc(mean(ROI_matrix(:,:,2,:),4)); title('2B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
subplot(1,4,3); imagesc(mean(ROI_matrix(:,:,3,:),4)); title('8B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
subplot(1,4,4); imagesc(mean(ROI_matrix(:,:,4,:),4)); title('I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
print(gcf, '-dtiff', ['../figures/Functional Connectivity (Audio-Motor group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);


%Note: I was going to do a similar thing but load all voxels for each ROI,
%but that's WAY too large (FC matrix has 3e9 elements)
for s = 1:nSubs
    subject = subjects(s);
    
    load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
    n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
    
    %Note: # voxels per ROI should be the same for each subject
    voxel_nums = zeros(nROIs,1);
    for ROI = 1:nROIs
        voxel_nums(ROI) = size(data_scramble{ROI},1);
    end
    
    for cond = 1:n_scramble_cond
    %Create an empty matrix that's voxel (all ROI voxels) x TR
    voxel_matrix = zeros(nROIs*sum(voxel_nums),size(data_ROIavg_scramble,2)-2*(n_cropped_TRs));
    
    for ROI = 1:nROIs
        %Extract the data for the current ROI (according to reasonable order above);
        %"data" is voxel x TR x cond x rep
        data = data_scramble{ROI_order(ROI)};
        
        %Crop N TRs from beginning and end and average across reps
        data = mean(data(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:),4);
        
        %Load data from this ROI (voxel x TR) into each scramble condition matrix
        cond_data_1B = vertcat(data(:,:,1),cond_data_1B);
        
    end
end
