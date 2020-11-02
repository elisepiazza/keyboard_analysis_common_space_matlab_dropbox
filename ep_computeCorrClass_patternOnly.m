%ep_computeCorrClass_patternOnly
%For each subject, for each ROI, for each condition, for each rep:
%1) Create the training pattern (column of voxels, avg'd across TRs, avg'd across N-1 reps).  
%This gives you 4 training patterns (1 for each condition). The 4 test patterns are the held-out patterns (rep 1).
%For each of the 4 test patterns, compute the corr between that and all 4 training patterns. E.g., for 1B, [.8 .2 .3 .3]. 
%Since the highest correlation is the first, you predict condition 1B. Do this for all 4 conditions, and you get 4 predictions 
%(e.g., [1B 2B 2B I]. Your accuracy (prop. correct) for this fold is 75%.
%Do this for all 3 reps (folds) and you get a total of 3 values. Average across them to get a single value.


%TO ADD: in which areas is I_A more correlated w/ I than the other 3
%scrambles?

clear;

group = 'M';
data_type = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2'; %v1_original_regressors, v2_jamals_regressors, v3_jamals_regressors_smoothing=1, v4_jamals_regressors_smoothing=1_defaultGMmask, v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3, v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2, v7_15_regressors_no_smoothing_defaultGMmask_polort=2
n_cropped_TRs = 10;

if strcmp(group, 'AM')
    subjects = [103 115 120 123];
elseif strcmp(group, 'M')
    subjects = [105 108 117 121 122];
end

for s = 1:length(subjects)
    subject = subjects(s);
    
    load(['../../common_space_AFNI/reshaped_by_conditions/' data_type '/sub-' num2str(subject) '.mat']);
    n_scramble_reps = size(data_ROIavg_scramble,4);
    
    nROIs = length(ROIs);
    
    for ROI = 1:nROIs
        
        data_scramble_thisROI = data_scramble{ROI}; %The extracted 4D matrix should be V x T x cond x rep
        data_scramble_thisROI = data_scramble_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:); %Crop TRs
        data_scramble_thisROI = squeeze(mean(data_scramble_thisROI,2)); %Average across TRs, resulting in V x cond x rep        
        
        for i = 1:n_scramble_reps
            
            %The current held-out rep
            test_rep = i;
            
            %The current included reps
            train_reps = setdiff([1:n_scramble_reps],test_rep);
            
            %Create 1B training pattern (Vx1)
            train_1B = mean(data_scramble_thisROI(:,1,train_reps),3);
            
            %Create 2B training pattern (Vx1)
            train_2B = mean(data_scramble_thisROI(:,2,train_reps),3);
            
            %Create 8B training pattern (Vx1)
            train_8B = mean(data_scramble_thisROI(:,3,train_reps),3);
            
            %Create I training pattern (Vx1)
            train_I = mean(data_scramble_thisROI(:,4,train_reps),3);
            
            %Create 1B test (held-out) pattern (Vx1)
            test_1B = data_scramble_thisROI(:,1,test_rep);
            
            %Create 2B test (held-out) pattern (VxT)
            test_2B = data_scramble_thisROI(:,2,test_rep);
            
            %Create 8B test (held-out) pattern (Vx1)
            test_8B = data_scramble_thisROI(:,3,test_rep);
            
            %Create I test (held-out) pattern (Vx1)
            test_I = data_scramble_thisROI(:,4,test_rep);
            
            %Compute the correlations between test_1B and all 4 training patterns
            R1 = corrcoef(test_1B,train_1B); R2 = corrcoef(test_1B,train_2B); R3 = corrcoef(test_1B,train_8B); R4 = corrcoef(test_1B,train_I);
            accs = [R1(2,1) R2(2,1) R3(2,1) R4(2,1)]; 
            %Check whether the max corr is w/ 1B's training pattern
            [m, el] = max(accs); acc_1B = el==1;
            
            %Compute the correlations between test_2B and all 4 training patterns
            R1 = corrcoef(test_2B,train_1B); R2 = corrcoef(test_2B,train_2B); R3 = corrcoef(test_2B,train_8B); R4 = corrcoef(test_2B,train_I);
            accs = [R1(2,1) R2(2,1) R3(2,1) R4(2,1)]; 
            %Check whether the max corr is w/ 2B's training pattern
            [m, el] = max(accs); acc_2B = el==2;

            %Compute the correlations between test_8B and all 4 training patterns
            R1 = corrcoef(test_8B,train_1B); R2 = corrcoef(test_8B,train_2B); R3 = corrcoef(test_8B,train_8B); R4 = corrcoef(test_8B,train_I);
            accs = [R1(2,1) R2(2,1) R3(2,1) R4(2,1)]; 
            %Check whether max corr is w/ 8B's training pattern
            [m, el] = max(accs); acc_8B = el==3;

            %Compute the correlations between test_I and all 4 training patterns
            R1 = corrcoef(test_I,train_1B); R2 = corrcoef(test_I,train_2B); R3 = corrcoef(test_I,train_8B); R4 = corrcoef(test_I,train_I);
            accs = [R1(2,1) R2(2,1) R3(2,1) R4(2,1)]; 
            %Check whether max corr is w/ I's training pattern
            [m, el] = max(accs); acc_I = el==4;

            %Store the average classification acc (across conditions) for
            %this held-out rep
            avg_acc(i) = mean([acc_1B acc_2B acc_8B acc_I]);
        end
        
        %Store the average classification acc (across reps) for this ROI
        ROI_acc(ROI) = mean(avg_acc);
    end
    
%     figsize = [100 100 400 500];
%     figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
%     print(gcf, '-dtiff', ['../figures/sub-' num2str(subject) '/Corr classifier_' data_type '.tif']);
    
    ROI_acc_allSubs(:,s) = ROI_acc;
    
end

%Plot each subject's acc across ROIs
figsize = [100 100 400 500];
figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc_allSubs); xlabel('Subject'); ylabel('ROI'); set(gca, 'XTickLabel', subjects, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
print(gcf, '-dtiff', ['../figures/Correlation classifier/Corr classifier_patternOnly_' group 'group_' data_type(1:2) '_n_cropped_TRs=' num2str(n_cropped_TRs)]);

%Plot each subject's above-chance acc across ROIs
figsize = [100 100 400 500];
figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc_allSubs>.25); xlabel('Subject'); ylabel('ROI'); set(gca, 'XTickLabel', subjects, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;
print(gcf, '-dtiff', ['../figures/Correlation classifier/Corr classifier_patternOnly_' group 'group_' data_type(1:2) '_n_cropped_TRs=' num2str(n_cropped_TRs) '_above_chance.tif']);


