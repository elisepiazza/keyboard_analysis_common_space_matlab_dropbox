%ep_computeCorrClass
%For each subject, for each ROI, for each of the 4 scrambled conditions,
%make an averaged VxT matrix (loaf) across nRep-1 reps and hold out the remaining
%one. For each held-out loaf, which of the 4 other conditions is it most correlated with? (Chance = .25)
%Repeat for nReps.

clear;
subjects = [103 105 108 115 117 120 121 122 123]; nSubs = length(subjects);

preproc_type = 'AFNI';
preproc_param = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

n_cropped_TRs = 10;

ROI_acc_control_vs_scramble = zeros(10, 2, nSubs); %ROI x nConds x nSubs

for s = 1:nSubs
    subject = subjects(s);
    
    load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
    n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
    n_control_cond = size(data_ROIavg_control,3); n_control_reps = size(data_ROIavg_control,4);
    
    nROIs = length(ROIs);
        
    for ROI = 1:nROIs
        
        data_scramble_thisROI = data_scramble{ROI}; %The extracted 4D matrix should be V x T x cond x rep
        data_scramble_thisROI = data_scramble_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);
        
        data_control_thisROI = data_control{ROI}; %The extracted 4D matrix should be V x T x cond x rep
        data_control_thisROI = data_control_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);
        held_out_reps = randperm(n_control_reps);
        
        
        %Run classifier for control (I_A, I_I) vs. scramble conditions
        for i = 1:n_control_reps
            
            test_rep = held_out_reps(i);
            
            %Average of I_A training reps (VxT)
            train_I_A = mean(data_control_thisROI(:,:,2,setdiff([1:n_control_reps],test_rep)),4);
            
            %Average of I_I training reps (VxT)
            train_I_I = mean(data_control_thisROI(:,:,3,setdiff([1:n_control_reps],test_rep)),4);
            
            %Average 1B data across reps (VxT)
            test_1B = mean(data_scramble_thisROI(:,:,1,:),4);
            
            %Average 2B data across reps (VxT)
            test_2B = mean(data_scramble_thisROI(:,:,2,:),4);
            
            %Average 8B data across reps (VxT)
            test_8B = mean(data_scramble_thisROI(:,:,3,:),4);
            
            %Average I data across reps (VxT)
            test_I = mean(data_scramble_thisROI(:,:,4,:),4);
            
            %Is train_I_A most strongly correlated with test_I than test_1B, test_2B, or test_8B?
            R1 = corrcoef(train_I_A(:,:),test_I(:,:)); R2 = corrcoef(train_I_A(:,:),test_1B(:,:)); R3 = corrcoef(train_I_A(:,:),test_2B(:,:)); R4 = corrcoef(train_I_A(:,:),test_8B(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_I_A_vs_I(i) = sum(acc1) == 3;
            
            %Is train_I_I most strongly correlated with test_I than test_1B, test_2B, or test_8B?
            R1 = corrcoef(train_I_I(:,:),test_I(:,:)); R2 = corrcoef(train_I_I(:,:),test_1B(:,:)); R3 = corrcoef(train_I_I(:,:),test_2B(:,:)); R4 = corrcoef(train_I_I(:,:),test_8B(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_I_I_vs_I(i) = sum(acc1) == 3;
        end
        
        ROI_acc_control_vs_scramble(ROI,1,s) = mean(acc_I_A_vs_I);
        ROI_acc_control_vs_scramble(ROI,2,s) = mean(acc_I_I_vs_I);
        
    end
    
    if strcmp(preproc_type, 'AFNI')
        figsize = [100 100 300 500];
    elseif strcmp(preproc_type, 'Python')
        figsize = [100 100 300 350];
    end
    
    
    %Plot classification accuracy (ROI x cond) for this subject (I_A, I_I vs. scramble conditions)
    figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc_control_vs_scramble(:,:,s)); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', {'I_A', 'I_I'}, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
    print(gcf, '-dtiff', ['../figures/Correlation classifier/sub-' num2str(subject) '_controlvscramble_' preproc_type '_' preproc_param '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
        
end


%Plot classification accuracy (ROI x condition (I_A, I_I)
figsize = [100 100 300 400];
figure('Units', 'pixels', 'Position', figsize);
imagesc(mean(ROI_acc_control_vs_scramble,3));
title('Control v. Scramble conditions'); xlabel('Condition'); ylabel('ROI'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica', 'XTickLabel', {'I_A', 'I_I'}, 'YTickLabel', ROIs); colorbar; caxis([0 .5]);
print(gcf, '-dtiff', ['../figures/Correlation classifier/Summary by ROI_controlvsscramble_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);


% %Plot avg classification accuracy across subjects, for each preproc combo (control v. scramble)
% y = mean(avg_acc_control_vs_scramble);
% errors = std(avg_acc_control_vs_scramble)/sqrt(nSubs);
% x = 1:length(preproc_types);
% 
% figsize = [100 100 600 400]; barwidth = .5;
% figure('Units', 'pixels', 'Position', figsize);
% bar(x,y,barwidth,'facecolor', [.2 .8 .9]); hold on;
% errorbar(x,y,errors,'k.', 'LineWidth', 1)
% 
% xlabel('Preprocessing Param Type'); ylabel('Correlation classifier accuracy (% correct)'); xlim([.3 10.7]); ylim([0 1]); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
% print(gcf, '-dtiff', ['../figures/CC/Summary stats_Control v. Scramble_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
