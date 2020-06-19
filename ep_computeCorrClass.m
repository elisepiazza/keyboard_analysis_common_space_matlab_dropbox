%ep_computeCorrClass
%For each subject, for each ROI, for each of the 4 scrambled conditions,
%make an averaged VxT matrix (loaf) across nRep-1 reps and hold out the remaining
%one. For each held-out loaf, which of the 4 other conditions is it most correlated with? (Chance = .25)
%Repeat for nReps

%TO ADD: in which areas is I_A more correlated w/ I than the other 3
%scrambles?

clear;
subject = 123;

data_type = 'v2_jamals_regressors'; %v1_original_regressors, v2_jamals_regressors, v3_jamals_regressors_smoothing=1, v4_jamals_regressors_smoothing=1_defaultGMmask, v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3, v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2, v7_15_regressors_no_smoothing_defaultGMmask_polort=2 
load(['../reshaped_by_conditions/' data_type '/sub-' num2str(subject) '.mat']);
n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
choose = @(samples) samples(randi(numel(samples)));

nROIs = length(ROIs);

%Crop N TRs from beginning and end
n_cropped_TRs = 0; 
data_ROIavg_scramble = data_ROIavg_scramble(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);
data_ROIavg_control = data_ROIavg_control(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);

for ROI = 1:nROIs
    
    data_scramble_thisROI = data_scramble{ROI}; %The extracted 4D matrix should be V x T x cond x rep
    
    held_out_reps = randperm(n_scramble_reps);
    
    for i = 1:n_scramble_reps
        
            test_rep = held_out_reps(i);

            %Average of 1B training reps (VxT)
            train_1B = mean(data_scramble_thisROI(:,:,1,setdiff([1:n_scramble_reps],test_rep)),4);
            
            %Average of 2B training reps (VxT)
            train_2B = mean(data_scramble_thisROI(:,:,2,setdiff([1:n_scramble_reps],test_rep)),4);

            %Average of 8B training reps (VxT)
            train_8B = mean(data_scramble_thisROI(:,:,3,setdiff([1:n_scramble_reps],test_rep)),4);

            %Average of I training reps (VxT)
            train_I = mean(data_scramble_thisROI(:,:,4,setdiff([1:n_scramble_reps],test_rep)),4);
            
            %Held-out 1B run (VxT)
            test_1B = data_scramble_thisROI(:,:,1,test_rep);
            
            %Held-out 2B run (VxT)
            test_2B = data_scramble_thisROI(:,:,2,test_rep);
            
            %Held-out 8B run (VxT)
            test_8B = data_scramble_thisROI(:,:,3,test_rep);

            %Held-out I run (VxT)
            test_I = data_scramble_thisROI(:,:,4,test_rep);
            
            %Is train_1B most strongly correlated with its own held-out (test) loaf than the other 3? 
            R1 = corrcoef(train_1B(:,:),test_1B(:,:)); R2 = corrcoef(train_1B(:,:),test_I(:,:)); R3 = corrcoef(train_1B(:,:),test_8B(:,:)); R4 = corrcoef(train_1B(:,:),test_2B(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_1B(i) = sum(acc1) == 3;
          
            %Is train_2B most strongly correlated with its own held-out (test) loaf than the other 3? 
            R1 = corrcoef(train_2B(:,:),test_2B(:,:)); R2 = corrcoef(train_2B(:,:),test_I(:,:)); R3 = corrcoef(train_2B(:,:),test_8B(:,:)); R4 = corrcoef(train_2B(:,:),test_1B(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_2B(i) = sum(acc1) == 3;
         
            %Is train_8B most strongly correlated with its own held-out (test) loaf than the other 3? 
            R1 = corrcoef(train_8B(:,:),test_8B(:,:)); R2 = corrcoef(train_8B(:,:),test_I(:,:)); R3 = corrcoef(train_8B(:,:),test_2B(:,:)); R4 = corrcoef(train_8B(:,:),test_1B(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_8B(i) = sum(acc1) == 3;

            %Is train_1 most strongly correlated with its own held-out (test) loaf than the other 3? 
            R1 = corrcoef(train_I(:,:),test_I(:,:)); R2 = corrcoef(train_I(:,:),test_8B(:,:)); R3 = corrcoef(train_I(:,:),test_2B(:,:)); R4 = corrcoef(train_I(:,:),test_1B(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_I(i) = sum(acc1) == 3;
    end
    
    ROI_acc(ROI,1) = mean(acc_1B);
    ROI_acc(ROI,2) = mean(acc_2B);
    ROI_acc(ROI,3) = mean(acc_8B);
    ROI_acc(ROI,4) = mean(acc_I);

end

figsize = [100 100 400 500];
figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
print(gcf, '-dtiff', ['../figures/sub-' num2str(subject) '/Corr classifier_' data_type '.tif']);

%
% %Plot corr classifier results across runs
% N = length(held_out_runs);
% x = 1;
% y = mean(mean_acc);
% errors = std(mean_acc)/sqrt(N);
%
% figsize = [100 100 300 375]; barwidth = .6; barcolor = [.5 0 .9];
% figure('Units', 'pixels', 'Position', figsize);
% bar(x,y,barwidth,'facecolor',barcolor); hold on;
% errorbar(x,y,errors,'k.', 'LineWidth', 1)
%
% xlabel('Melody 1 vs. Melody 2'); ylabel('Mean accuracy across runs'); title(['Corr classifier (' ROI_names{whichROI} smooth_tags{smoothOrNot} ')']); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
% % print(gcf, '-dtiff', ['../figures/Corr classifier (' ROI_names{whichROI} smooth_tags{smoothOrNot} ').tif']);

