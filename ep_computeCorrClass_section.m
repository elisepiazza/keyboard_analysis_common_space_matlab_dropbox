%ep_computeCorrClass_section
%How well can you classify which of the 3 main sections of the piece
%they're in? (In the Intact, I_A, and I_I conditions)
%For each subject, for each ROI, for each of the 3 sections, make an averaged VxT matrix (loaf)
%across nRep-1 reps and hold out the remaining one. For each held-out loaf,
%which of the 3 sections of the Intact piece is it most correlated with? (Chance = .25)
%Repeat for nReps.

clear;
subjects = [103 105 108 115 117 120 121 122 123]; nSubs = length(subjects);

preproc_type = 'AFNI';
preproc_param = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

%section_starts_measures = [1 34 64 88];
%The following are in TRs (not measures!)
starts = [1 41 78 107] + 2;
lengths = 30;
n_sections = 4;

%Initialize empty matrices of accuracies
ROI_acc_sections = zeros(10,4,nSubs); %nROIs x nSections x nSubs


for s = 1:nSubs
    subject = subjects(s);
    
    load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
    n_scramble_reps = size(data_ROIavg_scramble,4);
    
    nROIs = length(ROIs);
    
    for ROI = 1:nROIs
        
        %Run classifier for 4 sections of Intact condition
        data_scramble_thisROI = data_scramble{ROI}; %The extracted 4D matrix should be V x T x cond x rep
        held_out_reps = randperm(n_scramble_reps);
        
        for i = 1:n_scramble_reps
            
            test_rep = held_out_reps(i);
            
            %Average of section 1 training reps (VxT)
            train_1 = mean(data_scramble_thisROI(:,starts(1):starts(1)+lengths-1,4,setdiff([1:n_scramble_reps],test_rep)),4);
            
            %Average of section 2 training reps (VxT)
            train_2 = mean(data_scramble_thisROI(:,starts(2):starts(2)+lengths-1,4,setdiff([1:n_scramble_reps],test_rep)),4);
            
            %Average of section 3 training reps (VxT)
            train_3 = mean(data_scramble_thisROI(:,starts(3):starts(3)+lengths-1,4,setdiff([1:n_scramble_reps],test_rep)),4);
            
            %Average of section 4 training reps (VxT)
            train_4 = mean(data_scramble_thisROI(:,starts(4):starts(4)+lengths-1,4,setdiff([1:n_scramble_reps],test_rep)),4);
            
            %Held-out section 1 rep (VxT)
            test_1 = data_scramble_thisROI(:,starts(1):starts(1)+lengths-1,4,test_rep);
            
            %Held-out section 2 rep (VxT)
            test_2 = data_scramble_thisROI(:,starts(2):starts(2)+lengths-1,4,test_rep);
            
            %Held-out section 3 rep (VxT)
            test_3 = data_scramble_thisROI(:,starts(3):starts(3)+lengths-1,4,test_rep);
            
            %Held-out section 4 rep (VxT)
            test_4 = data_scramble_thisROI(:,starts(4):starts(4)+lengths-1,4,test_rep);
            
            %Is train_1 most strongly correlated with its own held-out (test) loaf than the other 3?
            R1 = corrcoef(train_1(:,:),test_1(:,:)); R2 = corrcoef(train_1(:,:),test_2(:,:)); R3 = corrcoef(train_1(:,:),test_3(:,:)); R4 = corrcoef(train_1(:,:),test_4(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_1(i) = sum(acc1) == n_sections-1;
            
            %Is train_2 most strongly correlated with its own held-out (test) loaf than the other 3?
            R1 = corrcoef(train_2(:,:),test_2(:,:)); R2 = corrcoef(train_2(:,:),test_3(:,:)); R3 = corrcoef(train_2(:,:),test_4(:,:)); R4 = corrcoef(train_2(:,:),test_1(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_2(i) = sum(acc1) == n_sections-1;
            
            %Is train_3 most strongly correlated with its own held-out (test) loaf than the other 3?
            R1 = corrcoef(train_3(:,:),test_3(:,:)); R2 = corrcoef(train_3(:,:),test_4(:,:)); R3 = corrcoef(train_3(:,:),test_1(:,:)); R4 = corrcoef(train_3(:,:),test_2(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_3(i) = sum(acc1) == n_sections-1;
            
            %Is train_4 most strongly correlated with its own held-out (test) loaf than the other 3?
            R1 = corrcoef(train_4(:,:),test_4(:,:)); R2 = corrcoef(train_4(:,:),test_1(:,:)); R3 = corrcoef(train_4(:,:),test_2(:,:)); R4 = corrcoef(train_4(:,:),test_3(:,:));
            acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_4(i) = sum(acc1) == n_sections-1;
            
        end
        
        ROI_acc_sections(ROI,1,s) = mean(acc_1);
        ROI_acc_sections(ROI,2,s) = mean(acc_2);
        ROI_acc_sections(ROI,3,s) = mean(acc_3);
        ROI_acc_sections(ROI,4,s) = mean(acc_4);
    end
    
    %Plot classification accuracy (ROI x cond) for this subject (all 4 sections)
    figsize = [100 100 400 500];
    figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc_sections(:,:,s)); xlabel('Section'); ylabel('ROI'); set(gca, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
    print(gcf, '-dtiff', ['../figures/Correlation classifier/sub-' num2str(subject) '_sections_' preproc_type '_' preproc_param '.tif']);
    
end

%Plot classification accuracy, averaged across subjects (ROI x section)
figsize = [100 100 300 400];
figure('Units', 'pixels', 'Position', figsize);
imagesc(mean(ROI_acc_sections,3));
title('Classify sections'); xlabel('Section'); ylabel('ROI'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica', 'YTickLabel', ROIs); colorbar; caxis([0 1]);
print(gcf, '-dtiff', ['../figures/Correlation classifier/Summary by ROI_section.tif']);

