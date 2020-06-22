%ep_computeCorrClass
%For each subject, for each ROI, for each of the 4 scrambled conditions,
%make an averaged VxT matrix (loaf) across nRep-1 reps and hold out the remaining
%one. For each held-out loaf, which of the 4 other conditions is it most correlated with? (Chance = .25)
%Repeat for nReps.

%Note: since there is some randomness in which reps are held-out, this will
%generate slightly different results across runs of the script. Add a random seed. 

%TO ADD: in which areas is I_A more correlated w/ I than the other 3
%scrambles?

clear;
subjects = [103 105 108 115 117 120 121 122 123]; nSubs = length(subjects);

preproc_types = {'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'Python', 'Python', 'Python'};
preproc_params = {'v1_original_regressors', 'v2_jamals_regressors', 'v3_jamals_regressors_smoothing=1', ...
    'v4_jamals_regressors_smoothing=1_defaultGMmask', 'v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3', ...
    'v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2', 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2', ...
    'HPF=.01Hz', 'HPF=.03Hz', 'HPF=.06Hz'};

n_cropped_TRs = 10;

avg_acc_scramble = zeros(nSubs, length(preproc_types));
avg_acc_control = zeros(nSubs, length(preproc_types));

for p = 1:length(preproc_types)
    
    preproc_type = preproc_types{p};
    preproc_param = preproc_params{p};
    
    for s = 1:nSubs
        subject = subjects(s);
        
        load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
        n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
        n_control_cond = size(data_ROIavg_control,3); n_control_reps = size(data_ROIavg_control,4);
        
        nROIs = length(ROIs);
        
        %Initialize empty matrices of accuracies
        ROI_acc_scramble = zeros(nROIs,n_scramble_cond);
        ROI_acc_control = zeros(nROIs,n_control_cond);
                
        for ROI = 1:nROIs
            
            %Run classifier for scramble conditions
            data_scramble_thisROI = data_scramble{ROI}; %The extracted 4D matrix should be V x T x cond x rep
            data_scramble_thisROI = data_scramble_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);
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
                R1 = corrcoef(train_1B(:,:),test_1B(:,:)); R2 = corrcoef(train_1B(:,:),test_2B(:,:)); R3 = corrcoef(train_1B(:,:),test_8B(:,:)); R4 = corrcoef(train_1B(:,:),test_I(:,:));
                acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_1B(i) = sum(acc1) == n_scramble_reps;
                
                %Is train_2B most strongly correlated with its own held-out (test) loaf than the other 3?
                R1 = corrcoef(train_2B(:,:),test_2B(:,:)); R2 = corrcoef(train_2B(:,:),test_8B(:,:)); R3 = corrcoef(train_2B(:,:),test_I(:,:)); R4 = corrcoef(train_2B(:,:),test_1B(:,:));
                acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_2B(i) = sum(acc1) == n_scramble_reps;
                
                %Is train_8B most strongly correlated with its own held-out (test) loaf than the other 3?
                R1 = corrcoef(train_8B(:,:),test_8B(:,:)); R2 = corrcoef(train_8B(:,:),test_I(:,:)); R3 = corrcoef(train_8B(:,:),test_1B(:,:)); R4 = corrcoef(train_8B(:,:),test_2B(:,:));
                acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_8B(i) = sum(acc1) == n_scramble_reps;
                
                %Is train_1 most strongly correlated with its own held-out (test) loaf than the other 3?
                R1 = corrcoef(train_I(:,:),test_I(:,:)); R2 = corrcoef(train_I(:,:),test_1B(:,:)); R3 = corrcoef(train_I(:,:),test_2B(:,:)); R4 = corrcoef(train_I(:,:),test_8B(:,:));
                acc1 = R1(2,1) > [R2(2,1) R3(2,1) R4(2,1)]; acc_I(i) = sum(acc1) == n_scramble_reps;
            end
            
            ROI_acc_scramble(ROI,1) = mean(acc_1B);
            ROI_acc_scramble(ROI,2) = mean(acc_2B);
            ROI_acc_scramble(ROI,3) = mean(acc_8B);
            ROI_acc_scramble(ROI,4) = mean(acc_I);
            
            %Run classifier for control conditions
            data_control_thisROI = data_control{ROI}; %The extracted 4D matrix should be V x T x cond x rep
            data_control_thisROI = data_control_thisROI(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);
            held_out_reps = randperm(n_control_reps);
            
            for i = 1:n_control_reps
                
                test_rep = held_out_reps(i);
                
                %Average of I_N training reps (VxT)
                train_I_N = mean(data_control_thisROI(:,:,1,setdiff([1:n_control_reps],test_rep)),4);
                
                %Average of I_A training reps (VxT)
                train_I_A = mean(data_control_thisROI(:,:,2,setdiff([1:n_control_reps],test_rep)),4);
                
                %Average of I_I training reps (VxT)
                train_I_I = mean(data_control_thisROI(:,:,3,setdiff([1:n_control_reps],test_rep)),4);
                                
                %Held-out I_N run (VxT)
                test_I_N = data_control_thisROI(:,:,1,test_rep);
                
                %Held-out I_A run (VxT)
                test_I_A = data_control_thisROI(:,:,2,test_rep);
                
                %Held-out I_I run (VxT)
                test_I_I = data_control_thisROI(:,:,3,test_rep);
                                
                %Is train_I_N most strongly correlated with its own held-out (test) loaf than the other 2?
                R1 = corrcoef(train_I_N(:,:),test_I_N(:,:)); R2 = corrcoef(train_I_N(:,:),test_I_A(:,:)); R3 = corrcoef(train_I_N(:,:),test_I_I(:,:)); 
                acc1 = R1(2,1) > [R2(2,1) R3(2,1)]; acc_I_N(i) = sum(acc1) == n_control_reps;
                
                %Is train_I_A most strongly correlated with its own held-out (test) loaf than the other 2?
                R1 = corrcoef(train_I_A(:,:),test_I_A(:,:)); R2 = corrcoef(train_I_A(:,:),test_I_I(:,:)); R3 = corrcoef(train_I_A(:,:),test_I_N(:,:)); 
                acc1 = R1(2,1) > [R2(2,1) R3(2,1)]; acc_I_A(i) = sum(acc1) == n_control_reps;
                
                %Is train_I_I most strongly correlated with its own held-out (test) loaf than the other 2?
                R1 = corrcoef(train_I_I(:,:),test_I_I(:,:)); R2 = corrcoef(train_I_I(:,:),test_I_N(:,:)); R3 = corrcoef(train_I_I(:,:),test_I_A(:,:)); 
                acc1 = R1(2,1) > [R2(2,1) R3(2,1)]; acc_I_I(i) = sum(acc1) == n_control_reps;
                
            end
            
            ROI_acc_control(ROI,1) = mean(acc_I_N);
            ROI_acc_control(ROI,2) = mean(acc_I_A);
            ROI_acc_control(ROI,3) = mean(acc_I_I);
            
        end
        
        if strcmp(preproc_type, 'AFNI')
            figsize = [100 100 400 500];
        elseif strcmp(preproc_type, 'Python')
            figsize = [100 100 400 350];
        end
        
        %Plot classification accuracy (ROI x cond) for this subject (scramble conditions)
        figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc_scramble); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
        print(gcf, '-dtiff', ['../figures/Correlation classifier/sub-' num2str(subject) '_scramble_' preproc_type '_' preproc_param '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
        

        %Plot classification accuracy (ROI x cond) for this subject (control conditions)
        figure('Units', 'pixels', 'Position', figsize); imagesc(ROI_acc_control); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', control_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([0 1]);
        print(gcf, '-dtiff', ['../figures/Correlation classifier/sub-' num2str(subject) '_control_' preproc_type '_' preproc_param '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
        
        %Average classification accuracy across ROIs/conditions (subject x preproc combo)
        avg_acc_scramble(s,p) = mean(mean(ROI_acc_scramble));
        avg_acc_control(s,p) = mean(mean(ROI_acc_control));

    end
end

%Plot avg classification accuracy across subjects, for each preproc combo (scramble condition)
y = mean(avg_acc_scramble);
errors = std(avg_acc_scramble)/sqrt(nSubs);
x = 1:length(preproc_types);

figsize = [100 100 600 400]; barwidth = .5;
figure('Units', 'pixels', 'Position', figsize);
bar(x,y,barwidth,'facecolor', [.2 .8 .9]); hold on;
errorbar(x,y,errors,'k.', 'LineWidth', 1)

xticklab = preproc_params;
xlabel('Preprocessing Param Type'); ylabel('Correlation classifier accuracy (% correct)'); xlim([.3 10.7]); ylim([0 1]); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
print(gcf, '-dtiff', ['../figures/Correlation classifier/Summary stats_Scramble conditions_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);


%Plot avg classification accuracy across subjects, for each preproc combo (control condition)
y = mean(avg_acc_control);
errors = std(avg_acc_control)/sqrt(nSubs);
x = 1:length(preproc_types);

figsize = [100 100 600 400]; barwidth = .5;
figure('Units', 'pixels', 'Position', figsize);
bar(x,y,barwidth,'facecolor', [.2 .8 .9]); hold on;
errorbar(x,y,errors,'k.', 'LineWidth', 1)

xticklab = preproc_params;
xlabel('Preprocessing Param Type'); ylabel('Correlation classifier accuracy (% correct)'); xlim([.3 10.7]); ylim([0 1]); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
print(gcf, '-dtiff', ['../figures/Correlation classifier/Summary stats_Control conditions_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
