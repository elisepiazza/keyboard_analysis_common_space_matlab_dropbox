%ep_computeISFC
%For each subject, for each ROI, for each of the 4 scrambled conditions + 3 control conditions,
%compute inter-subject correlation between the subjects' ROI-avg'd time
%series and all ROIs in each other subject
%1B (avg across reps), s1, correlate s1's 10 ROIs (ROI x TR data) with s2's 10 ROIs

clear;
group = 'M';
n_cropped_TRs = 0;

%The exact reps you want to include
scramble_reps_to_include = [1 2 3]; control_reps_to_include = [1 2];

preproc_type = 'AFNI'; preproc_params = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

all_subjects = [103 105 108 115 117 120 121 122 123];
groups = {'AM', 'M', 'M', 'AM', 'M', 'AM', 'M', 'M', 'AM'};

subjects = all_subjects(find(strcmp(groups,group))); nSubs = length(subjects);

ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', 'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC'}; nROIs = length(ROIs);
ROI_order = [9 3 4 2 5 6 1 7 8 10];

filepath = ['../../common_space_AFNI/reshaped_by_conditions/' preproc_params '/sub-'];
barcolor = [.9 .5 0];

nTRs = 148; 

%Total # of conditions and reps
n_scramble_cond = 4; n_scramble_reps = 3;
n_control_cond = 3; n_control_reps = 2;

%Initialize empty giant data matrices (ROI x TR x cond x rep x sub)
data_ROIavg_scramble_allSubs = zeros(nROIs,nTRs,n_scramble_cond,n_scramble_reps,nSubs);
data_ROIavg_control_allSubs = zeros(nROIs,nTRs,n_control_cond,n_control_reps,nSubs);

%Load data from all subs into giant matrices
for s = 1:nSubs
    load([filepath num2str(subjects(s)) '.mat']);
    
    data_ROIavg_scramble_allSubs(:,:,:,:,s) = data_ROIavg_scramble;
    data_ROIavg_control_allSubs(:,:,:,:,s) = data_ROIavg_control;
end

%Initialize empty ISFC matrices
ISFC_mat_scramble = zeros(nROIs,nROIs,n_scramble_cond,nSubs);
ISFC_mat_control = zeros(nROIs,nROIs,n_control_cond,nSubs);

%Compute ISC
%1B (avg across reps), s1, correlate s1's ROI x TR data with avg of others' ROI x TR data

%For scramble conditions
for cond = 1:n_scramble_cond
    for s = 1:nSubs
        otherSubs = setdiff(1:nSubs,s);
        
        %For this subject, extract the rep-averaged (ROI x TR) data for this condition
        currSubData = mean(data_ROIavg_scramble_allSubs(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,cond,scramble_reps_to_include,s),4);
        
        %Average the equivalent (ROI x TR) data across the other N subjects
        otherSubsData = mean(data_ROIavg_scramble_allSubs(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,cond,scramble_reps_to_include,otherSubs),4);
        avg_otherSubsData = mean(otherSubsData,5);
        
        ISFC_mat_scramble(:,:,cond,s) = corr(currSubData',avg_otherSubsData');                
    end    
end

%For control conditions
for cond = 1:n_control_cond
    for s = 1:nSubs
        otherSubs = setdiff(1:nSubs,s);
        
        %For this subject, extract the rep-averaged (ROI x TR) data for this condition
        currSubData = mean(data_ROIavg_control_allSubs(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,cond,control_reps_to_include,s),4);
        
        %Average the equivalent (ROI x TR) data across the other N subjects
        otherSubsData = mean(data_ROIavg_control_allSubs(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,cond,control_reps_to_include,otherSubs),4);
        avg_otherSubsData = mean(otherSubsData,5);
        
        ISFC_mat_control(:,:,cond,s) = corr(currSubData',avg_otherSubsData');
                
    end    
end

%For each scramble condition, plot the group-averaged ISFC matrix
figsize = [100 100 2000 300]; 
figure('Units', 'pixels', 'Position', figsize);
subplot(1,4,1); imagesc(mean(ISFC_mat_scramble(:,:,1,:),4)); title('1B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
subplot(1,4,2); imagesc(mean(ISFC_mat_scramble(:,:,2,:),4)); title('2B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
subplot(1,4,3); imagesc(mean(ISFC_mat_scramble(:,:,3,:),4)); title('8B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
subplot(1,4,4); imagesc(mean(ISFC_mat_scramble(:,:,4,:),4)); title('I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar;  caxis([-.1 .4]);
print(gcf, '-dtiff', ['../figures/ISFC/ISFC (scramble, ' group ' group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);

%For each control condition, plot the group-averaged ISFC matrix
figsize = [100 100 1000 300]; 
figure('Units', 'pixels', 'Position', figsize);
subplot(1,3,1); imagesc(mean(ISFC_mat_control(:,:,1,:),4)); title('I_N'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
subplot(1,3,2); imagesc(mean(ISFC_mat_scramble(:,:,2,:),4)); title('I_A'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
subplot(1,3,3); imagesc(mean(ISFC_mat_scramble(:,:,3,:),4)); title('I_I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
print(gcf, '-dtiff', ['../figures/ISFC/ISFC (control, ' group ' group)_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
