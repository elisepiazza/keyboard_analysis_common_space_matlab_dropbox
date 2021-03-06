%ep_computeISFC
%For each subject, for each ROI, for each of the 4 scrambled conditions + 3 control conditions,
%compute inter-subject correlation between the subjects' ROI-avg'd time series and all ROIs in each other subject
%For each condition, for each subject, extract their rep-averaged ROI x TR data and correlate with avg of others' ROI x TR data 

clear;
group = 'AM';
nROIs = 10;

n_cropped_TRs = 10;

%The exact reps you want to include
scramble_reps_to_include = [1]; control_reps_to_include = [1 2];

preproc_type = 'AFNI'; preproc_params = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

all_subjects = [103 105 108 115 117 120 121 122 123];
groups = {'AM', 'M', 'M', 'AM', 'M', 'AM', 'M', 'M', 'AM'};

subjects = all_subjects(find(strcmp(groups,group))); nSubs = length(subjects);

all_ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', 'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC', 'Hipp', 'lTPJ', 'rTPJ', 'PMC', 'V1'}; 
ROIs = all_ROIs(1:nROIs);

if nROIs == 10
    ROI_order = [9 3 4 5 6 1 7 8 10 2];
elseif nROIs == 15
    ROI_order = [9 3 4 5 12 13 6 1 7 14 8 10 11 15 2];
end

filepath = ['../../common_space_AFNI/reshaped_by_conditions/' preproc_params '/nROIs=' num2str(nROIs) '/sub-'];
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

%Compute ISFC
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
        
        %Correlate the current subject's ROI x TR data with the average ROI
        %x TR data across the other subjects
        ISFC_mat_scramble(:,:,cond,s) = corr(currSubData',avg_otherSubsData'); 
        
        %save a null distribution of r values for this channel and subject
        nPerm = 100;
        parfor p = 1:nPerm
            thisdata_phasescrambled = phase_scramble_time_series(currSubData',0); %make P (1000) phase-scrambled time series
            ISFC_mat_scramble_all(:,:,cond,s,p) = corr(thisdata_phasescrambled, avg_otherSubsData); %save the r value between this subject and channel's scrambled time series and avg of others (non-scrambled)
        end
    
    end    
end

% %For control conditions
% for cond = 1:n_control_cond
%     for s = 1:nSubs
%         otherSubs = setdiff(1:nSubs,s);
%         
%         %For this subject, extract the rep-averaged (ROI x TR) data for this condition
%         currSubData = mean(data_ROIavg_control_allSubs(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,cond,control_reps_to_include,s),4);
%         
%         %Average the equivalent (ROI x TR) data across the other N subjects
%         otherSubsData = mean(data_ROIavg_control_allSubs(ROI_order,n_cropped_TRs+1:end-n_cropped_TRs,cond,control_reps_to_include,otherSubs),4);
%         avg_otherSubsData = mean(otherSubsData,5);
%         
%         ISFC_mat_control(:,:,cond,s) = corr(currSubData',avg_otherSubsData');
%                 
%     end    
% end

% %For each scramble condition, plot the group-averaged ISFC matrix
% figsize = [100 100 1100 250]; 
% figure('Units', 'pixels', 'Position', figsize);
% subplot(1,4,1); imagesc(mean(ISFC_mat_scramble(:,:,1,:),4)); title('1B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]);
% subplot(1,4,2); imagesc(mean(ISFC_mat_scramble(:,:,2,:),4)); title('2B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]);
% subplot(1,4,3); imagesc(mean(ISFC_mat_scramble(:,:,3,:),4)); title('8B'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]);
% subplot(1,4,4); imagesc(mean(ISFC_mat_scramble(:,:,4,:),4)); title('I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]);  
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC (scramble, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(scramble_reps_to_include) '.tif']);
% 
% %For each control condition, plot the group-averaged ISFC matrix
% figsize = [100 100 800 250]; 
% figure('Units', 'pixels', 'Position', figsize);
% subplot(1,3,1); imagesc(mean(ISFC_mat_control(:,:,1,:),4)); title('I_N'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]);
% subplot(1,3,2); imagesc(mean(ISFC_mat_control(:,:,2,:),4)); title('I_A'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]);
% subplot(1,3,3); imagesc(mean(ISFC_mat_control(:,:,3,:),4)); title('I_I'); xlabel('ROIs'); ylabel('ROIs'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); caxis([-.1 .4]); 
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC (control, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(control_reps_to_include) '.tif']);
% 

%For each scramble condition, plot schemaball figure
% lineColor = [0 1 0];
% nodeColor = [0 0 1];
% 
% h = schemaball(mean(ISFC_mat_scramble(:,:,1,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (scramble, 1B, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(scramble_reps_to_include) '.tif']);
% 
% h = schemaball(mean(ISFC_mat_scramble(:,:,2,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (scramble, 2B, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(scramble_reps_to_include) '.tif']);
% 
% h = schemaball(mean(ISFC_mat_scramble(:,:,3,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (scramble, 8B, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(scramble_reps_to_include) '.tif']);
% 
% h = schemaball(mean(ISFC_mat_scramble(:,:,4,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (scramble, I, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(scramble_reps_to_include) '.tif']);
% 
% 
% %For each control condition, plot schemaball figure
% h = schemaball(mean(ISFC_mat_control(:,:,1,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (control, I_N, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(control_reps_to_include) '.tif']);
% 
% h = schemaball(mean(ISFC_mat_control(:,:,2,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (control, I_A, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(control_reps_to_include) '.tif']);
% 
% h = schemaball(mean(ISFC_mat_control(:,:,3,:),4), ROIs(ROI_order), lineColor, nodeColor); set(gca, 'FontSize', 16, 'FontName', 'Helvetica'); 
% set(h.l(~isnan(h.l)), 'LineWidth', 3); set(h.s, 'MarkerEdgeColor', 'blue', 'LineWidth', 2, 'SizeData', 100); set(gcf, 'InvertHardcopy', 'off');
% print(gcf, '-dtiff', ['../figures/ISFC/ISFC_schemaball (control, I_I, ' group ' group)_nROIs=' num2str(nROIs) '_nTRs_cropped=' num2str(n_cropped_TRs) '_rep' num2str(control_reps_to_include) '.tif']);
% 


% %Cluster the matrices
% group_avg_1B = mean(ISFC_mat_scramble(:,:,1,:),4);
% z = linkage(group_avg_1B,'ward');
% c = cluster(z,'Maxclust',4);
% scatter3(group_avg_1B(:,1),group_avg_1B(:,2),group_avg_1B(:,3),10,c)
% 
% load fisheriris
% Z = linkage(meas,'average','chebychev');
% T = cluster(Z,'maxclust',3);
% cutoff = median([Z(end-2,3) Z(end-1,3)]);
% dendrogram(Z,'ColorThreshold',cutoff)

