%ep_computeISC
%For each subject, for each ROI, for each of the 4 scrambled conditions,
%compute inter-subject correlation between the subjects' ROI-avg'd time series

clear;
group = 'AM';

preproc_types = {'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'Python', 'Python', 'Python'};
preproc_params = {'v1_original_regressors', 'v2_jamals_regressors', 'v3_jamals_regressors_smoothing=1', ...
    'v4_jamals_regressors_smoothing=1_defaultGMmask', 'v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3', ...
    'v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2', 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2', ...
    'HPF=.01Hz', 'HPF=.03Hz', 'HPF=.06Hz'};

% preproc_type = 'Python'; %'AFNI', 'Python'
% preproc_params = 'HPF=.01Hz'; 
% %AFNI parameter choices:
% %v1_original_regressors
% %v2_jamals_regressors
% %v3_jamals_regressors_smoothing=1
% %v4_jamals_regressors_smoothing=1_defaultGMmask
% %v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3
% %v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2
% %v7_15_regressors_no_smoothing_defaultGMmask_polort=2 
% %Python parameter choices:
% %HPF=.01Hz, HPF=.03Hz, HPF=.06Hz
n_cropped_TRs = 0;

all_subjects = [103 105 115 117 120 121 122 123]; 
groups = {'AM', 'M', 'AM', 'M', 'AM', 'M', 'M', 'AM'}; subjects = all_subjects(find(strcmp(groups,group))); nSubs = length(subjects);

conditions = {'1B', '2B', '8B', 'I'}; nCond = length(conditions);

if strcmp(preproc_type, 'AFNI')
    ROIs = {'AngularG', 'Cerebellum', 'HeschlsG', 'STG', 'MotorCortex', ...
    'TPJ', 'PCC', 'Precuneus', 'A1', 'mPFC'};
    filepath = ['../../common_space_AFNI/reshaped_by_conditions/' preproc_params '/sub-'];
    barcolor = [.9 .5 0];

elseif strcmp(preproc_type, 'Python')
    ROIs = {'A1', 'AngularGyrus', 'Erez-DMN', 'Precuneus', 'vmPFC'};
    filepath = ['../../common_space_Python/reshaped_by_conditions/' preproc_params '/sub-'];
    barcolor = [.3 .5 .9];
end

nTRs = 148; nReps = 3; nROIs = length(ROIs);
ISC_mat = zeros(nROIs,nCond,nSubs);
data_ROIavg_scramble_allSubs = zeros(nROIs,nTRs,nCond,nReps,nSubs);

%Load data from all subs into an ROI x TR x cond x rep x sub matrix
for s = 1:nSubs
    if strcmp(preproc_type,'AFNI')
        load([filepath num2str(subjects(s)) '.mat']);
        
    elseif strcmp(preproc_type,'Python')
        load([filepath num2str(subjects(s)) '.mat']);
    end
        
    data_ROIavg_scramble_allSubs(:,:,:,:,s) = data_ROIavg_scramble; 
end

%Compute ISC
for ROI = 1:nROIs
    for cond = 1:nCond
        for s = 1:nSubs
            otherSubs = setdiff(1:nSubs,s);
            %For this subject, extract the z-scored, rep-averaged time series for
            %this ROI, this condition
            currSubData = zscore(mean(data_ROIavg_scramble_allSubs(ROI,n_cropped_TRs+1:end-n_cropped_TRs,cond,:,s),4));
            
            %Average the equivalent time series across the other N subjects
            otherSubsData = zscore(mean(data_ROIavg_scramble_allSubs(ROI,n_cropped_TRs+1:end-n_cropped_TRs,cond,:,otherSubs),4));
            avg_otherSubsData = mean(otherSubsData,5);
            
            [ISC_r, ISC_p] = corrcoef(currSubData, avg_otherSubsData);
            ISC_mat(ROI,cond,s) = ISC_r(2,1);
        end
        
%         figure('Units', 'pixels', 'Position', [100 100 1000 375]);
%         plot(currSubData,'LineWidth',2); hold on; plot(squeeze(otherSubsData),'LineWidth',2);
%         xlabel('TR'); ylabel('BOLD'); title([ROIs{ROI} ' time series (' conditions{cond} '), ' group ' group, crop=' num2str(n_cropped_TRs)]); ylim([-5 5]); set(gca,'FontSize',16);
%         print(gcf, '-dtiff', ['../figures/Time series (' conditions{cond} '_' ROIs{ROI} ')_' group ' group, crop=' num2str(n_cropped_TRs) '_' AFNI_data_type '.tif']);       
    end
    
    %Plot ISC for scramble conditions
    x = 1:nCond; 
    %Reshape into sub x cond data
    data_to_plot = squeeze(ISC_mat(ROI,:,:))'; 
    y = mean(data_to_plot);
    errors = std(data_to_plot)/sqrt(nSubs);
%     y = nanmean(ISC_mat(ROI,:,:),3);
%     errors = nanstd(ISC_mat(ROI,:,:),3)/sqrt(N);
    
    figsize = [100 100 400 375]; barwidth = .5; 
    figure('Units', 'pixels', 'Position', figsize);
    bar(x,y,barwidth,'facecolor', barcolor); hold on;
    errorbar(x,y,errors,'k.', 'LineWidth', 1)
    
    xticklab = conditions;
    xlabel('Condition'); ylabel('ISC by condition (r)'); title([ROIs{ROI}]); xlim([.3 4.7]); ylim([0 .6]); set(gca, 'XTickLabel', conditions, 'FontSize', 16, 'FontName', 'Helvetica');
    print(gcf, '-dtiff', ['../figures/ISC/ISC by condition, ' group ' group, (' ROIs{ROI} ')_' preproc_type '_' preproc_params '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);

end
