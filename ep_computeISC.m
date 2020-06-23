%ep_computeISC
%For each subject, for each ROI, for each of the 4 scrambled conditions + 3 control conditions,
%compute inter-subject correlation between the subjects' ROI-avg'd time series

clear;
group = 'AM';

%The exact reps you want to include
scramble_reps_to_include = [3]; control_reps_to_include = [2]; 

preproc_type = 'AFNI'; preproc_params = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';

n_cropped_TRs = 10;

all_subjects = [103 105 108 115 117 120 121 122 123]; 
groups = {'AM', 'M', 'M', 'AM', 'M', 'AM', 'M', 'M', 'AM'}; 

subjects = all_subjects(find(strcmp(groups,group))); nSubs = length(subjects);

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

nTRs = 148; nROIs = length(ROIs);

%Total # of conditions and reps
n_scramble_cond = 4; n_scramble_reps = 3;
n_control_cond = 3; n_control_reps = 2;

%Initialize empty giant data matrices (ROI x TR x cond x rep x sub)
data_ROIavg_scramble_allSubs = zeros(nROIs,nTRs,n_scramble_cond,n_scramble_reps,nSubs);
data_ROIavg_control_allSubs = zeros(nROIs,nTRs,n_control_cond,n_control_reps,nSubs);

%Load data from all subs into giant matrices
for s = 1:nSubs
    if strcmp(preproc_type,'AFNI')
        load([filepath num2str(subjects(s)) '.mat']);
        
    elseif strcmp(preproc_type,'Python')
        load([filepath num2str(subjects(s)) '.mat']);
    end
            
    data_ROIavg_scramble_allSubs(:,:,:,:,s) = data_ROIavg_scramble; 
    data_ROIavg_control_allSubs(:,:,:,:,s) = data_ROIavg_control; 
end

%Initialize empty ISC matrices
ISC_mat_scramble = zeros(nROIs,n_scramble_cond,nSubs);
ISC_mat_control = zeros(nROIs,n_control_cond,nSubs);

%Compute ISC
for ROI = 1:nROIs
    
    %For scramble conditions
    for cond = 1:n_scramble_cond
        for s = 1:nSubs
            otherSubs = setdiff(1:nSubs,s);
            %For this subject, extract the rep-averaged time series for
            %this ROI, this condition
            currSubData = mean(data_ROIavg_scramble_allSubs(ROI,n_cropped_TRs+1:end-n_cropped_TRs,cond,scramble_reps_to_include,s),4);
            
            %Average the equivalent time series across the other N subjects
            otherSubsData = mean(data_ROIavg_scramble_allSubs(ROI,n_cropped_TRs+1:end-n_cropped_TRs,cond,scramble_reps_to_include,otherSubs),4);
            avg_otherSubsData = mean(otherSubsData,5);
            
            [ISC_r, ISC_p] = corrcoef(currSubData, avg_otherSubsData);
            ISC_mat_scramble(ROI,cond,s) = ISC_r(2,1);
        end
        
%         figure('Units', 'pixels', 'Position', [100 100 1000 375]);
%         plot(currSubData,'LineWidth',2); hold on; plot(squeeze(otherSubsData),'LineWidth',2);
%         xlabel('TR'); ylabel('BOLD'); title([ROIs{ROI} ' time series (' conditions{cond} '), ' group ' group, crop=' num2str(n_cropped_TRs)]); ylim([-5 5]); set(gca,'FontSize',16);
%         print(gcf, '-dtiff', ['../figures/Time series (' conditions{cond} '_' ROIs{ROI} ')_' group ' group, crop=' num2str(n_cropped_TRs) '_' AFNI_data_type '.tif']);       
    end
    
     %For control conditions
    for cond = 1:n_control_cond
        for s = 1:nSubs
            otherSubs = setdiff(1:nSubs,s);
            %For this subject, extract the rep-averaged time series for
            %this ROI, this condition
            currSubData = mean(data_ROIavg_control_allSubs(ROI,n_cropped_TRs+1:end-n_cropped_TRs,cond,control_reps_to_include,s),4);
            
            %Average the equivalent time series across the other N subjects
            otherSubsData = mean(data_ROIavg_control_allSubs(ROI,n_cropped_TRs+1:end-n_cropped_TRs,cond,control_reps_to_include,otherSubs),4);
            avg_otherSubsData = mean(otherSubsData,5);
            
            [ISC_r, ISC_p] = corrcoef(currSubData, avg_otherSubsData);
            ISC_mat_control(ROI,cond,s) = ISC_r(2,1);
        end
        
%         figure('Units', 'pixels', 'Position', [100 100 1000 375]);
%         plot(currSubData,'LineWidth',2); hold on; plot(squeeze(otherSubsData),'LineWidth',2);
%         xlabel('TR'); ylabel('BOLD'); title([ROIs{ROI} ' time series (' conditions{cond} '), ' group ' group, crop=' num2str(n_cropped_TRs)]); ylim([-5 5]); set(gca,'FontSize',16);
%         print(gcf, '-dtiff', ['../figures/Time series (' conditions{cond} '_' ROIs{ROI} ')_' group ' group, crop=' num2str(n_cropped_TRs) '_' AFNI_data_type '.tif']);       
    end
    
    
    %Plot ISC for scramble conditions
    x = 1:n_scramble_cond; 
    %Reshape into sub x cond data
    data_to_plot = squeeze(ISC_mat_scramble(ROI,:,:))'; 
    y = mean(data_to_plot);
    errors = std(data_to_plot)/sqrt(nSubs);
    
    figsize = [100 100 400 375]; barwidth = .5; 
    figure('Units', 'pixels', 'Position', figsize);
    bar(x,y,barwidth,'facecolor', barcolor); hold on;
    errorbar(x,y,errors,'k.', 'LineWidth', 1)
    
    xlabel('Condition'); ylabel('ISC by condition (r)'); title([ROIs{ROI}]); xlim([.3 4.7]); ylim([0 .6]); set(gca, 'XTickLabel', scramble_conditions, 'FontSize', 16, 'FontName', 'Helvetica');
    print(gcf, '-dtiff', ['../figures/ISC/ISC by condition (Scramble conds), rep1, ' group ' group, (' ROIs{ROI} ')_' preproc_type '_' preproc_params '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);

    
    %Plot ISC for control conditions
    x = 1:n_control_cond; 
    %Reshape into sub x cond data
    data_to_plot = squeeze(ISC_mat_control(ROI,:,:))'; 
    y = mean(data_to_plot);
    errors = std(data_to_plot)/sqrt(nSubs);
%     y = nanmean(ISC_mat(ROI,:,:),3);
%     errors = nanstd(ISC_mat(ROI,:,:),3)/sqrt(N);
    
    figsize = [100 100 400 375]; barwidth = .5; 
    figure('Units', 'pixels', 'Position', figsize);
    bar(x,y,barwidth,'facecolor', barcolor); hold on;
    errorbar(x,y,errors,'k.', 'LineWidth', 1)
    
    xlabel('Condition'); ylabel('ISC by condition (r)'); title([ROIs{ROI}]); xlim([.3 4.7]); ylim([0 .6]); set(gca, 'XTickLabel', control_conditions, 'FontSize', 16, 'FontName', 'Helvetica');
    print(gcf, '-dtiff', ['../figures/ISC/ISC by condition (Control conds), rep1, ' group ' group, (' ROIs{ROI} ')_' preproc_type '_' preproc_params '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);

    
end
