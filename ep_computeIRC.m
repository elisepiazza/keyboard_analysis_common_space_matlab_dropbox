%ep_computeIRC
%For each subject, for each ROI, for each of the 4 scrambled conditions, compute inter-rep
%correlation between reps of the same condition vs. reps of different
%conditions

clear;
subject = 103;

data_type = 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2';
load(['../reshaped_by_conditions/' data_type '/sub-' num2str(subject) '.mat']);
n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
choose = @(samples) samples(randi(numel(samples)));

nROIs = length(ROIs);

%Crop N TRs from beginning and end
n_cropped_TRs = 0; 
data_ROIavg_scramble = data_ROIavg_scramble(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);
data_ROIavg_control = data_ROIavg_control(:,n_cropped_TRs+1:end-n_cropped_TRs,:,:);

%Initialize empty matrices of IRCs
IRC_real_mat_scramble = zeros(nROIs,n_scramble_cond); 
IRC_rand_mat_scramble = zeros(nROIs,n_scramble_cond);

for ROI = 1:nROIs
    
    %Scrambles only
    for cond = 1:n_scramble_cond
        
%         %Plot ROI time course for each rep
%         figure('Units', 'pixels', 'Position', [100 100 1000 375]);
%         colors = {[1 0 0], [1 .5 0], [0 1 0], [0 .5 1], [0 0 1]};
%         for rep = 1:n_scramble_reps
%             plot(data_ROIavg_scramble(ROI,:,cond,rep), 'color', colors{rep}); hold on; plot(nanmean(data_ROIavg_scramble(ROI,:,cond,:),4),'k','LineWidth',2);
%             xlabel('TR'); ylabel('BOLD'); title([ROIs{ROI} ' time series (' scramble_conditions{cond} ')']); ylim([-10 10]); set(gca, 'FontSize', 16);
%         end
% %             print(gcf, '-dtiff', ['../figures/s' num2str(subject) '/Time series by rep (' conditions{cond} '_' ROIs{ROI} smooth_tags{smoothOrNot} ').tif']);
        
        %Compute IRC (for each rep of a given condition, correlate the average ROI time series for that rep w/ the average across the other reps
        for rep = 1:n_scramble_reps
            others = nanmean(data_ROIavg_scramble(ROI,:,cond,setdiff([1:n_scramble_reps],rep)),4);
            rs_and_ps = corrcoef(data_ROIavg_scramble(ROI,:,cond,rep),others);
            IRC_real(rep,cond) = rs_and_ps(2,1);
            
            %Pick a random condition that's not this one
            rand_cond = choose(setdiff(1:length(scramble_conditions),cond)); 
            %Average the reps for another condition (and the other reps) that's not this one
            others_rand_cond = nanmean(data_ROIavg_scramble(ROI,:,rand_cond,setdiff([1:n_scramble_reps],rep)),4); 
            rs_and_ps_rand = corrcoef(data_ROIavg_scramble(ROI,:,cond,rep),others_rand_cond);
            IRC_rand(rep,cond) = rs_and_ps_rand(2,1);
        end
        
    end
    
    IRC_real_mat_scramble(ROI,:) = mean(IRC_real);
    IRC_rand_mat_scramble(ROI,:) = mean(IRC_rand);
    
    
%     %Plot ISC for scramble conditions
%     N = 3;
%     x = 1:4;
%     y = nanmean(ISC(:,1:4));
%     errors = nanstd(ISC(:,1:4))/sqrt(N);
%     
%     figsize = [100 100 400 375]; barwidth = .5; barcolor = [.9 .5 0];
%     figure('Units', 'pixels', 'Position', figsize);
%     bar(x,y,barwidth,'facecolor', barcolor); hold on;
%     errorbar(x,y,errors,'k.', 'LineWidth', 1)
%     
%     xticklab = conditions;
%     xlabel('Condition'); ylabel('Cross-rep ISC by condition (r)'); title([ROIs{ROI}]); xlim([.3 7.7]); ylim([0 .6]); set(gca, 'XTickLabel', xticklab, 'FontSize', 16, 'FontName', 'Helvetica');
%     % print(gcf, '-dtiff', ['../figures/s' num2str(subject) '/Cross-rep ISC by condition (' ROIs{ROI} smooth_tags{smoothOrNot} ').tif']);
%     
%     %Plot ISC for all conditions
%     N = nReps;
%     x = 1:4;
%     y = nanmean(ISC_rand(:,1:4));
%     errors = nanstd(ISC_rand(:,1:4))/sqrt(N);
%     
%     figsize = [100 100 400 375]; barwidth = .5; barcolor = [.5 .9 0];
%     figure('Units', 'pixels', 'Position', figsize);
%     bar(x,y,barwidth,'facecolor', barcolor); hold on;
%     errorbar(x,y,errors,'k.', 'LineWidth', 1)
%     
%     xticklab = conditions;
%     xlabel('Condition'); ylabel('Cross-cond, cross-rep ISC by condition (r)'); title([ROIs{ROI}]); xlim([.3 7.7]); ylim([0 .6]); set(gca, 'XTickLabel', xticklab, 'FontSize', 16, 'FontName', 'Helvetica');
%     
end

figsize = [100 100 400 500]; figure('Units', 'pixels', 'Position', figsize); imagesc(IRC_real_mat_scramble); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
print(gcf, '-dtiff', ['../figures/sub-' num2str(subject) '/Inter-rep, within-condition correlation_' data_type '.tif']);


figsize = [100 100 400 500]; figure('Units', 'pixels', 'Position', figsize); imagesc(IRC_rand_mat_scramble); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
print(gcf, '-dtiff', ['../figures/sub-' num2str(subject) '/Inter-rep, between-condition correlation_' data_type '.tif']);


% %Plot ISC for each successive run for melody 1, melody 2
% figure('Units', 'pixels', 'Position', [100 100 900 375]);
% subplot(1,2,1); plot(corr_melody1, 'bo', 'MarkerSize', 16, 'MarkerFaceColor', 'b'); xlim([0 6]); title([ROI_names{whichROI} smooth_tags{smoothOrNot} ' (Melody 1)']); xlabel('Run'); ylabel('ISC (btwn this run and others)'); set(gca, 'FontSize', 16);
% subplot(1,2,2); plot(corr_melody2, 'bo', 'MarkerSize', 16, 'MarkerFaceColor', 'b'); xlim([0 6]); title([ROI_names{whichROI} smooth_tags{smoothOrNot} ' (Melody 2)']); xlabel('Run'); ylabel('ISC (btwn this run and others'); set(gca, 'FontSize', 16);
% print(gcf, '-dtiff', ['../figures/ISC by run (' ROI_names{whichROI} smooth_tags{smoothOrNot} ').tif']);


