%ep_computeIRC
%For each subject, for each ROI, for each of the 4 scrambled conditions, compute inter-rep
%correlation between reps of the same condition vs. reps of different
%conditions

clear;
subjects = [103 105 108 115 117 120 121 122 123]; nSubs = length(subjects);

preproc_types = {'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'AFNI', 'Python', 'Python', 'Python'};
preproc_params = {'v1_original_regressors', 'v2_jamals_regressors', 'v3_jamals_regressors_smoothing=1', ...
    'v4_jamals_regressors_smoothing=1_defaultGMmask', 'v5_jamals_regressors_smoothing=1_defaultGMmask_polort=3', ...
    'v6_jamals_regressors_smoothing=1_defaultGMmask_polort=2', 'v7_15_regressors_no_smoothing_defaultGMmask_polort=2', ...
    'HPF=.01Hz', 'HPF=.03Hz', 'HPF=.06Hz'};

n_cropped_TRs = 10;

avg_corr_within = zeros(nSubs, length(preproc_types));
avg_corr_between = zeros(nSubs, length(preproc_types));

for p = 1:length(preproc_types)
    
    preproc_type = preproc_types{p};
    preproc_param = preproc_params{p};
    
    for s = 1:nSubs
        subject = subjects(s);
        
        load(['../../common_space_' preproc_type '/reshaped_by_conditions/' preproc_param '/sub-' num2str(subject) '.mat']);
        n_scramble_cond = size(data_ROIavg_scramble,3); n_scramble_reps = size(data_ROIavg_scramble,4);
        choose = @(samples) samples(randi(numel(samples)));
        
        nROIs = length(ROIs);
        
        %Crop N TRs from beginning and end
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
            
        end
        
        
        if strcmp(preproc_type, 'AFNI')
            figsize = [100 100 400 500];
        elseif strcmp(preproc_type, 'Python')
            figsize = [100 100 400 350];
        end
        
        figure('Units', 'pixels', 'Position', figsize); imagesc(IRC_real_mat_scramble); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
        print(gcf, '-dtiff', ['../figures/IRCC/sub-' num2str(subject) '_Inter-rep, within-condition correlation_' preproc_type '_' preproc_param '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
        
        figure('Units', 'pixels', 'Position', figsize); imagesc(IRC_rand_mat_scramble); xlabel('Condition'); ylabel('ROI'); set(gca, 'XTickLabel', scramble_conditions, 'YTickLabel', ROIs, 'FontSize', 16, 'FontName', 'Helvetica'); colorbar; caxis([-.1 .4]);
        print(gcf, '-dtiff', ['../figures/IRCC/sub-' num2str(subject) '_Inter-rep, between-condition correlation_' '_' preproc_type '_' preproc_param '_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
        
        %Average IRC across ROIs/conditions (subject x preproc combo)        
        avg_corr_within(s,p) = mean(mean(IRC_real_mat_scramble));
        avg_corr_between(s,p) = mean(mean(IRC_rand_mat_scramble));

    end
end


%Plot avg IRC across subjects, for each preproc combo
y_within = mean(avg_corr_within);
errors_within = std(avg_corr_within)/sqrt(nSubs);

y_between = mean(avg_corr_between);
errors_between = std(avg_corr_between)/sqrt(nSubs);

y = horzcat(y_within', y_between');
errors = horzcat(errors_within', errors_between');

figsize = [100 100 800 400]; barwidth = 1;
figure('Units', 'pixels', 'Position', figsize);

b = bar(y, 'grouped', 'BarWidth', barwidth); set(b(1), 'FaceColor', [.7 .1 .5]); set(b(2), 'FaceColor', 'b');
hold on; nbars = size(y, 2);

% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',y,errors,'k.','linestyle','none', 'LineWidth' ,1);
hold off;

xlabel('Preprocessing Param Type'); ylabel('Inter-rep correlation'); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
print(gcf, '-dtiff', ['../figures/IRC/Summary stats_nTRs_cropped=' num2str(n_cropped_TRs) '.tif']);
legend({'Within-condition', 'Between-condition'});
