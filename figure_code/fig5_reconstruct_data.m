% Analyze the reconstructed original data
% Authored by Derek Doss and Jared Shless

clear all
close all
clc

restoredefaultpath

addpath Z:\shared_toolboxes\Violinplot-Matlab-master\

% load in the data_struct
load("Z:\PCA_project\data\brain_revision\pats_strict_5.mat");

% bands to analyze
band_names = ["delta", "theta", "alpha", "beta", "gamma_low", "gamma_high"];

% where to save the figures
out_dir = "Z:\PCA_project\code\brain_revision_code\freq_band_recreation\weight_separation_violin";

% save the p-values so we can correct for multiple comparisons
anova_pvals = nan(length(band_names),3);
anova_counter = 1;

% which components to reconstruct with
components_to_recon = ["ISH", "ExceptISH", "All"];

% cohorts to analyze, aka PC1, PC2+, or all
cohorts_to_analyze = ["All", "PC1", "PC2"];

% loop through bands
for kk = 1:length(band_names)

    % frequency band
    band_str = band_names(kk);
    for mm = 1:length(components_to_recon)

        % in - out strength
        soz_inout = [];
        pz_inout =[]; 
        niz_inout = [];
    

        % ISH component position for each patient
        ish_component_position = [];

        % loop through each patient
        for T = 1:size(pats,2) 
            
            % store the ish component position for each patient
            ish_component_position = [ish_component_position, pats(T).ish_component_position.(band_str)];

            % All indices
            all_indices = 1:length(pats(T).SOZ)-1;

            % ISH index
            ish_pos = pats(T).ish_component_position.(band_str);
            
            % All except ISH indices
            non_ish_ind = setdiff(all_indices, ish_pos);

            % what indicies should be used
            inds_to_use = [];

            switch components_to_recon(mm)
                case "ISH"
                    inds_to_use = ish_pos;
                case "ExceptISH"
                    inds_to_use = non_ish_ind;
                case "All"
                    inds_to_use = all_indices;
            end
            
            % All weights and coefficients
            weight_data = pats(T).weights.(band_str);
            coeff_data = pats(T).coeff_matrix.(band_str);
    
            % Weights and coefficients to use
            weights_to_analyze = weight_data(:,inds_to_use);
            coeffs_to_analyze = coeff_data(:,inds_to_use);
        
            % Reconstructed matrices
            recon_data = coeffs_to_analyze * weights_to_analyze'; % 2n x n
            recon_data_IN = recon_data(1:end/2, :); % inward coeffs for ish coeff
            recon_data_OUT = recon_data(end/2 + 1 :end,:); % outward coeffs for ish coeff
            
            % Nan the diagonals so we can do Z-scoring by the row!
            recon_data_IN(logical(eye(size(recon_data_IN)))) = nan;
            recon_data_OUT(logical(eye(size(recon_data_OUT)))) = nan; 

            % Z-Score by the row so we can mean the columns for strength! 
            recon_IN_Z = (recon_data_IN - repmat(mean(recon_data_IN, 2, 'omitnan'), [1, size(recon_data_IN, 2)])) ./ repmat(std(recon_data_IN, 0, 2, 'omitnan'), ...
                [1, size(recon_data_IN, 2)]);
            recon_OUT_Z = (recon_data_OUT - repmat(mean(recon_data_OUT, 2, 'omitnan'), [1, size(recon_data_OUT, 2)])) ./ repmat(std(recon_data_OUT, 0, 2, 'omitnan'), ...
                [1, size(recon_data_OUT, 2)]);
        
            % PZ, SOZ, and NIZ indices for this patient 
            ind_SOZ = find(pats(T).SOZ == 1); % SOZ indices
            ind_NIZ = find(pats(T).SOZ == 0); % working NIZ indices
            ind_NIZ = sort([ind_NIZ;find(pats(T).SOZ == 3)]); % NIZ indices
            ind_PZ = find(pats(T).SOZ == 2);
        
            % In - Out Strength
            soz_inout  = [soz_inout, mean(recon_IN_Z(:,ind_SOZ), 'all', 'omitnan') - mean(recon_OUT_Z(:,ind_SOZ), 'all', 'omitnan')];
            pz_inout = [pz_inout, mean(recon_IN_Z(:,ind_PZ), 'all', 'omitnan') - mean(recon_OUT_Z(:,ind_PZ), 'all', 'omitnan')];
            niz_inout = [niz_inout, mean(recon_IN_Z(:,ind_NIZ), 'all', 'omitnan') - mean(recon_OUT_Z(:,ind_NIZ), 'all', 'omitnan')];
        end

    
        % plot for each cohort
        for ii = 1:length(cohorts_to_analyze)
            
            % variable to store the cohort indices
            tmp_cohort_idxs = [];

            % figure out which patients to use
            switch cohorts_to_analyze(ii)
                case "All"
                    tmp_cohort_idxs = 1:length(ish_component_position);
                case "PC1"
                    tmp_cohort_idxs = find(ish_component_position==1);
                case "PC2"
                    tmp_cohort_idxs = find(ish_component_position~=1);
            end
            [anova_pvals(anova_counter, :)] = plot_reconstruction(soz_inout, pz_inout, niz_inout, cohorts_to_analyze(ii), tmp_cohort_idxs);
            anova_counter = anova_counter + 1;
        end

    end
   
end




function [P_anova] = plot_reconstruction(soz_inout, pz_inout, niz_inout, cohort_str, cohort_idxs)

    
    metric_str = "PDC";
    
    norm_style = "Row/Col Z-Score";
    ymin = -1.5;
    ymax = 5;
    yticks([0:1:5])
    
    tmp_soz_data = soz_inout(cohort_idxs);
    tmp_pz_data = pz_inout(cohort_idxs);
    tmp_non_data = niz_inout(cohort_idxs);
    
    % ORIGINAL INWARD - OUTWARD
    f = figure('Position', [100 100 540 900]);
    boxplot([tmp_soz_data,tmp_pz_data,tmp_non_data],'Labels',{'SOZ','PZ','Non'},'whisker',10);

    title(sprintf("%s %s\nINWARD-OUTWARD %s",...
        metric_str, cohort_str, band_str));
    ylim([ymin ymax])
    set(gcf, 'color', 'w');
    set(gca, 'TickDir', 'out');
    ylabel(sprintf('%s Inward-Outward %s Strength',norm_style, metric_str))
    hold off

    % Stats on IN-OUT bar chart
    [P_anova,ANOVATAB,STATS] = anova1([tmp_soz_data, tmp_pz_data, tmp_non_data],[],'off');
    annotation('textbox', [0.15, 0.07, 0.1, 0.1], 'String', sprintf("One-way ANOVA p=%0.2e",P_anova))
    [c,m,h,nms] = multcompare(STATS,'display','off');
    annotation('textbox', [0.60, 0.10, 0.1, 0.1], 'String', sprintf("SOZ-PZ p=%0.2e\nSOZ-Non p=%0.2e\nPZ-Non p=%0.2e",c(1,end),c(2,end),c(3,end)))

    
    savename = sprintf("%s_%s_InwardMOutward_%s", metric_str, band_str, cohort_str);
    saveas(f,fullfile(out_dir,savename),'eps')
    saveas(f,fullfile(out_dir,savename),'fig')
    saveas(f,fullfile(out_dir,savename),'svg')
    saveas(f,fullfile(out_dir,savename),'png')

end
