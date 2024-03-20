% Analyze the explained variance and weight differences for the ISH component across node designations (SOZ,
% PZ, NIZ)
% Authored by Derek Doss and Jared Shless

clear all
close all
clc

addpath Z:\shared_toolboxes\Violinplot-Matlab-master\

% load the struct that contains the functional connectivity and PCA data
load("Z:\PCA_project\data\brain_revision\with_delta\pats_strict_5.mat");

% bands to analyze
band_names = ["delta", "theta", "alpha", "beta", "gamma_low", "gamma_high"];

% where to store the figures
out_dir = "Z:\PCA_project\code\brain_revision_code\freq_band_recreation\fig3_weight_separation_delta";

% xlsx to save all the stats
xlsx_fpath = sprintf("%s\\all_stats.xlsx", out_dir);
xlsx_table = table();
row_labels = [""];

counter = 1;
for ii = 1:length(band_names)
    row_labels(counter, 1) = sprintf("%s", band_names(ii));
    counter = counter + 1;
end

xlsx_table.("Group") = row_labels;

% array to store the statistical data
stats_table_all = nan(length(band_names), 6);

% loop through each band
for kk = 1:length(band_names)

    % SOZ - NIZ weight difference for the ISH component
    soz_niz_weight_diff_pc1 = [];
    soz_niz_weight_diff_pc2 = [];
    
    % explained variance of PC1 cohort
    var_pc1 = [];
    % explained variance of PC2+ cohort
    var_pc2 = [];

    % loop through patients
    for D = 1:size(pats,2) 
        
        % indices of SOZ nodes
        index_SOZ = find(pats(D).SOZ == 1);

        % indices of NIZ nodes
        index_NIZ = find(pats(D).SOZ == 0);
        index_NIZ = sort([index_NIZ;find(pats(D).SOZ == 3)]);
        
        % position of the ISH component
        ish_pos = pats(D).ish_component_position.(band_names(kk));
        
        % weights 
        weights = pats(D).weights.(band_names(kk));
        
        % temporary variable
        explained_tmp_var = pats(D).explained.(band_names(kk));
        
        % PC 1 cohort
        if pats(D).ish_component_position.(band_names(kk)) == 1 
            % SOZ - NIZ weight differences for the ISH component
            soz_niz_weight_diff_pc1 = [soz_niz_weight_diff_pc1, abs(mean(weights(index_SOZ, ish_pos)) - mean(weights(index_NIZ, ish_pos)))];
            
            % Explained variance
            var_pc1 = [var_pc1, explained_tmp_var(ish_pos)];
        % PC 2 cohort
        elseif pats(D).ish_component_position.(band_names(kk)) ~= 1 
            % SOZ - NIZ weight differences for the ISH component
            soz_niz_weight_diff_pc2 = [soz_niz_weight_diff_pc2, abs(mean(weights(index_SOZ, ish_pos)) - mean(weights(index_NIZ, ish_pos)))];
            
            % explained variance
            var_pc2 = [var_pc2, explained_tmp_var(ish_pos)];
        end 
    end 
    
    
    %% Weight Differences Plots
    
    A = [soz_niz_weight_diff_pc1]'; A = A(~isnan(A));
    B = [soz_niz_weight_diff_pc2]'; B = B(~isnan(B));
    
    [~, p, ci, stats] = ttest2(A, B);
    
    % save the p-values so we can correct for multiple comparisons
    stats_table_all(kk, 1) = p;
    % lower CI
    stats_table_all(kk, 2) = ci(1);
    % higher CI
    stats_table_all(kk, 3) = ci(2);
    
    % plot the weight seperation
    violin_struct = struct;
    violin_struct.ISHA = A; violin_struct.ISHB= B;
    
    f = figure('Position',[0 0 600 800]);
    ylim([-12 30]);
    violins = violinplot(violin_struct);
    set(gcf,'color','w');

    title(sprintf("SOZ-NIZ Weights of ISH Component, %s", band_names(kk)));
    subtitle(num2str(p));
    ylabel('Weight Difference')
    ax = gca; 
    ax.FontSize = 18; 
    ax.LineWidth = 3;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.04];
       
    violins(1,1).ViolinColor = [0.251 0.192 0.459];
    violins(1,2).ViolinColor =  [0.165 0.486 0.286];

    set(gcf,'color','w');
    saveas(f, sprintf("%s\\weight_diffs_%s.png", out_dir, band_names(kk)));
    saveas(f, sprintf("%s\\weight_diffs_%s.svg", out_dir, band_names(kk)));
    close(f);

    %% violin plot of explained variance
    
    
    A = [var_pc1]';
    B = [var_pc2]';
    
    [~, p, ci, stats] = ttest2(A, B);

    % save the p-values so we can correct for multiple comparisons
    stats_table_all(kk, 4) = p;
    % lower CI
    stats_table_all(kk, 5) = ci(1);
    % higher CI
    stats_table_all(kk, 6) = ci(2);
    
    violin_struct = struct;
    violin_struct.ISHA = A; 
    violin_struct.ISHB= B;
    
    f = figure('Position',[0 0 600 800]);
    ylim([-10 60]);
    violins = violinplot(violin_struct);
    set(gcf,'color','w');
    % 
    title(sprintf("Explained Variance for PC1 & PC2, %s band", band_names(kk)));
    subtitle(p);
    ylabel('Explained Variance')
    ax = gca; 
    ax.FontSize = 18; 
    ax.LineWidth = 3;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.04];
    
    violins(1,1).ViolinColor = [0.251 0.192 0.459];
    violins(1,2).ViolinColor =  [0.165 0.486 0.286];
    
    saveas(f, sprintf("%s\\explained_var_%s.png", out_dir, band_names(kk)));
    saveas(f, sprintf("%s\\explained_var_%s.svg", out_dir, band_names(kk)));
    close(f)

end

% save the stats
xlsx_table.("WeightDiffs") = stats_table_all(:, 1);
xlsx_table.("WeightDiffsLowCI") = stats_table_all(:, 2);
xlsx_table.("WeightDiffsHighCI") = stats_table_all(:, 3);
xlsx_table.("ExplainedVar") = stats_table_all(:, 4);
xlsx_table.("ExplainedVarLowCI") = stats_table_all(:, 5);
xlsx_table.("ExplainedVarHighCI") = stats_table_all(:, 6);

writetable(xlsx_table, "Z:\PCA_project\code\brain_revision_code\freq_band_recreation\fig3_weight_separation_delta\fig3_stats.xlsx");