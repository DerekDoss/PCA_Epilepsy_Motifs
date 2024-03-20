% Analyze the weights of the ISH component for SOZ, PZ, and NIZ nodes
% Authored by Derek Doss and Jared Shless

clear all
close all
clc

restoredefaultpath
addpath Z:\shared_toolboxes\Derek_functions\
addpath Z:\shared_toolboxes\Violinplot-Matlab-master\

% load in the data struct
load("Z:\PCA_project\data\brain_revision\with_delta\pats_strict_5.mat");

% where to save the plots
out_dir = "Z:\PCA_project\code\brain_revision_code\freq_band_recreation\fig4_weights_violin";
 
% what frequency bands to analyze
band_names = ["delta","theta", "alpha", "beta", "gamma_low", "gamma_high"];

% xlsx to save all the stats
xlsx_fpath = sprintf("%s\\all_stats.xlsx", out_dir);
xlsx_table = table();
row_labels = [""];

counter = 1;
for ii = 1:length(band_names)
    row_labels(counter, 1) = sprintf("%s_PC1", band_names(ii));
    counter = counter + 1;
    row_labels(counter, 1) = sprintf("%s_PC2", band_names(ii));
    counter = counter + 1;
end

xlsx_table.("Group") = row_labels;

% table to store the stats
stats_table_all = [];

% counter for stats table
counter = 1;

% loop through bands
for kk = 1:length(band_names)   

    SOZ_components = [];
    NIZ_components =[];
    PZ_components = [];
    
    SOZ_components_PC2 = [];
    PZ_components_PC2 = [];
    NIZ_components_PC2 = [];
    
    for T = 1:size(pats,2) % for all patients in the cohort
    
        component_position = pats(T).ish_component_position.(band_names(kk));

        weights_vector = pats(T).weights.(band_names(kk));

        % indices for SOZ, PZ, and NIZ nodes
        index_SOZ = find(pats(T).SOZ == 1); % SOZ indices
        index_NIZ = find(pats(T).SOZ == 0); % working NIZ indices
        index_NIZ = sort([index_NIZ;find(pats(T).SOZ == 3)]); % NIZ indices
        index_PZ = find(pats(T).SOZ == 2);

        % PC1 cohort
        if pats(T).ish_component_position.(band_names(kk)) == 1
               
            SOZ_components = [SOZ_components, mean(weights_vector(index_SOZ,component_position))];
            NIZ_components = [NIZ_components, mean(weights_vector(index_NIZ,component_position))];
            PZ_components = [PZ_components, mean(weights_vector(index_PZ,component_position))];
        % PC2 cohort
        else
        
            SOZ_components_PC2 = [SOZ_components_PC2, mean(weights_vector(index_SOZ,component_position))];
            NIZ_components_PC2 = [NIZ_components_PC2, mean(weights_vector(index_NIZ,component_position))];
            PZ_components_PC2 = [PZ_components_PC2, mean(weights_vector(index_PZ,component_position))];
        end
    end 
    
    SOZ_components = [SOZ_components]'; 
    PZ_components = [PZ_components]'; 
    NIZ_components = [NIZ_components]'; 

    
    SOZ_components_PC2 = [SOZ_components_PC2]'; 
    PZ_components_PC2 = [PZ_components_PC2]';
    NIZ_components_PC2 = [NIZ_components_PC2]'; 
    
    
    % ANOVA with multcompare
    [p_vals, ~, STATS] = anova1([SOZ_components, PZ_components, NIZ_components],[],'off');
    [c,m,h,nms] = multcompare(STATS,'display','off');
    % 1 -> SOZ vs PZ
    % 2 -> SOZ vs NIZ
    % 3 -> PZ vs NIZ
    multcompare_pvals = c(:, end);
    
    % save the p-values so we can correct for multiple comparisons
    stats_table_all(counter, 1) = p_vals;
    stats_table_all(counter, 2:4) = multcompare_pvals';
    counter = counter + 1;
    
    soz_rgb_color = [0.6350 0.0780 0.1840];
    pz_rgb_color = [0.8500 0.3250 0.0980];
    niz_rgb_color = [0 0.4470 0.7410];
    
    f = figure('visible','on','Position',[0 0 600 800]);

    v = struct();
    v.SOZ = SOZ_components;
    v.PZ = PZ_components;
    v.NIZ = NIZ_components;
    violins = violinplot(v, {'SOZ', 'PZ', 'NIZ'}, 'ShowMean', true);

    violins(1,1).ViolinColor = soz_rgb_color;
    violins(1,2).ViolinColor =  pz_rgb_color;
    violins(1,3).ViolinColor =  niz_rgb_color;
    
    ylim([-10, 20])
    yticks([-10:5:20])
    ylabel("Average Weight")

    xlim([0, 4])
    xticklabels(["SOZ", "PZ", "NIZ"]);
    ax = gca;
    ax.Box = 'off';
    ax.LineWidth = 2;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.04];
    ax.FontSize = 20;
    
    title(sprintf("PC1 %s band", band_names(kk)));

    subtitle(sprintf("P=%d, SOZ-PZ=%d, SOZ-NIZ=%d, PZ-NIZ=%d", p_vals, multcompare_pvals), 'FontSize', 9);
    
    saveas(f, sprintf("%s\\%s_pc1_cohort_errorbar.png", out_dir, band_names(kk)));
    saveas(f, sprintf("%s\\%s_pc1_cohort_errorbar.svg", out_dir, band_names(kk)));
    saveas(f, sprintf("%s\\%s_pc1_cohort_errorbar.eps", out_dir, band_names(kk)));
    close(f);
    
    %% PC2+
    
    % ANOVA with multcompare
    [p_vals, ~, STATS] = anova1([SOZ_components_PC2, PZ_components_PC2, NIZ_components_PC2],[],'off');
    [c,m,h,nms] = multcompare(STATS,'display','off');
    % 1 -> SOZ vs PZ
    % 2 -> SOZ vs NIZ
    % 3 -> PZ vs NIZ
    multcompare_pvals = c(:, end);

    % save the p-values so we can correct for multiple comparisons
    stats_table_all(counter, 1) = p_vals;
    stats_table_all(counter, 2:4) = multcompare_pvals';
    counter = counter + 1;
    
    % colors for the plot
    soz_rgb_color = [0.6350 0.0780 0.1840];
    pz_rgb_color = [0.8500 0.3250 0.0980];
    niz_rgb_color = [0 0.4470 0.7410];
    
    line_width_points = 3;
    
    f = figure('visible','on','Position',[0 0 600 800]);

    v = struct();
    v.SOZ = SOZ_components_PC2;
    v.PZ = PZ_components_PC2;
    v.NIZ = NIZ_components_PC2;
    violins = violinplot(v, {'SOZ', 'PZ', 'NIZ'}, 'ShowMean', true);

    violins(1,1).ViolinColor = soz_rgb_color;
    violins(1,2).ViolinColor =  pz_rgb_color;
    violins(1,3).ViolinColor =  niz_rgb_color;
    
    ylim([-10, 20])
    yticks([-10:5:20])
    ylabel("Average Weight")
    
    xlim([0, 4])
    xticklabels(["SOZ", "PZ", "NIZ"]);
    ax = gca;
    ax.Box = 'off';
    ax.LineWidth = 2;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.04];
    ax.FontSize = 20;
    
    title(sprintf("PC2+ %s band", band_names(kk)));
    subtitle(sprintf("P=%d, SOZ-PZ=%d, SOZ-NIZ=%d, PZ-NIZ=%d", p_vals, multcompare_pvals), 'FontSize', 9);
    
    saveas(f, sprintf("%s\\%s_pc2_cohort_errorbar.png", out_dir, band_names(kk)));
    saveas(f, sprintf("%s\\%s_pc2_cohort_errorbar.svg", out_dir, band_names(kk)));
    saveas(f, sprintf("%s\\%s_pc2_cohort_errorbar.eps", out_dir, band_names(kk)));
    close(f);

end

% save the stats table
xlsx_table.("ANOVA_PVAL") = stats_table_all(:, 1);
xlsx_table.("SOZ-PZ_PVAL") = stats_table_all(:, 2);
xlsx_table.("SOZ-NIZ_PVAL") = stats_table_all(:, 3);
xlsx_table.("PZ-NIZ_PVAL") = stats_table_all(:, 4);

writetable(xlsx_table, xlsx_fpath);
