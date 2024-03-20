% Compute ISH component position
% Authored by Derek Doss and Jared Shless

clear all
close all
clc

% load in the PCA struct
load("Z:\PCA_project\data\brain_revision\with_delta\pat_pca_struct.mat")

bands_to_analyze = ["delta","theta", "alpha", "beta", "gamma_low", "gamma_high"];

for ii = 1:length(pats)
    
    % concatenates the SOZ designations two times over
    node_type_indicies = [pats(ii).soz_designation; pats(ii).soz_designation];

    % inward FC = 0, outward = 1
    in_out_indicides = [zeros(length(pats(ii).soz_designation), 1); ones(length(pats(ii).soz_designation), 1)];

    soz_indicies = find(node_type_indicies == 1);

    % inward FC indicies
    soz_in_indices = find(in_out_indicides == 0 & node_type_indicies == 1);

    % outward FC indicies
    soz_out_indices = find(in_out_indicides == 1 & node_type_indicies == 1);

    for jj = 1:length(bands_to_analyze)
        
        % variables (in/out connectivity) x number of components
        tmp_components = pats(ii).coeff_matrix.(bands_to_analyze(jj));

        soz_in_coeffs = mean(tmp_components(soz_in_indices, 1:10), 1, 'omitnan');
        soz_out_coeffs = mean(tmp_components(soz_out_indices, 1:10), 1, 'omitnan');

        soz_in_minus_out = soz_in_coeffs - soz_out_coeffs;

        pats(ii).soz_in_minus_out_coeffs.(bands_to_analyze(jj)) = soz_in_minus_out;

        [~, pats(ii).ish_component_position.(bands_to_analyze(jj))] = max(soz_in_minus_out);

        pats(ii).ISH_degree_val.(bands_to_analyze(jj)) = max(soz_in_minus_out);
        pats(ii).first_degree_val.(bands_to_analyze(jj)) = soz_in_coeffs(1) - soz_out_coeffs(1);

    end

end

save("Z:\PCA_project\data\brain_revision\with_delta\pats_strict_5.mat", "pats");