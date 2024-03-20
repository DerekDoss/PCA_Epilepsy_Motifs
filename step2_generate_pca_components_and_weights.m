% runs a PCA algorithm on the PDC data and saves it to a struct
% Importing all 81 pats 
% Code authored by Derek Doss and Jared Shless


clear all
close all
restoredefaultpath

% save path for the final struct
struct_save_fpath = "Z:\PCA_project\data\brain_revision\with_delta\pat_pca_struct.mat";

file_keyword = 'FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz';

% Loads data for all 81 patients; called pats.mat. 
load("Z:\PCA_project\data\brain_revision\with_delta\pat_structs_with_dkregion_and_soz.mat")

pdc_mats = {};

% classic frequency bands from 'hdrs'
band_names = ["delta","theta", "alpha", "beta", "gamma_low", "gamma_high"];

% ordered by [patient, band, row, column]
% outputs a struct array with pdc matrices for each band
% each row is a pt., each col is a double of data in the freq. field
tmp = vertcat(pats.long);

% all of the thetas for us to index into, and...:
pdc_mats(:,1) = {tmp(:).delta};
pdc_mats(:,2) = {tmp(:).theta};
% all of the alphas...
pdc_mats(:,3) = {tmp(:).alpha};
% all of the betas...
pdc_mats(:,4) = {tmp(:).beta};
% all of the low gammas...
pdc_mats(:,5) = {tmp(:).gamma_low};
% ...and, all of the high gammas
pdc_mats(:,6) = {tmp(:).gamma_high};


input_matrix = {};
weight_matrix = {}; % pre-allocate for weights
explained_matrix = {}; % pre-allocate explained matrix 

for B = 1:length(band_names) % for each band; keep alpha 
   
    for p = 1:size(pats,2) % for all 81 patients in  cohort 

        % FC matrix
        FC_matrix = pats(p).long.(band_names{B});

        % set the diagonal of the FC matrix to zero for z-scoring
        FC_matrix(logical(eye(size(FC_matrix)))) = nan;

        % FC outward connectivity
        mat_out = FC_matrix;
        % FC inward connectivity
        mat_in = FC_matrix';

       % variables x observations aka 2*nodes x nodes   
       prelim_matrix = cat(1, mat_out, mat_in);
       
       % diagonal in combined matrix space
       concatenated_diagonal = cat(1, eye(size(mat_out)), eye(size(mat_in)) );

       % variables x observations
       input_matrix = (prelim_matrix - repmat(nanmean(prelim_matrix,2),[1, size(prelim_matrix,2)])) ./ repmat(nanstd(prelim_matrix,[],2),[1, size(prelim_matrix,2)]);

       % set the diagonal of the input matrix to zero
       input_matrix(logical(concatenated_diagonal)) = 0;

       
       % flip it such that it is observations x variables
       input_matrix = input_matrix';

       pats(p).input_matrix.(band_names{B}) = input_matrix; 

       [coeff,score,~,~,explained,mu] = pca(input_matrix);

       pats(p).coeff_matrix.(band_names{B}) = coeff;
       pats(p).weights.(band_names{B}) = score;
       pats(p).explained.(band_names{B}) = explained;

    end 
end

save(struct_save_fpath, 'pats');
 