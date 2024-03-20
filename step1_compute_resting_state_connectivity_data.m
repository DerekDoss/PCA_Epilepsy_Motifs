
%% Connectivity calculation for 2-minute epochs (epochs already generated)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original code by Graham Johnson
% Modified by Derek Doss

% Hard coded to do Delta [1-3 Hz], Theta[4-8 Hz], Alpha[8-12 Hz], Beta[13-30 Hz], Low
% Gamma[31-80 Hz], High Gamma[81-150 Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated on 2/20/2024 by Derek Doss to include delta


close all
clear all
restoredefaultpath

source_dir = 'Z:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data';
sub_folder_path = 'Unfiltered_Chunks\First_Collection\All_2minChunks_Bipole_UnusedChannelsDeleted_UnFiltered';
edf_keyword = 'Raw_bipole_unusedChannelsDeleted_Unfiltered';
soz_label_file = 'Z:\000_Data\SEEG\SEEG_EyesClosed_RestingState\labels\all_pats_bipole_soz_labels.csv';
metric_loc_in_struct = 3;

%Path to violin plot folder
addpath('Z:\shared_toolboxes\Violinplot-Matlab-master')
% Path to EDFread
addpath('Z:\shared_toolboxes\edfRead')
% Path to fieldtrip toolbox
addpath('Z:\shared_toolboxes\fieldtrip-20190819');

metric = 'PDC';
mvar_order = 10; % Will need to do 1 for short epochs
long_epoch_seg_length = 5; % seconds in which to define new trials


pat_dir_list = ["Epat26_dpat05", "Epat27_dpat06", "Epat30",...
    "Epat31_apat111", "Epat34_dpat11", "Epat35_apat115", "Epat37_dpat13",...
    "Epat38", "Epat39", "Spat28", "Spat30",...
    "Spat31", "Spat34", "Spat36", "Spat37",...
    "Epat43_apat117_vpat04", "Spat41", "Spat42", "Spat44",...
    "Spat39", "Spat40", "Epat02", "Epat03",...
    "Epat04_dpat07", "Epat05", "Epat06_apat102", "Epat08_dpat08",...
    "Epat09", "Epat10", "Epat11_apat113", "Epat13",...
    "Epat14", "Epat15", "Epat17_apat112", "Epat18_dpat03",...
    "Epat19_dpat04", "Epat20", "Epat21_pat107",...
    "Epat22", "Epat23_spat21_dpat02", "Epat24",...
    "Epat25_apat108_dpat10", "Epat28_apat106", "pat11", "pat33_dpat09",...
    "Spat12", "Spat13", "Spat14", "Spat17",...
    "Spat18", "Spat19", "Spat20", "Spat22",...
    "Spat23", "Spat24", "Spat26", "Spat29",...
    "Spat32", "Spat33", "Spat48", "Spat49",...
    "Spat50", "Spat51", "Spat52", "Spat53",...
    "Spat47", "Epat33", "Spat02", "Spat03",...
    "Spat05", "Spat06", "Spat07", "Spat08",...
    "Spat09", "Spat10", "Spat11", "Spat25",...,
    "Spat27", "Spat45", "Spat46", "Epat41"];



% Initialize the output
% Save into .mat
pats = struct;


% Iterate through the subject directories
for i = 1:length(pat_dir_list)% 1 to the amount of selected patients
    
    curr_dir = fullfile(source_dir, pat_dir_list(i),sub_folder_path);
    
    pats(i).subID = pat_dir_list(i);
    
    % Get directory contents
    contents = dir(curr_dir);
    edf_idxs = find(contains({contents.name}, edf_keyword));
    if isempty(edf_idxs); error("%d EDFs found in directory (expected > 0)\n%s",length(edf_idxs), curr_dir); end
    fprintf("\n%d EDFs found in directory:\n",length(edf_idxs))
    edf_names = string({contents(edf_idxs).name}');
    full_edf_paths = fullfile(curr_dir,edf_names);
    
    data = cell(length(edf_idxs),1);
    hdrs = cell(length(edf_idxs),1);
    
    % Import all of the EDFs
    for j = 1:length(edf_idxs)
        disp("Importing EDF...")
        [hdrs{j}, data{j}] = edfread(full_edf_paths(j));
        disp("EDF import complete")
    end
    
    %% Obtain connectivity metric for all 2-minute epochs for this subject
    % Initialize the 2-minute epoch results
    all_D = nan(length(edf_idxs),length(hdrs{1}.label),length(hdrs{1}.label));
    all_T = all_D;
    all_A = all_D;
    all_B = all_D;
    all_G_low = all_D;
    all_G_high = all_D;
    
    for j = 1:length(edf_idxs)
        
        fprintf("%d/%d Sub ID: %s, epoch %d\n",i,length(pat_dir_list),pat_dir_list(i),j)
        
        if j ==1
            % Save a copy of the labels for this patient
            pats(i).labels = string(hdrs{j}.label)';
        end
        
        % Filter the data
        data_curr = data{j};
        filt_data = FT_filt_A(data_curr,hdrs{j}.frequency(1));
        
        % Calculate connectivity
        % bands: Delta [1-3 Hz], Theta[4-8 Hz], Alpha[8-12 Hz], Beta[13-30 Hz], Low Gamma[31-80 Hz], High Gamma[81-150 Hz]
        [D,T,A,B,G_low,G_high] = pdc_calc(filt_data,hdrs{j},mvar_order,long_epoch_seg_length);
               
        all_D(j,:,:) = D;
        all_T(j,:,:) = T;
        all_A(j,:,:) = A;
        all_B(j,:,:) = B;
        all_G_low(j,:,:) = G_low;
        all_G_high(j,:,:) = G_high;
        
    end
    
    % Average the PDC across 2-minute epochs
    avg_D = squeeze(nanmean(all_D,1));
    avg_T = squeeze(nanmean(all_T,1));
    avg_A = squeeze(nanmean(all_A,1));
    avg_B = squeeze(nanmean(all_B,1));
    avg_G_low = squeeze(nanmean(all_G_low,1));
    avg_G_high = squeeze(nanmean(all_G_high,1));
    
    % Assign to master patient variable
    pats(i).long.delta = avg_D;
    pats(i).long.theta = avg_T;
    pats(i).long.alpha = avg_A;
    pats(i).long.beta = avg_B;
    pats(i).long.gamma_low = avg_G_low;
    pats(i).long.gamma_high = avg_G_high;
    
    
end


% save to struct
save('Z:\000_Data\SEEG\SEEG_EyesClosed_RestingState\code\connectivity_calculations\PDC_pats_with_delta.mat', 'pats');


