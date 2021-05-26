function ephysData = loadEphysDataJF(ephys_path)
load_parts.ephys=true; 
%AP_load_experimentJF;

spike_templates_0idx = readNPY([ephys_path filesep 'spike_templates.npy']);
spike_times = readNPY([ephys_path filesep  'spike_times.npy']);

pc_features = readNPY([ephys_path filesep  'pc_features.npy']);
pc_features_ind = readNPY([ephys_path filesep  'pc_feature_ind.npy']) + 1;

spike_times_full = spike_times;
spike_templates_full = spike_templates_0idx + 1;

%% obtain struct fields
channel_map = readNPY([ephys_path, filesep, 'channel_map.npy']);
channel_positions = readNPY([ephys_path, filesep, 'channel_positions.npy']);
%cluster_groups =  
ephys_sample_rate = 30000
%spike_depths = 
winv = readNPY([ephys_path, filesep, 'whitening_mat_inv.npy']);

cluster_filename = [ephys_path, filesep, 'cluster_group', '.tsv']
fid = fopen(cluster_filename)
cluster_groups = textscan(fid, '%d%s', 'HeaderLines', 1);
disp(cluster_groups)
fclose(fid)


templates_whitened = readNPY([ephys_path, filesep, 'templates.npy']);
channel_positions = readNPY([ephys_path, filesep, 'channel_positions.npy']);
channel_map = readNPY([ephys_path, filesep, 'channel_map.npy']);

template_amplitudes = readNPY([ephys_path, filesep, 'amplitudes.npy']);

% Default channel map/positions are from end: make from surface
% (hardcode this: kilosort2 drops channels)
max_depth = 3840;
channel_positions(:, 2) = max_depth - channel_positions(:, 2);

% Unwhiten templates
templates = zeros(size(templates_whitened));
for t = 1:size(templates_whitened, 1)
    templates(t, :, :) = squeeze(templates_whitened(t, :, :)) * winv;
end

% Get the waveform of all templates (channel with largest amplitude)
[~, max_site] = max(max(abs(templates), [], 2), [], 3);
templates_max = nan(size(templates, 1), size(templates, 2));
for curr_template = 1:size(templates, 1)
    templates_max(curr_template, :) = ...
        templates(curr_template, :, max_site(curr_template));
end
waveforms = templates_max;
% Get depth of each template
% (get min-max range for each channel)
template_chan_amp = squeeze(range(templates, 2));
% (zero-out low amplitude channels)
template_chan_amp_thresh = max(template_chan_amp, [], 2) * 0.5;
template_chan_amp_overthresh = template_chan_amp .* (template_chan_amp >= template_chan_amp_thresh);
% (get center-of-mass on thresholded channel amplitudes)
template_depths = sum(template_chan_amp_overthresh.*channel_positions(:, 2)', 2) ./ sum(template_chan_amp_overthresh, 2);

% Get the depth of each spike (templates are zero-indexed)
spike_depths = template_depths(spike_templates_0idx+1);

% Get trough-to-peak time for each template
templates_max_signfix = bsxfun(@times, templates_max, ...
    sign(abs(min(templates_max, [], 2))-abs(max(templates_max, [], 2))));

[~, waveform_trough] = min(templates_max, [], 2);
[~, waveform_peak_rel] = arrayfun(@(x) ...
    max(templates_max(x, waveform_trough(x):end), [], 2), ...
    transpose(1:size(templates_max, 1)));
waveform_peak = waveform_peak_rel + waveform_trough;

templateDuration = waveform_peak - waveform_trough;
templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

%% addendum
spike_times_timeline = spike_times_full
%% getting the 'good templates'

% Check that all used spike templates have a label
spike_templates_0idx_unique = unique(spike_templates_0idx);
if ~all(ismember(spike_templates_0idx_unique, uint32(cluster_groups{1}))) || ...
        ~all(ismember(cluster_groups{2}, {'good', 'mua', 'noise'}))
    warning([animal, ' ', day, ': not all templates labeled']);
end

% Define good units from labels
good_templates_idx = uint32(cluster_groups{1}( ...
    strcmp(cluster_groups{2}, 'good') | strcmp(cluster_groups{2}, 'mua')));
good_templates = ismember(0:size(templates, 1)-1, good_templates_idx);
% Throw out all non-good template data
templates = templates(good_templates, :, :);
template_depths = template_depths(good_templates);
waveforms = waveforms(good_templates, :);
templateDuration = templateDuration(good_templates);
templateDuration_us = templateDuration_us(good_templates);

% Throw out all non-good spike data
good_spike_idx = ismember(spike_templates_0idx, good_templates_idx);
spike_times = spike_times(good_spike_idx);
spike_times_full = spike_times_timeline;
spike_templates_full = spike_templates_0idx+1;
spike_templates_0idx = spike_templates_0idx(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_depths = spike_depths(good_spike_idx);
spike_times_timeline = spike_times_timeline(good_spike_idx);

% Rename the spike templates according to the remaining templates
% (and make 1-indexed from 0-indexed)
new_spike_idx = nan(max(spike_templates_0idx)+1, 1);
new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
spike_templates = new_spike_idx(spike_templates_0idx+1);
%% save in a structure
ephysData = struct;
% in progress: the lone '%%' after a line means I've re-implemented 
% obtaining that field in the 'obtain struct fields' section or it has
% already been obtained previously in this file
ephysData.channel_map = channel_map;%%
ephysData.channel_positions = channel_positions;%%
ephysData.cluster_groups = cluster_groups;%%
ephysData.ephys_sample_rate = ephys_sample_rate;%%
ephysData.spike_depths = spike_depths; %only good ones %%
ephysData.spike_templates = spike_templates; %only good ones %%
ephysData.spike_times_full = spike_times_full;%%
ephysData.spike_templates_full = spike_templates_full;%%
ephysData.spike_times_timeline = spike_times_timeline %% % spike_times_timeline; %only good ones
ephysData.template_amplitudes = template_amplitudes; %% %only good ones
ephysData.template_depths = template_depths;%%
ephysData.templates = templates; %%
ephysData.template_waveforms = waveforms; %%
ephysData.winv = winv;%%
ephysData.good_templates = good_templates;
ephysData.pc_features = pc_features(good_spike_idx,:,:);
ephysData.pc_features_full = pc_features;%%
ephysData.pc_features_ind_full = pc_features_ind;%%
ephysData.pc_features_ind = pc_features_ind(good_templates_idx+1,:);
%ephysData.str_templates = str_templates;
%ephysData.spike_rate = spike_rate;
ephysData.new_spike_idx = new_spike_idx;


end