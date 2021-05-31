function ephysData = loadEphysDataJF(ephys_path)
load_parts.ephys=true; 
isSpikeGlx=0;
loadClusters = 0;
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

% %% taken from ap_load
% 
% if ephys_exists && load_parts.ephys
%     if verbose;
%         disp('Estimating striatum boundaries on probe...');
%     end
% 
%     % str_align = alignment method ('none', 'depth', or 'kernel')
% 
%     % requires n_aligned_depths for alignment, set default
%     if ~exist('n_aligned_depths', 'var')
%         n_aligned_depths = 3;
%     end
% 
%     % if no alignment specified, default kernel
%     if ~exist('str_align', 'var')
%         str_align = 'kernel';
%     end
%     try
%     [str_depth, aligned_str_depth_group] = AP_align_striatum_ephysJF;
%     catch
%     end
% 
% end
% %% Classify spikes
% 
% if ephys_exists && load_parts.ephys && exist('str_depth')
%     if verbose;
%         disp('Classifying spikes...');
%     end
% 
%     str_templates = template_depths >= str_depth(1) & template_depths <= str_depth(2);
%     non_str_templates = ~str_templates;
% 
%     % Define the window to look for spiking statistics in (spikes go in and
%     % out, so take the bin with the largest firing rate for each cell and work
%     % with that one)
%     % spiking_stat_window = 60*5; % seconds
%     % spiking_stat_bins = min(spike_times_timeline):spiking_stat_window: ...
%     %     max(spike_times_timeline);
% 
%     % % (for whole session)
%     spiking_stat_window = max(spike_times_timeline) - min(spike_times_timeline);
%     spiking_stat_bins = [min(spike_times_timeline), max(spike_times_timeline)];
% 
%     % Get firing rate across the session
%     bin_spikes = nan(size(templates, 1), ...
%         length(spiking_stat_bins)-1);
%     for curr_template = unique(spike_templates)'
%         bin_spikes(curr_template, :) = ...
%             histcounts(spike_times_timeline(spike_templates == curr_template), ...
%             spiking_stat_bins);
%     end
%     min_spikes = 10;
%     use_spiking_stat_bins = bsxfun(@ge, bin_spikes, prctile(bin_spikes, 80, 2)) & bin_spikes > min_spikes;
%     spike_rate = sum(bin_spikes.*use_spiking_stat_bins, 2) ./ ...
%         (sum(use_spiking_stat_bins, 2) * spiking_stat_window);
% 
%     % Get proportion of ISI > 2s (Yamin/Cohen 2013) and CV2 (Stalnaker/Schoenbaum 2016)
%     prop_long_isi = nan(size(templates, 1), 1);
%     cv2 = nan(size(templates, 1), 1);
%     for curr_template = unique(spike_templates)'
% 
%         long_isi_total = 0;
%         isi_ratios = [];
%         for curr_bin = find(use_spiking_stat_bins(curr_template, :))
%             curr_spike_times = spike_times_timeline( ...
%                 spike_times_timeline > spiking_stat_bins(curr_bin) & ...
%                 spike_times_timeline < spiking_stat_bins(curr_bin+1) & ...
%                 spike_templates == curr_template);
%             curr_isi = diff(curr_spike_times);
% 
%             long_isi_total = long_isi_total + sum(curr_isi(curr_isi > 2));
% 
%             isi_ratios = [isi_ratios; (2 * abs(curr_isi(2:end)-curr_isi(1:end-1))) ./ ...
%                 (curr_isi(2:end) + curr_isi(1:end-1))];
%         end
% 
%         prop_long_isi(curr_template) = long_isi_total / ...
%             (sum(use_spiking_stat_bins(curr_template, :)) * spiking_stat_window);
%         cv2(curr_template) = nanmean(isi_ratios);
% 
%     end
% 
%     % Cortical classification (like Bartho JNeurophys 2004)
%     waveform_duration_cutoff = 400;
%     narrow = non_str_templates & templateDuration_us <= waveform_duration_cutoff;
%     wide = non_str_templates & templateDuration_us > waveform_duration_cutoff;
% 
%     % Striatum classification
%     prop_long_isi_cutoff = 0.35;
%     cv2_cutoff = 0.8;
% 
%     msn = str_templates & ...
%         templateDuration_us > waveform_duration_cutoff & ...
%         prop_long_isi >= prop_long_isi_cutoff;
% 
%     fsi = str_templates & ...
%         templateDuration_us <= waveform_duration_cutoff & ...
%         prop_long_isi < prop_long_isi_cutoff;
% 
%     tan = str_templates & ...
%         templateDuration_us > waveform_duration_cutoff & ...
%         prop_long_isi < prop_long_isi_cutoff;
% 
%     uin = str_templates & ~msn & ~fsi & ~tan;
% 
%     waveform_t = 1e3 * ((0:size(templates, 2) - 1) / ephys_sample_rate);
% 
%     if verbose
% 
%         % Plot the waveforms and spike statistics
%         figure;
% 
%         if any(non_str_templates)
%             subplot(2, 2, 1);
%             hold on;
%             p = plot(waveform_t, waveforms(non_str_templates, :)');
%             set(p(wide(non_str_templates)), 'color', 'k')
%             set(p(narrow(non_str_templates)), 'color', 'r')
%             xlabel('Time (ms)')
%             title('Not striatum');
%             legend([p(find(wide(non_str_templates), 1)), p(find(narrow(non_str_templates), 1))], {'Wide', 'Narrow'})
%         end
% 
%         subplot(2, 2, 2);
%         hold on;
%         p = plot(waveform_t, waveforms(str_templates, :)');
%         set(p(msn(str_templates)), 'color', 'm')
%         set(p(fsi(str_templates)), 'color', 'b')
%         set(p(tan(str_templates)), 'color', 'g')
%         set(p(uin(str_templates)), 'color', 'c')
%         xlabel('Time (ms)')
%         title('Striatum');
%         legend([p(find(msn(str_templates), 1)), p(find(fsi(str_templates), 1)), ...
%             p(find(tan(str_templates), 1)), p(find(uin(str_templates), 1))], {'MSN', 'FSI', 'TAN', 'UIN'});
% 
%         subplot(2, 2, 3);
%         hold on;
% 
%         stem3( ...
%             templateDuration_us(wide)/1000, ...
%             prop_long_isi(wide), ...
%             spike_rate(wide), 'k');
% 
%         stem3( ...
%             templateDuration_us(narrow)/1000, ...
%             prop_long_isi(narrow), ...
%             spike_rate(narrow), 'r');
% 
%         xlabel('waveform duration (ms)')
%         ylabel('frac long ISI')
%         zlabel('spike rate')
% 
%         set(gca, 'YDir', 'reverse')
%         set(gca, 'XDir', 'reverse')
%         view(3);
%         grid on;
%         axis vis3d;
% 
%         subplot(2, 2, 4);
%         hold on;
%         stem3( ...
%             templateDuration_us(msn)/1000, ...
%             prop_long_isi(msn), ...
%             spike_rate(msn), 'm');
% 
%         stem3( ...
%             templateDuration_us(fsi)/1000, ...
%             prop_long_isi(fsi), ...
%             spike_rate(fsi), 'b');
% 
%         stem3( ...
%             templateDuration_us(tan)/1000, ...
%             prop_long_isi(tan), ...
%             spike_rate(tan), 'g');
% 
%         stem3( ...
%             templateDuration_us(uin)/1000, ...
%             prop_long_isi(uin), ...
%             spike_rate(uin), 'c');
% 
%         xlabel('waveform duration (ms)')
%         ylabel('frac long ISI')
%         zlabel('spike rate')
% 
%         set(gca, 'YDir', 'reverse')
%         set(gca, 'XDir', 'reverse')
%         view(3);
%         grid on;
%         axis vis3d;
% 
%         % Plot depth vs. firing rate colored by cell type
%         celltype_labels = {'Wide', 'Narrow', 'MSN', 'FSI', 'TAN', 'UIN'};
%         celltypes = wide .* 1 + narrow .* 2 + msn .* 3 + fsi .* 4 + tan .* 5 + uin .* 6;
%         use_colors = ...
%             [0, 0, 0; ...
%             1, 0, 0; ...
%             1, 0, 1; ...
%             0, 0, 1; ...
%             0, 1, 0; ...
%             0, 1, 1];
% 
%         plot_celltypes = any([wide, narrow, msn, fsi, tan, uin], 1);
% 
%         norm_spike_n = mat2gray(log10(accumarray(spike_templates, 1)+1));
% 
%         figure('Position', [94, 122, 230, 820]);
%         gscatter(norm_spike_n, template_depths, celltypes, use_colors, [], 10);
%         xlim([0, 1])
%         set(gca, 'YDir', 'reverse');
%         xlabel('Norm log_{10} spike rate');
%         ylabel('Depth (\mum)');
%         legend(celltype_labels(plot_celltypes), 'location', 'NW');
%         ylim([0, max(channel_positions(:, 2))])
% 
%         drawnow;
% 
%     end
% 
% end
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
ephysData.templateDuration = templateDuration

end