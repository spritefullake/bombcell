%% Pipeline 
%% 0. Preliminary setup
% matlab path
addpath(genpath([pwd, filesep, 'helpers']))
%% 1. define parameters for quality metrics, ephys properties and
%% classification 
%define parameters
param = struct;
param.dontdo = 0; %re-calulate metrics and ephysparams if they already exist
% for calulating qMetrics
param.plotThis = 0; %plot metrics/params for each unit
param.dist = 0; %calculate distance metrics or not (this takes >3 timeslonger with)
param.driftdo = 1; %calculate slow drift, and metrics for chunks of time with most spikes present
param.chunkBychunk = 0; %calulate metrics for each chunk
param.tauR = 0.0010; %refractory period time (s)
param.tauC = 0.0002; %censored period time (s)
param.nChannelsIsoDist = 4; %like tetrodes
param.chanDisMax = 300; %maximum distance
param.raw = 1; %calculate metrics also for raw data
param.strOnly = 0; %only use str_templates
% for calulating eParams
param.ACGbinSize = 0.001; %bin size to calc. ACG
param.ACGduration = 1; %ACG full duration
param.maxFRbin = 10; %
param.histBins = 1000;
% to choose good units
param.minNumSpikes = 300;
param.minIsoDist = 0;
param.minLratio = 0;
param.minSScore = 0.01;
param.minSpatDeKlowbound = 1.5;

param.maxNumPeak = 3;
param.minAmpli = 77;
param.maxRPV = 2;
param.somaCluster = 1;
param.plotMetricsCtypes = 0;
% for burst merging - WIP, not implemented yet
param.maxPercMissing = 30;
param.maxChanDistance = 40;
param.waveformMinSim = 0.8;
param.spikeMaxLab = 0.15;
param.minPeakRatio = 0.7;
param.maxdt = 10;
% for cell-type classification
param.cellTypeDuration = 400;
param.cellTypePostS = 40;

% %% load experiment 
% % you can use any other script to load your data
ephys_path = uigetdir('H:\Neuropixels Data')

%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_16-20-06_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_16-20-06_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'

%'H:\Neuropixels Data\20201030-WT\20201030-WT-baseline+PTZ\2020-10-30_12-27-33_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'

%'H:\Neuropixels Data\20200915-WT\20200915-WT-baseline2+PTZ\2020-09-15_15-23-58_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_16-20-06_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20200915-WT\20200915-WT-baseline2+PTZ\2020-09-15_15-23-58_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_16-20-06_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_14-02-37_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20201023-WT\20201023-WT-baseline+PTZ\2020-10-23_10-44-16_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_14-02-37_extracted\continuous\Neuropix-PXI-slot4-probe2-AP' 
%'H:\Neuropixels Data\20200915-WT\20200915-WT-baseline2+PTZ\2020-09-15_15-23-58_extracted\continuous\Neuropix-PXI-slot4-probe2-AP';
ephysData = loadEphysDataJF(ephys_path);
%% run quality metrics 
getQualityMetrics;

%% run ephys properties
getEphysProperties;

keep qMetric ephysParams ephysData param ephys_path

%% classify cells 
classifyStriatum; 

%% show cell types percentages
figure();
pie([sum(msn), sum(tan), sum(fsi), sum(uin)]);
lgd = legend({'msn','tan','fsi','uin'});
lgd.Location = 'south';
title("Cell Type Percentages");
%save qualityMetrics, ephysProperties and classification
%% Plot histograms of spike/s counts per cell type
figure();
cellTypeNames = fields(cellTypes);
for iCellType=1:length(cellTypeNames)
    subplot(1,length(cellTypeNames),iCellType)
    cellTypeIndices = cellTypes.(cellTypeNames{iCellType});
    spike_data = ephysParams.spike_rateAP(cellTypeIndices);
    [~,edges] = histcounts(log10(spike_data));
    %new_edges = linspace(min(edges),max(edges),20);
    histogram(spike_data,10.^edges);
    set(gca, 'xscale','log');
    title(cellTypeNames{iCellType});
    ylabel('Count')
    xlabel('Spikes/s')
end
%% JF's plot checks; seems to work only if nothing is 'unsorted'
%very quick plotting-just to check 
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];
waveform_t = 1e3*((0:size(ephysData.template_waveforms,2)-1)/30000);
figure();
hold on;
for iCellType=1:4
    plot(waveform_t, nanmean(ephysData.template_waveforms(cellTypesClassif(iCellType).cells,:)),'Color',celltype_col(iCellType,:))
    plotshaded(waveform_t, [ nanmean(ephysData.template_waveforms(cellTypesClassif(iCellType).cells,:)) - ...
        nanstd(ephysData.template_waveforms(cellTypesClassif(iCellType).cells,:)); ...
        nanmean(ephysData.template_waveforms(cellTypesClassif(iCellType).cells,:)) + ...
        nanstd(ephysData.template_waveforms(cellTypesClassif(iCellType).cells,:))],celltype_col(iCellType,:))
end
makepretty;

figure();
hold on;
for iCellType=1:4
    subplot(1,4,iCellType)
    area(0:0.001:1, nanmean(ephysParams.ACG(cellTypesClassif(iCellType).cells,:)),'FaceColor',celltype_col(iCellType,:))
    plotshaded(0:0.001:1,[ nanmean(ephysParams.ACG(cellTypesClassif(iCellType).cells,:)) - ...
        nanstd(ephysParams.ACG(cellTypesClassif(iCellType).cells,:)); ...
        nanmean(ephysParams.ACG(cellTypesClassif(iCellType).cells,:)) + ...
        nanstd(ephysParams.ACG(cellTypesClassif(iCellType).cells,:))],celltype_col(iCellType,:))
end