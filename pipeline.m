%% Pipeline 
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
% % you can use any other script to load your data, you need to end up with:
% % xxxxx
% animals={'AP024'};
% curr_animal = 1; % (set which animal to use)
% animal = animals{curr_animal};
% protocol = 'vanillaChoiceworld'; % (this is the name of the Signals protocol)
% experiments = AP_find_experimentsJF(animal, protocol, true);
% experiments = experiments([experiments.imaging] & [experiments.ephys]);
% curr_day = 1; % (set which day to use)
% day = experiments(curr_day).day; % date
% thisDay = experiments(curr_day).day; % date
% thisDate = thisDay;
% experiment = experiments(curr_day).experiment; % experiment number
% load_parts.cam=false;
% load_parts.imaging=true;
% load_parts.ephys=true;

%loading
% [ephys_path, ephys_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys',[],[]);
%ephys_path = strcat(experiments(curr_day).location, '\ephys\kilosort2\');
% corona = 0;
% ephysData = loadEphysDataJF(ephys_path, animal, day, experiment); %load and format ephysData to use later 
ephys_path = 'H:\Neuropixels Data\20201023-WT\20201023-WT-baseline+PTZ\2020-10-23_10-44-16_extracted\continuous\Neuropix-PXI-slot4-probe2-AP'
%'H:\Neuropixels Data\20210216 striatum-938 HOM-20200607 CamKII-Gi-PrL-20210121 headplate\2021-02-16_14-02-37_extracted\continuous\Neuropix-PXI-slot4-probe2-AP' 
%'H:\Neuropixels Data\20200915-WT\20200915-WT-baseline2+PTZ\2020-09-15_15-23-58_extracted\continuous\Neuropix-PXI-slot4-probe2-AP';
ephysData = loadEphysDataJF(ephys_path);
%% run quality metrics 
getQualityMetrics;

%% run ephys properties
getEphysProperties;

keep qMetric ephysParams ephysData param  

%% classify cells 
classifyStriatum; 

%plot cells 

%save qualityMetrics, ephysProperties and classification