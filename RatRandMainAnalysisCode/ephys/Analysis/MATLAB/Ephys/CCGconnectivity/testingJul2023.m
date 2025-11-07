%% 45 v 20 

tarID = 20;
refID = 34;

tarShank = 2;

spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';

[data_dir, name, ~] = fileparts(spike_data_fullpath);

load(spike_data_fullpath, 'spikes')
load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
load(fullfile(data_dir, 'wake-basics.mat'),   'basics');

temp = spikes.RoyMaze1(refID).time;
temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
refSpikeTimes = temp/1e6; % convert to ms

temp = spikes.RoyMaze1(tarID).time;
temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
tarSpikeTimes = temp/1e6; % convert to ms

onsetTime = behavior.RoyMaze1.time(2,1);
% spikeTimeIndex = round((spikeTimes-onsetTime/1e6)*fs) + 1;

%% load union clusters

tarShNoiseST      = load([dataPath 'RoyMaze_shank' num2str(tarShank) '_spikeTimesNoise']);
tarShValidUnitsST = load([dataPath 'RoyMaze_shank' num2str(tarShank) '_spikeTimesValidUnits']);
% tarShValidUnitsST = setdiff(tarShValidUnitsST.validUnits,refSpikeTimes);
tarShValidUnitsST = setdiff(tarShValidUnitsST.validUnits,tarSpikeTimes);
tarShankAllST     = [tarShNoiseST.noiseCluster; tarShValidUnitsST];

%% CCG screen pairwise

figure(102)

jscale         = 1;
alpha_name     = 5;
duration       = 0.01;
fpass          = 300;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
filterFlag     = false;
fs        = 30000;
binSize   = 1/fs;

[GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(refSpikeTimes,tarSpikeTimes,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant,  'plot_pointwiseBands',false);
    
[SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(refSpikeTimes', tarSpikeTimes', ones(length(GSPExc),1),duration);

refSpikeTimesSync = [];
tarSpikeTimesSync = [];

for k = 1:length(GSPExc)

    if (GSPExc(k) == 1) 
        refSpikeTimesSync = [refSpikeTimesSync; SyncSpBinAll{k}(:,1)];
        tarSpikeTimesSync = [tarSpikeTimesSync; SyncSpBinAll{k}(:,2)];
    end

end

%% CCG screen union

figure(102)

jscale         = 1;
alpha_name     = 5;
duration       = 0.01;
fpass          = 300;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
filterFlag     = false;
fs        = 30000;
binSize   = 1/fs;

[GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(refSpikeTimes,tarShankAllST,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant,  'plot_pointwiseBands',false);
    
[SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(refSpikeTimes', tarShankAllST, ones(length(GSPExc),1),duration);

refSpikeTimesSync = [];
tarSpikeTimesSync = [];

for k = 1:length(GSPExc)

    if (GSPExc(k) == 1) 
        refSpikeTimesSync = [refSpikeTimesSync; SyncSpBinAll{k}(:,1)];
        tarSpikeTimesSync = [tarSpikeTimesSync; SyncSpBinAll{k}(:,2)];
    end

end

%%

tarCh                       = 4;
tarUnitSpikeCountAllSess    = length(spikes.RoyMaze1(tarID).time);
plotParam                   = 'ind';
synchSpikesFlag             = true; 
refUnitSpikeTimes           = refSpikeTimesSync;

curateClustersForPlots(tarShank, ...
                       tarCh, ...
                       tarUnitSpikeCountAllSess, ...
                       plotParam, ...
                       synchSpikesFlag, ...
                       refUnitSpikeTimes)
                   
                   
%%

tiledlayout(1,3)

nexttile
scatter(PC1noSynch{21}, PC2noSynch{21})
hold on
scatter(PC1synch{2}, PC2synch{2})
hold off
xlabel('PC1')
ylabel('PC2')
legend('noSynch','synch','location','southoutside')

nexttile
scatter(PC1noSynch{21}, PC3noSynch{21})
hold on
scatter(PC1synch{2}, PC3synch{2})
hold off
xlabel('PC1')
ylabel('PC3')
legend('noSynch','synch','location','southoutside')
title(["ch: " num2str(tarCh)])

nexttile
scatter(PC2noSynch{21}, PC3noSynch{21})
hold on
scatter(PC2synch{2}, PC3synch{2})
hold off
xlabel('PC2')
ylabel('PC3')
legend('noSynch','synch','location','southoutside')

%%

% refShNoiseST = load([dataPath 'RoyMaze_shank' num2str(1) '_spikeTimesNoise']);
% refShNoiseST = refShNoiseST.noiseCluster;
% 
% figure(102)
% 
% jscale         = 1;
% alpha_name     = 5;
% duration       = 0.01;
% fpass          = 300;
% fig_use        = 102;
% njitter        = 500;
% alpha          = 0.05;
% for_grant      = false;
% filterFlag     = false;
% fs        = 30000;
% binSize   = 1/fs;
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
%       CCG_jitter(refShNoiseST,tarShNoiseST,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_flag', true, ...
%                     'plot_output', get(fig_use, 'Number'), ...
%                     'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant,  'plot_pointwiseBands',false);

%%

% figure(102)
% 
% jscale         = 1;
% alpha_name     = 5;
% duration       = 0.01;
% fpass          = 300;
% fig_use        = 102;
% njitter        = 500;
% alpha          = 0.05;
% for_grant      = false;
% filterFlag     = false;
% fs        = 30000;
% binSize   = 1/fs;
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
%       CCG_jitter(refSpikeTimes,tarShankAllST,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_flag', true, ...
%                     'plot_output', get(fig_use, 'Number'), ...
%                     'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant,  'plot_pointwiseBands',false);
