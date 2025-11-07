%% case 79i chans 1,2,3

load('dataSpatialCCGs_neuron_79i_RoyMaze1.mat')

% constants
session_name   = 'RoyMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1;
alpha_name     = 5;
duration       = 0.007;
fs             = 30000;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
numWaves       = 100;
resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

peak = 122;
tWave = (-26:27)*(1/30);

CCGstart = 76;
CCGend   = 136;

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

% channel 1
res1 = dataNetPlot(3).cell1spikeTimes;
res2 = dataNetPlot(3).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan1 = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan1,SyncCCG_chan1,SyncSpBinAll_chan1,SyncSpIdxBin_chan1] = SyncSpikes(res1, res2, ones(211,1));

% channel 2
res1 = dataNetPlot(1).cell1spikeTimes;
res2 = dataNetPlot(1).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan2 = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan2,SyncCCG_chan2,SyncSpBinAll_chan2,SyncSpIdxBin_chan2] = SyncSpikes(res1, res2, ones(211,1));

% channel 3
res1 = dataNetPlot(2).cell1spikeTimes;
res2 = dataNetPlot(2).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan3 = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan3,SyncCCG_chan3,SyncSpBinAll_chan3,SyncSpIdxBin_chan3] = SyncSpikes(res1, res2, ones(211,1));

tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
nexttile

plot(tR(CCGstart:CCGend)*1000,CCG_chan1,'LineWidth',2)
hold on
plot(tR(CCGstart:CCGend)*1000,CCG_chan2,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan3,'LineWidth',2)
xline(tR(peak)*1000,'LineWidth',2)
hold off
xlabel('[ms]')
ylabel('prob.')
title('CCG')

legend({[num2str(dataNetPlot(3).cell2ID) dataNetPlot(3).cell2type ...
        ' (shank: ' num2str(dataNetPlot(3).cell2shank) ', ch: ' num2str(dataNetPlot(3).chCell2) ')'], ...
        [num2str(dataNetPlot(1).cell2ID) dataNetPlot(1).cell2type, ...  
        '   (shank: ' num2str(dataNetPlot(1).cell2shank) ', ch: ' num2str(dataNetPlot(1).chCell2) ')'], ...
        [num2str(dataNetPlot(2).cell2ID) dataNetPlot(2).cell2type, ...
        '   (shank: ' num2str(dataNetPlot(2).cell2shank) ', ch: ' num2str(dataNetPlot(2).chCell2) ')']}, ...
        'Location','northwest')
    
waveIdx_chan1 = SyncSpIdxBin_chan1{peak};
waveIdx_chan3 = SyncSpIdxBin_chan3{peak};

wave_ref_chan1 = dataNetPlot(3).waveDataCell1(waveIdx_chan1(:,1),:);
wave_tar_chan1 = dataNetPlot(3).waveDataCell2(waveIdx_chan1(:,2),:);

wave_ref_chan2 = dataNetPlot(1).waveDataCell1;
wave_tar_chan2 = dataNetPlot(1).waveDataCell2;

wave_ref_chan3 = dataNetPlot(2).waveDataCell1(waveIdx_chan3(:,1),:);
wave_tar_chan3 = dataNetPlot(2).waveDataCell2(waveIdx_chan3(:,2),:);

wave_ref_mean_chan1 = mean(wave_ref_chan1)/1000;
wave_ref_mean_chan2 = mean(wave_ref_chan2)/1000;
wave_ref_mean_chan3 = mean(wave_ref_chan3)/1000;

wave_tar_mean_chan1 = mean(wave_tar_chan1)/1000;
wave_tar_mean_chan2 = mean(wave_tar_chan2)/1000;
wave_tar_mean_chan3 = mean(wave_tar_chan3)/1000;

nexttile
plot(tWave,wave_ref_mean_chan1,'LineWidth',2)
hold on
plot(tWave,wave_ref_mean_chan2,'LineWidth',2)
plot(tWave,wave_ref_mean_chan3,'LineWidth',2)
hold off
ylabel('[mV]')
xlabel('[ms]')
title('reference')
set(gca, 'YDir','reverse')

legend({[num2str(dataNetPlot(3).cell2ID) dataNetPlot(3).cell2type '' ...
        ' (shank: ' num2str(dataNetPlot(3).cell1shank) ', ch: ' num2str(dataNetPlot(3).chCell1) ')'], ...
        [num2str(dataNetPlot(1).cell2ID) dataNetPlot(2).cell2type, ...  
        '   (shank: ' num2str(dataNetPlot(1).cell1shank) ', ch: ' num2str(dataNetPlot(1).chCell1) ')'], ...
        [num2str(dataNetPlot(2).cell2ID) dataNetPlot(1).cell2type, ...
        '   (shank: ' num2str(dataNetPlot(2).cell1shank) ', ch: ' num2str(dataNetPlot(2).chCell1) ')']}, ...
        'Location','northwest')
    

nexttile
plot(tWave,wave_tar_mean_chan1,'LineWidth',2)
hold on
plot(tWave,wave_tar_mean_chan2,'LineWidth',2)
plot(tWave,wave_tar_mean_chan3,'LineWidth',2)
hold off
ylabel('[mV]')
xlabel('[ms]')
title('target')
set(gca, 'YDir','reverse')

legend({[num2str(dataNetPlot(3).cell2ID) dataNetPlot(3).cell2type ...
        ' (shank: ' num2str(dataNetPlot(3).cell2shank) ', ch: ' num2str(dataNetPlot(3).chCell2) ')'], ...
        [num2str(dataNetPlot(1).cell2ID) dataNetPlot(2).cell2type, ...  
        '   (shank: ' num2str(dataNetPlot(1).cell2shank) ', ch: ' num2str(dataNetPlot(1).chCell2) ')'], ...
        [num2str(dataNetPlot(2).cell2ID) dataNetPlot(1).cell2type, ...
        '   (shank: ' num2str(dataNetPlot(2).cell2shank) ', ch: ' num2str(dataNetPlot(2).chCell2) ')']}, ...
        'Location','northwest')
    



