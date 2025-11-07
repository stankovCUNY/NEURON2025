%% case 45i chans 1,2,3

load('dataSpatialCCGs_neuron_45i_RoyMaze1.mat')

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

peak = [106 110];
tWave = (-26:27)*(1/30);

CCGstart = 76;
CCGend   = 136;

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

% channel 1
res1 = dataNetPlot(6).cell1spikeTimes;
res2 = dataNetPlot(6).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan1_18p = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan1_18p,SyncCCG_chan1_18p,SyncSpBinAll_chan1_18p,SyncSpIdxBin_chan1_18p] = SyncSpikes(res1, res2, ones(211,1));

res1 = dataNetPlot(8).cell1spikeTimes;
res2 = dataNetPlot(8).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan1_20i = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan1_20i,SyncCCG_chan1_20i,SyncSpBinAll_chan1_20i,SyncSpIdxBin_chan1_20i] = SyncSpikes(res1, res2, ones(211,1));

% channel 2
res1 = dataNetPlot(1).cell1spikeTimes;
res2 = dataNetPlot(1).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan2_3i = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan2_3i,SyncCCG_chan2_3i,SyncSpBinAll_chan2_3i,SyncSpIdxBin_chan2_3i] = SyncSpikes(res1, res2, ones(211,1));

res1 = dataNetPlot(3).cell1spikeTimes;
res2 = dataNetPlot(3).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan2_12p = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan2_12p,SyncCCG_chan2_12p,SyncSpBinAll_chan2_12p,SyncSpIdxBin_chan2_12p] = SyncSpikes(res1, res2, ones(211,1));

res1 = dataNetPlot(4).cell1spikeTimes;
res2 = dataNetPlot(4).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan2_15p = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan2_15p,SyncCCG_chan2_15p,SyncSpBinAll_chan2_15p,SyncSpIdxBin_chan2_15p] = SyncSpikes(res1, res2, ones(211,1));

res1 = dataNetPlot(5).cell1spikeTimes;
res2 = dataNetPlot(5).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan2_17p = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan2_17p,SyncCCG_chan2_17p,SyncSpBinAll_chan2_17p,SyncSpIdxBin_chan2_17p] = SyncSpikes(res1, res2, ones(211,1));

% channel 3
res1 = dataNetPlot(2).cell1spikeTimes;
res2 = dataNetPlot(2).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan3_4i = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan3_4i,SyncCCG_chan3_4i,SyncSpBinAll_chan3_4i,SyncSpIdxBin_chan3_4i] = SyncSpikes(res1, res2, ones(211,1));

res1 = dataNetPlot(7).cell1spikeTimes;
res2 = dataNetPlot(7).cell2spikeTimes;
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), ...
                'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
CCG_chan3_19p = ccgR(CCGstart:CCGend,1,2)/sum(ccgR(CCGstart:CCGend,1,2));
[SyncSp_chan3_19p,SyncCCG_chan3_19p,SyncSpBinAll_chan3_19p,SyncSpIdxBin_chan3_19p] = SyncSpikes(res1, res2, ones(211,1));


tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
nexttile

plot(tR(CCGstart:CCGend)*1000,CCG_chan1_18p,'LineWidth',2)
hold on
plot(tR(CCGstart:CCGend)*1000,CCG_chan1_20i,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan2_3i ,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan2_12p,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan2_15p,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan2_17p,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan3_4i ,'LineWidth',2)
plot(tR(CCGstart:CCGend)*1000,CCG_chan3_19p,'LineWidth',2)
% xline(tR(peak)*1000,'LineWidth',2)
hold off
xlabel('[ms]')
ylabel('prob.')
title('CCG')

legend({[num2str(dataNetPlot(6).cell2ID) dataNetPlot(6).cell2type ...
        ' (shank: ' num2str(dataNetPlot(6).cell2shank) ', ch: ' num2str(dataNetPlot(6).chCell2) ')'], ...
        [num2str(dataNetPlot(8).cell2ID) dataNetPlot(8).cell2type, ...  
        '  (shank: ' num2str(dataNetPlot(8).cell2shank) ', ch: ' num2str(dataNetPlot(8).chCell2) ')'], ...
        [num2str(dataNetPlot(1).cell2ID) dataNetPlot(1).cell2type, ...
        '    (shank: ' num2str(dataNetPlot(1).cell2shank) ', ch: ' num2str(dataNetPlot(1).chCell2) ')'], ...
        [num2str(dataNetPlot(3).cell2ID) dataNetPlot(3).cell2type ...
        ' (shank: ' num2str(dataNetPlot(3).cell2shank) ', ch: ' num2str(dataNetPlot(3).chCell2) ')'], ...
        [num2str(dataNetPlot(4).cell2ID) dataNetPlot(4).cell2type, ...  
        ' (shank: ' num2str(dataNetPlot(4).cell2shank) ', ch: ' num2str(dataNetPlot(4).chCell2) ')'], ...
        [num2str(dataNetPlot(5).cell2ID) dataNetPlot(5).cell2type, ...
        ' (shank: ' num2str(dataNetPlot(5).cell2shank) ', ch: ' num2str(dataNetPlot(5).chCell2) ')'], ...
        [num2str(dataNetPlot(2).cell2ID) dataNetPlot(2).cell2type, ...  
        '    (shank: ' num2str(dataNetPlot(2).cell2shank) ', ch: ' num2str(dataNetPlot(2).chCell2) ')'], ...
        [num2str(dataNetPlot(7).cell2ID) dataNetPlot(7).cell2type, ...
        ' (shank: ' num2str(dataNetPlot(7).cell2shank) ', ch: ' num2str(dataNetPlot(7).chCell2) ')']}, ...
        'Location','northwest')

waveIdx_chan1_18p = SyncSpIdxBin_chan1_18p{peak(1)};
waveIdx_chan1_20i = SyncSpIdxBin_chan1_20i{peak(1)};
waveIdx_chan2_3i  = SyncSpIdxBin_chan2_3i{peak(1)};
waveIdx_chan2_12p = SyncSpIdxBin_chan2_12p{peak(1)};
waveIdx_chan2_15p = SyncSpIdxBin_chan2_15p{peak(1)};
waveIdx_chan2_17p = SyncSpIdxBin_chan2_17p{peak(1)};
waveIdx_chan3_4i  = SyncSpIdxBin_chan3_4i{peak(1)};
waveIdx_chan3_19p = SyncSpIdxBin_chan3_19p{peak(1)};
    
wave_ref_chan1_18p
wave_ref_chan1_20i
wave_ref_chan2_3i
wave_ref_chan2_12p
wave_ref_chan2_15p
wave_ref_chan2_17p
wave_ref_chan3_4i
wave_ref_chan3_19p

wave_tar_chan1_18p
wave_tar_chan1_20i
wave_tar_chan2_3i
wave_tar_chan2_12p
wave_tar_chan2_15p
wave_tar_chan2_17p
wave_tar_chan3_4i
wave_tar_chan3_19p

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