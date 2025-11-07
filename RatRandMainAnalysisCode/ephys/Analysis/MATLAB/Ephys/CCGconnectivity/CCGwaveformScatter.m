close all; clear all; clc; 

twoTrough = [45  19;
             79  38;
             3   45;
             18  45;
             32  79;
             45  4 ];
         
% pairsAna = twoTrough;

oneTrough = [3   44;
             20  45;
             119 122;
             44  79;
             4   79;
             20  34;
             20  44;
             3   38;
             91  101;
             45  12;
             82  91; 
             3   25;
             70  91];
                  
pairsAna = [twoTrough; oneTrough];

CCGleftTroughs  = [98  105 97  98  95  100 86 98  104 102 86 103 104 108 101 106 95  99  95];
CCGpeaks        = [106 114 102 106 104 106 90 102 106 104 90 106 106 111 105 111 99  106 99];
CCGrightTroughs = [112 119 109 114 108 115 94 106 108 106 94 109 108 114 109 116 103 113 103];
% 
% RefLeftTrough  = [19 18 20 20 22 19 23 19 22 21 23 19 20 22 23 18 23 19 21];
% RefRightTrough = [34 35 33 35 36 34 34 36 35 35 32 37 37 34 36 34 37 33 35];
% TarLeftTrough  = [23 22 19 19 20 23 21 20 19 20 20 21 21 23 21 23 23 19 23];
% TarRightTrough = [37 35 34 34 35 34 35 34 41 35 35 35 36 35 33 36 36 34 37];
% 
% CCGpeak2trough = mean([CCGpeaks-CCGleftTroughs;CCGrightTroughs-CCGpeaks]*(1/30));
% REFpeak2trough = mean([27-RefLeftTrough;RefRightTrough-27]*(1/30));
% TARpeak2trough = mean([27-TarLeftTrough;TarRightTrough-27]*(1/30));
% 
% pair1 = [98  106 112];
% pair2 = [105 114 119];
% pair3 = [97  102 109];
% pair4 = [98  106 114];
% pair5 = [95  104 108];
% pair6 = [100 106 115];
% 
% REF1 = [19 27 34];
% REF2 = [18 27 35];
% REF3 = [20 27 33];
% REF4 = [20 27 35];
% REF5 = [22 27 36];
% REF6 = [19 27 34];
% 
% TAR1 = [];
% TAR2 = [];
% TAR3 = [];
% TAR4 = [];
% TAR5 = [];
% TAR6 = [];
% 
% scatter([REFpeak2trough TARpeak2trough],[CCGpeak2trough CCGpeak2trough])

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

CCGparams = struct; 
CCGparams.fs          = fs;
CCGparams.binSize     = binSize;
CCGparams.duration    = duration;
CCGparams.jscale      = jscale;
CCGparams.plot_output = true; 
CCGparams.alpha       = alpha;
CCGparams.njitter     = njitter;
CCGparams.for_grant   = for_grant;

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

% pair data 
datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

% get all pairs
cell_pairs = {jit_var_out.cell_pair}';

% convert to matrix
cell_pairs_mat = [];
for i = 1:size(cell_pairs,1)
    cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
end

%% get session maze time markers
spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
[data_dir, name, ~] = fileparts(spike_data_fullpath);
load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');

%% get waveforms
extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
RoySession1 = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);

tWave = (-26:27)*(1/30);

%%

waveDataSpCounts = [];
for i = 1:size(RoySession1.times,2)
    waveDataSpCounts(i) = size(RoySession1.times{i},1);
end

for j = 3:size(pairsAna,1)
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = pairsAna(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;

    n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
    n2shank = jit_var_out(current_pair_indices(2)).cell2shank;

    numSpCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                 length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                 length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

    numSpCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                 length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                 length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

    waveIdxCell1 = find(numSpCell1 == waveDataSpCounts);
    waveIdxCell2 = find(numSpCell2 == waveDataSpCounts);

    waveDataCell1 = RoySession1.allSpikeRawWaveforms{waveIdxCell1};
    waveDataCell2 = RoySession1.allSpikeRawWaveforms{waveIdxCell2};

    % filter waveDataSpCounts to maze sessions 
    waveDataCell1 = waveDataCell1(length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + 1: ...
                                  length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                                  length(jit_var_out(current_pair_indices(2)).cell1_spike_times),:,:);

    waveDataCell2 = waveDataCell2(length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + 1: ...
                                  length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                                  length(jit_var_out(current_pair_indices(2)).cell2_spike_times),:,:);
                 
    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    
    tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact');
    nexttile([1 2])
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant);
           
    ylims = get(gca,'ylim');
                
    if any(GSPExc)
        hold on;
        plot(tR(GSPExc == 1)*1000, 0.99*ylims(2),'k^');
    end
    if any(GSPInh)
        hold on;
        plot(tR(GSPInh == 1)*1000, 0.99*ylims(2),'kv');
    end
    ylim(ylims)
%     xlim([-1 1])
    
    title(sprintf('%s - %d%s (sh %d) v %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                n1shank, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                n2shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
%     set(gca,'FontSize',16)
    
     % plot waveforms 
    [SyncSp,SyncCCG,SyncSpBinAll,SyncSpIdxBin] = SyncSpikes(res1, res2, ones(211,1));
    
    SyncSpIdxCell1 = [];
    SyncSpIdxCell2 = [];
    
    % all sync spike idx
    for k = 1:211
        temp = SyncSpIdxBin{k};
        SyncSpIdxCell1 = [SyncSpIdxCell1; temp(:,1)];
        SyncSpIdxCell2 = [SyncSpIdxCell2; temp(:,2)];
    end
    
    temp = SyncSpIdxBin{CCGpeaks(j)};
    SyncSpIdxCell1Bin1 = temp(:,1);
    SyncSpIdxCell2Bin1 = temp(:,2);
    
     % remove redundant entries
    SyncSpIdxCell1 = unique(SyncSpIdxCell1);
    SyncSpIdxCell2 = unique(SyncSpIdxCell2);
    
    % remove peak spikes from all spikes
    SyncSpIdxCell1 = setdiff(SyncSpIdxCell1,[SyncSpIdxCell1Bin1]);
    SyncSpIdxCell2 = setdiff(SyncSpIdxCell2,[SyncSpIdxCell2Bin1]);
     
    waveDataCell1Bin1   = double(waveDataCell1(SyncSpIdxCell1Bin1,:));
    waveDataCell2Bin1   = double(waveDataCell2(SyncSpIdxCell2Bin1,:));
    
    waveDataCell1NoSync = double(waveDataCell1(SyncSpIdxCell1,:));
    waveDataCell2NoSync = double(waveDataCell2(SyncSpIdxCell2,:));
    
    tWaveLat1 = tWave + 1000*tR(CCGpeaks(j));
    
    yyaxis right
    hold on
    
    meanRefSyncWaveform   = mean(waveDataCell1Bin1,1)/1000;
    meanTarSyncWaveform   = mean(waveDataCell2Bin1,1)/1000;
    
    meanRefNoSyncWaveform = mean(waveDataCell1NoSync,1)/1000;
    meanTarNoSyncWaveform = mean(waveDataCell2NoSync,1)/1000;
    
    stanDevRefSyncWaveform   = std(waveDataCell1Bin1,0,1)/1000;
    stanDevTarSyncWaveform   = std(waveDataCell2Bin1,0,1)/1000;
    
    stanDevRefNoSyncWaveform = std(waveDataCell1NoSync,0,1)/1000;
    stanDevTarNoSyncWaveform = std(waveDataCell2NoSync,0,1)/1000;
    
    plot(tWaveLat1,meanRefSyncWaveform,      'color',[0      0.4470 0.7410],'LineWidth',2)
    hold on
    plot(tWaveLat1,meanTarSyncWaveform, '-', 'color',[0.8500 0.3250 0.0980],'LineWidth',2)
    hold off
    
    
    l = legend('point-wise acc. bands', ...
           'simultaneous acc. bands', ...
           'jitter surr. mean', ...
           'raw CCG', ...
           'Location','NorthEast');
%     l.FontSize = 8;
    
%     sprintf('ref. %d%s sync wave', ...
%            jit_var_out(current_pair_indices(2)).cell_pair(1),...
%            jit_var_out(current_pair_indices(2)).cell1type), ...
%            sprintf('ref. %d%s sync wave', ...
%            jit_var_out(current_pair_indices(2)).cell_pair(2),...
%            jit_var_out(current_pair_indices(2)).cell1type), ...
%            'Location','NorthEast'
       
    xlabel('[ms]')
    ylabel('[mV]')

%     set(gca,'FontSize',16)
    set(gca, 'YDir','reverse')
    
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    stanDevRefSyncWaveform   = std(waveDataCell1Bin1,0,1)/1000;
    stanDevTarSyncWaveform   = std(waveDataCell2Bin1,0,1)/1000;

    stanDevRefNoSyncWaveform = std(waveDataCell1NoSync,0,1)/1000;
    stanDevTarNoSyncWaveform = std(waveDataCell2NoSync,0,1)/1000;

    %%
    
    nexttile
    
    Async   = meanRefSyncWaveform   + stanDevRefSyncWaveform;
    Bsync   = meanRefSyncWaveform   - stanDevRefSyncWaveform;

    AnoSync = meanRefNoSyncWaveform + stanDevRefNoSyncWaveform;
    BnoSync = meanRefNoSyncWaveform - stanDevRefNoSyncWaveform;

    plot(tWave,meanRefSyncWaveform,  '-',  'color',[0      0.4470 0.7410],'LineWidth',2)
    hold on
    plot(tWave,meanRefNoSyncWaveform,'-', 'color','k','LineWidth',2)
    patch([tWave'; flipud(tWave')],[Async';   flipud(Bsync')],   [0      0.4470 0.7410],'EdgeAlpha',0, 'FaceAlpha', 0.1); 
    patch([tWave'; flipud(tWave')],[AnoSync'; flipud(BnoSync')], [125 125 125]/255,'EdgeAlpha',0, 'FaceAlpha', 0.1); 
    hold off
    
    xlim([tWave(1) tWave(end)])
    xlabel('[ms]')
    ylabel('[mV]')
    title(sprintf('reference - %d%s eSpike wave', ...
        jit_var_out(current_pair_indices(2)).cell_pair(1),...
        jit_var_out(current_pair_indices(2)).cell1type))
    l = legend('sync.','no-sync','Location','SouthEast');
%     set(gca,'FontSize',16)
%     l.FontSize = 6;

    %%
    
    nexttile
    
    Async   = meanTarSyncWaveform   + stanDevTarSyncWaveform;
    Bsync   = meanTarSyncWaveform   - stanDevTarSyncWaveform;

    AnoSync = meanTarNoSyncWaveform + stanDevTarNoSyncWaveform;
    BnoSync = meanTarNoSyncWaveform - stanDevTarNoSyncWaveform;

    plot(tWave,meanTarSyncWaveform,  '-',  'color',[0.8500 0.3250 0.0980],'LineWidth',2)
    hold on
    plot(tWave,meanTarNoSyncWaveform,'-', 'color','k','LineWidth',2)
    patch([tWave'; flipud(tWave')],[Async';   flipud(Bsync')],   [0.8500 0.3250 0.0980],'EdgeAlpha',0, 'FaceAlpha', 0.1); 
    patch([tWave'; flipud(tWave')],[AnoSync'; flipud(BnoSync')], [125 125 125]/255,'EdgeAlpha',0, 'FaceAlpha', 0.1); 
    hold off
    
    xlim([tWave(1) tWave(end)])
    xlabel('[ms]')
    ylabel('[mV]')
    title(sprintf('target - %d%s eSpike wave', ...
        jit_var_out(current_pair_indices(2)).cell_pair(2),...
        jit_var_out(current_pair_indices(2)).cell1type))
    l = legend('sync.','no-sync','Location','SouthEast');
%     set(gca,'FontSize',16)
%     l.FontSize = 6;
    
    %%
    
    save_file = fullfile(data_dir, ['waveformAna_all' session_name '_jscale' num2str(jscale)]);
    set(fig_use,'PaperOrientation','landscape');
    print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
            
%     printNK(['waveformAna_' session_name '_jscale' num2str(jscale)],...
%                     'SuFiSyn', 'hfig', fig_use, 'append', true);
            
    % reset plot
    
    close 102

    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

close 102