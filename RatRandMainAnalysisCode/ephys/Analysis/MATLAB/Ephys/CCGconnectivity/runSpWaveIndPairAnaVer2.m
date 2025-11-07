close all; clear all; clc; 

%% hardcoding pairs of interest for waveform analysis

pairAna = [20  45];
binGroup1 = 88:93;
binGroup2 = 101:107;
       
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

plotIdx = [1,3,5,7,9,11,13,15];


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
% extraDrivePath = '/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/'; % backup computer
extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/'; % regular computer
RoySession1 = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);

tWave = (-26:27)*(1/30);

%%

waveDataSpCounts = [];
for i = 1:size(RoySession1.times,2)
    waveDataSpCounts(i) = size(RoySession1.times{i},1);
end

pair_idx = find((pairAna(1) == cell_pairs_mat(:,1)) & ...
                    (pairAna(2) == cell_pairs_mat(:,2)));
    
n1maze  = jit_var_out(pair_idx(2)).cell1_spike_times;
n2maze  = jit_var_out(pair_idx(2)).cell2_spike_times;

n1shank = jit_var_out(pair_idx(2)).cell1shank;
n2shank = jit_var_out(pair_idx(2)).cell2shank;
    
numSpCell1 = length(jit_var_out(pair_idx(1)).cell1_spike_times) + ...
             length(jit_var_out(pair_idx(2)).cell1_spike_times) + ...
             length(jit_var_out(pair_idx(3)).cell1_spike_times);

numSpCell2 = length(jit_var_out(pair_idx(1)).cell2_spike_times) + ...
             length(jit_var_out(pair_idx(2)).cell2_spike_times) + ...
             length(jit_var_out(pair_idx(3)).cell2_spike_times);
                
waveIdxCell1 = find(numSpCell1 == waveDataSpCounts);
waveIdxCell2 = find(numSpCell2 == waveDataSpCounts);
    
waveDataCell1 = RoySession1.allSpikeRawWaveformsAllChans{waveIdxCell1};
waveDataCell2 = RoySession1.allSpikeRawWaveformsAllChans{waveIdxCell2};
    
% filter waveDataSpCounts to maze sessions 
waveDataCell1 = waveDataCell1(length(jit_var_out(pair_idx(1)).cell1_spike_times) + 1: ...
                              length(jit_var_out(pair_idx(1)).cell1_spike_times) + ...
                              length(jit_var_out(pair_idx(2)).cell1_spike_times),:,:);

waveDataCell2 = waveDataCell2(length(jit_var_out(pair_idx(1)).cell2_spike_times) + 1: ...
                              length(jit_var_out(pair_idx(1)).cell2_spike_times) + ...
                              length(jit_var_out(pair_idx(2)).cell2_spike_times),:,:);
    
% plot CCG 
res1 = jit_var_out(pair_idx(2)).cell1_spike_times;
res2 = jit_var_out(pair_idx(2)).cell2_spike_times;

for i = 1:8

    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                'plot_output', get(fig_use, 'Number'), 'subfig', 1, 'plot_flag', false, ...
                'subplot_size', [3, 1], 'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);

%     ylims = get(gca,'ylim');
% 
%     hold on
%     for k = 1:length(binGroup1)
%         line([tR(binGroup1(k))*1000,tR(binGroup1(k))*1000],ylims,'Color',[0 0.4470 0.7410],'LineWidth',2)
%     end
%     for k = 1:length(binGroup2)
%         line([tR(binGroup2(k))*1000,tR(binGroup2(k))*1000],ylims,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
%     end
%     hold off
% 
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1000, 0.95*ylims(2),'k^');
%     end
% 
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'kv');
%     end
% 
%     title(sprintf('%s - %d%s (sh %d) v %d%s (sh %d) - %s', ...
%                 jit_var_out(pair_idx(2)).epoch, ...
%                 jit_var_out(pair_idx(2)).cell_pair(1),jit_var_out(pair_idx(2)).cell1type, ...
%                 n1shank, ...
%                 jit_var_out(pair_idx(2)).cell_pair(2),jit_var_out(pair_idx(2)).cell2type, ...
%                 n2shank, ...
%                 conn_type{conn_type_idx(pair_idx(2))}))

    % plot waveforms 
    [SyncSp,SyncCCG,SyncSpBinAll,SyncSpIdxBin] = SyncSpikes(res1, res2, ones(211,1));
    
    SyncSpIdxCell1Bin1 = [];
    SyncSpIdxCell2Bin1 = [];
    SyncSpIdxCell1Bin2 = [];
    SyncSpIdxCell2Bin2 = [];
    
    BINS = binGroup1;
    for k = 1:length(BINS)
        temp = SyncSpIdxBin{BINS(k)};
        SyncSpIdxCell1Bin1 = [SyncSpIdxCell1Bin1; temp(:,1)];
        SyncSpIdxCell2Bin1 = [SyncSpIdxCell2Bin1; temp(:,2)];
    end
    
    BINS = binGroup2;
    for k = 1:length(BINS)
        temp = SyncSpIdxBin{BINS(k)};
        SyncSpIdxCell1Bin2 = [SyncSpIdxCell1Bin2; temp(:,1)];
        SyncSpIdxCell2Bin2 = [SyncSpIdxCell2Bin2; temp(:,2)];
    end

    waveDataCell1Bin1 = squeeze(double(waveDataCell1(SyncSpIdxCell1Bin1,i,:)));
    waveDataCell1Bin2 = squeeze(double(waveDataCell1(SyncSpIdxCell1Bin2,i,:)));
    
    waveDataCell2Bin1 = squeeze(double(waveDataCell2(SyncSpIdxCell2Bin1,i,:)));
    waveDataCell2Bin2 = squeeze(double(waveDataCell2(SyncSpIdxCell2Bin2,i,:)));
    
    subplot(8,2,plotIdx(i))
    hold on
    if size(waveDataCell1Bin1,1) <= numWaves
        for k = 1:size(waveDataCell1Bin1,1)
            patchline(tWave,waveDataCell1Bin1(k,:)/1000,'EdgeColor',[0 0.4470 0.7410],'LineWidth',1,'EdgeAlpha',0.25);
        end
    elseif size(waveDataCell1Bin1,1) > numWaves
        temp = randperm(size(waveDataCell1Bin1,1));
        temp = temp(1:numWaves);
        for k = 1:numWaves
            patchline(tWave,waveDataCell1Bin1(temp(k),:)/1000,'EdgeColor',[0 0.4470 0.7410],'LineWidth',1,'EdgeAlpha',0.25);
        end
    end
    
    if size(waveDataCell1Bin2,1) <= numWaves
        for k = 1:size(waveDataCell1Bin2,1)
            patchline(tWave,waveDataCell1Bin2(k,:)/1000,'EdgeColor',[0.8500 0.3250 0.0980],'LineWidth',1,'EdgeAlpha',0.25);
        end
    elseif size(waveDataCell1Bin2,1) > numWaves
        temp = randperm(size(waveDataCell1Bin2,1));
        temp = temp(1:numWaves);
        for k = 1:numWaves
            patchline(tWave,waveDataCell1Bin2(temp(k),:)/1000,'EdgeColor',[0.8500 0.3250 0.0980],'LineWidth',1,'EdgeAlpha',0.25);
        end
    end
    
    plot(tWave,mean(waveDataCell1Bin1,1)/1000,'color',[0 0.4470 0.7410],'LineWidth',2)
    plot(tWave,mean(waveDataCell1Bin2,1)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
    
%     patch([tWave'; fliplr(tWave)'], ...
%           [mean(waveDataCell1Bin1,1)'/1000 - std(waveDataCell1Bin1,0,1)'/1000;...
%           flipud(mean(waveDataCell1Bin1,1)'/1000 + std(waveDataCell1Bin1,0,1)'/1000)], ...
%           [0 0.4470 0.7410],'FaceAlpha',0.5,'EdgeAlpha',0);
%     patch([tWave'; fliplr(tWave)'], ...
%           [mean(waveDataCell1Bin2,1)'/1000 - std(waveDataCell1Bin2,0,1)'/1000;...
%           flipud(mean(waveDataCell1Bin2,1)'/1000 + std(waveDataCell1Bin2,0,1)'/1000)], ...
%           [0.8500 0.3250 0.0980],'FaceAlpha',0.5,'EdgeAlpha',0);    
    hold off
%     xlabel('[ms]')
    ylabel('[mV]')
    xlim([tWave(1),tWave(end)])
    ylim([-3 2])
    title(sprintf('%d%s channel %d (Max chan: 2)', ...
                jit_var_out(pair_idx(2)).cell_pair(1), ...
                jit_var_out(pair_idx(2)).cell1type,i))
    
    subplot(8,2,plotIdx(i)+1)
    
    hold on
    if size(waveDataCell2Bin1,1) <= numWaves
        for k = 1:size(waveDataCell2Bin1,1)
            patchline(tWave,waveDataCell2Bin1(k,:)/1000,'EdgeColor',[0 0.4470 0.7410],'LineWidth',1,'EdgeAlpha',0.25);
        end
    elseif size(waveDataCell2Bin1,1) > numWaves
        temp = randperm(size(waveDataCell2Bin1,1));
        temp = temp(1:numWaves);
        for k = 1:numWaves
            patchline(tWave,waveDataCell2Bin1(temp(k),:)/1000,'EdgeColor',[0 0.4470 0.7410],'LineWidth',1,'EdgeAlpha',0.25);
        end
    end
        
    if size(waveDataCell2Bin2,1) <= numWaves
        for k = 1:size(waveDataCell2Bin2,1)
            patchline(tWave,waveDataCell2Bin2(k,:)/1000,'EdgeColor',[0.8500 0.3250 0.0980],'LineWidth',1,'EdgeAlpha',0.25);
        end
    elseif size(waveDataCell2Bin2,1) > numWaves
        temp = randperm(size(waveDataCell2Bin2,1));
        temp = temp(1:numWaves);
        for k = 1:numWaves
            patchline(tWave,waveDataCell2Bin2(temp(k),:)/1000,'EdgeColor',[0.8500 0.3250 0.0980],'LineWidth',1,'EdgeAlpha',0.25);
        end
    end
        
    plot(tWave,mean(waveDataCell2Bin1,1)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2)
    plot(tWave,mean(waveDataCell2Bin2,1)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
    
%     patch([tWave'; fliplr(tWave)'], ...
%           [mean(waveDataCell2Bin1,1)'/1000 - std(waveDataCell2Bin1,0,1)'/1000;...
%           flipud(mean(waveDataCell2Bin1,1)'/1000 + std(waveDataCell2Bin1,0,1)'/1000)], ...
%           [0 0.4470 0.7410],'FaceAlpha',0.5,'EdgeAlpha',0);
    
%     patch([tWave'; fliplr(tWave)'], ...
%           [mean(waveDataCell2Bin2,1)'/1000 - std(waveDataCell2Bin2,0,1)'/1000;...
%           flipud(mean(waveDataCell2Bin2,1)'/1000 + std(waveDataCell2Bin2,0,1)'/1000)], ...
%           [0.8500 0.3250 0.0980],'FaceAlpha',0.5,'EdgeAlpha',0);
    hold off
%     xlabel('[ms]')
    ylabel('[mV]')
    xlim([tWave(1),tWave(end)])
    ylim([-5 2])
    title(sprintf('%d%s channel %d (Max chan: 3)', ...
                jit_var_out(pair_idx(2)).cell_pair(2), ...
                jit_var_out(pair_idx(2)).cell2type,i))
        
%     printNK(['waveformAna_' session_name '_jscale' num2str(jscale)],...
%                     'SuFiSyn', 'hfig', fig_use, 'append', true);
            
    % reset plot
%     close 102
% 
%     hcomb = figure(102);
%     arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

end

save_file = fullfile(data_dir, ['waveformAna_20v45_allChans_Ver2' session_name '_jscale' num2str(jscale)]);
set(fig_use,'PaperOrientation','landscape');
print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');

close all
                         
                         