close all; clear all; clc; 

%% hardcoding pairs of interest for waveform analysis

pairsAna = [3   44;             20  45;             20  79;             34  79;
            79  91;             79  118;            91  118;            44  79;
            91  109;            73  91;             119 91;             45  4;
            45  19;             79  38;             3   45;             4   79;
            18  45;             32  79;             79  38;             3   38;
            91  101;            45  12;             82  91;             3   25;
            70  91;             20  34;             20  44;             40  79;
            79  58;             3   32;             4   25;             20  38;
            46  79;             79  87;             87  119;            98  119;
            119 122;            45  15;             17  45;             79  19;
            91  68;             79  95;             91  110];

binGroup1 = {90;                88:93;              90;                 89;
             91;                106;                112;                81:93;
             119;               89;                 93:98;              105:106;
             105:107;           111:115;            101:104;            88:91;
             104:108;           101:105;            111:115;            111:114;
             105;               110:111;            98:100;             105:107;
             99:100;            106:107;            106;                102:104;
             126;               110;                105:106;            106;
             96:97;             111:112;            105:106;            105:106;
             106;               110;                102;                122;
             106;               112;                119};
       
binGroup2 = {106;               101:107;            104;                97:102;
             105:106;           117:119;            115:116;            103:104;
             133;               99:100;             106:108;            119:121;
             121;               125;                [];                 [];
             [];                [];                 [];                 [];
             [];                [];                 [];                 [];
             [];                [];                 [];                 [];
             [];                [];                 [];                 [];
             [];                [];                 [];                 [];
             [];                [];                 [];                 [];
             [];                 [];                []};
       
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

for j = 1:size(pairsAna,1)
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = pairsAna(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
   
    n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
    n2shank = jit_var_out(current_pair_indices(2)).cell2shank;
    
    numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                     length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                     length(jit_var_out(current_pair_indices(3)).cell1_spike_times);
                
    numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                     length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                     length(jit_var_out(current_pair_indices(3)).cell2_spike_times);
                
    waveIdxCell1 = find(numSpikesCell1 == waveDataSpCounts);
    waveIdxCell2 = find(numSpikesCell2 == waveDataSpCounts);
    
    waveDataCell1 = RoySession1.allSpikeRawWaveforms{waveIdxCell1};
    waveDataCell2 = RoySession1.allSpikeRawWaveforms{waveIdxCell2};
    
    % filter waveDataSpCounts to maze sessions 
    waveDataCell1 = waveDataCell1(length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + 1: ...
                                  length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                                  length(jit_var_out(current_pair_indices(2)).cell1_spike_times),:);
    
    waveDataCell2 = waveDataCell2(length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + 1: ...
                                  length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                                  length(jit_var_out(current_pair_indices(2)).cell2_spike_times),:);
    
    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    
%     tiledlayout(3,1, 'Padding', 'none', 'TileSpacing', 'compact');
%     nexttile
    
    
    yyaxis left
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant);
            
    ylims = get(gca,'ylim');

    hold on
    for k = 1:length(binGroup1{idx})
        line([tR(binGroup1{idx}(k))*1000,tR(binGroup1{idx}(k))*1000],ylims,'Color',[0 0.4470 0.7410],'LineWidth',1)
    end
    for k = 1:length(binGroup2{idx})
        line([tR(binGroup2{idx}(k))*1000,tR(binGroup2{idx}(k))*1000],ylims,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
    end
    hold off
    
    if any(GSPExc)
        hold on;
        plot(tR(GSPExc == 1)*1000, 0.95*ylims(2),'k^');
    end

    if any(GSPInh)
        hold on;
        plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'kv');
    end
    
    ylim(ylims)
    
    title(sprintf('%s - %d%s (sh %d) v %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                n1shank, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                n2shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
    set(gca,'FontSize',16)
    
    % plot waveforms 
    [SyncSp,SyncCCG,SyncSpBinAll,SyncSpIdxBin] = SyncSpikes(res1, res2, ones(211,1));
    
    SyncSpIdxCell1 = [];
    SyncSpIdxCell2 = [];
    SyncSpIdxCell1Bin1 = [];
    SyncSpIdxCell2Bin1 = [];
    SyncSpIdxCell1Bin2 = [];
    SyncSpIdxCell2Bin2 = [];
    
    % all sync spike idx
    for k = 1:211
        temp = SyncSpIdxBin{k};
        SyncSpIdxCell1 = [SyncSpIdxCell1; temp(:,1)];
        SyncSpIdxCell2 = [SyncSpIdxCell2; temp(:,2)];
    end
    
    BINS = binGroup1{idx};
    for k = 1:length(BINS)
        temp = SyncSpIdxBin{BINS(k)};
        SyncSpIdxCell1Bin1 = [SyncSpIdxCell1Bin1; temp(:,1)];
        SyncSpIdxCell2Bin1 = [SyncSpIdxCell2Bin1; temp(:,2)];
    end
    
    BINS = binGroup2{idx};
    for k = 1:length(BINS)
        temp = SyncSpIdxBin{BINS(k)};
        SyncSpIdxCell1Bin2 = [SyncSpIdxCell1Bin2; temp(:,1)];
        SyncSpIdxCell2Bin2 = [SyncSpIdxCell2Bin2; temp(:,2)];
    end
    
    % remove redundant entries
    SyncSpIdxCell1 = unique(SyncSpIdxCell1);
    SyncSpIdxCell2 = unique(SyncSpIdxCell2);
    
    % remove peak spikes from all spikes
    SyncSpIdxCell1 = setdiff(SyncSpIdxCell1,[SyncSpIdxCell1Bin1; SyncSpIdxCell1Bin2]);
    SyncSpIdxCell2 = setdiff(SyncSpIdxCell2,[SyncSpIdxCell2Bin1; SyncSpIdxCell2Bin2]);
    
    waveDataCell1NoSync = double(waveDataCell1(SyncSpIdxCell1,:));
    waveDataCell2NoSync = double(waveDataCell2(SyncSpIdxCell2,:));
    
    waveDataCell1Bin1   = double(waveDataCell1(SyncSpIdxCell1Bin1,:));
    waveDataCell1Bin2   = double(waveDataCell1(SyncSpIdxCell1Bin2,:));
    
    waveDataCell2Bin1   = double(waveDataCell2(SyncSpIdxCell2Bin1,:));
    waveDataCell2Bin2   = double(waveDataCell2(SyncSpIdxCell2Bin2,:));
    
    tWaveLat1 = tWave + 1000*mean(tR(binGroup1{j}));
    tWaveLat2 = tWave + 1000*mean(tR(binGroup2{j}));

    yyaxis right
    hold on
    
    plot(tWaveLat1,mean(waveDataCell1Bin1,1)/1000,'--', 'color',[0      0.4470 0.7410],'LineWidth',2)
    plot(tWaveLat1,mean(waveDataCell2Bin1,1)/1000,'-.', 'color',[0      0.4470 0.7410],'LineWidth',2)
    
    hold on

    plot(tWaveLat2,mean(waveDataCell1Bin2,1)/1000,'--', 'color',[0.8500 0.3250 0.0980],'LineWidth',2)
    plot(tWaveLat2,mean(waveDataCell2Bin2,1)/1000,'-.', 'color',[0.8500 0.3250 0.0980],'LineWidth',2)

    hold off
    xlabel('[ms]')
    ylabel('[mV]')

    set(gca,'FontSize',16)
    set(gca, 'YDir','reverse')
    
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    save_file = fullfile(data_dir, ['waveformAna_' session_name '_jscale' num2str(jscale)]);
    set(fig_use,'PaperOrientation','landscape');
    print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
            
    % reset plot
    close 102

    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

close 102


                         
                         
                         