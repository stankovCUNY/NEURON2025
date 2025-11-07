%% 

clear all; clc

% constants
session_name   = 'RoyMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1;
alpha_name     = 5;
duration       = 0.001;
fs             = 30000;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
numWaves       = 100;
resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

chSimpleLayout = ...
    [1 8 9  16 17 24 25 32 33 40 41 48 49 56 57 64 ...
     2 7 10 15 18 23 26 31 34 39 42 47 50 55 58 63 ...
     3 6 11 14 19 22 27 30 35 38 43 46 51 54 59 62 ...
     4 5 12 13 20 21 28 29 36 37 44 45 52 53 60 61];

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));


% pair data 
datapath   = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";
figurepath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/figures/";

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

%% get waveforms
extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
RoySession1 = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);

% get all pairs
cell_pairs = {jit_var_out.cell_pair}';

% convert to matrix
cell_pairs_mat = [];
for i = 1:size(cell_pairs,1)
    cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
end

%% (NOTE: the indexing work here is a mess and redundant but it works!)

%maze only
idxMaze = cellfun(@(a) strmatch(a,{jit_var_out.epoch}),{'Maze'},'uniform',false);

jit_var_outMaze = jit_var_out(idxMaze{1});

matIdx = [];
for i = 1:size(jit_var_outMaze,2)
    matIdx = [matIdx; jit_var_outMaze(i).cell_pair];
end

% unique pairs only
[~,~,IC] = unique(sort(matIdx,2),'row','first');

matIdxNew = [];
for i = 1:max(IC)
    tmp = find(IC == i);
    matIdxNew = [matIdxNew; matIdx(tmp(1),:)];
end
    
matIdx = matIdxNew;

excHist = zeros(211,1);
inhHist = zeros(211,1);

pairCCHstack        = [];
pairCCHstackIDs     = {};
pairCCHstackIDschar = {};

tiledlayout(4,16, 'Padding', 'none', 'TileSpacing', 'compact'); 

for j = 1:size(matIdx,1)
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = matIdx(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    excHist = excHist + jit_var_out(current_pair_indices(2)).GSPExc;
    inhHist = inhHist + jit_var_out(current_pair_indices(2)).GSPInh;
    
    stack = zeros(211,1);
    stack(find(jit_var_out(current_pair_indices(2)).GSPExc)) =  1;
    stack(find(jit_var_out(current_pair_indices(2)).GSPInh)) = -1;
    
    if sum(jit_var_out(current_pair_indices(2)).GSPExc + ...
           jit_var_out(current_pair_indices(2)).GSPInh) ~= 0
        pairCCHstack = [pairCCHstack stack];
        pairCCHstackIDs{j} = current_pair;
    end
end
pairCCHstackIDs = pairCCHstackIDs(~cellfun('isempty',pairCCHstackIDs));

for i = 1:size(pairCCHstackIDs,2)
    matTemp = pairCCHstackIDs{i};
    pairCCHstackIDschar{i} = sprintf("%d %d",matTemp(1),matTemp(2));
end

waveDataSpCounts = [];
for i = 1:size(RoySession1.times,2)
    waveDataSpCounts(i) = size(RoySession1.times{i},1);
end

%% sig pairs
 
pairsAna = [];
for i = 1:length(pairCCHstackIDs)
    pairsAna = [pairsAna; pairCCHstackIDs{i}];
end

% unique neurons
refNeuron = unique(pairsAna(:));

tiledlayout(4,16, 'Padding', 'none', 'TileSpacing', 'compact'); 
tWave = (-26:27)*(1/30);

for i = 1:length(refNeuron)
    
    colOne  = (refNeuron(i) == pairsAna(:,1));
    colTwo  = (refNeuron(i) == pairsAna(:,2));
    colBoth = (colOne | colTwo);
    
    colIdxFlipped = find(colTwo);
    
    copyNetwork = pairsAna;
    
    for j = 1:length(colIdxFlipped)
        copyNetwork(colIdxFlipped(j),:) = fliplr(pairsAna(colIdxFlipped(j),:));
    end
    
    refNetwork = copyNetwork(find(colBoth),:);
    
    dataNetPlot = struct;
    for j = 1:size(refNetwork,1)

        flipped = false;

        idx = j; % pair to run, right now we are working with pre-maze time period

        current_pair = refNetwork(idx,:);

        current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));

        if isempty(current_pair_indices)
            current_pair = fliplr(current_pair);
            current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));
            flipped = true;
        end

        numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

        numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

        waveIdxCell1 = find(numSpikesCell1 == waveDataSpCounts);
        waveIdxCell2 = find(numSpikesCell2 == waveDataSpCounts);

        if ~flipped
            
            waveDataCell1 = RoySession1.allSpikeRawWaveforms{waveIdxCell1};
            waveDataCell2 = RoySession1.allSpikeRawWaveforms{waveIdxCell2};

            % filter waveDataSpCounts to maze sessions 
            dataNetPlot(j).waveDataCell1 = waveDataCell1(length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + 1: ...
                                                         length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                                                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times),:);

            dataNetPlot(j).waveDataCell2 = waveDataCell2(length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + 1: ...
                                                         length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                                                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times),:);

            dataNetPlot(j).chCell1 = RoySession1.maxWaveformCh(waveIdxCell1);
            dataNetPlot(j).chCell2 = RoySession1.maxWaveformCh(waveIdxCell2);

            dataNetPlot(j).cell1ID = jit_var_out(current_pair_indices(2)).cell_pair(1);
            dataNetPlot(j).cell2ID = jit_var_out(current_pair_indices(2)).cell_pair(2);

            dataNetPlot(j).cell1shank = jit_var_out(current_pair_indices(2)).cell1shank;
            dataNetPlot(j).cell2shank = jit_var_out(current_pair_indices(2)).cell2shank;

            dataNetPlot(j).cell1type = jit_var_out(current_pair_indices(2)).cell1type;
            dataNetPlot(j).cell2type = jit_var_out(current_pair_indices(2)).cell2type;

            res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
            res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
            
            dataNetPlot(j).cell1spikeTimes = res1;
            dataNetPlot(j).cell2spikeTimes = res2;

        elseif flipped
            
            waveDataCell2 = RoySession1.allSpikeRawWaveforms{waveIdxCell1};
            waveDataCell1 = RoySession1.allSpikeRawWaveforms{waveIdxCell2};

            % filter waveDataSpCounts to maze sessions 
            dataNetPlot(j).waveDataCell2 = waveDataCell2(length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + 1: ...
                                                         length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                                                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times),:);

            dataNetPlot(j).waveDataCell1 = waveDataCell1(length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + 1: ...
                                                         length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                                                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times),:);

            dataNetPlot(j).chCell1 = RoySession1.maxWaveformCh(waveIdxCell2);
            dataNetPlot(j).chCell2 = RoySession1.maxWaveformCh(waveIdxCell1);

            dataNetPlot(j).cell1ID = jit_var_out(current_pair_indices(2)).cell_pair(2);
            dataNetPlot(j).cell2ID = jit_var_out(current_pair_indices(2)).cell_pair(1);

            dataNetPlot(j).cell1shank = jit_var_out(current_pair_indices(2)).cell2shank;
            dataNetPlot(j).cell2shank = jit_var_out(current_pair_indices(2)).cell1shank;

            dataNetPlot(j).cell1type = jit_var_out(current_pair_indices(2)).cell2type;
            dataNetPlot(j).cell2type = jit_var_out(current_pair_indices(2)).cell1type;

            res1 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
            res2 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
            
            dataNetPlot(j).cell1spikeTimes = res1;
            dataNetPlot(j).cell2spikeTimes = res2;

        end

        %% make CCGs

        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
              CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);

        CCG = ccgR(:,1,2)/sum(ccgR(:,1,2));
        dataNetPlot(j).CCG = CCG;
        nSigPeaks   = sum(GSPExc);
        nSigTroughs = sum(GSPInh);
        
        if (nSigPeaks + nSigTroughs) == 0
            dataNetPlot(j) = [];
        end
%             ylims = get(gca,'ylim');
%             ylim(ylims)


%             title(sprintf('%d%s v %d%s', ...
%                         ,, ...
%                         ,jit_var_out(current_pair_indices(2)).cell2type))

%             if (j == 1 || j == 2); xlabel(''); end
%             if (j == 2 || j == 4); ylabel(''); end

%             make_figure_pretty(102)
    %     set(gca,'FontSize',16)

    end
    
    emptyRowFields = [];
    for j = 1:size(dataNetPlot,2)
        emptyRowFields = [emptyRowFields isempty(dataNetPlot(j).CCG)];
    end
    dataNetPlot(find(emptyRowFields)) = [];
    
    targetCh =[dataNetPlot.chCell2];
    uniqueCh = unique(targetCh);
    
    for j = 1:length(uniqueCh)
        clear legendStr
        nexttile(find(chSimpleLayout == uniqueCh(j)))

        occurances = find(targetCh == uniqueCh(j));
        plot(tR*1000,dataNetPlot(occurances(1)).CCG)
        if length(occurances) > 1
            hold on
            for k = 2:length(occurances)
                 plot(tR*1000,dataNetPlot(occurances(k)).CCG)
            end
            hold off
            
            legendStr{1} = ['tar'  num2str(dataNetPlot(occurances(1)).cell2ID) dataNetPlot(occurances(1)).cell2type];
            for k = 2:length(occurances)
                legendStr{k} = ['tar'  num2str(dataNetPlot(occurances(k)).cell2ID) dataNetPlot(occurances(1)).cell2type];
            end
        else 
            legendStr = ['tar'  num2str(dataNetPlot(occurances(1)).cell2ID) dataNetPlot(occurances(1)).cell2type];
        end
        
        
        legend(legendStr,'Location','southoutside')
        
        ylabel('prob.')
        xlabel('[ms]')
        title(sprintf('sh:%d ch:%d',dataNetPlot(occurances(1)).cell2shank,uniqueCh(j)))
    end
    
    if isempty(dataNetPlot)
        continue
    end
    
    refChPlot = find(chSimpleLayout == dataNetPlot(1).chCell1);
    nexttile(refChPlot)
    plot(tWave,mean(dataNetPlot(occurances(1)).waveDataCell1)/1000)
    xlim([-0.5 0.5])
    ylabel('[mV]')
    xlabel('[ms]')
    title(['\color{red}' sprintf('sh:%d ch:%d',dataNetPlot(1).cell1shank,chSimpleLayout(refChPlot))])
    set(gca, 'YDir','reverse')
    legendStr = ['ref' num2str(dataNetPlot(1).cell1ID) dataNetPlot(1).cell1type];
    legend(legendStr,'Location','southoutside')
    
    save_file = fullfile(datapath, ['spatialCCGs_neuron_' num2str(dataNetPlot(1).cell1ID) dataNetPlot(1).cell1type '_' session_name]);
%     set(fig_use,'PaperOrientation','landscape');
    print(fig_use, save_file,'-djpeg',resolution_use);
    
    close all
    
    save([char(datapath) 'dataSpatialCCGs_neuron_' ...
        num2str(dataNetPlot(1).cell1ID) dataNetPlot(1).cell1type ... 
        '_' session_name '.mat'],'dataNetPlot')
    
    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    tiledlayout(4,16, 'Padding', 'none', 'TileSpacing', 'compact'); 

end

