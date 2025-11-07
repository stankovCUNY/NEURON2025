
jscale         = 1;
alpha_name     = 5;
duration       = 0.007;
fs             = 30000;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
plotFlag       = true;
resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

%% Roy 82 v 91 (flipped)

traceType = 'raw'; % 'pca' or 'raw'

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
    % pos = [1720 2562 2*560*0.4 2*420*0.4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
pos = [1720 2562 2*560*0.4 420*0.35]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

%    tiledlayout(2,2,'Padding','none','TileSpacing','compact')
tiledlayout(1,3,'Padding','none','TileSpacing','compact')

load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable

ref = pairGroupStatsTable.tarSpikeTimes{23,1};
tar = pairGroupStatsTable.refSpikeTimes{23,1}; 
refCelltype = pairGroupStatsTable.tarCellExplorerType{23,1};
tarCelltype = pairGroupStatsTable.refCellExplorerType{23,1};
d = pairGroupStatsTable.pairDistance(23);
region = pairGroupStatsTable.brainRegion{23,1};

tarShank = 3;
tarCh    = 23;

nexttile(1)
[GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
          CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                'norm_flag', true);
xlim([-1,1])
ylims = get(gca,'ylim');
ylabel('Spike Probability')
xlabel('[ms]')
    title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

hold on
scatter(tR(109)*1e3,ccgR(109,1,2)/length(ref),100,'bo','LineWidth',1)
scatter(tR(113)*1e3,ccgR(113,1,2)/length(ref),100,'ro','LineWidth',1)
hold off
set(gcf, 'Renderer', 'Painters');
box off

lagLimit = duration/2;

% removing the synchronous spikes                
[CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(ref',tar',binSize,lagLimit);

% peak
refPeak = spikeTimesInBin{109}(:,1);
tarPeak = spikeTimesInBin{109}(:,2);

% trough
refTrough = spikeTimesInBin{113}(:,1);
tarTrough = spikeTimesInBin{113}(:,2);

    tar = setdiff(tar,[tarPeak; tarTrough]);

% % inset
%     hAxes2 = axes('Position',[.2 .7 .1 .2]);
%     plot(tR*1e3,ccgR(:,1,2),'k')
%     hold on
%     scatter(tR(110)*1e3,ccgR(110,1,2),20,'ro','LineWidth',1)
%     scatter(tR(114)*1e3,ccgR(114,1,2),20,'bo','LineWidth',1)
%     hold off
%     set(hAxes2,'xtick',[])
%     set(hAxes2,'ytick',[])
%     set(hAxes2,'xlim',[-1 1])

% PCA traces

if strcmp(traceType,'pca')

    fpass      = 300;

    preLength  = 36;
    postLength = 36;

    nPCsRecovery = 3;

    chanDataDir  = ['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/'];

    chanData = load([ chanDataDir 'ch' num2str(tarCh) 'highpass300hz.mat'],'data');

    noDemeanFlag = true;

    spikeTimeIndxCell = round((tar-chanData.data.onsetTime/1e6)*fs) + 1;
    [waveAvg,waveforms] = waveformAvg(double(chanData.data.channel), ... 
                                              spikeTimeIndxCell, ... 
                                              preLength,postLength,fpass,fs,false,noDemeanFlag);
    [coeff, score] = pca(waveforms,'NumComponents',nPCsRecovery);

    [~,idxPeak,~  ] = intersect(tar,tarPeak);
    [~,idxTrough,~] = intersect(tar,tarTrough);

    filterSpon   = score * coeff(randperm(size(coeff,1),10000),:)';
    filterPeak   = score * coeff(idxPeak,:)';
    filterTrough = score * coeff(idxTrough,:)';

    tiledlayout(1,2)

    nexttile; plot(filterSpon(:,1:100000),'k'); ylim([-6 2])
    nexttile; plot(filterPeak,'k');             ylim([-6 2])
    nexttile; plot(filterTrough,'k');           ylim([-6 2])

    nexttile
    patchline(pairGroupStatsTable.tWave{1,1}*1000, ... 
                 filterSpon(:,1), ...
                 'linewidth',2,'edgealpha',0.1)
    hold on 
    for i = 2:length(filterSpon)
        patchline(pairGroupStatsTable.tWave{1,1}*1000, ... 
                 filterSpon(:,i), ...
                 'linewidth',2,'edgealpha',0.1)
    end
    hold off
    ylim([-6 2])

        nexttile(2)
    nexttile
    patchline(pairGroupStatsTable.tWave{1,1}, ... 
                 filterSpon(:,1), ...
                 'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
    hold on
    for i = 2:length(filterSpon)
        patchline(pairGroupStatsTable.tWave{1,1}, ... 
                 filterSpon(:,i), ...
                 'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
    end
    for i = 1:length(filterPeak)
        patchline(pairGroupStatsTable.tWave{1,1}, ... 
                 filterPeak(:,i), ...
                 'linewidth',1,'edgealpha',0.1,'edgecolor','r')
    end
    hold off
    ylim([-6 2])
    xlim([-1,1])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    set(gca, 'YDir','reverse')
    box off

        nexttile(4)
    nexttile
    patchline(pairGroupStatsTable.tWave{1,1}, ... 
                 filterSpon(:,1), ...
                 'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
    hold on
    for i = 2:length(filterSpon)
        patchline(pairGroupStatsTable.tWave{1,1}, ... 
                 filterSpon(:,i), ...
                 'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
    end
    for i = 1:length(filterTrough)
        patchline(pairGroupStatsTable.tWave{1,1}, ... 
                 filterTrough(:,i), ...
                 'linewidth',1,'edgealpha',0.2,'edgecolor','b') 
    end
    hold off
    ylim([-6 2])
    xlim([-1,1])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    set(gca, 'YDir','reverse')
    set(gcf, 'Renderer', 'Painters');
    box off

    legend('spon','peak','trough')

% snippets 
elseif strcmp(traceType,'raw')

    % weird timeshift in Roy's data
    timeShift      = 514;

    chanData       = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarCh) 'highpass300hz.mat']);

    % chanData       = chanData.data.channel;

    filterFlag     = false; 
    noDemeanFlag   = false;
    fpass          = 300;
    preLength      = 5*30;
    postLength     = 5*30;

    bx = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-behavior.mat');

    % tarLat       = round(tar*fs)'       - bx.behavior.RoyMaze1.datFrame(2,1) - timeShift;
    % tarPeakLat   = round(tarPeak*fs)'   - bx.behavior.RoyMaze1.datFrame(2,1) - timeShift;
    % tarTroughLat = round(tarTrough*fs)' - bx.behavior.RoyMaze1.datFrame(2,1) - timeShift;

    tarLat       = round((tar      -chanData.data.onsetTime/1e6)*fs) + 1;
    tarPeakLat   = round((tarPeak  -chanData.data.onsetTime/1e6)*fs) + 1;
    tarTroughLat = round((tarTrough-chanData.data.onsetTime/1e6)*fs) + 1;

    tarLatRandPerm = sort(tarLat(randperm(length(tarLat),10000)));

    [tarSponWaveMean,     tarSponWaveforms]      = waveformAvg(chanData.data.channel,tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
    [tarPeakLatWaveMean,  tarPeakLatWaveforms]   = waveformAvg(chanData.data.channel,tarPeakLat,    preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
    [tarTroughLatWaveMean,tarTroughLatWaveforms] = waveformAvg(chanData.data.channel,tarTroughLat,  preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

        % nexttile(2)
    nexttile

    patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,1), ...
                     'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
    hold on
    for i = 2:size(tarSponWaveforms,2)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                 tarSponWaveforms(:,i), ...
                 'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
    end
    for i = 1:size(tarPeakLatWaveforms,2)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                 tarPeakLatWaveforms(:,i), ...
                 'linewidth',0.01,'edgealpha',0.1,'edgecolor','b')
    end
    hold off
    ylim([-8 3])
    xlim([-1,1])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    set(gca, 'YDir','reverse')
    box off
    set(gcf, 'Renderer', 'Painters');

        % nexttile(4)
    nexttile

    patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,1), ...
                     'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
    hold on
    for i = 2:size(tarSponWaveforms,2)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                 tarSponWaveforms(:,i), ...
                 'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
    end
    for i = 1:100 %size(tarTroughLatWaveforms,2)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                 tarTroughLatWaveforms(:,i), ...
                 'linewidth',0.01,'edgealpha',0.1,'edgecolor','r')
    end
    hold off
    ylim([-8 3])
    xlim([-1,1])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    set(gca, 'YDir','reverse')
    box off
    set(gcf, 'Renderer', 'Painters');

end

% cluster plots

% Monitor specific plot settings.
screensize = get(0,'screensize');
%initiate figure
hcomb = figure(103);

res_type = 'QHD';
pos = [1720 2562 3*560*0.33 2*420*0.5]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

tiledlayout(2,3,'Padding','none')

% source for tarUnitSpikeCountAllSess
% spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
% [data_dir, name, ~] = fileparts(spike_data_fullpath);
% load(spike_data_fullpath, 'spikes')

% peak
tarUnitSpikeCountAllSess = 96594; % unit 82
plotParam = 'ind';
synchSpikesFlag = true;
refUnitSpikeTimes = refPeak;

curateClustersForPlots(tarShank, ...
                       tarCh, ...
                       tarUnitSpikeCountAllSess, ...
                       plotParam, ...
                       synchSpikesFlag, ...
                       refUnitSpikeTimes, ...
                       'b', ...
                       109)

% trough                   
% tarShank = 1;
% tarCh    = 1;
tarUnitSpikeCountAllSess = 96594; % unit 82
plotParam = 'ind';
synchSpikesFlag = true;
refUnitSpikeTimes = refTrough;

curateClustersForPlots(tarShank, ...
                       tarCh, ...
                       tarUnitSpikeCountAllSess, ...
                       plotParam, ...
                       synchSpikesFlag, ...
                       refUnitSpikeTimes,...
                       'r', ...
                       113)

set(gcf, 'Renderer', 'Painters');

%%

% filterFlag     = false; 
% noDemeanFlag   = false;
% fpass          = 300;
% fs             = 3e4;
% preLength      = 5*30;
% postLength     = 5*30;
% 
% figPath = '/media/nasko/WD_BLACK31/ClusterAnaTemp/figures/';
% 
% resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
% fig_use        = 102;
% 
% animalNames    = {'N','S','U','AG1','AG2','Roy'};
% 
% NspikesPlot    = 1000;
% 
% ratListItem    = [5,  ...  % Rat AG
%                   3];      % Rat U
% pairListItem   = [11,  ... % unit 33 vs. unit 138 
%                   13];     % unit 52 vs. unit 129
% 
% itemNo = 2;
% 
% loopRat   = ratListItem(itemNo);
% loopPairs = pairListItem(itemNo);
% 
% if loopRat == 1
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat')  
%     neuronsN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/unitsGlobalSnapshots.mat')
% 
%     channelSet = {1:16, ...
%                   17:32, ... 
%                   33:48, ...
%                   49:64, ...
%                   65:80, ...
%                   81:96, ...
%                   97:112, ...
%                   113:128};
% 
% elseif loopRat == 2
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStatsClusterAna.mat') 
%     neuronsS = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/unitsGlobalSnapshots.mat')
% 
%     channelSet = {1:10, ...
%                   11:21, ...
%                   22:32, ... 
%                   33:43, ...
%                   44:54, ...
%                   55:63, ...
%                   64:79, ...
%                   80:95, ...
%                   96:111, ...
%                   112:127, ...
%                   128:143, ...
%                   144:159, ...
%                   160:175, ...
%                   176:191};
% 
% elseif loopRat == 3
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAna.mat') 
%     neuronsU = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
% elseif loopRat == 4
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStatsClusterAna.mat')
%     UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
%     cell_metrics = load([UtkuPath 'AG_2019-12-23_NSD' '/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);
% elseif loopRat == 5
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStatsClusterAna.mat')
%     UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
%     cell_metrics = load([UtkuPath 'AG_2019-12-27_NSD' '/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);
% elseif loopRat == 6
%     load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
% 
%     spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
%     load(spike_data_fullpath, 'spikes')
% 
%     bx = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-behavior.mat');
% 
%     cellIDs         = 1:size(spikes.RoyMaze1,2);
%     neuronTemp      = cell2num({spikes.RoyMaze1.id}');
% 
%     cellTypeList = strings(size(spikes.RoyMaze1,2),1);
%     cellTypeList([spikes.RoyMaze1.quality] == 8) = "i";
%     cellTypeList(neuronTemp(:,2) == 1)           = "noise";
%     cellTypeList(strcmp(cellTypeList,""))        = "p"; % double check!!
% end
% 
% refSpikeTimes               = pairGroupStatTable.tarSpikeTimes{loopPairs};
% refShank                    = pairGroupStatTable.tarShank(loopPairs);
% refChan                     = pairGroupStatTable.tarChannel(loopPairs);
% refID                       = pairGroupStatTable.tarNeuronID(loopPairs);
% % refSpikeTimesSynch          = pairGroupStatTable.tarSynch{loopPairs};
% refCellType                 = pairGroupStatTable.tarCellExplorerType{loopPairs};
% 
% tarSpikeTimes               = pairGroupStatTable.refSpikeTimes{loopPairs};
% tarShank                    = pairGroupStatTable.refShank(loopPairs);
% tarChan                     = pairGroupStatTable.refChannel(loopPairs);
% tarID                       = pairGroupStatTable.refNeuronID(loopPairs);
% % tarSpikeTimesSynch          = pairGroupStatTable.refSynch{loopPairs};
% tarCellType                 = pairGroupStatTable.refCellExplorerType{loopPairs};
% 
% d                           = pairGroupStatTable.pairDistance(loopPairs);
% 
% tWave = pairGroupStatTable.tWave{loopPairs};
% 
% lagLimit = duration/2;
% 
% % removing the synchronous spikes                
% [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(refSpikeTimes',tarSpikeTimes',binSize,lagLimit);
% 
% % peak (trough too)
% if itemNo == 1
%     refSpikeTimesSynch = spikeTimesInBin{119}(:,1);
%     tarSpikeTimesSynch = spikeTimesInBin{119}(:,2);
% elseif itemNo == 2
%     refSpikeTimesSynch = spikeTimesInBin{113}(:,1);
%     tarSpikeTimesSynch = spikeTimesInBin{113}(:,2);
% end
% 
% [~,refSpikeTimesSynchIdx,~] = intersect(refSpikeTimes,refSpikeTimesSynch);
% [~,tarSpikeTimesSynchIdx,~] = intersect(tarSpikeTimes,tarSpikeTimesSynch);
% 
%     %% 
% 
% titleStrRef = ['rat ' animalNames{loopRat} ': '...
%                 num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')' ' - ' ...
%                 num2str(refID) refCellType ' (sh: ' num2str(refShank) ')' ' pair ' num2str(d) '\mum '];
% 
% titleStrTar = ['rat ' animalNames{loopRat} ': '...
%                 num2str(refID) refCellType ' (sh: ' num2str(refShank) ')' ' - ' ...
%                 num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')' ' pair ' num2str(d) '\mum '];
% 
% saveStr  = ['rat ' animalNames{loopRat} ': '...
%             num2str(refID) refCellType ' (sh ' num2str(refShank) ')' ' - ' ...
%             num2str(tarID) tarCellType ' (sh ' num2str(tarShank) ')' ' pair ' num2str(d) ' microns.jpeg'];
% 
% % [~,tarChan] = max(max(abs(unitsGlobalSnapshots.waveforms{tarID,1}(channelSet{tarShank},:)),[],2));
% 
% if loopRat == 1
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'    num2str(refShank) '.mat']) 
%     clusterValNeuronsRef = clusterValNeuron;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'    num2str(tarShank) '.mat']) 
%     clusterValNeuronsTar = clusterValNeuron;
% 
%     loopIdxTar = find((neuronsN.neuron_ids + 1) == tarID);
%     loopIdxRef = find((neuronsN.neuron_ids + 1) == refID);
% 
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
%     clusterValNeuronsPCAvalsRef = clusterValAllChans;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
%     clusterValNeuronsPCAvalsTar = clusterValAllChans;
% 
%     clear clusterValNeuron
%     clear clusterValAllChans
% elseif loopRat == 2
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'    num2str(refShank) '.mat']);
%     clusterValNeuronsRef = clusterValNeuron;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'    num2str(tarShank) '.mat']);
%     clusterValNeuronsTar = clusterValNeuron;
% 
%     loopIdxTar = find((neuronsS.neuron_ids + 1) == tarID);
%     loopIdxRef = find((neuronsS.neuron_ids + 1) == refID);
% 
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
%     clusterValNeuronsPCAvalsRef = clusterValAllChans;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
%     clusterValNeuronsPCAvalsTar = clusterValAllChans;
% 
%     clear clusterValNeuron
%     clear clusterValAllChans
% elseif loopRat == 3
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'    num2str(refShank) '.mat'])
%     clusterValNeuronsRef = clusterValNeuron;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'    num2str(tarShank) '.mat'])
%     clusterValNeuronsTar = clusterValNeuron;
% 
%     loopIdxTar = find((neuronsU.neuron_ids + 1) == tarID);
%     loopIdxRef = find((neuronsU.neuron_ids + 1) == refID);
% 
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
%     clusterValNeuronsPCAvalsRef = clusterValAllChans;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
%     clusterValNeuronsPCAvalsTar = clusterValAllChans;
% 
%     clear clusterValNeuron
%     clear clusterValAllChans
% elseif loopRat == 4
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(refShank) '.mat'])
%     clusterValNeuronsRef = clusterValNeuron;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(tarShank) '.mat'])
%     clusterValNeuronsTar = clusterValNeuron;
% 
%     loopIdxTar = find(cell_metrics.cell_metrics.UID == tarID);
%     loopIdxRef = find(cell_metrics.cell_metrics.UID == refID);
% 
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
%     clusterValNeuronsPCAvalsRef = clusterValAllChans;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
%     clusterValNeuronsPCAvalsTar = clusterValAllChans;
% 
%     clear clusterValNeuron
%     clear clusterValAllChans
% elseif loopRat == 5
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(refShank) '.mat'])
%     clusterValNeuronsRef = clusterValNeuron;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(tarShank) '.mat'])
%     clusterValNeuronsTar = clusterValNeuron;
% 
%     loopIdxTar = find(cell_metrics.cell_metrics.UID == tarID);
%     loopIdxRef = find(cell_metrics.cell_metrics.UID == refID);
% 
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
%     clusterValNeuronsPCAvalsRef = clusterValAllChans;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
%     clusterValNeuronsPCAvalsTar = clusterValAllChans;
% 
%     clear clusterValNeuron
%     clear clusterValAllChans
% elseif loopRat == 6
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(refShank) '.mat'])
%     clusterValNeuronsRef = clusterValNeuron;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(tarShank) '.mat'])
%     clusterValNeuronsTar = clusterValNeuron;
% 
%     loopIdxTar = find(cellIDs == tarID);
%     loopIdxRef = find(cellIDs == refID);
% 
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
%     clusterValNeuronsPCAvalsRef = clusterValAllChans;
%     load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
%     clusterValNeuronsPCAvalsTar = clusterValAllChans;
% 
%     clear clusterValNeuron
%     clear clusterValAllChans
% end
% 
% %% gather PC values 
% 
% % PCAchans = clusterValNeuron.templateChannels{clusterValNeuron.cellID == tarID};
% % PCAchansIdx = ~isnan(PCAchans);
% 
% if loopRat == 6
%     cluIdxRef  = 1;
%     cluIdxTar  = 1;
% else
%     cluIdxRef  = clusterValNeuronsRef.clusIdx(find(clusterValNeuronsRef.cellID == refID));
%     cluIdxTar  = clusterValNeuronsTar.clusIdx(find(clusterValNeuronsTar.cellID == tarID));
% end
% 
% PCAvalsRef = clusterValNeuronsPCAvalsRef;
% PCAvalsTar = clusterValNeuronsPCAvalsTar;
% % PCAvals = {PCAvals{1,PCAchansIdx}};
% 
% PCAref      = PCAvalsRef{cluIdxRef};
% PCAtar      = PCAvalsTar{cluIdxTar};
% PCAsynchRef = PCAref(refSpikeTimesSynchIdx(refSpikeTimesSynchIdx < size(PCAref,1)),:); % needed cutoff
% PCAsynchTar = PCAtar(tarSpikeTimesSynchIdx(tarSpikeTimesSynchIdx < size(PCAtar,1)),:); % needed cutoff
% 
% %% template indice split 
% 
% spikeTemplatesTar             = clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)};
% spikeTemplatesSynchTar        = clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)}(tarSpikeTimesSynchIdx(tarSpikeTimesSynchIdx < size(PCAtar,1))); % needed cutoff 
% 
% spikeTemplatesRef             = clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)};
% spikeTemplatesSynchRef        = clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)}(refSpikeTimesSynchIdx(refSpikeTimesSynchIdx < size(PCAref,1)));
% 
% %% getting the waveform snippets 
% 
% % timeShift      = 514;
% % chanData       = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarCh) 'highpass300hz.mat']);
% 
% if loopRat == 1
%     chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(tarChan) 'highpass300hz.mat']);
%     chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(refChan) 'highpass300hz.mat']);
% elseif loopRat == 2
%     chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(tarChan) 'highpass300hz.mat']);
%     chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(refChan) 'highpass300hz.mat']);
% elseif loopRat == 3
%     chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(tarChan) 'highpass300hz.mat']);
%     chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(refChan) 'highpass300hz.mat']);
% elseif loopRat == 4
%     chanDataRef = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(refShank) '/ch' num2str(refChan) 'highpass300hz.mat']);
%     chanDataTar = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(tarShank) '/ch' num2str(tarChan) 'highpass300hz.mat']);
% elseif loopRat == 5
%     chanDataRef = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(refShank) '/ch' num2str(refChan) 'highpass300hz.mat']);
%     chanDataTar = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(tarShank) '/ch' num2str(tarChan) 'highpass300hz.mat']);
% elseif loopRat == 6
%     chanDataRef = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(refChan) 'highpass300hz.mat']);
%     chanDataTar = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarChan) 'highpass300hz.mat']);
% end
% 
% if loopRat == 6
%     refLat       = round((refSpikeTimes       - chanDataRef.data.onsetTime/1e6)*fs) + 1;
%     refSynchLat  = round((refSpikeTimesSynch  - chanDataRef.data.onsetTime/1e6)*fs) + 1;
% 
%     tarLat       = round((tarSpikeTimes       - chanDataTar.data.onsetTime/1e6)*fs) + 1;
%     tarSynchLat  = round((tarSpikeTimesSynch  - chanDataTar.data.onsetTime/1e6)*fs) + 1;
% else
%     refLat       = round(refSpikeTimes*fs)';
%     refSynchLat  = round(refSpikeTimesSynch*fs)';
% 
%     tarLat       = round(tarSpikeTimes*fs)';
%     tarSynchLat  = round(tarSpikeTimesSynch*fs)';
% end
% 
% %%
% 
% if length(refLat) > NspikesPlot
%     refLatRandPerm = sort(refLat(randperm(length(refLat),NspikesPlot)));
% else
%     refLatRandPerm = refLat;
% end
% 
% if length(refSynchLat) > NspikesPlot
%     refSynchLatRandPerm = sort(refSynchLat(randperm(length(refSynchLat),NspikesPlot)));
% else
%     refSynchLatRandPerm = refSynchLat;
% end
% 
% if loopRat == 6
%     [refSponWaveMean,      refSponWaveforms]                                       = waveformAvg(chanDataRef.data.channel,refLatRandPerm,     preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     [refSynchLatWaveMean,  refPeakLatWaveforms]                                    = waveformAvg(chanDataRef.data.channel,refSynchLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
% else
%     [refSponWaveMean,      refSponWaveforms]                                       = waveformAvg(chanDataRef.data,        refLatRandPerm,      preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     [refSynchLatWaveMean,  refPeakLatWaveforms]                                    = waveformAvg(chanDataRef.data,        refSynchLatRandPerm, preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
% end
% 
% if length(tarLat) > NspikesPlot
%     tarLatRandPerm = sort(tarLat(randperm(length(tarLat),NspikesPlot)));
% else
%     tarLatRandPerm = tarLat;
% end
% 
% if length(tarSynchLat) > NspikesPlot
%     tarSynchLatRandPerm = sort(tarSynchLat(randperm(length(tarSynchLat),NspikesPlot)));
% else
%     tarSynchLatRandPerm = tarSynchLat;
% end
% 
% 
% if loopRat == 6
%     [tarSponWaveMean,      tarSponWaveforms]     = waveformAvg(chanDataTar.data.channel, tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     [tarSynchLatWaveMean,  tarPeakLatWaveforms]  = waveformAvg(chanDataTar.data.channel, tarSynchLatRandPerm,    preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
% else
%     [tarSponWaveMean,      tarSponWaveforms]     = waveformAvg(chanDataTar.data,         tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     [tarSynchLatWaveMean,  tarPeakLatWaveforms]  = waveformAvg(chanDataTar.data,         tarSynchLatRandPerm,    preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
% end
% 
% %%    
% 
% if itemNo == 1
%     synchPlotColor = 'r';
% elseif itemNo == 2
%     synchPlotColor = 'b';
% end
% 
% % Monitor specific plot settings.
% screensize = get(0,'screensize');
% % initiate figure
% hcomb = figure(102);
% 
% res_type = 'QHD';
% %     pos = [1720 2562 2*560*0.4 2*420*0.4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
% pos = [1720 2562 2*560*0.4 2*420*0.35]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
% arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% 
% %     tiledlayout(2,2,'Padding','none','TileSpacing','compact')
%     tiledlayout(2,3,'Padding','none','TileSpacing','compact')
% nexttile(1)
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%           CCG_jitter(refSpikeTimes,tarSpikeTimes,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
%                                 'plot_output', get(fig_use, 'Number'), ...
%                                 'njitter', njitter, 'alpha', alpha,...
%                                 'for_grant', for_grant, 'plot_pointwiseBands', false, ...
%                                 'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
%                                 'norm_flag', true);
% xlim([-1,1])
% ylims = get(gca,'ylim');
% ylabel('Spike Probability')
% xlabel('[ms]')
% %     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% hold on
% if itemNo == 1
%     scatter(tR(119)*1e3,ccgR(119,1,2)/length(refSpikeTimes),100,'ro','LineWidth',1)
% elseif itemNo == 2
%     scatter(tR(113)*1e3,ccgR(113,1,2)/length(refSpikeTimes),100,'bo','LineWidth',1)
% end
% 
% hold off
% set(gcf, 'Renderer', 'Painters');
% box off
% 
% % tR     = pairGroupStatTable.CCGbinLagTimes{loopPairs}*1000;
% % ccgR   = pairGroupStatTable.pairRawCCG{loopPairs};
% % GSPExc = pairGroupStatTable.GSPExc{loopPairs};
% % 
% % nexttile(1)
% % plot(tR,ccgR,'k','LineWidth',1)
% % hold on
% % scatter(tR(find(GSPExc)),ccgR(find(GSPExc)),100,'ko','LineWidth',1)
% % hold off
% % 
% % xlim([-1,1])
% % ylims = get(gca,'ylim');
% % ylabel('Spike Probability')
% % xlabel('[ms]')
% % title(titleStrTar)
% % set(gca,'FontSize',5)
% % set(gca,'FontName','Arial')
% % box off
% 
% %%
% 
% nexttile(2)
% patchline((-preLength:postLength-1)*(30/1000), ... 
%                  tarSponWaveforms(:,1), ...
%                  'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
% hold on
% for i = 2:size(tarSponWaveforms,2)
%     patchline((-preLength:postLength-1)*(30/1000), ... 
%              tarSponWaveforms(:,i), ...
%              'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
% end
% if itemNo == 1
%     for i = 1:size(tarPeakLatWaveforms,2)
%         patchline((-preLength:postLength-1)*(30/1000), ... 
%                  tarPeakLatWaveforms(:,i), ...
%                  'linewidth',1,'edgealpha',0.1,'edgecolor','r')
%     end
% elseif itemNo == 2
%     for i = 1:size(tarPeakLatWaveforms,2)
%         patchline((-preLength:postLength-1)*(30/1000), ... 
%                  tarPeakLatWaveforms(:,i), ...
%                  'linewidth',1,'edgealpha',0.1,'edgecolor','b')
%     end
% end
% 
% hold off
% % ylim([-6 2])
% xlim([-1,1])
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% xlabel('[ms]');
% ylabel('[mV]');
% set(gca, 'YDir','reverse')
% title('waveform snippets')
% box off
% 
% %%
% 
% sz = [];
% 
% nexttile(4)
% scatter(PCAtar(:,1),        PCAtar(:,2),        sz,'+','MarkerEdgeColor',  [.7 .7 .7],    'LineWidth',0.25)
% hold on
% scatter(PCAsynchTar(:,1),   PCAsynchTar(:,2),   sz,'+','MarkerEdgeColor',  synchPlotColor,'LineWidth',0.25)
% hold off
% xlabel('PC1')
% ylabel('PC2')
% % xlim([-2500 2500])
% % ylim([-2500 2500])
% daspect([1 1 1])
% title('PC1 v PC2')
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% nexttile(5)
% scatter(PCAtar(:,2),        PCAtar(:,3),        sz,'+','MarkerEdgeColor',  [.7 .7 .7],    'LineWidth',0.25)
% hold on
% scatter(PCAsynchTar(:,2),   PCAsynchTar(:,3),   sz,'+','MarkerEdgeColor',  synchPlotColor,'LineWidth',0.25)
% hold off
% xlabel('PC1')
% ylabel('PC3')
% % xlim([-2500 2500])
% % ylim([-2500 2500])
% daspect([1 1 1])
% title('PC1 v PC3')
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% nexttile(6)
% scatter(PCAtar(:,1),     PCAtar(:,3),        sz,'+','MarkerEdgeColor',  [.7 .7 .7],    'LineWidth',0.25)
% hold on
% scatter(PCAsynchTar(:,1),PCAsynchTar(:,3),   sz,'+','MarkerEdgeColor',  synchPlotColor,'LineWidth',0.25)
% hold off
% xlabel('PC2')
% ylabel('PC3')
% % xlim([-2500 2500])
% % ylim([-2500 2500])
% daspect([1 1 1])
% title('PC2 v PC3')
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% % legend('spon largest template','synch largest template', ...
% %        'spon all other templates','synch all other templates', ...
% %        'Location','eastoutside')


%%

