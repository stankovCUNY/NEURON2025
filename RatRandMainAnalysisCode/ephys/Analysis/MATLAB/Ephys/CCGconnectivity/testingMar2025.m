
filtarFlag     = false; 
noDemeanFlag   = false;
fpass          = 300;
fs             = 3e4;
preLength      = 5*30;
postLength     = 5*30;

figPath = '/media/nasko/WD_BLACK31/ClusterAnaTemp/figures/';

resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
fig_use        = 102;

animalNames = {'N','S','U','AG1','AG2','Roy'};

NspikesPlot = 1000;
NsortSet    = 500;

for loopRat = 1:6

    if loopRat == 1
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat')  
        neuronsN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/unitsGlobalSnapshots.mat')

        channelSet = {1:16, ...
                      17:32, ... 
                      33:48, ...
                      49:64, ...
                      65:80, ...
                      81:96, ...
                      97:112, ...
                      113:128};

    elseif loopRat == 2
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStatsClusterAna.mat') 
        neuronsS = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/unitsGlobalSnapshots.mat')

        channelSet = {1:10, ...
                      11:21, ...
                      22:32, ... 
                      33:43, ...
                      44:54, ...
                      55:63, ...
                      64:79, ...
                      80:95, ...
                      96:111, ...
                      112:127, ...
                      128:143, ...
                      144:159, ...
                      160:175, ...
                      176:191};

    elseif loopRat == 3
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAna.mat') 
        neuronsU = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
    elseif loopRat == 4
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStatsClusterAna.mat')
        UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
        cell_metrics = load([UtkuPath 'AG_2019-12-23_NSD' '/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);
    elseif loopRat == 5
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStatsClusterAna.mat')
        UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
        cell_metrics = load([UtkuPath 'AG_2019-12-27_NSD' '/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);
    elseif loopRat == 6
        load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')

        spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
        load(spike_data_fullpath, 'spikes')

        bx = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-behavior.mat');

        cellIDs         = 1:size(spikes.RoyMaze1,2);
        neuronTemp      = cell2num({spikes.RoyMaze1.id}');

        cellTypeList = strings(size(spikes.RoyMaze1,2),1);
        cellTypeList([spikes.RoyMaze1.quality] == 8) = "i";
        cellTypeList(neuronTemp(:,2) == 1)           = "noise";
        cellTypeList(strcmp(cellTypeList,""))        = "p"; % double check!!
    end

    for loopPairs = 1:size(pairGroupStatTable,1)

        refShank           = pairGroupStatTable.refShank(loopPairs);
        refChan            = pairGroupStatTable.refChannel(loopPairs);
        refID              = pairGroupStatTable.refNeuronID(loopPairs);
        refSpikeTimes      = pairGroupStatTable.refSpikeTimes{loopPairs};
        refSpikeTimesSynch = pairGroupStatTable.refSynch{loopPairs};
        refCellType        = pairGroupStatTable.refCellExplorerType{loopPairs};
        [~,refSpikeTimesSynchIdx,~] = intersect(refSpikeTimes,refSpikeTimesSynch);

        tarShank           = pairGroupStatTable.tarShank(loopPairs);
        tarChan            = pairGroupStatTable.tarChannel(loopPairs);
        tarID              = pairGroupStatTable.tarNeuronID(loopPairs);
        tarSpikeTimes      = pairGroupStatTable.tarSpikeTimes{loopPairs};
        tarSpikeTimesSynch = pairGroupStatTable.tarSynch{loopPairs};
        tarCellType        = pairGroupStatTable.tarCellExplorerType{loopPairs};
        [~,tarSpikeTimesSynchIdx,~] = intersect(tarSpikeTimes,tarSpikeTimesSynch);

        d                  = pairGroupStatTable.pairDistance(loopPairs);

        tWave = pairGroupStatTable.tWave{loopPairs};

        %% 

        titleStrRef = ['rat ' animalNames{loopRat} ': '...
                        num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')' ' - ' ...
                        num2str(refID) refCellType ' (sh: ' num2str(refShank) ')' ' pair ' num2str(d) '\mum '];

        titleStrTar = ['rat ' animalNames{loopRat} ': '...
                        num2str(refID) refCellType ' (sh: ' num2str(refShank) ')' ' - ' ...
                        num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')' ' pair ' num2str(d) '\mum '];
        
        titleStrRefOnly = ['rat ' animalNames{loopRat} ': ' num2str(refID) refCellType ' (sh: ' num2str(refShank) ')'];
        titleStrTarOnly = ['rat ' animalNames{loopRat} ': ' num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')'];

        saveStr  = ['rat ' animalNames{loopRat} ': '...
                    num2str(refID) refCellType ' (sh ' num2str(refShank) ')' ' - ' ...
                    num2str(tarID) tarCellType ' (sh ' num2str(tarShank) ')' ' pair ' num2str(d) ' microns.jpeg'];

        % [~,tarChan] = max(max(abs(unitsGlobalSnapshots.waveforms{tarID,1}(channelSet{tarShank},:)),[],2));

        if loopRat == 1
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'    num2str(refShank) '.mat']) 
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'    num2str(tarShank) '.mat']) 
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find((neuronsN.neuron_ids + 1) == tarID);
            loopIdxRef = find((neuronsN.neuron_ids + 1) == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 2
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'    num2str(refShank) '.mat']);
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'    num2str(tarShank) '.mat']);
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find((neuronsS.neuron_ids + 1) == tarID);
            loopIdxRef = find((neuronsS.neuron_ids + 1) == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 3
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'    num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'    num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find((neuronsU.neuron_ids + 1) == tarID);
            loopIdxRef = find((neuronsU.neuron_ids + 1) == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 4
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find(cell_metrics.cell_metrics.UID == tarID);
            loopIdxRef = find(cell_metrics.cell_metrics.UID == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 5
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find(cell_metrics.cell_metrics.UID == tarID);
            loopIdxRef = find(cell_metrics.cell_metrics.UID == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 6
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find(cellIDs == tarID);
            loopIdxRef = find(cellIDs == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        end

        %%

        % PCAchans = clusterValNeuron.templateChannels{clusterValNeuron.cellID == tarID};
        % PCAchansIdx = ~isnan(PCAchans);

        if loopRat == 6
            cluIdxRef  = 1;
            cluIdxTar  = 1;
        else
            cluIdxRef  = clusterValNeuronsRef.clusIdx(find(clusterValNeuronsRef.cellID == refID));
            cluIdxTar  = clusterValNeuronsTar.clusIdx(find(clusterValNeuronsTar.cellID == tarID));
        end

        PCAvalsRef = clusterValNeuronsPCAvalsRef;
        PCAvalsTar = clusterValNeuronsPCAvalsTar;
        % PCAvals = {PCAvals{1,PCAchansIdx}};

        PCAref      = PCAvalsRef{cluIdxRef};
        PCAtar      = PCAvalsTar{cluIdxTar};
        PCAsynchRef = PCAref(refSpikeTimesSynchIdx(refSpikeTimesSynchIdx < size(PCAref,1)),:); % needed cutoff
        PCAsynchTar = PCAtar(tarSpikeTimesSynchIdx(tarSpikeTimesSynchIdx < size(PCAtar,1)),:); % needed cutoff
        
        %%
        
        [refSortedPC1,refSortedPC1idx] = sort(PCAref(:,1),'ascend');
        [refSortedPC2,refSortedPC2idx] = sort(PCAref(:,2),'ascend');
        [refSortedPC3,refSortedPC3idx] = sort(PCAref(:,3),'ascend');

        refSortedPC1idxMinSet = refSortedPC1idx(1:NsortSet);
        refSortedPC1idxMaxSet = refSortedPC1idx(end-NsortSet:end);
        
        refSortedPC2idxMinSet = refSortedPC2idx(1:NsortSet);
        refSortedPC2idxMaxSet = refSortedPC2idx(end-NsortSet:end);

        refSortedPC3idxMinSet = refSortedPC3idx(1:NsortSet);
        refSortedPC3idxMaxSet = refSortedPC3idx(end-NsortSet:end);
        
        [tarSortedPC1,tarSortedPC1idx] = sort(PCAtar(:,1),'ascend');
        [tarSortedPC2,tarSortedPC2idx] = sort(PCAtar(:,2),'ascend');
        [tarSortedPC3,tarSortedPC3idx] = sort(PCAtar(:,3),'ascend');

        tarSortedPC1idxMinSet = tarSortedPC1idx(1:NsortSet);
        tarSortedPC1idxMaxSet = tarSortedPC1idx(end-NsortSet:end);
        
        tarSortedPC2idxMinSet = tarSortedPC2idx(1:NsortSet);
        tarSortedPC2idxMaxSet = tarSortedPC2idx(end-NsortSet:end);

        tarSortedPC3idxMinSet = tarSortedPC3idx(1:NsortSet);
        tarSortedPC3idxMaxSet = tarSortedPC3idx(end-NsortSet:end);

        %% getting snippets

        % timeShift      = 514;
        % chanData       = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarCh) 'highpass300hz.mat']);

        if loopRat == 1
            chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(tarChan) 'highpass300hz.mat']);
            chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(refChan) 'highpass300hz.mat']);
        elseif loopRat == 2
            chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(tarChan) 'highpass300hz.mat']);
            chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(refChan) 'highpass300hz.mat']);
        elseif loopRat == 3
            chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(tarChan) 'highpass300hz.mat']);
            chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(refChan) 'highpass300hz.mat']);
        elseif loopRat == 4
            chanDataRef = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(refShank) '/ch' num2str(refChan) 'highpass300hz.mat']);
            chanDataTar = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(tarShank) '/ch' num2str(tarChan) 'highpass300hz.mat']);
        elseif loopRat == 5
            chanDataRef = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(refShank) '/ch' num2str(refChan) 'highpass300hz.mat']);
            chanDataTar = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(tarShank) '/ch' num2str(tarChan) 'highpass300hz.mat']);
        elseif loopRat == 6
            chanDataRef = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(refChan) 'highpass300hz.mat']);
            chanDataTar = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarChan) 'highpass300hz.mat']);
        end

        if loopRat == 6
            refLat       = round((refSpikeTimes       - chanDataRef.data.onsetTime/1e6)*fs) + 1;
            refSynchLat  = round((refSpikeTimesSynch  - chanDataRef.data.onsetTime/1e6)*fs) + 1;

            tarLat       = round((tarSpikeTimes       - chanDataTar.data.onsetTime/1e6)*fs) + 1;
            tarSynchLat  = round((tarSpikeTimesSynch  - chanDataTar.data.onsetTime/1e6)*fs) + 1;
        else
            refLat       = round(refSpikeTimes*fs)';
            refSynchLat  = round(refSpikeTimesSynch*fs)';

            tarLat       = round(tarSpikeTimes*fs)';
            tarSynchLat  = round(tarSpikeTimesSynch*fs)';
        end

        if length(refLat) > NspikesPlot
            refLatRandPerm = sort(refLat(randperm(length(refLat),NspikesPlot)));
        else
            refLatRandPerm = refLat;
        end

        if length(refSynchLat) > NspikesPlot
            refSynchLatRandPerm = sort(refSynchLat(randperm(length(refSynchLat),NspikesPlot)));
        else
            refSynchLatRandPerm = refSynchLat;
        end

        if loopRat == 6
            [refSponWaveMean,      refSponWaveforms]      = waveformAvg(chanDataRef.data.channel,refLatRandPerm,     preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refSynchLatWaveMean,  refPeakLatWaveforms]   = waveformAvg(chanDataRef.data.channel,refSynchLatRandPerm,preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [refPC1minSetWaveMean, refPC1minSetWaveforms] = waveformAvg(chanDataRef.data.channel,refLat(refSortedPC1idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refPC1maxSetWaveMean, refPC1maxSetWaveforms] = waveformAvg(chanDataRef.data.channel,refLat(refSortedPC1idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [refPC2minSetWaveMean, refPC2minSetWaveforms] = waveformAvg(chanDataRef.data.channel,refLat(refSortedPC2idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refPC2maxSetWaveMean, refPC2maxSetWaveforms] = waveformAvg(chanDataRef.data.channel,refLat(refSortedPC2idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [refPC3minSetWaveMean, refPC3minSetWaveforms] = waveformAvg(chanDataRef.data.channel,refLat(refSortedPC3idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refPC3maxSetWaveMean, refPC3maxSetWaveforms] = waveformAvg(chanDataRef.data.channel,refLat(refSortedPC3idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
        else
            [refSponWaveMean,      refSponWaveforms]      = waveformAvg(chanDataRef.data,       refLatRandPerm,     preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refSynchLatWaveMean,  refPeakLatWaveforms]   = waveformAvg(chanDataRef.data,       refSynchLatRandPerm,preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [refPC1minSetWaveMean, refPC1minSetWaveforms] = waveformAvg(chanDataRef.data,refLat(refSortedPC1idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refPC1maxSetWaveMean, refPC1maxSetWaveforms] = waveformAvg(chanDataRef.data,refLat(refSortedPC1idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [refPC2minSetWaveMean, refPC2minSetWaveforms] = waveformAvg(chanDataRef.data,refLat(refSortedPC2idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refPC2maxSetWaveMean, refPC2maxSetWaveforms] = waveformAvg(chanDataRef.data,refLat(refSortedPC2idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [refPC3minSetWaveMean, refPC3minSetWaveforms] = waveformAvg(chanDataRef.data,refLat(refSortedPC3idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [refPC3maxSetWaveMean, refPC3maxSetWaveforms] = waveformAvg(chanDataRef.data,refLat(refSortedPC3idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
        end

        if length(tarLat) > NspikesPlot
            tarLatRandPerm = sort(tarLat(randperm(length(tarLat),NspikesPlot)));
        else
            tarLatRandPerm = tarLat;
        end

        if length(tarSynchLat) > NspikesPlot
            tarSynchLatRandPerm = sort(tarSynchLat(randperm(length(tarSynchLat),NspikesPlot)));
        else
            tarSynchLatRandPerm = tarSynchLat;
        end

        if loopRat == 6
            [tarSponWaveMean,      tarSponWaveforms]     = waveformAvg(chanDataTar.data.channel, tarLatRandPerm,preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarSynchLatWaveMean,  tarPeakLatWaveforms]  = waveformAvg(chanDataTar.data.channel, tarSynchLatRandPerm,    preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [tarPC1minSetWaveMean, tarPC1minSetWaveforms] = waveformAvg(chanDataTar.data.channel,tarLat(tarSortedPC1idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarPC1maxSetWaveMean, tarPC1maxSetWaveforms] = waveformAvg(chanDataTar.data.channel,tarLat(tarSortedPC1idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [tarPC2minSetWaveMean, tarPC2minSetWaveforms] = waveformAvg(chanDataTar.data.channel,tarLat(tarSortedPC2idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarPC2maxSetWaveMean, tarPC2maxSetWaveforms] = waveformAvg(chanDataTar.data.channel,tarLat(tarSortedPC2idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [tarPC3minSetWaveMean, tarPC3minSetWaveforms] = waveformAvg(chanDataTar.data.channel,tarLat(tarSortedPC3idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarPC3maxSetWaveMean, tarPC3maxSetWaveforms] = waveformAvg(chanDataTar.data.channel,tarLat(tarSortedPC3idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
        else
            [tarSponWaveMean,      tarSponWaveforms]     = waveformAvg(chanDataTar.data,         tarLatRandPerm,preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarSynchLatWaveMean,  tarPeakLatWaveforms]  = waveformAvg(chanDataTar.data,         tarSynchLatRandPerm,    preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [tarPC1minSetWaveMean, tarPC1minSetWaveforms] = waveformAvg(chanDataTar.data,tarLat(tarSortedPC1idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarPC1maxSetWaveMean, tarPC1maxSetWaveforms] = waveformAvg(chanDataTar.data,tarLat(tarSortedPC1idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [tarPC2minSetWaveMean, tarPC2minSetWaveforms] = waveformAvg(chanDataTar.data,tarLat(tarSortedPC2idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarPC2maxSetWaveMean, tarPC2maxSetWaveforms] = waveformAvg(chanDataTar.data,tarLat(tarSortedPC2idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);

            [tarPC3minSetWaveMean, tarPC3minSetWaveforms] = waveformAvg(chanDataTar.data,tarLat(tarSortedPC3idxMinSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
            [tarPC3maxSetWaveMean, tarPC3maxSetWaveforms] = waveformAvg(chanDataTar.data,tarLat(tarSortedPC3idxMaxSet),preLength,postLength,fpass,fs,filtarFlag,noDemeanFlag);
        end
        % tarRandPermIdx = sort(randperm(length(tarSpikeTimes),10000));

        % tarSponWaveforms  = snippetsData.waveforms(:,tarRandPermIdx);
        % tarSynchWaveforms = snippetsData.waveforms(:,tarSpikeTimesSynchIdx);

        %% ref PC1 split
        
        figure(102)
        tiledlayout(3,5)

        nexttile(1)
        scatter(PCAref(:,1),                 PCAref(:,2),                 [],[.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC1idxMinSet,1),PCAref(refSortedPC1idxMinSet,2),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC1idxMaxSet,1),PCAref(refSortedPC1idxMaxSet,2),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        legend('all spikes','min PC1','max PC1','Location','north')
        title(titleStrRefOnly)
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(2)
        scatter(PCAref(:,1),                 PCAref(:,3),                 [], [.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC1idxMinSet,1),PCAref(refSortedPC1idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC1idxMaxSet,1),PCAref(refSortedPC1idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(3)
        scatter(PCAref(:,2),                 PCAref(:,3),                 [], [.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC1idxMinSet,2),PCAref(refSortedPC1idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC1idxMaxSet,2),PCAref(refSortedPC1idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %% ref PC2 split

        nexttile(6)
        scatter(PCAref(:,1),                 PCAref(:,2),                 [],[.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC2idxMinSet,1),PCAref(refSortedPC2idxMinSet,2),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC2idxMaxSet,1),PCAref(refSortedPC2idxMaxSet,2),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        legend('all spikes','min PC2','max PC2','Location','north')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(7)
        scatter(PCAref(:,1),                 PCAref(:,3),                 [], [.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC2idxMinSet,1),PCAref(refSortedPC2idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC2idxMaxSet,1),PCAref(refSortedPC2idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(8)
        scatter(PCAref(:,2),                 PCAref(:,3),                 [], [.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC2idxMinSet,2),PCAref(refSortedPC2idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC2idxMaxSet,2),PCAref(refSortedPC2idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %% ref PC3 split

        nexttile(11)
        scatter(PCAref(:,1),                 PCAref(:,2),                 [],[.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC3idxMinSet,1),PCAref(refSortedPC3idxMinSet,2),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC3idxMaxSet,1),PCAref(refSortedPC3idxMaxSet,2),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        legend('all spikes','min PC3','max PC3','Location','north')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(12)
        scatter(PCAref(:,1),                 PCAref(:,3),                 [], [.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC3idxMinSet,1),PCAref(refSortedPC3idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC3idxMaxSet,1),PCAref(refSortedPC3idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(13)
        scatter(PCAref(:,2),                 PCAref(:,3),                 [], [.7 .7 .7])
        hold on
        scatter(PCAref(refSortedPC3idxMinSet,2),PCAref(refSortedPC3idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAref(refSortedPC3idxMaxSet,2),PCAref(refSortedPC3idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %% ref snippet plots 
        
        nexttile(4)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         refSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(refSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(refPC1minSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPC1minSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0 0.4470 0.7410])
        end
        for i = 1:size(refPC1maxSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPC1maxSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0.8500 0.3250 0.0980])
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(5)
        plot((-preLength:postLength-1)*(30/1000),refSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),refPC1minSetWaveMean,'Color',[0 0.4470 0.7410],'LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),refPC1maxSetWaveMean,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(9)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         refSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(refSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(refPC2minSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPC2minSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0 0.4470 0.7410])
        end
        for i = 1:size(refPC2maxSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPC2maxSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0.8500 0.3250 0.0980])
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(10)
        plot((-preLength:postLength-1)*(30/1000),refSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),refPC2minSetWaveMean,'Color',[0 0.4470 0.7410],'LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),refPC2maxSetWaveMean,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(14)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         refSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(refSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(refPC3minSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPC3minSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0 0.4470 0.7410])
        end
        for i = 1:size(refPC3maxSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPC3maxSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0.8500 0.3250 0.0980])
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(15)
        plot((-preLength:postLength-1)*(30/1000),refSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),refPC3minSetWaveMean,'Color',[0 0.4470 0.7410],'LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),refPC3maxSetWaveMean,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        save_file = fullfile(figPath, [titleStrRefOnly + " PC split"]);
        print(fig_use, save_file,'-djpeg',resolution_use);

        %% tar PC1 split
        
        figure(102)
        tiledlayout(3,5)

        nexttile(1)
        scatter(PCAtar(:,1),                    PCAtar(:,2),                    [],[.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC1idxMinSet,1),PCAtar(tarSortedPC1idxMinSet,2),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC1idxMaxSet,1),PCAtar(tarSortedPC1idxMaxSet,2),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        legend('all spikes','min PC1','max PC1','Location','north')
        title(titleStrTarOnly)
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(2)
        scatter(PCAtar(:,1),                    PCAtar(:,3),                    [], [.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC1idxMinSet,1),PCAtar(tarSortedPC1idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC1idxMaxSet,1),PCAtar(tarSortedPC1idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(3)
        scatter(PCAtar(:,2),                    PCAtar(:,3),                    [], [.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC1idxMinSet,2),PCAtar(tarSortedPC1idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC1idxMaxSet,2),PCAtar(tarSortedPC1idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %% tar PC2 split

        nexttile(6)
        scatter(PCAtar(:,1),                    PCAtar(:,2),                    [],[.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC2idxMinSet,1),PCAtar(tarSortedPC2idxMinSet,2),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC2idxMaxSet,1),PCAtar(tarSortedPC2idxMaxSet,2),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        legend('all spikes','min PC2','max PC2','Location','north')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(7)
        scatter(PCAtar(:,1),                    PCAtar(:,3),                    [], [.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC2idxMinSet,1),PCAtar(tarSortedPC2idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC2idxMaxSet,1),PCAtar(tarSortedPC2idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(8)
        scatter(PCAtar(:,2),                    PCAtar(:,3),                    [], [.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC2idxMinSet,2),PCAtar(tarSortedPC2idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC2idxMaxSet,2),PCAtar(tarSortedPC2idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %% tar PC3 split

        nexttile(11)
        scatter(PCAtar(:,1),                    PCAtar(:,2),                    [],[.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC3idxMinSet,1),PCAtar(tarSortedPC3idxMinSet,2),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC3idxMaxSet,1),PCAtar(tarSortedPC3idxMaxSet,2),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        legend('all spikes','min PC3','max PC3','Location','north')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(12)
        scatter(PCAtar(:,1),                    PCAtar(:,3),                    [], [.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC3idxMinSet,1),PCAtar(tarSortedPC3idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC3idxMaxSet,1),PCAtar(tarSortedPC3idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(13)
        scatter(PCAtar(:,2),                    PCAtar(:,3),                    [], [.7 .7 .7])
        hold on
        scatter(PCAtar(tarSortedPC3idxMinSet,2),PCAtar(tarSortedPC3idxMinSet,3),[],[0 0.4470 0.7410])
        scatter(PCAtar(tarSortedPC3idxMaxSet,2),PCAtar(tarSortedPC3idxMaxSet,3),[],[0.8500 0.3250 0.0980])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %% tar snippet plots 
        
        nexttile(4)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         tarSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(tarSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(tarPC1minSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPC1minSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0 0.4470 0.7410])
        end
        for i = 1:size(tarPC1maxSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPC1maxSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0.8500 0.3250 0.0980])
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(5)
        plot((-preLength:postLength-1)*(30/1000),tarSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),tarPC1minSetWaveMean,'Color',[0 0.4470 0.7410],'LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),tarPC1maxSetWaveMean,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold off
        
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(9)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         tarSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(tarSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(tarPC2minSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPC2minSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0 0.4470 0.7410])
        end
        for i = 1:size(tarPC2maxSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPC2maxSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0.8500 0.3250 0.0980])
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(10)
        plot((-preLength:postLength-1)*(30/1000),tarSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),tarPC2minSetWaveMean,'Color',[0 0.4470 0.7410],'LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),tarPC2maxSetWaveMean,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(14)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         tarSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(tarSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(tarPC3minSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPC3minSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0 0.4470 0.7410])
        end
        for i = 1:size(tarPC3maxSetWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPC3maxSetWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor',[0.8500 0.3250 0.0980])
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        nexttile(15)
        plot((-preLength:postLength-1)*(30/1000),tarSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),tarPC3minSetWaveMean,'Color',[0 0.4470 0.7410],'LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),tarPC3maxSetWaveMean,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

        save_file = fullfile(figPath, [titleStrTarOnly + " PC split"]);
        print(fig_use, save_file,'-djpeg',resolution_use);


        %%    

%         synchPlotColor = 'r';
% 
%         figure(102)
%         tiledlayout(4,3)
% 
%         tR     = pairGroupStatTable.CCGbinLagTimes{loopPairs}*1000;
%         ccgR   = pairGroupStatTable.pairRawCCG{loopPairs};
%         GSPExc = pairGroupStatTable.GSPExc{loopPairs};
% 
%         nexttile(1)
%         plot(tR,ccgR,'k','LineWidth',1)
%         hold on
%         scatter(tR(find(GSPExc)),ccgR(find(GSPExc)),100,'ro','LineWidth',1)
%         hold off
% 
%         xlim([-1,1])
%         ylims = get(gca,'ylim');
%         ylabel('Spike Probability')
%         xlabel('[ms]')
%         title(titleStrTar)
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
%         box off
% 
%         nexttile(7)
%         plot(tR,flipud(ccgR),'k','LineWidth',1)
%         hold on
%         scatter(tR(find(flipud(GSPExc))),flipud(ccgR(find(GSPExc))),100,'ro','LineWidth',1)
%         hold off
% 
%         xlim([-1,1])
%         ylims = get(gca,'ylim');
%         ylabel('Spike Probability')
%         xlabel('[ms]')
%         title(titleStrRef)
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
%         box off
% 
%         %%
% 
%         nexttile(2)
%         patchline((-preLength:postLength-1)*(30/1000), ... 
%                          tarSponWaveforms(:,1), ...
%                          'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
%         hold on
%         for i = 2:size(tarSponWaveforms,2)
%             patchline((-preLength:postLength-1)*(30/1000), ... 
%                      tarSponWaveforms(:,i), ...
%                      'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
%         end
%         for i = 1:size(tarPeakLatWaveforms,2)
%             patchline((-preLength:postLength-1)*(30/1000), ... 
%                      tarPeakLatWaveforms(:,i), ...
%                      'linewidth',1,'edgealpha',0.1,'edgecolor','r')
%         end
%         hold off
%         % ylim([-6 2])
%         xlim([-1,1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
%         xlabel('[ms]');
%         ylabel('[mV]');
%         set(gca, 'YDir','reverse')
%         box off
% 
%         nexttile(3)
%         plot((-preLength:postLength-1)*(30/1000),tarSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
%         hold on
%         plot((-preLength:postLength-1)*(30/1000),tarSynchLatWaveMean,'r','LineWidth',1)
%         hold off
% 
%         xlim([-1,1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
%         xlabel('[ms]');
%         ylabel('[mV]');
%         set(gca, 'YDir','reverse')
%         box off
% 
%         nexttile(8)
%         patchline((-preLength:postLength-1)*(30/1000), ... 
%                          refSponWaveforms(:,1), ...
%                          'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
%         hold on
%         for i = 2:size(refSponWaveforms,2)
%             patchline((-preLength:postLength-1)*(30/1000), ... 
%                      refSponWaveforms(:,i), ...
%                      'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
%         end
%         for i = 1:size(refPeakLatWaveforms,2)
%             patchline((-preLength:postLength-1)*(30/1000), ... 
%                      refPeakLatWaveforms(:,i), ...
%                      'linewidth',1,'edgealpha',0.1,'edgecolor','r')
%         end
%         hold off
%         % ylim([-6 2])
%         xlim([-1,1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
%         xlabel('[ms]');
%         ylabel('[mV]');
%         set(gca, 'YDir','reverse')
%         box off
% 
%         nexttile(9)
%         plot((-preLength:postLength-1)*(30/1000),refSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
%         hold on
%         plot((-preLength:postLength-1)*(30/1000),refSynchLatWaveMean,'r','LineWidth',1)
%         hold off
% 
%         xlim([-1,1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
%         xlabel('[ms]');
%         ylabel('[mV]');
%         set(gca, 'YDir','reverse')
%         box off
% 
%         %%
% 
%         nexttile(4)
%         scatter(PCAtar(:,1), PCAtar(:,2), [], [.7 .7 .7])
%         hold on
%         scatter(PCAsynchTar(:,1), PCAsynchTar(:,2), [], synchPlotColor)
%         hold off
%         xlabel('PC1')
%         ylabel('PC2')
%         % xlim([-2500 2500])
%         % ylim([-2500 2500])
%         daspect([1 1 1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
% %             legend('noSynch','synch','location','southoutside')
% 
%         nexttile(5)
%         scatter(PCAtar(:,2), PCAtar(:,3), [], [.7 .7 .7])
%         hold on
%         scatter(PCAsynchTar(:,2), PCAsynchTar(:,3), [], synchPlotColor)
%         hold off
%         xlabel('PC1')
%         ylabel('PC3')
%         % xlim([-2500 2500])
%         % ylim([-2500 2500])
%         daspect([1 1 1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
% %             legend('noSynch','synch','location','southoutside')
% %             title(["ch: " num2str(tarCh)])
% 
%         nexttile(6)
%         scatter(PCAtar(:,1), PCAtar(:,3), [], [.7 .7 .7])
%         hold on
%         scatter(PCAsynchTar(:,1), PCAsynchTar(:,3), [], synchPlotColor)
%         hold off
%         xlabel('PC2')
%         ylabel('PC3')
%         % xlim([-2500 2500])
%         % ylim([-2500 2500])
%         daspect([1 1 1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
% %             legend('noSynch','synch','location','southoutside')
% 
%         nexttile(10)
%         scatter(PCAref(:,1), PCAref(:,2), [], [.7 .7 .7])
%         hold on
%         scatter(PCAsynchRef(:,1), PCAsynchRef(:,2), [], synchPlotColor)
%         hold off
%         xlabel('PC1')
%         ylabel('PC2')
%         % xlim([-2500 2500])
%         % ylim([-2500 2500])
%         daspect([1 1 1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
% %             legend('noSynch','synch','location','southoutside')
% 
%         nexttile(11)
%         scatter(PCAref(:,2), PCAref(:,3), [], [.7 .7 .7])
%         hold on
%         scatter(PCAsynchRef(:,2), PCAsynchRef(:,3), [], synchPlotColor)
%         hold off
%         xlabel('PC1')
%         ylabel('PC3')
%         % xlim([-2500 2500])
%         % ylim([-2500 2500])
%         daspect([1 1 1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
% %             legend('noSynch','synch','location','southoutside')
% %             title(["ch: " num2str(tarCh)])
% 
%         nexttile(12)
%         scatter(PCAref(:,1), PCAref(:,3), [], [.7 .7 .7])
%         hold on
%         scatter(PCAsynchRef(:,1), PCAsynchRef(:,3), [], synchPlotColor)
%         hold off
%         xlabel('PC2')
%         ylabel('PC3')
%         % xlim([-2500 2500])
%         % ylim([-2500 2500])
%         daspect([1 1 1])
%         set(gca,'FontSize',5)
%         set(gca,'FontName','Arial')
% %             legend('noSynch','synch','location','southoutside')

        % save_file = fullfile(figPath, saveStr);
        % print(fig_use, save_file,'-djpeg',resolution_use);


        %%


        % PCA = [];
        % for loopPCAchans = sum(PCAchansIdx)
        %     PCA = [PCA double(PCAvals{1,loopPCAchans})];
        % end
        % 
        % COM = mean(PCA);
        % 
        % for loopSpikeTimes = 1:size(PCA,1)
        %     D(loopSpikeTimes)   = sqrt(sum((COM - PCA(loopSpikeTimes,:)) .^ 2));
        % end
        % 
        % %%
        % [~,~,idx] = intersect(pairGroupStatTable.tarSynch{loopPairs},pairGroupStatTable.tarSpikeTimes{loopPairs});
        % for loopSpikeTimes = 1:size(idx,1)
        %     Dsynch(loopSpikeTimes)   = sqrt(sum((COM - PCA(idx(loopSpikeTimes),:)) .^ 2));
        % end
        % 
        % histogram(D,     'Normalization','probability','BinWidth',1)
        % hold on
        % histogram(Dsynch,'Normalization','probability','BinWidth',1)
        % hold off

    end

end

%% exquisite graphs

load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable
G = graph( string(pairGroupStatsTable.refNeuronID), string(pairGroupStatsTable.tarNeuronID)); 
plot(G)

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/groupStatsRatN')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable
G = graph( string(pairGroupStatsTable.refNeuronID), string(pairGroupStatsTable.tarNeuronID)); 
plot(G)

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/groupStatsRatS')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable
G = graph( string(pairGroupStatsTable.refNeuronID), string(pairGroupStatsTable.tarNeuronID)); 
plot(G)

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/groupStatsRatU')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable
G = graph( string(pairGroupStatsTable.refNeuronID), string(pairGroupStatsTable.tarNeuronID)); 
plot(G)
    
load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/groupStatsRatAGday1.mat')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable
G = graph( string(pairGroupStatsTable.refNeuronID), string(pairGroupStatsTable.tarNeuronID)); 
plot(G)

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/groupStatsRatAGday2.mat')
pairGroupStatsTable = pairGroupStatTable;
clear pairGroupStatTable
G = graph( string(pairGroupStatsTable.refNeuronID), string(pairGroupStatsTable.tarNeuronID)); 
plot(G)


