function steinmetzLoadNPYs

    savePath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';
    
    fs = 30000;
    
    spikeTimesAll = readNPY('spikes.times.npy');
    clusters      = readNPY('spikes.clusters.npy');
    waveforms     = readNPY('clusters.waveforms.npy');

    % readNPY('clusters.CCF_APDVLR.npy');
    % readNPY('clusters.acronym.npy');

    brainArea = steinmetzBrainArea;
    
    brainAreaList = unique(brainArea);

%     neuronIdx = find(strcmp(brainArea,'CA1') | strcmp(brainArea,'CA3'));
%     neuronIdx = find(strcmp(brainArea,'VPL'));
%     neuronIdx = find(strcmp(brainArea,'VPM'));
%     neuronIdx = find(strcmp(brainArea,'MOs6a'));
    
    for j = 1:length(brainAreaList)
        neuronIdx = find(strcmp(brainArea,brainAreaList{j}));
        
        display([num2str((j/length(brainAreaList))*100) '% done'])
        
        for i = 1:length(neuronIdx)
            
            [~,peakChan] = max(max(abs(squeeze(waveforms(neuronIdx(i),:,:)))));
            avgWaveforms = squeeze(waveforms(neuronIdx(i),:,:));
            [~,troughTime] = min(avgWaveforms(:,peakChan));
            [~,troughToPeakLength] = max(avgWaveforms(troughTime:end,peakChan));
            troughToPeakLength = (troughToPeakLength/fs)*1000; % convert from samples to ms
            spikeTimes = spikeTimesAll(neuronIdx(i) == clusters);

            neurons.putativeCellType{i}    = cellClass(spikeTimes,troughToPeakLength,fs);
            neurons.troughToPeakLength(i)  = troughToPeakLength;
            neurons.unitID(i)              = neuronIdx(i);
            neurons.peakChan(i)            = peakChan;
            neurons.probeID{i}             = 1;
            neurons.brainRegion{i}         = brainArea{neuronIdx(i)};
            neurons.spikeTimes{i}          = spikeTimes;
            neurons.avgWaveforms{i}        = avgWaveforms;
            
        end

        neurons.probeX = readNPY('channels.xcoords.npy');
        neurons.probeY = readNPY('channels.ycoords.npy');
        
        if strcmp(brainAreaList{j},'MOs2/3')
            brainRegionSaveStr = 'MOs2and3';
        elseif strcmp(brainAreaList{j},'ORBl2/3')
            brainRegionSaveStr = 'ORBl2and3';
        elseif strcmp(brainAreaList{j},'VISam2/3')
            brainRegionSaveStr = 'VISam2and3';
        elseif strcmp(brainAreaList{j},'VISpm2/3')
            brainRegionSaveStr = 'VISpm2and3';
        else
            brainRegionSaveStr = brainAreaList{j};
        end
        
        saveName = ['neuronsSteinmetzDataArea' brainRegionSaveStr '.mat'];
        save([savePath saveName],'neurons')
        
        clear neurons
    end

    % scatter(x,y)
    % axis equal
    
end