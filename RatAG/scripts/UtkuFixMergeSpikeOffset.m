function UtkuFixMergeSpikeOffset(dayNo)
    
    fs = 3e4;
    t  = (-35:36)*(1000/fs);

    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    if strcmp(dayNo,'1')
        data_dir    = [UtkuPath 'AG_2019-12-23_NSD' '/'];
        snippetsDir = '/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/snippetsPeakChan/';
        saveName    = 'AG_2019-12_23_NSD_adjustedSpikeTimes.mat';
    elseif strcmp(dayNo,'2')
        data_dir    = [UtkuPath 'AG_2019-12-27_NSD' '/'];
        snippetsDir = '/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/snippetsPeakChan/';
        saveName    = 'AG_2019-12-27_NSD_adjustedSpikeTimes.mat';
    elseif strcmp(dayNo,'B1')
        data_dir = [UtkuPath 'B1_2022-06-24_NSD_CA1_24hrs/2022-06-24_NSD_CA1_24hrs/2022-06-24_NSD_CA1_24hrscrs-merged_cleaned.GUI' '/'];
    end    

    UtkuData = loadUtkuV2(dayNo);
    
    adjustedSpikeTimes = table; 
    
    for loopUnits = 1:length(UtkuData.cell_metrics.UID) 
        
        cellID     = UtkuData.cell_metrics.UID(loopUnits);
        peakChan   = UtkuData.cell_metrics.maxWaveformCh1(loopUnits);
        spikeTimes = UtkuData.cell_metrics.spikes.times{loopUnits};
        
        if exist([snippetsDir 'unit' num2str(cellID) 'peakCh' num2str(peakChan) 'snippets.mat'],'file') > 0
            load([snippetsDir 'unit' num2str(cellID) 'peakCh' num2str(peakChan) 'snippets.mat']);
            
            %% determine positive or negative spike
            if sum(abs(snippetsData.waveformAvg(34:38))) > 0
                positiveWav = true; 
            end
            
            [~,idx]         = max(abs(snippetsData.waveforms));
            modeWavePeakIdx = mode(idx);
            
            offset             = idx - modeWavePeakIdx;
            
            adjustedSpikeTimes.cellID(loopUnits)     = cellID;
            adjustedSpikeTimes.spikeTimes{loopUnits} = spikeTimes - offset'*(1e3/fs);
            
        end
        
    end
    
    save([data_dir saveName],'adjustedSpikeTimes')
    
end

