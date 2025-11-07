function [waveAvgAll,waveforms,chans] = extractWaveformsWholeShank(spikeTimes,shank,channelSet,cellID,ratID,fs,preLength,postLength)

    loopShanks = channelSet{shank};
    
    for i = 1:length(loopShanks)
        
        chanData = load(['/media/nasko/WD_BLACK3/Rat' ratID '_temp/rawDataCh' num2str(loopShanks(i)) '.mat'],'data');
        spikeTimeIndx = round(spikeTimes*fs) + 1;
        [waveAvg,waves] = waveformAvg(double(chanData.data),spikeTimeIndx,preLength,postLength,300,fs,false);
        
        waveAvgAll(i,:) = waveAvg;
        waveforms{i}    = waves;
        
    end
    
    chans        = loopShanks;
    
end