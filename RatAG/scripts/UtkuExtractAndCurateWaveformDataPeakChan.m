function UtkuExtractAndCurateWaveformDataPeakChan
    
    flagAdjusted = true;

    fs         = 30000;
    fpass      = 300;
    
    preLength  = 36;
    postLength = 36;

    channelShanks = [ones(32,1); 2*ones(32,1); 3*ones(32,1); 4*ones(32,1); 5*ones(32,1); 6*ones(32,1);];

    for loopDays = 2:2

        disp(['looping day: ' num2str(loopDays)])

        if loopDays == 1
            chanDataDir  = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/'];
            if flagAdjusted
                load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/AG_2019-12_23_NSD_adjustedSpikeTimes.mat');
            end
        elseif loopDays == 2
            chanDataDir  = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/'];
             if flagAdjusted
                load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/AG_2019-12-27_NSD_adjustedSpikeTimes.mat');
            end
        end

        snippetsDir = [chanDataDir 'snippetsPeakChan/'];

        neurons = loadUtkuV2(num2str(loopDays));
             
        % loop units
        for loopUnits = 1:length(neurons.cell_metrics.UID)
            
            neuronID = neurons.cell_metrics.UID(loopUnits);
            
            if flagAdjusted            
                if isempty(find(adjustedSpikeTimes.cellID == neuronID))
                    continue 
                else 
                    spikeTimes = adjustedSpikeTimes.spikeTimes{find(adjustedSpikeTimes.cellID == neuronID),1};
                end
            else
                spikeTimes =  neurons.cell_metrics.spikes.times{1,loopUnits};
            end
                
            spikeTimeIndxCell = round(spikeTimes*fs) + 1;

            peakChan = neurons.cell_metrics.maxWaveformCh1(loopUnits);

            chanData = load([ chanDataDir 'shank' num2str(channelShanks(peakChan)) '/ch' num2str(peakChan) 'highpass300hz.mat'],'data');

            [waveAvg,waveforms] = waveformAvg(double(chanData.data), ... 
                                                         spikeTimeIndxCell, ... 
                                                         preLength,postLength,fpass,fs,false);

            snippetsData.waveforms   = waveforms;
            snippetsData.waveformAvg = waveAvg;
            
            if flagAdjusted
                saveName = ['unit' num2str(neuronID) 'peakCh' num2str(peakChan) 'adjustedSnippets.mat' ];
            else
                saveName = ['unit' num2str(neuronID) 'peakCh' num2str(peakChan) 'snippets.mat' ];
            end
            save([snippetsDir saveName],'snippetsData')

        end
    end
end
 