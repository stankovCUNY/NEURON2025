function UtkuExtractAndCurateWaveformData
    
%     parpool(2)

    flagAdjusted = true;

    % how many random spikes
    N          = 1000; 

    fs         = 30000;
    fpass      = 300;
    
    preLength  = 36;
    postLength = 36;

    channelShanks = [ones(32,1); 2*ones(32,1); 3*ones(32,1); 4*ones(32,1); 5*ones(32,1); 6*ones(32,1);];
    
    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';   
    
    for loopDays = 1:2

        disp(['looping day: ' num2str(loopDays)])

        if loopDays == 1
            data_dir    = [UtkuPath 'AG_2019-12-23_NSD' '/'];
            chanDataDir = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/'];
            snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        elseif loopDays == 2
            data_dir    = [UtkuPath 'AG_2019-12-27_NSD' '/'];
            chanDataDir  = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/'];
            snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        end
                
        if flagAdjusted
            saveName = ['ratAGday' num2str(loopDays) '_1000randomAdjustedSpikeTimeData.mat' ];
        else
            saveName = ['ratAGday' num2str(loopDays) '_1000randomSpikeTimeData.mat' ];
        end

        if ~(exist([snippetsDir saveName],'file') == 2)
            
            neurons = loadUtkuV2(num2str(loopDays));
            if flagAdjusted
               if loopDays == 1
                   load([data_dir 'AG_2019-12_23_NSD_adjustedSpikeTimes.mat']);
               elseif loopDays == 2
                   load([data_dir 'AG_2019-12-27_NSD_adjustedSpikeTimes.mat']);
               end
            end

            spikeTimeData = table;
            for loopUnits = 1:length(neurons.cell_metrics.UID) 
                
                neuronID   = neurons.cell_metrics.UID(loopUnits);
                if flagAdjusted            
                    if isempty(find(adjustedSpikeTimes.cellID == neuronID))
                        continue 
                    else 
                        spikeTimes = adjustedSpikeTimes.spikeTimes{find(adjustedSpikeTimes.cellID == neuronID),1};
                    end
                else
                    spikeTimes = neurons.cell_metrics.spikes.times{1,loopUnits};
                end

                if length(spikeTimes) > N
                    spikeTimes = spikeTimes(randperm(length(spikeTimes),N));
                end

                spikeTimeData.neuronID(loopUnits)    = neuronID;
                spikeTimeData.spikeTimes{loopUnits}  = spikeTimes;

            end

            save([snippetsDir saveName],'spikeTimeData')

        else
            disp('1000 random spike time data file exists')
        end
    end
    
    for loopDays = 1:2

        disp(['looping day: ' num2str(loopDays)])

        if loopDays == 1
            chanDataDir  = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/'];
        elseif loopDays == 2
            chanDataDir  = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD//'];
        end

        snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];

        saveName = ['ratAGday' num2str(loopDays) '_1000randomSpikeTimeData.mat' ];
        
        neurons = loadUtkuV2(num2str(loopDays));
        
        load([snippetsDir '/' saveName]);
        
        % loop channels
        for loopChannel = 1:192
            
            disp(['looping channel: ' num2str(loopChannel)])
            
            snippetsData            = table;
            snippetsData.neuronID   = spikeTimeData.neuronID;
            snippetsData.spikeTimes = spikeTimeData.spikeTimes;
            
            chanData = load([ chanDataDir 'shank' num2str(channelShanks(loopChannel)) '/ch' num2str(loopChannel) 'highpass300hz.mat'],'data');
            
            % loop units
            for loopUnits = 1:length(neurons.cell_metrics.UID)
                
%                 disp([num2str((loopUnits/length(neurons.neuron_ids)*100)) '%'])
                
                spikeTimes = spikeTimeData.spikeTimes{loopUnits};


                % load time series data for that channel 
                spikeTimeIdxCell = round(spikeTimes*fs) + 1;

                [waveAvg,waveforms] = waveformAvg(double(chanData.data), ... 
                                                             spikeTimeIdxCell, ... 
                                                             preLength,postLength,fpass,fs,false);
                
                snippetsData.waveforms{loopUnits}   = waveforms;
                snippetsData.waveformAvg{loopUnits} = waveAvg;

            end
            
            if flagAdjusted
                saveName = ['ch' num2str(loopChannel) 'adjustedSnippets.mat' ];
            else
                saveName = ['ch' num2str(loopChannel) 'snippets.mat' ];
            end
            
            save([snippetsDir 'shank' num2str(channelShanks(loopChannel)) '/' saveName],'snippetsData')
            
        end
    end
        
end
 