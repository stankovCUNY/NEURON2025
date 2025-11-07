function BapunExtractAndCurateWaveformData
    
    flagAdjusted = true;

    % how many random spikes
    N          = 1000; 

    fs         = 30000;
    fpass      = 300;
    
    preLength  = 36;
    postLength = 36;

    animalIDs  = {'N','S','U'};
    numOfChans = [128,195,192];

    for loopAnimal = 1:length(animalIDs)
        
        RatID = animalIDs{loopAnimal};
        
        disp(['looping rat: ' RatID])
        
        chanDataDir  = ['/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/Rat' RatID '/'];
        snippertsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        
        if flagAdjusted
            saveName = ['rat' RatID '1000randomAdjustedSpikeTimeData.mat' ];
        else
            saveName = ['rat' RatID '1000randomSpikeTimeData.mat' ];
        end
        
        if ~(exist([snippertsDir saveName],'file') == 2)
            if strcmp(RatID,'N')
                neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
                if flagAdjusted
                    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatN_adjustedSpikeTimes.mat');
                end
            elseif strcmp(RatID,'S')
                neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
                if flagAdjusted
                    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatS_adjustedSpikeTimes.mat');
                end
            elseif strcmp(RatID,'U')
                neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
                 if flagAdjusted
                     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatU_adjustedSpikeTimes.mat');
                 end
            end

            spikeTimeData = table;
            for loopUnits = 1:length(neurons.neuron_ids) 
                
                neuronID = double(neurons.neuron_ids(loopUnits) + 1);
                if flagAdjusted            
                    if isempty(find(adjustedSpikeTimes.cellID == neuronID))
                        continue 
                    else 
                        spikeTimes = adjustedSpikeTimes.spikeTimes{find(adjustedSpikeTimes.cellID == neuronID),1};
                    end
                else
                    spikeTimes = neurons.spiketrains{1,loopUnits};
                end

                if length(spikeTimes) > N
                    spikeTimes = spikeTimes(randperm(length(spikeTimes),N));
                end

                neuronID   = neurons.neuron_ids(loopUnits) + 1;

                spikeTimeData.neuronID(loopUnits)    = neuronID;
                spikeTimeData.spikeTimes{loopUnits}  = spikeTimes;

            end
            
            save([snippertsDir saveName],'spikeTimeData')
            
        else
            disp('1000 random spike time data file exists')
        end
        
    end
    
    for loopAnimal = 1:length(animalIDs)
        
        RatID = animalIDs{loopAnimal};
        
        disp(['looping rat: ' RatID])
        
        chanDataDir  = ['/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/Rat' RatID '/'];
        snippertsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        
        if flagAdjusted
            saveName = ['rat' RatID '1000randomAdjustedSpikeTimeData.mat' ];
        else
            saveName = ['rat' RatID '1000randomSpikeTimeData.mat' ];
        end
        
        if strcmp(RatID,'N')
            neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        elseif strcmp(RatID,'S')
            neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        elseif strcmp(RatID,'U')
            neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
        end
        
        load([snippertsDir '/' saveName]);
        
        loopLenght = numOfChans(strcmp(animalIDs,RatID));
        % loop channels
        for loopChannel = 1:loopLenght
            
            if strcmp(animalIDs{loopAnimal},'S')
                if (loopChannel == 52) || (loopChannel >= 193)
                    continue
                end
            end
            
            disp(['looping channel: ' num2str(loopChannel)])
            
            snippetsData            = table;
            snippetsData.neuronID   = spikeTimeData.neuronID;
            snippetsData.spikeTimes = spikeTimeData.spikeTimes;
            
            chanData = load([ chanDataDir 'ch' num2str(loopChannel) 'highpass300hz.mat'],'data');
            
            % loop units
            for loopUnits = 1:length(neurons.neuron_ids) 
                
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
            
            save([snippertsDir saveName],'snippetsData')
            
        end
    end
        
end
 