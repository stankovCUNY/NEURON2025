function BapunGlobalSnapshot 
    
    preLength  = 36;
    postLength = 36;
    
    animalIDs  = {'N','S','U'};
    numOfChans = [128,195,192];
    
    for loopAnimal = 1:length(animalIDs)
        
        RatID = animalIDs{loopAnimal};
        
        disp(['looping rat: ' RatID])
        
        chanDataDir  = ['/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/Rat' RatID '/'];
        snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        
        if strcmp(RatID,'N')
            dataPath   = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/';
            neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
            numOfChans = 128;
        elseif strcmp(RatID,'S')
            dataPath   = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/';
            neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
            numOfChans = 195;
        elseif strcmp(RatID,'U')
            dataPath   = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/';
            neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
            numOfChans = 192;
        end
    
        unitsGlobalSnapshots = table; 
        
        for loopUnits = 1:length(neurons.neuron_ids) 
            
            disp(['unit ' num2str(loopUnits) '/' num2str(length(neurons.neuron_ids))])
            
            neuronID = neurons.neuron_ids(loopUnits) + 1;
            
            peakChan   = neurons.peak_channels(loopUnits) + 1;
            if isempty(peakChan)
               continue 
            end

            globalSnapshot = zeros(numOfChans,preLength + postLength);
            for loopChannels = 1:numOfChans
                
                if strcmp(animalIDs{loopAnimal},'S')
                    if (loopChannels == 52) || (loopChannels >= 193)
                        continue
                    end
                end
            
                
                load([snippetsDir '/ch' num2str(loopChannels) 'snippets.mat'])
                globalSnapshot(loopChannels,:) = snippetsData.waveformAvg{loopUnits}';
            end

            unitsGlobalSnapshots.neuronID(loopUnits)  = neuronID;
            unitsGlobalSnapshots.waveforms{loopUnits} = globalSnapshot;

        end

        unitsGlobalSnapshots(find(unitsGlobalSnapshots.neuronID == 0),:) = [];

        saveName = 'unitsGlobalSnapshots.mat';
        save([dataPath saveName],'unitsGlobalSnapshots')
        
    end
    
end