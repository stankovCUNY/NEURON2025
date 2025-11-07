function UtkuGlobalSnapshot 
    
    preLength  = 36;
    postLength = 36;
    
    numOfChans = 192;
    
    channelShanks = [ones(32,1); 2*ones(32,1); 3*ones(32,1); 4*ones(32,1); 5*ones(32,1); 6*ones(32,1);];
    
    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/'; 
    
    for loopDays = 1:2
                
        disp(['looping day: ' num2str(loopDays)])
        
        if loopDays == 1
            dataPath    = [UtkuPath 'AG_2019-12-23_NSD' '/'];
            chanDataDir = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/'];
            snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        elseif loopDays == 2
            dataPath    = [UtkuPath 'AG_2019-12-27_NSD' '/'];
            chanDataDir = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/'];
            snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        end
        
        neurons = loadUtkuV2(num2str(loopDays));
        
        unitsGlobalSnapshots = table; 
        
        for loopUnits = 1:length(neurons.cell_metrics.UID) 
            
            disp(['unit ' num2str(loopUnits) '/' num2str(length(neurons.cell_metrics.UID))])
            
            neuronID   = neurons.cell_metrics.UID(loopUnits);

            globalSnapshot = zeros(numOfChans,preLength + postLength);
            parfor loopChannels = 1:numOfChans
                
                snippetsData =  load_func([snippetsDir 'shank' num2str(channelShanks(loopChannels)) '/ch' num2str(loopChannels) 'snippets.mat']);
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

function snippetsData = load_func(file)

    load( file );
    
end