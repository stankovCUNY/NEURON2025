function HiroGlobalSnapshot 

    numOfChans = 64;
    
    preLength  = 36;
    postLength = 36;
    
    spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
    
    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')
    
    extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
    Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);
    
    waveDataSpCounts = [];
    for i = 1:size(Session.times,2)
        waveDataSpCounts(i) = size(Session.times{i},1);
    end
    
    %%
    
    chanDataDir = ['/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/'];
    snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
    
    unitsGlobalSnapshots = table; 
    
    for loopUnits = 1:length(spikes.RoyMaze1) 
        
        tic
        peakChan = Session.maxWaveformCh(length(spikes.RoyMaze1(loopUnits).time) == waveDataSpCounts) + 1;
        if isempty(peakChan)
           continue 
        end
        
        neuronID = loopUnits;
        
        globalSnapshot = zeros(numOfChans,preLength + postLength);
        for loopChannels = 1:numOfChans
            load([snippetsDir '/ch' num2str(loopChannels) 'snippets.mat'])
            globalSnapshot(loopChannels,:) = snippetsData.waveformAvg{loopUnits}';
        end
        
        unitsGlobalSnapshots.neuronID(loopUnits)  = loopUnits;
        unitsGlobalSnapshots.waveforms{loopUnits} = globalSnapshot;
        toc 
        
    end
    
    unitsGlobalSnapshots(find(unitsGlobalSnapshots.neuronID == 0),:) = [];
    
    saveName = 'unitsGlobalSnapshots.mat';
    save([dataPath saveName],'unitsGlobalSnapshots')
    
end