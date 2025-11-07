function HiroExtractAndCurateWaveformDataPeakChan
    
    % how many random spikes
    
    fs         = 30000;
    fpass      = 300;
    
    preLength  = 36;
    postLength = 36;
    
    chanDataDir  = ['/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/'];
    snippetsDir = [chanDataDir 'snippetsPeakChan/'];
    
    %%
    
    spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
    
    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-basics.mat'),   'basics');
    
    extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
    Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);
    
    waveDataSpCounts = [];
    for i = 1:size(Session.times,2)
        waveDataSpCounts(i) = size(Session.times{i},1);
    end
    
    %%
    
    spikeTimeData = table;
    for loopUnits = 1:length(spikes.RoyMaze1) 

        temp = spikes.RoyMaze1(loopUnits).time;
        peakChan = Session.maxWaveformCh(length(spikes.RoyMaze1(loopUnits).time) == waveDataSpCounts) + 1;
        if isempty(peakChan)
           continue 
        end
        temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
        spikeTimes = temp/1e6; % convert to ms
        
        chanData = load([ chanDataDir 'ch' num2str(peakChan) 'highpass300hz.mat'],'data');
        
        spikeTimeIndxCell = round((spikeTimes-chanData.data.onsetTime/1e6)*fs) + 1;

        [waveAvg,waveforms] = waveformAvg(double(chanData.data.channel), ... 
                                          spikeTimeIndxCell, ... 
                                          preLength,postLength,fpass,fs,false);

        neuronID = loopUnits;
                                      
        snippetsData.waveforms   = waveforms;
        snippetsData.waveformAvg = waveAvg;

        saveName = ['unit' num2str(neuronID) 'peakCh' num2str(peakChan) 'snippets.mat' ];
        save([snippetsDir saveName],'snippetsData')

    end
        
end
 