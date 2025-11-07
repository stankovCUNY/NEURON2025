function HiroExtractAndCurateWaveformData
    
    % how many random spikes
    N          = 1000; 

    fs         = 30000;
    fpass      = 300;
    
    preLength  = 36;
    postLength = 36;

    numOfChans = 64;
     
    chanDataDir  = ['/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/'];
    snippertsDir = [chanDataDir 'snippetsPerChanPerUnit/'];

    saveName = 'Roymaze1000randomSpikeTimeData.mat';
    
    %%
    
    spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
    
    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-basics.mat'),   'basics');
    
    %%
    
    if ~(exist([snippertsDir saveName],'file') == 2)

        spikeTimeData = table;
        for loopUnits = 1:length(spikes.RoyMaze1) 
            
            temp = spikes.RoyMaze1(loopUnits).time;
            temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
            temp = temp/1e6; % convert to ms

            if length(temp) > N
                spikeTimes = temp(randperm(length(temp),N));
            end
            spikeTimes = sort(spikeTimes);
            
            neuronID = loopUnits;

            spikeTimeData.neuronID(loopUnits)    = neuronID;
            spikeTimeData.spikeTimes{loopUnits}  = spikeTimes;

        end

        save([snippertsDir saveName],'spikeTimeData')

    else
        disp('1000 random spike time data file exists')
    end
        
    %%
                
    load([snippertsDir '/' saveName]);
        
    % loop channels
    for loopChannel = 1:numOfChans

        disp(['looping channel: ' num2str(loopChannel)])

        snippetsData            = table;
        snippetsData.neuronID   = spikeTimeData.neuronID;
        snippetsData.spikeTimes = spikeTimeData.spikeTimes;

        chanData = load([ chanDataDir 'ch' num2str(loopChannel) 'highpass300hz.mat'],'data');

        % loop units
        for loopUnits = 1:length(spikes.RoyMaze1)

%             spikeTimes = spikeTimeData.spikeTimes{loopUnits}/1e6 - (chanData.data.onsetTime/1e6)*fs;
            spikeTimes = spikeTimeData.spikeTimes{loopUnits};
            
            % load time series data for that channel 
            spikeTimeIndxCell = round((spikeTimes-chanData.data.onsetTime/1e6)*fs) + 1;

            [waveAvg,waveforms] = waveformAvg(double(chanData.data.channel), ... 
                                              spikeTimeIndxCell, ... 
                                              preLength,postLength,fpass,fs,false);

            snippetsData.waveforms{loopUnits}   = waveforms;
            snippetsData.waveformAvg{loopUnits} = waveAvg;

        end

        saveName = ['ch' num2str(loopChannel) 'snippets.mat' ];
        save([snippertsDir saveName],'snippetsData')

    end
        
end
 