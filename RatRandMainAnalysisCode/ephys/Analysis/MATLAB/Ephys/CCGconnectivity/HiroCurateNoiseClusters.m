function HiroCurateNoiseClusters 
    
    spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
    
    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-basics.mat'),   'basics');
            
    % loop shanks
    for loopShank = 1:8
        
        % renew cell array
        unitStack = {};
        
        for i = 1:size(spikes.RoyMaze1,2)
            if spikes.RoyMaze1(i).id(1) == loopShank
                temp = spikes.RoyMaze1(i).time;
                temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
                unitStack{i} = temp;
            else
                unitStack{i} = [];
            end
            unitStack = unitStack(~cellfun('isempty',unitStack));
            
        end
        
        % noise cluster
        noiseCluster = unitStack{1};
        noiseCluster = noiseCluster'/1e6; % convert to seconds 
        saveNameNoise = ['RoyMaze_shank' num2str(loopShank) '_spikeTimesNoise'];
        save([dataPath saveNameNoise],'noiseCluster')

        % valid units cluster
        validUnits = [];
        for loopCluster = 2:size(unitStack,2)
            validUnits = [validUnits unitStack{loopCluster}];
        end
        validUnits = validUnits'/1e6; % convert to seconds 
        saveNameValidUnits = ['RoyMaze_shank' num2str(loopShank) '_spikeTimesValidUnits'];
        save([dataPath saveNameValidUnits],'validUnits')
                
    end

end