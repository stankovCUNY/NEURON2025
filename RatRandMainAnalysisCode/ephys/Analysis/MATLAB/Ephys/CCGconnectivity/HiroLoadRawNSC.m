function [chanData, onsetTime] = HiroLoadRawNSC(session_name,channel)

    if strcmp(session_name,'RoyMaze1')
        filepath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy_ncs_data/maze/';
    elseif strcmp(session_name,'KevinMaze1')
        filepath = '/media/nasko/WD_BLACK3/Kevin_ncs_data/maze/';
    end
    
    [ncs] = read_neuralynx_ncs([filepath 'CSC' num2str(channel) '.ncs']); % channel 

    chanData = ncs.dat; % 512 samples per one time period which spans (1/30) * 512 ms
    chanData = reshape(chanData,[],1);
    
    onsetTime = double(ncs.TimeStamp(1));
end