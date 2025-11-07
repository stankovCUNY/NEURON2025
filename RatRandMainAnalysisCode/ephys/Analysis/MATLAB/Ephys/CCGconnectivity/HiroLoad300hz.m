function [chanData, onsetTime] = H(session_name,channel)

    if strcmp(session_name,'RoyMaze1')
        filepath = '/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/';
    elseif strcmp(session_name,'KevinMaze1')
        filepath = '/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Kevin/maze/highpass300hz/';
    end
    
    load([filepath 'ch' num2str(channel) 'highpass300hz.mat']);
    
    chanData  = data.channel;
    onsetTime = data.onsetTime;
end