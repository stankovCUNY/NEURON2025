function filtSpikeTimes = filterSpikes(spikeTimes,condMarkers)
    
    filtSpikeTimes = [];

    for i = 1:size(condMarkers,1)
        filtSpikeTimes = [filtSpikeTimes; spikeTimes((spikeTimes > condMarkers(i,1)) & (spikeTimes < condMarkers(i,2)))];
    end

end