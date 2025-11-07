function artifactSpikeTimes = zeroLagArtifactFix
    neurons = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');

    allSpikeTimes  = [];
    allShankLabels = [];

    for i = 1:length(neurons.neuron_type)

        if strcmp(neurons.neuron_type(i,:),'mua  ')
            continue
        end

        allSpikeTimes  = [allSpikeTimes neurons.spiketrains{i}];
        allShankLabels = [allShankLabels (double(neurons.shank_ids(i))+1)*ones(1,length(neurons.spiketrains{i}))];

    end

    %list of elements
    [uniqueSpikeTimes,idx] = unique(allSpikeTimes);    
    
    %provides a count of each element's occurrence
    count            = hist(allSpikeTimes,uniqueSpikeTimes); 

    artifactSpikeTimes  = uniqueSpikeTimes(find(count == 2));
    
%     parpool(6)
    parfor i = 1:length(artifactSpikeTimes)
        
        
        spikePairShanks = allShankLabels(find(artifactSpikeTimes(i) == allSpikeTimes));
        
        if spikePairShanks(1) == spikePairShanks(2)
            sameShank(i) = 1;
        elseif spikePairShanks(1) ~= spikePairShanks(2)
            sameShank(i) = 0;
        end
        
        
    end
        
end
