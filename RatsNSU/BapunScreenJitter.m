function [filteredPairs,nullSpikeTimes] = BapunScreenJitter(RatID,pairType,nullFlag)
    
    dataDir = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    if strcmp(pairType,'all')
        idx = 1:211;
    elseif strcmp(pairType,'exquisite')
        idx = 91:121;
    end
    
    if nullFlag
        load([dataDir 'Rat' RatID '/Bapuns_Rat' RatID 'null_jitterScreen_pairs.mat'])
    else
        load([dataDir 'Rat' RatID '/Bapuns_Rat' RatID '_jitterScreen_pairs.mat'])
    end
    
    filter = zeros(size(pairCCGscreen.pair,1),1);

    for i = 1:size(pairCCGscreen.pair,1)
        
        if (pairCCGscreen.pair(i,1) == 0) && (pairCCGscreen.pair(i,2) == 0)
           continue 
        end
        
        exc = pairCCGscreen.GSPExc{i};
        inh = pairCCGscreen.GSPInh{i};
        
        if sum(exc(idx)) + sum(inh(idx)) > 0
            filter(i) = 1;
        end

    end

    filteredPairs = pairCCGscreen.pair(find(filter),:);

    % unique pairs only
    [~,~,IC] = unique(sort(filteredPairs,2),'row','first');

    filteredPairsNew = [];
    for i = 1:max(IC)
        tmp = find(IC == i);
        filteredPairsNew = [filteredPairsNew; filteredPairs(tmp(1),:)];
    end

    filteredPairs = filteredPairsNew;
    
    %% null spike Times
    if nullFlag
        nullSpikeTimes.spikeTimes = pairCCGscreen.spikeTimes;
        nullSpikeTimes.spikeInd   = pairCCGscreen.spikeInd;
    else
        nullSpikeTimes.spikeTimes = [];
        nullSpikeTimes.spikeInd   = [];
    end
    
end