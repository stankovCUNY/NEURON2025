function [filteredPairs,nullSpikeTimes] = UtkuScreenJitterV2(dayNo,pairType,nullFlag)
    
    if strcmp(pairType,'all')
        idx = 1:211;
    elseif strcmp(pairType,'exquisite')
        idx = 91:121;
    end
    
    if nullFlag
        load(['AG' dayNo 'null_jitterScreen_pairs.mat'])
    else
        load(['AG' dayNo '_jitterScreen_pairs.mat'])
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