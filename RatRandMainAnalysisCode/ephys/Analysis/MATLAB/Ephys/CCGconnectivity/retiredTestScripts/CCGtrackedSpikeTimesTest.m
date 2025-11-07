function [CCG,spikeTimesInBin] = CCGtrackedSpikeTimesTest(R,T,binSize,lagLimit)
    
    % INPUT:
    % n1 - neuron 1 spike times [seconds]
    % n2 - neuron 2 spike times [seconds]
    % binSize  [seconds]
    % lagLimit [seconds]
    %
    % OUTPUT:
    % CCG [raw]
    % spikeTimesInBin - spikes times from each neuron belong in that bin
    
    % Atanas Stankov, CUNY 2023
    
    %% variable initiation    
    buffer = 1; % buffer needed to avoid undercoundting at desired lag time bins at the edges
    intHalfWin  = floor(lagLimit/binSize)+buffer;
    intLagTimes = -intHalfWin:intHalfWin;
    halfWin     = intHalfWin*binSize;
    
    %% pre-filter spikes outside of CCG window
    [sortedUnion, idxUnion] = sort([R; T]);
%     diffUnion   = diff(sortedUnion);
%     filtUnion   = diffUnion(diffUnion < halfWin);
   
    %%
    
%     tic
%     for i = 1:length(R)
% 
%         % find nearest spike of n1(i) in n2
%         delta   = T-R(i);
% 
%         [~,idx] = min(abs(delta));
% 
%         allMinDelta(i) = delta(idx);
%         allIdx(i)      = idx;
% 
%     end
%     toc
    
%     allIdx = dsearchn(T,R);
%     allIdx = binarySearchNearest(T, R);
    allIdx = find(ismember(T, R));
    allMinDelta = T(allIdx) - R;
    
    %% remove n1 spike times that have nearest n2 neuron spike times outside of CCG window limit
    allIdxFilt      = allIdx(     abs(allMinDelta) <= halfWin);
    allMinDeltaFilt = allMinDelta(abs(allMinDelta) <= halfWin);
    
    %% find all n2 spike times around each n1 spike time within the CCG window
%     n2NearSpikes = {};
%     for i = 1:length(allIdxFilt)
% 
%         n1Marker = T(allIdxFilt(i)) - allMinDeltaFilt(i);
%         n2NearSpikes{i} = T((T >=  (n1Marker - halfWin)) & ... 
%                             (T <=  (n1Marker + halfWin)));
% 
%     end
    n1Marker = T(allIdxFilt) - allMinDeltaFilt;
    n2NearSpikes = cell(length(n1Marker), 1);
    
    for i = 1:length(n1Marker)
        % Find the lower and upper bounds of the range using binary search
        lowerIdx = findLowerBound(T, n1Marker(i) - halfWin);
        upperIdx = findUpperBound(T, n1Marker(i) + halfWin);

        if lowerIdx <= upperIdx
            % Extract the elements within the range
            n2NearSpikes{i} = T(lowerIdx:upperIdx);
        else
            % No elements within the range
            n2NearSpikes{i} = [];
        end
    end
    
    %% list spike times in each CCG bin
    spikeTimesInBin{length(intLagTimes)} = [];
    for i = 1:length(allIdxFilt)

        n1Marker = T(allIdxFilt(i)) - allMinDeltaFilt(i);
        delta    = n2NearSpikes{i}   - n1Marker;
        CCGlag   = round(delta/binSize);

        for j = 1:length(delta)
            spikeTimesInBin{find(intLagTimes == CCGlag(j))} = ... 
                [spikeTimesInBin{find(intLagTimes == CCGlag(j))}; ...
                 n1Marker n2NearSpikes{i}(j)];
        end   
    end
    
    %% compute CCG with edge bins removed 
    for i = 2:length(intLagTimes)-1
        CCG(i-1)    = size(spikeTimesInBin{i},1);
    end
    
    % remove extra edge bins for the tracked spike times
    spikeTimesInBin = {spikeTimesInBin{2:end-1}};
end

% Function to find the nearest points using binary search
function nearestIndices = binarySearchNearest(sortedArray, targets)
    n = numel(targets);
    nearestIndices = zeros(size(targets));

    left = 1;
    right = numel(sortedArray);

    for i = 1:n
        target = targets(i);
        minDistance = Inf;

        while left <= right
            mid = floor((left + right) / 2);
            currentDistance = abs(sortedArray(mid) - target);

            if currentDistance < minDistance
                minDistance = currentDistance;
                nearestIndices(i) = mid;
            end

            if sortedArray(mid) == target
                break;
            elseif sortedArray(mid) < target
                left = mid + 1;
            else
                right = mid - 1;
            end
        end

        % Reset left and right for the next target
        left = 1;
        right = numel(sortedArray);
    end
end

% Binary search for the lower bound
function lowerIdx = findLowerBound(arr, target)
    lower = 1;
    upper = length(arr);
    lowerIdx = [];

    while lower <= upper
        mid = floor((lower + upper) / 2);
        if arr(mid) >= target
            lowerIdx = mid;
            upper = mid - 1;
        else
            lower = mid + 1;
        end
    end
end

% Binary search for the upper bound
function upperIdx = findUpperBound(arr, target)
    lower = 1;
    upper = length(arr);
    upperIdx = [];

    while lower <= upper
        mid = floor((lower + upper) / 2);
        if arr(mid) <= target
            upperIdx = mid;
            lower = mid + 1;
        else
            upper = mid - 1;
        end
    end
end