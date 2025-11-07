function [syncSp,syncCCG] = SyncSpikes(n1,n2,LOI,binSize,lag)
    
    % All inputs are assumped in seconds except for lag and bin size which
    % are interger lengths of samples.
    % INPUT:
    % n1 - neuron 1 spike times
    % n2 - neuron 2 spike times
    % LOI - lags of interest, AKA significant bins 
    % binary: 1 - significant, 0 - ignore
    % lag: labels for lags, right now only using the lag size
    %
    % OUTPUT:
    % SyncSp  - informed spikes, AKA synchronous spikes
    % SyncCCG - bin counts only for significant bins. I use this for sanity
    % checks to compare with full CCG.
    
    % Atanas Stankov CUNY 2021
        
    if sum(LOI) == 0
        
        syncSp  = []; % informed spikes
        syncCCG = zeros(211,1);
        return;
    else

        fs = 30000;
        
%         n1  = round(n1*fs);
%         n2  = round(n2*fs);
        
        allIdx   = zeros(length(n1),1);
        allMinDelta = zeros(length(n1),1);

        halfWin     = 105*(1/fs);
        intLagTimes = -105:105;
%         part        = 50; % make this adaptive in the future

        parfor i = 1:length(n1)

            % find nearest spike of n1(i) in n2
            delta = n2-n1(i);
            [~,idx] = min(abs(delta));

            allMinDelta(i) = delta(idx);
            allIdx(i)      = idx;

        end
        
%         for i = 1:ceil(length(n1)/part)
%             if i == ceil(length(n1)/part)
%                 n1_part = n1(1+(i-1)*part:end);
%                 n1_part = repmat(n1_part,1,length(n2))';
%                 n2_part = repmat(n2,1,length(n1) - (i-1)*part);
%                 
% %                 find nearest spike of n1(i) in n2
%                 delta = n2_part-n1_part;
%                 [~,idx] = min(abs(delta));
%                 
% %                 tricky indices
%                 idxMin = idx + 1:length(n2):(length(n1) - (i-1)*part)*length(n2);
%                 allMinDelta(1+(i-1)*part:end) = delta(idxMin);
%                 allIdx(     1+(i-1)*part:end) = idx;
%             elseif i < ceil(length(n1)/part)
%                 n1_part = n1(1+(i-1)*part:i*part);
%                 n1_part = repmat(n1_part,1,length(n2))';
%                 n2_part = repmat(n2,1,part);
% 
% %                 find nearest spike of n1(i) in n2
%                 delta = n2_part-n1_part;
%                 [~,idx] = min(abs(delta));
%                 
% %                 tricky indices
%                 idxMin = idx + 1:length(n2):part*length(n2);
%                 allMinDelta(1+(i-1)*part:i*part) = delta(idxMin);
%                 allIdx(     1+(i-1)*part:i*part) = idx;
%             end
%         end
        
        % remove n1 neurons that have nearest n2 neuro further than CCG window
        n1Filt          = n1(         abs(allMinDelta) <= halfWin);
        allIdxFilt      = allIdx(     abs(allMinDelta) <= halfWin);
        allMinDeltaFilt = allMinDelta(abs(allMinDelta) <= halfWin);

        % find all n2 spikes around each n1 spike within CCG window
        n2NearSpikes = {};
        for i = 1:length(allIdxFilt)

            n1Marker = n2(allIdxFilt(i)) - allMinDeltaFilt(i);
            n2NearSpikes{i} = n2((n2 >=  (n1Marker - halfWin)) & ... 
                                 (n2 <=  (n1Marker + halfWin)));

        end

        % compute CCG
        syncSpAll{211} = [];
        for i = 1:length(allIdxFilt)

            n1Marker = n2(allIdxFilt(i)) - allMinDeltaFilt(i);
            delta    = n2NearSpikes{i}   - n1Marker;
            CCGlag   = round(delta*fs);

            for j = 1:length(delta)
                syncSpAll{find(intLagTimes == CCGlag(j))} = ... 
                    [syncSpAll{find(intLagTimes == CCGlag(j))}; ...
                     n1Marker n2NearSpikes{i}(j)];
            end   
        end
        
        
        for i = 1:211
            syncCCG(i)    = size(syncSpAll{i},1);
        end
        syncCCG(~LOI) = 0;
        
        syncSp = [];
        LOIidx = find(LOI);
        for i = 1:length(LOIidx)
            syncSp = [syncSp; syncSpAll{LOIidx(i)}];
        end
        syncSp = sort(min(syncSp,[],2),2);
        
    end
end