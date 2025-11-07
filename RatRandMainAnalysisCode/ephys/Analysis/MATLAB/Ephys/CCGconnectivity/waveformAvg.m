function [spikeAverage,spikeWaves] = waveformAvg(timeSeries,markerIdx,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag)
    
    % initialize 
    spikeAverage = zeros(preLength + postLength,1);
    spikeWaves   = [];
    
    % generate high-pass filter
    [b, a, k] = butter(9, fpass / fs * 2, 'high');
    [sos_high,g] = zp2sos(b, a, k);
    
    if filterFlag
        tic
%             temp     = highpass(temp,fpass,fs);
%         timeSeries = sosfilt(sos_high, timeSeries);
        timeSeries = filtfilt(sos_high,g,timeSeries);
        toc
    end
    
    normConstant = length(markerIdx);
    
    for i = 1:length(markerIdx) 
        
        idx = markerIdx(i);
        
        % avoiding indexing issues with times less than pre and post
        if (idx < preLength + 1) || (idx > length(timeSeries) - postLength)
            continue
        end
%         temp         = timeSeries(idx - preLength + 1:idx + postLength) - mean(timeSeries(idx - preLength:idx - preLength+14));  % base line length 0.5msec
        
        if noDemeanFlag
            temp     = timeSeries(idx - preLength + 1:idx + postLength);
        else
            temp     = timeSeries(idx - preLength + 1:idx + postLength) - mean(timeSeries(idx - preLength:idx));  % base line length 1msec
        end
        
        temp = double(temp);
        
        temp         = temp/1000; % converting units to mV
        spikeAverage = spikeAverage + temp/normConstant;
        spikeWaves   = [spikeWaves temp];
        
    end
end