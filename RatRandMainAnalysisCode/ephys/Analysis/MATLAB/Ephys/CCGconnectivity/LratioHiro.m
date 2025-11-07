function [Lratio] = LratioHiro(resClus,resNonClus,shank)
    
    %% constants
    
    session_name  = 'RoyMaze1';
    shankChanList = {1:8  , ...
                     9:16 , ...
                     17:24, ...
                     25:32, ...
                     33:40, ...
                     41:48, ...
                     49:56, ...
                     57:64};
    
    fpass      = 300;
    fs         = 30000;
    filterFlag = false;
    
    preLength  = 36;
    postLength = 36;

    %% cluster unit
    
    chanList = shankChanList{shank};
    
    tic
    for i = 1:8
        i
        [chanDataClu, onsetTime] = HiroLoad300hz(session_name,chanList(i));
        spikeTimeIdxClu = round((resClus-onsetTime/1e6)*fs) + 1; 
        [avgWaveClu(:,i), waveClu{i}] = waveformAvg(chanDataClu,spikeTimeIdxClu,preLength,postLength,fpass,fs,filterFlag);
    end
    
    for i = 1:8
        i
        [chanDataNonClu, onsetTime] = HiroLoad300hz(session_name,chanList(i));
        spikeTimeIdxNonClu = round((resNonClus-onsetTime/1e6)*fs) + 1; 
        [avgWaveNonClu(:,i), waveNonClu{i}] = waveformAvg(chanDataNonClu,spikeTimeIdxNonClu,preLength,postLength,fpass,fs,filterFlag);
    end
    toc
    
    [deltaTroughAmps,Idx] = sort(min(avgWaveClu-avgWaveNonClu,[],1));
    deltaTroughAmps       = deltaTroughAmps(1:4);
    Idx                   = Idx(1:4);
    
    waveCluSort           = waveClu(Idx);
    waveNonCluSort        = waveNonClu(Idx); 
    
    for i = 1:4
        troughAmpsClu(i,:)    = min(waveCluSort{i});
        troughAmpsNonClu(i,:) = min(waveNonCluSort{i});
        meanTroughAmpsClu(i)  = mean(troughAmpsClu(i,:));
    end
    
    sigma = troughAmpsClu*troughAmpsClu';
    meanTroughAmpsClu = meanTroughAmpsClu';
    
    %% non-cluster spikes
    
    for i = 1:size(troughAmpsNonClu,2)
        Dsq(i) = (troughAmpsNonClu(:,i) - meanTroughAmpsClu)'*(sigma^(-1))*(troughAmpsNonClu(:,i) - meanTroughAmpsClu);
    end
    
    for i = 1:size(troughAmpsNonClu,2)
        L(i) = 1-chi2cdf(Dsq(i),4);
    end
    L = sum(L);
    
    Lratio = L/length(waveClu{1});
end
