function Lratio
    
    %% cluster neuron
    
    [waveAvg1,waveforms1,chans1] = extractWaveformsWholeShank(res1,shank1,channelSet,cellID1,'U',fs,preLength,postLength);
    
    [meanTroughAmps,Idx] = sort(min(waveAvg1,[],2));
    meanTroughAmps       = meanTroughAmps(1:4);
    
    waveformSort         = waveforms1(Idx);
    
    for i = 1:4
        troughAmps(i,:) = min(waveforms1{i});
    end
    
    sigma = troughAmps*troughAmps';
    
    %% non-cluster neuron
    [waveAvg2,waveforms2,chans2] = extractWaveformsWholeShank(res2,shank2,channelSet,cellID2,'U',fs,preLength,postLength);
        
    for i = 1:4
        troughAmps(i,:) = min(waveforms2{i});
    end
    
    for i = 1:size(troughAmps,2)
        Dsq(i) = (troughAmps(:,i) - meanTroughAmps)'*(sigma^(-1))*(troughAmps(:,i) - meanTroughAmps);
    end
    
    for i = 1:size(troughAmps,2)
        L(i) = 1-chi2cdf(Dsq(i),4);
    end
    L = sum(L);
    
    Lratio = L/length(waveforms1{1});
end
