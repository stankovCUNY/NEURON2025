function [N,NperTrial] = PSTH(spikeTimes,trialTimes,trialLength,preTrialLength,binSize,fs,color)
    
    PSTHout   = [];
    NperTrial = zeros(size(trialTimes,1),(preTrialLength + trialLength)/binSize);
    
%     trialTimes(trialTimes < (max(spikeTimes) - trialLength/fs)) = [];
    
    for i = 1:size(trialTimes,1)
    
        PSTHtemp = [];
        PSTHtemp = spikeTimes((spikeTimes > trialTimes(i)) & (spikeTimes < (trialTimes(i) + (preTrialLength+trialLength)/fs)))';
        PSTHtemp = PSTHtemp - trialTimes(i);
        
        PSTHout   = [PSTHout PSTHtemp];
        NperTrial(i,:) = histcounts(PSTHtemp,(preTrialLength + trialLength)/binSize,'BinWidth',binSize/fs,'BinLimits',[0,(preTrialLength + trialLength)/fs]);
    end
    
    N = histcounts(PSTHout,(preTrialLength + trialLength)/binSize,'BinWidth',binSize/fs,'BinLimits',[0,(preTrialLength + trialLength)/fs]);
%     histogram(PSTHout,'BinWidth',binSize,'BinLimits',[0,preTrialLength + trialLength],'FaceColor',color,'EdgeColor',color)

end