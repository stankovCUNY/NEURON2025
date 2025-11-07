function [phaseData,phaseBinSpikeTimes] = phaseAnalysis(res,ch,freqRange,session_name,fs)

    winSize = 0;
    
    phase = phaseExtraction(ch,freqRange,session_name,fs);
    t     = (0:1/fs:(1/fs)*length(phase))';

    %% curate spike times
    
    startTime = res(1)*fs; % first sample time
    
    res  = res*fs;         % change units from seconds to samples
    res  = res-startTime;  % reset first sample to 0 by using first sample time
    res  = round(res) + 1; % correct for miniscule decimals bc of change of units
    
    %% estimating firing rates
    
    noSamples = length(t);

    spikeBinTar = zeros(noSamples,1);
    spikeBinTar(res) = 1;
    
    if winSize == 0
        spikeCount = spikeBinTar;
    else
        spikeCount = movsum(spikeBinTar,winSize); % how to normalize here?
    end

    %% calculating phase preference (yes the code below is redundant and lazy!!!)
    
    phaseBins   = -pi:pi/4:pi;
    phaseLabels = zeros(noSamples,1);

    for i = 1:length(phase)

        if ((phase(i) > phaseBins(1)) && (phase(i) < phaseBins(2)))
            phaseLabels(i) = 1;
        elseif ((phase(i) > phaseBins(2)) && (phase(i) < phaseBins(3)))
            phaseLabels(i) = 2;
        elseif ((phase(i) > phaseBins(3)) && (phase(i) < phaseBins(4)))
            phaseLabels(i) = 3;
        elseif ((phase(i) > phaseBins(4)) && (phase(i) < phaseBins(5)))
            phaseLabels(i) = 4;
        elseif ((phase(i) > phaseBins(5)) && (phase(i) < phaseBins(6)))
            phaseLabels(i) = 5;
        elseif ((phase(i) > phaseBins(6)) && (phase(i) < phaseBins(7)))
            phaseLabels(i) = 6;
        elseif ((phase(i) > phaseBins(7)) && (phase(i) < phaseBins(8)))
            phaseLabels(i) = 7;
        elseif ((phase(i) > phaseBins(8)) && (phase(i) < phaseBins(9)))
            phaseLabels(i) = 8;
        end

    end
    
    phaseLabelsGroup1     = zeros(noSamples,1);
    phaseLabelsGroup1(phaseLabels == 1) = 1;
    phaseLatGroup1        = find(phaseLabelsGroup1);
    phaseSpikeTimesGroup1 = (find(spikeCount & phaseLabelsGroup1) + startTime - 1)/fs;
    clear phaseLabelsGroup1
    groupPhase1           = ((mean(spikeCount(phaseLatGroup1)))/mean(spikeCount));

    phaseLabelsGroup2     = zeros(noSamples,1);
    phaseLabelsGroup2(phaseLabels == 2) = 1;
    phaseLatGroup2        = find(phaseLabelsGroup2);
    phaseSpikeTimesGroup2 = (find(spikeCount & phaseLabelsGroup2) + startTime - 1)/fs;
    clear phaseLabelsGroup2
    groupPhase2           = ((mean(spikeCount(phaseLatGroup2)))/mean(spikeCount));

    phaseLabelsGroup3     = zeros(noSamples,1);
    phaseLabelsGroup3(phaseLabels == 3) = 1;
    phaseLatGroup3        = find(phaseLabelsGroup3);
    phaseSpikeTimesGroup3 = (find(spikeCount & phaseLabelsGroup3) + startTime - 1)/fs;
    clear phaseLabelsGroup3
    groupPhase3           = ((mean(spikeCount(phaseLatGroup3)))/mean(spikeCount));

    phaseLabelsGroup4     = zeros(noSamples,1);
    phaseLabelsGroup4(phaseLabels == 4) = 1;
    phaseLatGroup4        = find(phaseLabelsGroup4);
    phaseSpikeTimesGroup4 = (find(spikeCount & phaseLabelsGroup4) + startTime - 1)/fs;
    clear phaseLabelsGroup4
    groupPhase4           = ((mean(spikeCount(phaseLatGroup4)))/mean(spikeCount));

    phaseLabelsGroup5     = zeros(noSamples,1);
    phaseLabelsGroup5(phaseLabels == 5) = 1;
    phaseLatGroup5        = find(phaseLabelsGroup5);
    phaseSpikeTimesGroup5 = (find(spikeCount & phaseLabelsGroup5) + startTime - 1)/fs;
    clear phaseLabelsGroup5
    groupPhase5           = ((mean(spikeCount(phaseLatGroup5)))/mean(spikeCount));

    phaseLabelsGroup6     = zeros(noSamples,1);
    phaseLabelsGroup6(phaseLabels == 6) = 1;
    phaseLatGroup6        = find(phaseLabelsGroup6);
    phaseSpikeTimesGroup6 = (find(spikeCount & phaseLabelsGroup6) + startTime - 1)/fs;
    clear phaseLabelsGroup6
    groupPhase6           = ((mean(spikeCount(phaseLatGroup6)))/mean(spikeCount));

    phaseLabelsGroup7     = zeros(noSamples,1);
    phaseLabelsGroup7(phaseLabels == 7) = 1;
    phaseLatGroup7        = find(phaseLabelsGroup7);
    phaseSpikeTimesGroup7 = (find(spikeCount & phaseLabelsGroup7) + startTime - 1)/fs;
    clear phaseLabelsGroup7
    groupPhase7           = ((mean(spikeCount(phaseLatGroup7)))/mean(spikeCount));

    phaseLabelsGroup8     = zeros(noSamples,1);
    phaseLabelsGroup8(phaseLabels == 8) = 1;
    phaseLatGroup8        = find(phaseLabelsGroup8);
    phaseSpikeTimesGroup8 = (find(spikeCount & phaseLabelsGroup8) + startTime - 1)/fs;
    clear phaseLabelsGroup8
    groupPhase8           = ((mean(spikeCount(phaseLatGroup8)))/mean(spikeCount));

    phaseData = [groupPhase1 groupPhase2 groupPhase3 groupPhase4 ...
                 groupPhase5 groupPhase6 groupPhase7 groupPhase8];
             
    phaseBinSpikeTimes = {{phaseSpikeTimesGroup1}, ...
                          {phaseSpikeTimesGroup2}, ...
                          {phaseSpikeTimesGroup3}, ...
                          {phaseSpikeTimesGroup4}, ...
                          {phaseSpikeTimesGroup5}, ... 
                          {phaseSpikeTimesGroup6}, ...
                          {phaseSpikeTimesGroup7}, ...
                          {phaseSpikeTimesGroup8}};
    
end