function BapunRunSpWaveAnaMedianSplit(RatID,metric,pairNo,testType)
    
    pairType      = 'exquisite';
    nullFlag      = false;
    filteredPairs = BapunScreenJitter(RatID,pairType,nullFlag);
    
    if strcmp(RatID,'N')
        
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
        
        chanID_xml = [13  2   14  1   15  0   12  3   11  4   10  5   9   6   8   7   ...
                      29  18  30  17  31  16  28  19  27  20  26  21  25  22  24  23  ...
                      45  34  46  33  47  32  44  35  43  36  42  37  41  38  40  39  ...
                      61  50  62  49  63  48  60  51  59  52  58  53  57  54  56  55  ...
                      77  66  78  65  79  64  76  67  75  68  74  69  73  70  72  71  ...
                      93  82  94  81  95  80  92  83  91  84  90  85  89  86  88  87  ...
                      109 98  110 97  111 96  108 99  107 100 106 101 105 102 104 103 ...
                      125 114 126 113 127 112 124 115 123 116 122 117 121 118 120 199
                      ];  
        
        % hard coding from probe file
        channelSet = {1:16, ...
                      17:32, ... 
                      33:48, ...
                      49:64, ...
                      65:80, ...
                      81:96, ...
                      97:112, ...
                      113:128};

    elseif strcmp(RatID,'S')
        
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');
        
        chanID_xml = [19  3   26  20  5   8   14  0   4   25     ...
                      2   17  16  6   24  21  10  30  22  9   1  ...
                      27  18  31  12  23  28  15  11  29  7   13 ...
                      36  48  58  56  45  32  53  47  33  40  50 ...
                      34  41  52  54  43  37  57  60  49  44  63 ...
                      46  61  59  39  35  62  42  38  55         ...
                      77  66  78  65  79  64  76  67  75  68  74  69  73  70  72  71 ...
                      93  82  94  81  95  80  92  83  91  84  90  85  89  86  88  87 ...
                      109 98  110 97  111 96  108 99  107 100 106 101 105 102 104 103 ...
                      125 114 126 113 127 112 124 115 123 116 122 117 121 118 120 119 ...
                      141 130 142 129 143 128 140 131 139 132 138 133 137 134 136 135 ...
                      157 146 158 145 159 144 156 147 155 148 154 149 153 150 152 151 ...
                      173 162 174 161 175 160 172 163 171 164 170 165 169 166 168 167 ...
                      189 178 190 177 191 176 188 179 187 180 186 181 185 182 184 183];
        
        % hard coding from probe file
        channelSet = {1:10, ...
                      11:21, ...
                      22:32, ... 
                      33:43, ...
                      44:54, ...
                      55:63, ...
                      64:79, ...
                      80:95, ...
                      96:111, ...
                      112:127, ...
                      128:143, ...
                      144:159, ...
                      160:175, ...
                      176:191};
        
    elseif strcmp(RatID,'U')
        
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
        probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');
        
        % hard coding from probe file
        channelSet = {1:16, ...
                      17:32, ... 
                      33:48, ...
                      49:64, ...
                      65:80, ...
                      81:96, ...
                      97:112, ...
                      113:128, ...
                      129:144, ...
                      145:160, ...
                      161:176, ...
                      177:192};
        
    end
    
    matWaveIdx = find(cell2mat(probegroup.connected));
    
    %%

    jscale         = 5;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    fpass          = 300;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    Nwaveforms     = 100;
    filterFlag     = true;
    
    figpath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    preLength  = 36;
    postLength = 36;
    
    lfpDownsampleFactor = 100;

    preLengthLFP  = fs/lfpDownsampleFactor;
    postLengthLFP = fs/lfpDownsampleFactor;

    tLFP = (-(preLengthLFP-1):postLengthLFP)/((fs)/lfpDownsampleFactor);
    
    waveTimeLong   = (-preLength+1:postLength)*(1/30);
    
    if isempty(pairNo)
        % Monitor specific plot settings.
        screensize = get(0,'screensize');
        % initiate figure
        hcomb = figure(102);
    
        res_type = 'QHD';
        pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    end

    %% save variables
    medSplit.pairs       = [];
    medSplit.metric      = metric;
    medSplit.refAboveMed = {};
    medSplit.refBelowMed = {};
    medSplit.tarAboveMed = {};
    medSplit.tarBelowMed = {};
    medSplit.refCh       = [];
    medSplit.tarCh       = [];

    %% using loop for all vs one pairs
    if isempty(pairNo)
        startLoop = 1;
        endLoop   = size(filteredPairs,1);
    else 
        startLoop = pairNo;
        endLoop   = pairNo;
    end 
    
    for i = startLoop:endLoop
        
        % tiledlayout(7,2)
        
        res1    = neurons.spiketrains{1,filteredPairs(i,1)};
        res2    = neurons.spiketrains{1,filteredPairs(i,2)};

        cellID1 = filteredPairs(i,1);
        cellID2 = filteredPairs(i,2);

        ch1     = double(neurons.peak_channels(filteredPairs(i,1))) + 1;
        ch2     = double(neurons.peak_channels(filteredPairs(i,2))) + 1;

        shank1  = double(neurons.shank_ids(filteredPairs(i,1))) + 1;
        shank2  = double(neurons.shank_ids(filteredPairs(i,2))) + 1;
        
        cell1type = neurons.neuron_type(filteredPairs(i,1),:);
        cell2type = neurons.neuron_type(filteredPairs(i,2),:);
        
        % reject MUAs at this stage
        if strcmp(cell1type,'mua  ') || strcmp(cell2type,'mua  ')
            continue
        elseif shank1 == shank2
            continue
        end
        
        if strcmp(cell1type,'pyr  ')
            cell1type = 'p';
        elseif strcmp(cell1type,'inter')
            cell1type = 'i';
        end
        
        if strcmp(cell2type,'pyr  ')
            cell2type = 'p';
        elseif strcmp(cell2type,'inter') 
            cell2type = 'i';
        end

        %%
        
        pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' ...
                   num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
        refStr  = [num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')'];
        tarStr  = [num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
        
        saveStr = ['Rat ' RatID '  (median split - ' metric ') - '  num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' ...
                                                                    num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];       
               
        chanDataCell1 = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/Rat' RatID '/ch' num2str(ch1) 'highpass300hz.mat'],'data');
        spikeTimeIndxCell1 = round(res1*fs) + 1;
        dataTemp = double(chanDataCell1.data);
        clear chanDataCell1
        [spikeAvgMaxChanCell1,waveformsMaxCell1] = waveformAvg(dataTemp,spikeTimeIndxCell1,preLength,postLength,300,fs,false,false);
        dataTemp = downsample(dataTemp,lfpDownsampleFactor);
        [lfpAvgMaxChanCell1,lfpMaxCell1]         = waveformAvg(dataTemp,round(spikeTimeIndxCell1/lfpDownsampleFactor),preLengthLFP,postLengthLFP,300,fs,false,false);
        clear dataTemp

        if strcmp(metric,'symmetry')
            
            [~,troughIdx]  = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            troughIdx      = troughIdx + 14;
            
            rightSide      = waveformsMaxCell1(troughIdx+1:troughIdx+14,:);
            leftSide       = flipud(waveformsMaxCell1(troughIdx-14:troughIdx-1,:));
            clear troughIdx

            squaredError   = sum((rightSide-leftSide).^2);
            medianForSplit = median(squaredError);
        
            aboveMedBin    = squaredError > medianForSplit;
            belowMedBin    = squaredError < medianForSplit;
            
            nexttile(1)
            histogram(squaredError)
            xline(medianForSplit,'r','median','LineWidth',1)
        
        elseif strcmp(metric,'troughAmplitude')

            troughMins     = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            % medianForSplit = median(troughMins);
            % 
            % aboveMedBin    = troughMins > medianForSplit;
            % belowMedBin    = troughMins < medianForSplit;

            % use quartiles for median split 
            Q = quantile(troughMins,[0.25, 0.5, 0.75]);

            aboveMedBin    = troughMins > Q(3);
            belowMedBin    = troughMins < Q(1);
            
            surrIdx = [];
            for loopSurr = 1:1000 
                surrIdx = [surrIdx; sort(randperm(length(troughMins),sum(aboveMedBin)))];
            end

            if isempty(pairNo)
                nexttile(1)
                histogram(troughMins)
                xline(medianForSplit,'r','median','LineWidth',1)
            end


        elseif strcmp(metric,'peakToTrough')

            troughMins     = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            peakMax        = max(waveformsMaxCell1(preLength:preLength+14,:));

            peakToTrough   = peakMax - troughMins;
            medianForSplit = median(peakToTrough);

            aboveMedBin    = peakToTrough > medianForSplit;
            belowMedBin    = peakToTrough < medianForSplit;
            
            nexttile(1)
            histogram(peakToTrough)
            xline(medianForSplit,'r','median','LineWidth',1)

        elseif strcmp(metric,'halfWidth')

            troughMins     = min(waveformsMaxCell1(preLength-14:preLength+14,:));

            peakToBaseline        = mean(waveformsMaxCell1(1:preLength,:)) - troughMins;
            peakToBaselineHalf    = peakToBaseline/2;
            peakToBaselineHalfAdj = peakToBaselineHalf + troughMins;

            [~,idxRight] = min(abs(waveformsMaxCell1(preLength:preLength+14,:) - peakToBaselineHalfAdj));
            [~,idxLeft]  = min(abs(waveformsMaxCell1(preLength-14-1:preLength-1,:) - peakToBaselineHalfAdj));

            idxRight = idxRight + preLength;
            idxLeft  = idxLeft  + (preLength-14);

            timeHalfWidthRight = waveTimeLong(idxRight);
            timeHalfWidthLeft  = waveTimeLong(idxLeft);

            halfWidth = timeHalfWidthRight - timeHalfWidthLeft;

            % medianForSplit = median(halfWidth);
            % 
            % aboveMedBin    = halfWidth > medianForSplit;
            % belowMedBin    = halfWidth < medianForSplit;

            Q = quantile(halfWidth,[0.25, 0.5, 0.75]);

            aboveMedBin    = halfWidth > Q(3);
            belowMedBin    = halfWidth < Q(1);

            surrIdx = [];
            for loopSurr = 1:1000 
                surrIdx = [surrIdx; sort(randperm(length(halfWidth),sum(aboveMedBin)))];
            end
            
            if isempty(pairNo)
                nexttile(1)
                histogram(halfWidth)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        elseif strcmp(metric,'spikeWidth')

            [~,idxTrough]     = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            [~,idxPeak]       = max(waveformsMaxCell1(preLength:preLength+14,:));

            timeTrough = waveTimeLong(idxTrough + preLength-14);
            timePeak   = waveTimeLong(idxPeak   + preLength);

            spikeWidth = timePeak - timeTrough;

            medianForSplit = median(spikeWidth);

            aboveMedBin    = spikeWidth > medianForSplit;
            belowMedBin    = spikeWidth < medianForSplit;
                       
            nexttile(1)
            histogram(spikeWidth)
            xline(medianForSplit,'r','median','LineWidth',1)

        end

        if isempty(pairNo)
            title(['eSpike: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ') metric: ' metric])
        end

        resAboveMed  = res1(aboveMedBin);
        resBelowMed  = res1(belowMedBin);
        
        medSplit.refAboveMed{i} = resAboveMed;
        medSplit.refBelowMed{i} = resBelowMed;
        medSplit.refCh          = [medSplit.refCh ch1];
        medSplit.tarCh          = [medSplit.tarCh ch2];
    
        waveAboveMed = waveformsMaxCell1(:,aboveMedBin);
        waveBelowMed = waveformsMaxCell1(:,belowMedBin);

        lfpAboveMed  = lfpMaxCell1(:,aboveMedBin(1:size(lfpMaxCell1,2)));
        lfpBelowMed  = lfpMaxCell1(:,belowMedBin(1:size(lfpMaxCell1,2)));
        
        if (size(waveAboveMed,2) < Nwaveforms) || (size(waveBelowMed,2) < Nwaveforms)
            continue
        else
            waveAboveMedRandIdx = randperm(size(waveAboveMed,2),Nwaveforms);
            waveBelowMedRandIdx = randperm(size(waveBelowMed,2),Nwaveforms);
        end
        
        % hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        % hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRabove,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resAboveMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        % hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRbelow,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resBelowMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);
        
        for loopSurr = 1:1000 
            [ccgRtemp, tR] = CCG([res1(surrIdx(loopSurr,:))';res2'],[ones(size(res1(surrIdx(loopSurr,:))))';2*ones(size(res2))'], ...
                    'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                    'norm', 'counts');

            ccgj(:,loopSurr) = ccgRtemp(:,1,2);
        end

        if strcmp(testType,'noMinus')

            [aPoint, bPoint, aSim, bSim, ~, ~] = ... 
                computeBandsV2(ccgR(:,1,2),ccgj,tR);

        elseif strcmp(testType,'minus')
        
            ccgRdiff = ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2))-ccgRabove(:,1,2)/sum(ccgRabove(:,1,2));
        
            ccgjDiff = [];
            for loopSurr = 1:500
    
                ccgjDiff = [ccgjDiff ccgj(:,loopSurr)/sum(ccgj(:,loopSurr)) - ccgj(:,loopSurr + 500)/sum(ccgj(:,loopSurr + 500))];
    
            end
        
            [aPoint, bPoint, aSim, bSim, ~, ~] = ... 
                computeBandsV2(ccgRdiff,ccgjDiff,tR);

        end 

        tRms = tR*1000;
        
        if isempty(pairNo)
            nexttile(3)
            patchline(waveTimeLong,waveAboveMed(:,waveAboveMedRandIdx(1)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(waveTimeLong,waveAboveMed(:,waveAboveMedRandIdx(k)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['above median (n = ' num2str(size(waveAboveMed,2)) ')'])
            
            nexttile(5)
            patchline(waveTimeLong,waveBelowMed(:,waveBelowMedRandIdx(1)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(waveTimeLong,waveBelowMed(:,waveBelowMedRandIdx(k)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
            
            nexttile(7)
            plot(waveTimeLong,spikeAvgMaxChanCell1,'LineWidth',1)
            hold on
            plot(waveTimeLong,mean(waveAboveMed,2),'LineWidth',1);
            plot(waveTimeLong,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            set(gca,'YDir','reverse')
            legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
                   ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
                   ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
    
            nexttile(9)
    %         plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
            plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'LineWidth',1)
            hold on
    %         plot(tR*1000,ccgRabove(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
    %         plot(tR*1000,ccgRbelow(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
            plot(tR*1000,ccgRabove(:,1,2)/sum(ccgRabove(:,1,2)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2)),'LineWidth',1)
            hold off
            title(['normalized CCG:' pairStr])
            legend('all','above median','below median')
            
            nexttile(11)
    %         plot(tR*1000,ccgR(:,1,1)/2,'LineWidth',1)
            plot(tR*1000,ccgR(:,1,1)/sum(ccgR(:,1,1)),'LineWidth',1)
            hold on
    %         plot(tR*1000,ccgRabove(:,1,1)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
    %         plot(tR*1000,ccgRbelow(:,1,1)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
            plot(tR*1000,ccgRabove(:,1,1)/sum(ccgRabove(:,1,1)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,1)/sum(ccgRbelow(:,1,1)),'LineWidth',1)
            hold off
            title(['normalized ACG: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')'])
            legend('all','above median','below median')
            
            nexttile(13)
            plot(tLFP,lfpAvgMaxChanCell1,'LineWidth',1)
            hold on
            plot(tLFP,mean(lfpAboveMed,2),'LineWidth',1);
            plot(tLFP,mean(lfpBelowMed,2),'LineWidth',1);
            hold off
            set(gca,'YDir','reverse')
            title('spike-triggered LFP')
            legend('all','above median','below median')

        else 
            if (pairNo == 86) & strcmp(RatID,'S')
                nexttile(1)
            elseif pairNo == 90 & strcmp(RatID,'S')
                nexttile(3)
            elseif pairNo == 18 & strcmp(RatID,'U')
                nexttile(9)
            elseif pairNo == 40 & strcmp(RatID,'U')
                nexttile(11)
            end
            plot(waveTimeLong,spikeAvgMaxChanCell1,'k','LineWidth',1)
            hold on
            plot(waveTimeLong,mean(waveAboveMed,2),'LineWidth',1);
            plot(waveTimeLong,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            xlim([-1 1])
            ylabel('[mV]')
            xlabel('[ms]')
            title({['Rat ' RatID ' waveform: '],refStr})
            set(gca,'YDir','reverse')
            % legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
            %        ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
            %        ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
            if pairNo == 20
                % legend('all','> med','< med')
                if strcmp(testType,'noMinus')
                    legend('all','Q_4','Q_1')
                elseif strcmp(testType,'minus')
                    legend(Q_4 - Q_1')
                end 
            end
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')
            
            if (pairNo == 86) & strcmp(RatID,'S')
                nexttile(5)
            elseif pairNo == 90 & strcmp(RatID,'S')
                nexttile(7)
            elseif pairNo == 18 & strcmp(RatID,'U')
                nexttile(13)
            elseif pairNo == 40 & strcmp(RatID,'U')
                nexttile(15)
            end

            if strcmp(testType,'noMinus')
                %         plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
                % patch([tRms; flipud(tRms)],[aSim/sum(mean(ccgj,2));   flipud(bSim/sum(mean(ccgj,2)))],   [200 200 200]/255,'EdgeAlpha',0); % simultaneous bands
                patch([tRms; flipud(tRms)],[aPoint/sum(mean(ccgj,2)); flipud(bPoint/sum(mean(ccgj,2)))], [200 200 200]/255,'EdgeAlpha',0); % point-wise bands
                
                hold on 
                plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'k','LineWidth',1)
        %         plot(tR*1000,ccgRabove(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
        %         plot(tR*1000,ccgRbelow(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
                plot(tR*1000,ccgRabove(:,1,2)/sum(ccgRabove(:,1,2)),'LineWidth',1)
                plot(tR*1000,ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2)),'LineWidth',1)
                hold off
            elseif strcmp(testType,'minus')
                % patch([tRms; flipud(tRms)],[aSim;   flipud(bSim)],   [200 200 200]/255,'EdgeAlpha',0); % simultaneous bands
                patch([tRms; flipud(tRms)],[aPoint; flipud(bPoint)], [200 200 200]/255,'EdgeAlpha',0); % point-wise bands

                hold on
                plot(tR*1000,ccgRdiff,'k','LineWidth',1)
                hold off
            end 

            xlim([-1 1])
            ylabel('Spike Probability')
            xlabel('[ms]')
            title({['CCG: ' refStr],tarStr,newline})
            % legend('all','above median','below median')
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')

        end

        clear waveformsMaxCell1
        clear lfpMaxCell1

        %%
        
        res2    = neurons.spiketrains{1,filteredPairs(i,1)};
        res1    = neurons.spiketrains{1,filteredPairs(i,2)};

        cellID2 = filteredPairs(i,1);
        cellID1 = filteredPairs(i,2);

        ch2     = double(neurons.peak_channels(filteredPairs(i,1))) + 1;
        ch1     = double(neurons.peak_channels(filteredPairs(i,2))) + 1;

        shank2  = double(neurons.shank_ids(filteredPairs(i,1))) + 1;
        shank1  = double(neurons.shank_ids(filteredPairs(i,2))) + 1;
        
        cell2type = neurons.neuron_type(filteredPairs(i,1),:);
        cell1type = neurons.neuron_type(filteredPairs(i,2),:);
        
       % reject MUAs at this stage
        if strcmp(cell1type,'mua  ') || strcmp(cell2type,'mua  ')
            continue
        elseif shank1 == shank2
            continue
        end
        
        if strcmp(cell1type,'pyr  ')
            cell1type = 'p';
        elseif strcmp(cell1type,'inter')
            cell1type = 'i';
        end
        
        if strcmp(cell2type,'pyr  ')
            cell2type = 'p';
        elseif strcmp(cell2type,'inter') 
            cell2type = 'i';
        end
        
        pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' ...
                   num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
        refStr  = [num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')'];
        tarStr  = [num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
        
        chanDataCell1 = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/Rat' RatID '/ch' num2str(ch1) 'highpass300hz.mat'],'data');
        spikeTimeIndxCell1 = round(res1*fs) + 1;
        dataTemp = double(chanDataCell1.data);
        clear chanDataCell1
        [spikeAvgMaxChanCell1,waveformsMaxCell1] = waveformAvg(dataTemp,spikeTimeIndxCell1,preLength,postLength,300,fs,false,false);
        dataTemp = downsample(dataTemp,lfpDownsampleFactor);
        [lfpAvgMaxChanCell1,lfpMaxCell1]         = waveformAvg(dataTemp,round(spikeTimeIndxCell1/lfpDownsampleFactor),preLengthLFP,postLengthLFP,300,fs,false,false);
        clear dataTemp

        if strcmp(metric,'symmetry')
            
            [~,troughIdx]  = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            troughIdx      = troughIdx + 14;
            
            rightSide      = waveformsMaxCell1(troughIdx+1:troughIdx+14,:);
            leftSide       = flipud(waveformsMaxCell1(troughIdx-14:troughIdx-1,:));
            clear troughIdx

            squaredError   = sum((rightSide-leftSide).^2);
            medianForSplit = median(squaredError);
        
            aboveMedBin    = squaredError > medianForSplit;
            belowMedBin    = squaredError < medianForSplit;
            
            nexttile(2)
            histogram(squaredError)
            xline(medianForSplit,'r','median','LineWidth',1)
            
        elseif strcmp(metric,'troughAmplitude')

            troughMins     = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            % medianForSplit = median(troughMins);
            % 
            % aboveMedBin    = troughMins > medianForSplit;
            % belowMedBin    = troughMins < medianForSplit;
            
            % use quartiles for median split 
            Q = quantile(troughMins,[0.25, 0.5, 0.75]);

            aboveMedBin    = troughMins > Q(3);
            belowMedBin    = troughMins < Q(1);

            surrIdx = [];
            for loopSurr = 1:1000 
                surrIdx = [surrIdx; sort(randperm(length(troughMins),sum(aboveMedBin)))];
            end

        elseif strcmp(metric,'peakToTrough')

            troughMins     = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            peakMax        = max(waveformsMaxCell1(preLength:preLength+14,:));

            peakToTrough   = peakMax - troughMins;
            medianForSplit = median(peakToTrough);

            aboveMedBin    = peakToTrough > medianForSplit;
            belowMedBin    = peakToTrough < medianForSplit;
            
            nexttile(2)
            histogram(peakToTrough)
            xline(medianForSplit,'r','median','LineWidth',1)

        elseif strcmp(metric,'halfWidth')

            troughMins     = min(waveformsMaxCell1(preLength-14:preLength+14,:));

            peakToBaseline        = mean(waveformsMaxCell1(1:preLength,:)) - troughMins;
            peakToBaselineHalf    = peakToBaseline/2;
            peakToBaselineHalfAdj = peakToBaselineHalf + troughMins;

            [~,idxRight] = min(abs(waveformsMaxCell1(preLength:preLength+14,:) - peakToBaselineHalfAdj));
            [~,idxLeft]  = min(abs(waveformsMaxCell1(preLength-14-1:preLength-1,:) - peakToBaselineHalfAdj));

            idxRight = idxRight + preLength;
            idxLeft  = idxLeft  + (preLength-14);

            timeHalfWidthRight = waveTimeLong(idxRight);
            timeHalfWidthLeft  = waveTimeLong(idxLeft);

            halfWidth = timeHalfWidthRight - timeHalfWidthLeft;

            % medianForSplit = median(halfWidth);
            % 
            % aboveMedBin    = halfWidth > medianForSplit;
            % belowMedBin    = halfWidth < medianForSplit;
    
            % use quartiles for median split 
            Q = quantile(halfWidth,[0.25, 0.5, 0.75]);

            aboveMedBin    = halfWidth > Q(3);
            belowMedBin    = halfWidth < Q(1);

            surrIdx = [];
            for loopSurr = 1:1000 
                surrIdx = [surrIdx; sort(randperm(length(halfWidth),sum(aboveMedBin)))];
            end
            
            if isempty(pairNo)
                nexttile(1)
                histogram(halfWidth)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        elseif strcmp(metric,'spikeWidth')

            [~,idxTrough]     = min(waveformsMaxCell1(preLength-14:preLength+14,:));
            [~,idxPeak]       = max(waveformsMaxCell1(preLength:preLength+14,:));

            timeTrough = waveTimeLong(idxTrough + preLength-14);
            timePeak   = waveTimeLong(idxPeak   + preLength);

            spikeWidth = timePeak - timeTrough;

            medianForSplit = median(spikeWidth);

            aboveMedBin    = spikeWidth > medianForSplit;
            belowMedBin    = spikeWidth < medianForSplit;
    
            nexttile(2)
            histogram(spikeWidth)
            xline(medianForSplit,'r','median','LineWidth',1)

        end
        
        if isempty(pairNo)
            title(['eSpike: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ') metric: ' metric])
        end

        resAboveMed  = res1(aboveMedBin);
        resBelowMed  = res1(belowMedBin);
        
        medSplit.tarAboveMed{i} = resAboveMed;
        medSplit.tarBelowMed{i} = resBelowMed;

        waveAboveMed = waveformsMaxCell1(:,aboveMedBin);
        waveBelowMed = waveformsMaxCell1(:,belowMedBin);

        lfpAboveMed  = lfpMaxCell1(:,aboveMedBin(1:size(lfpMaxCell1,2)));
        lfpBelowMed  = lfpMaxCell1(:,belowMedBin(1:size(lfpMaxCell1,2)));
        
        if (size(waveAboveMed,2) < Nwaveforms) || (size(waveBelowMed,2) < Nwaveforms)
            continue
        else
            waveAboveMedRandIdx = randperm(size(waveAboveMed,2),Nwaveforms);
            waveBelowMedRandIdx = randperm(size(waveBelowMed,2),Nwaveforms);
        end
        
        % hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        % hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRabove,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resAboveMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        % hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRbelow,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resBelowMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        for loopSurr = 1:1000 
            [ccgRtemp, tR] = CCG([res1(surrIdx(loopSurr,:))';res2'],[ones(size(res1(surrIdx(loopSurr,:))))';2*ones(size(res2))'], ...
                    'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                    'norm', 'counts');

            ccgj(:,loopSurr) = ccgRtemp(:,1,2);
        end
        
        if strcmp(testType,'noMinus')

            [aPoint, bPoint, aSim, bSim, ~, ~] = ... 
                computeBandsV2(ccgR(:,1,2),ccgj,tR);

        elseif strcmp(testType,'minus')
            
            ccgRdiff = ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2))-ccgRabove(:,1,2)/sum(ccgRabove(:,1,2));
            
            ccgjDiff = [];
            for loopSurr = 1:500
    
                ccgjDiff = [ccgjDiff ccgj(:,loopSurr)/sum(ccgj(:,loopSurr)) - ccgj(:,loopSurr + 500)/sum(ccgj(:,loopSurr + 500))];
    
            end
            
            [aPoint, bPoint, aSim, bSim, ~, ~] = ... 
                computeBandsV2(ccgRdiff,ccgjDiff,tR);

        end 

        tRms = tR*1000;
        
        if isempty(pairNo)

            nexttile(4)
            patchline(waveTimeLong,waveAboveMed(:,waveAboveMedRandIdx(1)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(waveTimeLong,waveAboveMed(:,waveAboveMedRandIdx(k)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['above median (n = ' num2str(size(waveAboveMed,2)) ')'])
            
            nexttile(6)
            patchline(waveTimeLong,waveBelowMed(:,waveBelowMedRandIdx(1)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(waveTimeLong,waveBelowMed(:,waveBelowMedRandIdx(k)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
                        
            nexttile(8)
            plot(waveTimeLong,spikeAvgMaxChanCell1,'LineWidth',1)
            hold on
            plot(waveTimeLong,mean(waveAboveMed,2),'LineWidth',1);
            plot(waveTimeLong,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            set(gca,'YDir','reverse')
            legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
                   ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
                   ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
    
            nexttile(10)
    %         plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
            plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'LineWidth',1)
            hold on
    %         plot(tR*1000,ccgRabove(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
    %         plot(tR*1000,ccgRbelow(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
            plot(tR*1000,ccgRabove(:,1,2)/sum(ccgRabove(:,1,2)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2)),'LineWidth',1)
            hold off
            title(['normalized CCG:' pairStr])
            legend('all','above median','below median')
    
            nexttile(12)
    %         plot(tR*1000,ccgR(:,1,1)/2,'LineWidth',1)
            plot(tR*1000,ccgR(:,1,1)/sum(ccgR(:,1,1)),'LineWidth',1)
            hold on
    %         plot(tR*1000,ccgRabove(:,1,1)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
    %         plot(tR*1000,ccgRbelow(:,1,1)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
            plot(tR*1000,ccgRabove(:,1,1)/sum(ccgRabove(:,1,1)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,1)/sum(ccgRbelow(:,1,1)),'LineWidth',1)
            hold off
            title(['normalized ACG: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')'])
            legend('all','above median','below median')
    
            nexttile(12)
            plot(tR*1000,ccgR(:,1,1)/sum(ccgR(:,1,1)),'LineWidth',1)
            hold on
            plot(tR*1000,ccgRabove(:,1,1)/sum(ccgRabove(:,1,1)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,1)/sum(ccgRbelow(:,1,1)),'LineWidth',1)
            hold off
    
            nexttile(14)
            plot(tLFP,lfpAvgMaxChanCell1,'LineWidth',1)
            hold on
            plot(tLFP,mean(lfpAboveMed,2),'LineWidth',1);
            plot(tLFP,mean(lfpBelowMed,2),'LineWidth',1);
            hold off
            set(gca,'YDir','reverse')
            title('spike-triggered LFP')
            legend('all','above median','below median')

        else 
            
            if (pairNo == 86) & strcmp(RatID,'S')
                nexttile(2)
            elseif pairNo == 90 & strcmp(RatID,'S')
                nexttile(4)
            elseif pairNo == 18 & strcmp(RatID,'U')
                nexttile(10)
            elseif pairNo == 40 & strcmp(RatID,'U')
                nexttile(12)
            end
            plot(waveTimeLong,spikeAvgMaxChanCell1,'k','LineWidth',1)
            hold on
            plot(waveTimeLong,mean(waveAboveMed,2),'LineWidth',1);
            plot(waveTimeLong,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            xlim([-1 1])
            ylabel('[mV]')
            xlabel('[ms]')
            title({['Rat ' RatID ' waveform: '],refStr})
            set(gca,'YDir','reverse')
            % legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
            %        ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
            %        ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')
            
            if (pairNo == 86) & strcmp(RatID,'S')
                nexttile(6)
            elseif pairNo == 90 & strcmp(RatID,'S')
                nexttile(8)
            elseif pairNo == 18 & strcmp(RatID,'U')
                nexttile(14)
            elseif pairNo == 40 & strcmp(RatID,'U')
                nexttile(16)
            end

            if strcmp(testType,'noMinus')
                % plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
                % patch([tRms; flipud(tRms)],[aSim/sum(mean(ccgj,2));   flipud(bSim/sum(mean(ccgj,2)))],   [200 200 200]/255,'EdgeAlpha',0); % simultaneous bands
                patch([tRms; flipud(tRms)],[aPoint/sum(mean(ccgj,2)); flipud(bPoint/sum(mean(ccgj,2)))], [200 200 200]/255,'EdgeAlpha',0); % point-wise bands
                
                hold on 
                plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'k','LineWidth',1)
        %         plot(tR*1000,ccgRabove(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
        %         plot(tR*1000,ccgRbelow(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
                plot(tR*1000,ccgRabove(:,1,2)/sum(ccgRabove(:,1,2)),'LineWidth',1)
                plot(tR*1000,ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2)),'LineWidth',1)
                hold off
            elseif strcmp(testType,'minus')
                % patch([tRms; flipud(tRms)],[aSim;   flipud(bSim)],   [200 200 200]/255,'EdgeAlpha',0); % simultaneous bands
                patch([tRms; flipud(tRms)],[aPoint; flipud(bPoint)], [200 200 200]/255,'EdgeAlpha',0); % point-wise bands
                hold on
                plot(tR*1000,ccgRdiff,'k','LineWidth',1)
                hold off
            end 

            xlim([-1 1])
            ylabel('Spike Probability')
            xlabel('[ms]')
            title({['CCG: ' refStr],tarStr,newline})
            % legend('all','above median','below median')
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')

        end

        clear waveformsMaxCell1
        clear lfpMaxCell1
        
        %%
        if isempty(pairNo)
            save_file = fullfile(figpath, saveStr);
            print(fig_use, save_file,'-dpng',resolution_use);
    
            % reset plot
            close 102
    
            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));        
        end
    end
    
end