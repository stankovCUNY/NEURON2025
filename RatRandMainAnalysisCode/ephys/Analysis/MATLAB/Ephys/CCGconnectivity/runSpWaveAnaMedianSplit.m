function runSpWaveAnaMedianSplit(metric,state,pairNo)

    % close all
    
    %% metric option:
    %
    % 'troughAmplitude'
    % 'peakToTrough'
    % 'halfWidth'
    % 'spikeWidth'

    %% state options:
    % 
    % 'allStates'
    % 'spindle'
    % 'ripple'
    % 'mobility'
    % 'theta'
    % 'gamma'
    % 'noTheta'
    % 'noGamma'
    % 'noThetaAndGamma'
    % 'noStates'

    % condition markers
    [spindleMaze,rippleMaze,mobileMaze,thetaMaze,gammaMaze] = makeBXandFreqFlags;


    %% hardcoding pairs of interest for waveform analysis

    % constants
    session_name   = 'RoyMaze1';
    % session_name   = 'KevinMaze1';
    conn_type      = {'GapPairs','ExcPairs','InhPairs'};
    jscale         = 1;
    alpha_name     = 5;
%     duration       = CCGwidth;
%     duration       = 0.002;
    duration       = 0.007;
    fpass          = 300;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    filterFlag     = false;
    Nwaveforms     = 100;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using
    fs             = 30000;
    binSize        = 1/fs;
    
    % -1.167 to 1.2ms
%     preLength  = 36;
%     postLength = 36;
    metricParam = 14;
    
    % -3.47 to 3.5ms
    preLength  = 150;
    postLength = 150;
    
    spikeTime = (1/30)*(-preLength+1:postLength);

    if strcmp(session_name,'RoyMaze1')
%         pairsAna = [16 38]; % only monosynapse
        pairsAna = [3   44;             20  45;             20  79;             34  79;
                    79  91;             79  118;            91  118;            44  79;
                    91  109;            73  91;             119 91;             45  4;
                    45  19;             79  38;             3   45;             4   79;
                    18  45;             32  79;             79  38;             3   38;
                    91  101;            45  12;             82  91;             3   25;
                    70  91;             20  34;             20  44;             40  79;
                    79  58;             3   32;             4   25;             20  38;
                    46  79;             79  87;             87  119;            98  119;
                    119 122;            45  15;             17  45;             79  19;
                    91  68;             79  95;             91  110;            79  45];
        % pairsAna = [pairsAna; fliplr(pairsAna)];

    elseif strcmp(session_name,'KevinMaze1')

        pairsAna = [70  22;             70  54;             70  50;             3   22];

    end
    
    maxSimilarityChs = [
                         8    10
                         8    16
                         8    24
                         9    24
                        24    25
                        24    40
                        25    39
                        10    24
                        25    35
                        17    25
                        34    25
                        16     8
                        16     8
                        24    11
                         7    16
                         8    24
                         6    16
                        16    24
                        24    11
                         8    10
                        25    40
                        16     8
                        17    25
                         8    15
                        17    25
                         7     9
                         8     9
                        11    24
                        24    15
                         7     9
                         8    14
                         8    11
                         9    24
                        24    25
                        25    34
                        25    34
                        34    48
                        16     6
                         7    16
                        24     8
                        25    24
                        24    32
                        25    33
                        24    16];

    shankChanList ={1:8  , ...
                    9:16 , ...
                    17:24, ...
                    25:32, ...
                    33:40, ...
                    41:48, ...
                    49:56, ...
                    57:64};
    
    if isempty(pairNo)
        % Monitor specific plot settings.
        screensize = get(0,'screensize');
        % initiate figure
        hcomb = figure(102);
    
        res_type = 'QHD';
        pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    end

    % pair data 
    datapath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
    figpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/figures/';

    [jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

    % get all pairs
    cell_pairs = {jit_var_out.cell_pair}';

    % convert to matrix
    cell_pairs_mat = [];
    for i = 1:size(cell_pairs,1)
        cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
    end

    %% get session maze time markers
    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    [data_dir, name, ~] = fileparts(spike_data_fullpath);
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-spikes.mat'));

    %% get channels from waveform files. Not computationally savvy but it's fastest solution
    extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';

    if strcmp(session_name,'RoyMaze1')
        Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);
    elseif strcmp(session_name,'KevinMaze1')
        Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Kevin-maze1'],'noPrompts',true);
    end
    
    %% save variables
    medSplit.pairs  = [];
    medSplit.metric = metric;
    medSplit.refAboveMed = {};
    medSplit.refBelowMed = {};
    medSplit.tarAboveMed = {};
    medSplit.tarBelowMed = {};
    medSplit.refCh       = [];
    medSplit.tarCh       = [];
    
    %%

    waveDataSpCounts = [];
    for i = 1:size(Session.times,2)
        waveDataSpCounts(i) = size(Session.times{i},1);
    end
    
    %% using loop for all vs one pairs
    if isempty(pairNo)
        startLoop = 1;
        endLoop   = size(filteredPairs,1);
    else 
        startLoop = pairNo;
        endLoop   = pairNo;
    end 

    for loopPairs = startLoop:endLoop
        
        if isempty(pairNo)
            tiledlayout(6,2)
        end

        idx = loopPairs; % pair to run, right now we are working with pre-maze time period

        current_pair = pairsAna(idx,:);

        current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));
        
        medSplit.pairs = [medSplit.pairs; current_pair];
        %%

        n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
        n2shank = jit_var_out(current_pair_indices(2)).cell2shank;

        numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

        numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

        waveIdxCell1 = find(numSpikesCell1 == waveDataSpCounts);
        waveIdxCell2 = find(numSpikesCell2 == waveDataSpCounts);

%         chCell1 = Session.maxWaveformCh(waveIdxCell1) + 1;
%         chCell2 = Session.maxWaveformCh(waveIdxCell2) + 1;
        
        cellID1 = current_pair(1);
        cellID2 = current_pair(2);

        load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')
        
        refIdx = find(unitsGlobalSnapshots.neuronID == cellID1);
        tarIdx = find(unitsGlobalSnapshots.neuronID == cellID2);

        refWave = unitsGlobalSnapshots.waveforms{refIdx};
        tarWave = unitsGlobalSnapshots.waveforms{tarIdx};

        [~,chCell1] = min(min(refWave'));
        [~,chCell2] = min(min(tarWave'));
        
%         chCell1 = maxSimilarityChs(loopPairs,1);
%         chCell2 = maxSimilarityChs(loopPairs,2);

        cell1type = jit_var_out(current_pair_indices(2)).cell1type;
        cell2type = jit_var_out(current_pair_indices(2)).cell2type;

        % plot CCG 
        res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        if strcmp(state,'spindle')

            res1 = filterSpikes(res1,spindleMaze);
            res2 = filterSpikes(res2,spindleMaze);

        elseif strcmp(state,'ripple')

            res1 = filterSpikes(res1,rippleMaze);
            res2 = filterSpikes(res2,rippleMaze);

        elseif strcmp(state,'mobility')

            res1 = filterSpikes(res1,mobileMaze);
            res2 = filterSpikes(res2,mobileMaze);

        elseif strcmp(state,'theta')

            res1 = filterSpikes(res1,thetaMaze);
            res2 = filterSpikes(res2,thetaMaze);

        elseif strcmp(state,'gamma')

            res1 = filterSpikes(res1,gammaMaze);
            res2 = filterSpikes(res2,gammaMaze);
        
        elseif strcmp(state,'noRipple')
                        
            res1Ripple = filterSpikes(res1,rippleMaze);
            res1       = setdiff(res1,res1Ripple);

            res2Ripple = filterSpikes(res2,rippleMaze);
            res2       = setdiff(res2,res2Ripple);
            
        elseif strcmp(state,'noTheta')

            res1Theta = filterSpikes(res1,thetaMaze);
            res1      = setdiff(res1,res1Theta);

            res2Theta = filterSpikes(res2,thetaMaze);
            res2      = setdiff(res2,res2Theta);

        elseif strcmp(state,'noGamma')

            res1Gamma = filterSpikes(res1,gammaMaze);
            res1      = setdiff(res1,res1Gamma);

            res2Gamma = filterSpikes(res2,gammaMaze);
            res2      = setdiff(res2,res2Gamma);
        
        elseif strcmp(state,'noThetaAndRipple')
            
            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Ripple   = filterSpikes(res1,rippleMaze);

            res1ThetaAndRipple = [res1Theta; res1Ripple];
            res1ThetaAndRipple = sort(unique(res1ThetaAndRipple));
            res1               = setdiff(res1,res1ThetaAndRipple);

            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Ripple   = filterSpikes(res2,rippleMaze);

            res2ThetaAndRipple = [res2Theta; res2Ripple];
            res2ThetaAndRipple = sort(unique(res2ThetaAndRipple));
            res2               = setdiff(res2,res2ThetaAndRipple);
            
        elseif strcmp(state,'noThetaAndGamma')

            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1ThetaAndRipple = [res1Theta; res1Gamma];
            res1ThetaAndRipple = sort(unique(res1ThetaAndRipple));
            res1               = setdiff(res1,res1ThetaAndRipple);

            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2ThetaAndRipple = [res2Theta; res2Gamma];
            res2ThetaAndRipple = sort(unique(res2ThetaAndRipple));
            res2               = setdiff(res2,res2ThetaAndRipple);

        elseif strcmp(state,'noStates')

            res1Spindle  = filterSpikes(res1,spindleMaze);
            res1Ripple   = filterSpikes(res1,rippleMaze);
            res1Mobility = filterSpikes(res1,mobileMaze);
            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1State = [res1Spindle; res1Ripple; res1Mobility; res1Theta; res1Gamma];
            res1State = sort(unique(res1State));
            res1      = setdiff(res1,res1State);

            res2Spindle  = filterSpikes(res2,spindleMaze);
            res2Ripple   = filterSpikes(res2,rippleMaze);
            res2Mobility = filterSpikes(res2,mobileMaze);
            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2State = [res2Spindle; res2Ripple; res2Mobility; res2Theta; res2Gamma];
            res2State = sort(unique(res2State));
            res2      = setdiff(res2,res2State);

        end
        
        [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,chCell1);
%         [chanDataAtCell1, onsetTime] = HiroLoadRawNSC(session_name,chCell1);
%         [chanDataAtCell2, onsetTime] = HiroLoadRawNSC(session_name,chCell2);

        pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                   num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];
        refStr  = [num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'];
        tarStr  = [num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];

        saveStr = [session_name ' (median split - ' metric ' - condition: ' state ') - '  num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                                                                   num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];


%         res1sync = [];
%         res2sync = [];

    %     for k = 1:length(GSPExc)
    %        
    %         if (GSPExc(k) == 1) || (GSPInh(k) == 1)
    %             res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
    %             res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
    %         end
    %             
    %     end

        % remove sync spikes for waveforms
%         res1nosync = setdiff(res1,res1sync);
%         res2nosync = setdiff(res2,res2sync);

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;

        % randomly selected spikes
%         spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;

        [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag,false);    
    %     [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);

        if strcmp(metric,'troughAmplitude')

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
            for loopSurr = 1:500 
                surrIdx = [surrIdx; sort(randperm(length(troughMins),sum(aboveMedBin)))];
            end

        elseif strcmp(metric,'peakToTrough')

            troughMins     = min(waveformsMaxCell1(preLength-metricParam:preLength+metricParam,:));
            peakMax        = max(waveformsMaxCell1(preLength:preLength+metricParam,:));

            peakToTrough   = peakMax - troughMins;
            medianForSplit = median(peakToTrough);

            aboveMedBin    = peakToTrough > medianForSplit;
            belowMedBin    = peakToTrough < medianForSplit;
            
            if isempty(pairNo)
                nexttile(1)
                histogram(peakToTrough)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        elseif strcmp(metric,'halfWidth')

            troughMins     = min(waveformsMaxCell1(preLength-metricParam:preLength+metricParam,:));
    %         peakMax        = max(waveformsMaxCell1(preLength:preLength+14,:));

            peakToBaseline        = mean(waveformsMaxCell1(1:preLength,:)) - troughMins;
            peakToBaselineHalf    = peakToBaseline/2;
            peakToBaselineHalfAdj = peakToBaselineHalf + troughMins;

            [~,idxRight] = min(abs(waveformsMaxCell1(preLength:preLength+metricParam,:) - peakToBaselineHalfAdj));
            [~,idxLeft]  = min(abs(waveformsMaxCell1(preLength-metricParam-1:preLength-1,:) - peakToBaselineHalfAdj));

            idxRight = idxRight + preLength;
            idxLeft  = idxLeft  + (preLength-metricParam);

            timeHalfWidthRight = spikeTime(idxRight);
            timeHalfWidthLeft  = spikeTime(idxLeft);

            halfWidth = timeHalfWidthRight - timeHalfWidthLeft;

            medianForSplit = median(halfWidth);
    %         topQuantileThr    = quantile(halfWidth,0.75);
    %         bottomQuantileThr = quantile(halfWidth,0.25);

            aboveMedBin    = halfWidth > medianForSplit;
            belowMedBin    = halfWidth < medianForSplit;

    %         aboveMedBin    = halfWidth > topQuantileThr;
    %         belowMedBin    = halfWidth < bottomQuantileThr;
    
            if isempty(pairNo)
                nexttile(1)
                histogram(halfWidth)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        elseif strcmp(metric,'spikeWidth')

            [~,idxTrough]     = min(waveformsMaxCell1(preLength-metricParam:preLength+metricParam,:));
            [~,idxPeak]       = max(waveformsMaxCell1(preLength:preLength+metricParam,:));

            timeTrough = spikeTime(idxTrough + preLength-metricParam);
            timePeak   = spikeTime(idxPeak   + preLength);

            spikeWidth = timePeak - timeTrough;

            medianForSplit = median(spikeWidth);
    %         topQuantileThr    = quantile(spikeWidth,0.75);
    %         bottomQuantileThr = quantile(spikeWidth,0.25);

            aboveMedBin    = spikeWidth > medianForSplit;
            belowMedBin    = spikeWidth < medianForSplit;
                        
    %         aboveMedBin    = spikeWidth > topQuantileThr;
    %         belowMedBin    = spikeWidth < bottomQuantileThr;
    
            if isempty(pairNo)
                nexttile(1)
                histogram(spikeWidth)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        end

        if isempty(pairNo)
            title(['eSpike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ') metric: ' metric ', condition: ' state])
        end

        resAboveMed  = res1(aboveMedBin);
        resBelowMed  = res1(belowMedBin);
        
        medSplit.refAboveMed{loopPairs} = resAboveMed;
        medSplit.refBelowMed{loopPairs} = resBelowMed;
        medSplit.refCh          = [medSplit.refCh chCell1];
        medSplit.tarCh          = [medSplit.tarCh chCell2];
    
        waveAboveMed = waveformsMaxCell1(:,aboveMedBin);
        waveBelowMed = waveformsMaxCell1(:,belowMedBin);
        
        waveAboveMedRandIdx = randperm(size(waveAboveMed,2),Nwaveforms);
        waveBelowMedRandIdx = randperm(size(waveBelowMed,2),Nwaveforms);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRabove,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resAboveMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRbelow,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resBelowMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        for loopSurr = 1:500 
            [ccgRtemp, tR] = CCG([res1(surrIdx(loopSurr,:));res2],[ones(size(res1(surrIdx(loopSurr,:))));2*ones(size(res2))], ...
                    'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                    'norm', 'counts');

            ccgj(:,loopSurr) = ccgRtemp(:,1,2);
        end

        [aPoint, bPoint, aSim, bSim, ~, ~] = ... 
            computeBandsV2(ccgR(:,2,1),ccgj,tR);

        tRms = tR*1000;
        
        if isempty(pairNo)

            nexttile(3)
            patchline(spikeTime,waveAboveMed(:,waveAboveMedRandIdx(1)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(spikeTime,waveAboveMed(:,waveAboveMedRandIdx(k)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['above median (n = ' num2str(size(waveAboveMed,2)) ')'])
            
            nexttile(5)
            patchline(spikeTime,waveBelowMed(:,waveBelowMedRandIdx(1)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(spikeTime,waveBelowMed(:,waveBelowMedRandIdx(k)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
            
            nexttile(7)
            plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',1)
            hold on
            plot(spikeTime,mean(waveAboveMed,2),'LineWidth',1);
            plot(spikeTime,mean(waveBelowMed,2),'LineWidth',1);
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
            title(['normalized ACG: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
            legend('all','above median','below median')

        else

            nexttile(9)
            plot(spikeTime,spikeAvgMaxChanCell1,'k','LineWidth',1)
            hold on
            plot(spikeTime,mean(waveAboveMed,2),'LineWidth',1);
            plot(spikeTime,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            xlim([-1 1])
            ylabel('[mV]')
            xlabel('[ms]')
            title(['avg. spike waveform: ' refStr])
            set(gca,'YDir','reverse')
            % legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
            %        ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
            %        ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')
    
            nexttile(11)
    %         plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
            % patch([tRms; flipud(tRms)],[aSim/sum(mean(ccgj,2));   flipud(bSim/sum(mean(ccgj,2)))],   [200 200 200]/255,'EdgeAlpha',0); % simultaneous bands
            patch([tRms; flipud(tRms)],[aPoint/sum(mean(ccgj,2)); flipud(bPoint/sum(mean(ccgj,2)))], [125 125 125]/255,'EdgeAlpha',0); % point-wise bands
            hold on
            plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'k','LineWidth',1)
    %         plot(tR*1000,ccgRabove(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
    %         plot(tR*1000,ccgRbelow(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
            plot(tR*1000,ccgRabove(:,1,2)/sum(ccgRabove(:,1,2)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2)),'LineWidth',1)
            hold off
            xlim([-1 1])
            ylabel('Spike Probability')
            xlabel('[ms]')
            title({['CCG: ' pairStr],newline})
            % legend('all','above median','below median')
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')

        end

        %%

        n2maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        n1maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        n2shank = jit_var_out(current_pair_indices(2)).cell1shank;
        n1shank = jit_var_out(current_pair_indices(2)).cell2shank;

        numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

        numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

        waveIdxCell2 = find(numSpikesCell1 == waveDataSpCounts);
        waveIdxCell1 = find(numSpikesCell2 == waveDataSpCounts);

%         chCell2 = Session.maxWaveformCh(waveIdxCell1);
%         chCell1 = Session.maxWaveformCh(waveIdxCell2);
        
        cellID2 = current_pair(1);
        cellID1 = current_pair(2);

        load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')
        
        refIdx = find(unitsGlobalSnapshots.neuronID == cellID2);
        tarIdx = find(unitsGlobalSnapshots.neuronID == cellID1);

        refWave = unitsGlobalSnapshots.waveforms{refIdx};
        tarWave = unitsGlobalSnapshots.waveforms{tarIdx};

        [~,chCell2] = min(min(refWave'));
        [~,chCell1] = min(min(tarWave'));

%         chCell1 = maxSimilarityChs(loopPairs,2);
%         chCell2 = maxSimilarityChs(loopPairs,1);

        cell2type = jit_var_out(current_pair_indices(2)).cell1type;
        cell1type = jit_var_out(current_pair_indices(2)).cell2type;

        % plot CCG 
        res2 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        res1 = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        if strcmp(state,'spindle')

            res1 = filterSpikes(res1,spindleMaze);
            res2 = filterSpikes(res2,spindleMaze);

        elseif strcmp(state,'ripple')

            res1 = filterSpikes(res1,rippleMaze);
            res2 = filterSpikes(res2,rippleMaze);

        elseif strcmp(state,'mobility')

            res1 = filterSpikes(res1,mobileMaze);
            res2 = filterSpikes(res2,mobileMaze);

        elseif strcmp(state,'theta')

            res1 = filterSpikes(res1,thetaMaze);
            res2 = filterSpikes(res2,thetaMaze);

        elseif strcmp(state,'gamma')

            res1 = filterSpikes(res1,gammaMaze);
            res2 = filterSpikes(res2,gammaMaze);

        elseif strcmp(state,'noTheta')

            res1Theta = filterSpikes(res1,thetaMaze);
            res1      = setdiff(res1,res1Theta);

            res2Theta = filterSpikes(res2,thetaMaze);
            res2      = setdiff(res2,res2Theta);

        elseif strcmp(state,'noGamma')

            res1Gamma = filterSpikes(res1,gammaMaze);
            res1      = setdiff(res1,res1Gamma);

            res2Gamma = filterSpikes(res2,gammaMaze);
            res2      = setdiff(res2,res2Gamma);

        elseif strcmp(state,'noThetaAndGamma')

            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1ThetaAndRipple = [res1Theta; res1Gamma];
            res1ThetaAndRipple = sort(unique(res1ThetaAndRipple));
            res1      = setdiff(res1,res1ThetaAndRipple);

            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2ThetaAndRipple = [res2Theta; res2Gamma];
            res2ThetaAndRipple = sort(unique(res2ThetaAndRipple));
            res2      = setdiff(res2,res2ThetaAndRipple);

        elseif strcmp(state,'noStates')

            res1Spindle  = filterSpikes(res1,spindleMaze);
            res1Ripple   = filterSpikes(res1,rippleMaze);
            res1Mobility = filterSpikes(res1,mobileMaze);
            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1State = [res1Spindle; res1Ripple; res1Mobility; res1Theta; res1Gamma];
            res1State = sort(unique(res1State));
            res1      = setdiff(res1,res1State);

            res2Spindle  = filterSpikes(res2,spindleMaze);
            res2Ripple   = filterSpikes(res2,rippleMaze);
            res2Mobility = filterSpikes(res2,mobileMaze);
            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2State = [res2Spindle; res2Ripple; res2Mobility; res2Theta; res2Gamma];
            res2State = sort(unique(res2State));
            res2      = setdiff(res2,res2State);

        end
        
        [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,chCell1);
%         [chanDataAtCell1, onsetTime] = HiroLoadRawNSC(session_name,chCell1);
%         [chanDataAtCell2, onsetTime] = HiroLoadRawNSC(session_name,chCell2);

        pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                   num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];
        refStr  = [num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'];
        tarStr  = [num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];

        res1sync = [];
        res2sync = [];

    %     for k = 1:length(GSPExc)
    %        
    %        if (GSPExc(k) == 1) || (GSPInh(k) == 1)
    %            res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
    %            res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
    %        end
    %             
    %     end

        % remove sync spikes for waveforms
%         res1nosync = setdiff(res1,res1sync);
%         res2nosync = setdiff(res2,res2sync);

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;

        % randomly selected spikes
%         spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;

        [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag,false);    
    %     [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);

        if strcmp(metric,'troughAmplitude')

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
            for loopSurr = 1:500 
                surrIdx = [surrIdx; sort(randperm(length(troughMins),sum(aboveMedBin)))];
            end

        elseif strcmp(metric,'peakToTrough')

            troughMins     = min(waveformsMaxCell1(preLength-metricParam:preLength+metricParam,:));
            peakMax        = max(waveformsMaxCell1(preLength:preLength+metricParam,:));

            peakToTrough   = peakMax - troughMins;
            medianForSplit = median(peakToTrough);

            aboveMedBin    = peakToTrough > medianForSplit;
            belowMedBin    = peakToTrough < medianForSplit;
            
            if isempty(pairNo)
                nexttile(2)
                histogram(peakToTrough)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        elseif strcmp(metric,'halfWidth')

            troughMins     = min(waveformsMaxCell1(preLength-metricParam:preLength+metricParam,:));
    %         peakMax        = max(waveformsMaxCell1(preLength:preLength+14,:));

            peakToBaseline        = mean(waveformsMaxCell1(1:preLength,:)) - troughMins;
            peakToBaselineHalf    = peakToBaseline/2;
            peakToBaselineHalfAdj = peakToBaselineHalf + troughMins;

            [~,idxRight] = min(abs(waveformsMaxCell1(preLength:preLength+metricParam,:) - peakToBaselineHalfAdj));
            [~,idxLeft]  = min(abs(waveformsMaxCell1(preLength-metricParam-1:preLength-1,:) - peakToBaselineHalfAdj));

            idxRight = idxRight + preLength;
            idxLeft  = idxLeft  + (preLength-metricParam);

            timeHalfWidthRight = spikeTime(idxRight);
            timeHalfWidthLeft  = spikeTime(idxLeft);

            halfWidth = timeHalfWidthRight - timeHalfWidthLeft;

            medianForSplit = median(halfWidth);
    %         topQuantileThr    = quantile(halfWidth,0.75);
    %         bottomQuantileThr = quantile(halfWidth,0.25);

            aboveMedBin    = halfWidth > medianForSplit;
            belowMedBin    = halfWidth < medianForSplit;

    %         aboveMedBin    = halfWidth > topQuantileThr;
    %         belowMedBin    = halfWidth < bottomQuantileThr;
        
            if isempty(pairNo)
                nexttile(2)
                histogram(halfWidth)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        elseif strcmp(metric,'spikeWidth')

            [~,idxTrough]     = min(waveformsMaxCell1(preLength-metricParam:preLength+metricParam,:));
            [~,idxPeak]       = max(waveformsMaxCell1(preLength:preLength+metricParam,:));

            timeTrough = spikeTime(idxTrough + preLength-metricParam);
            timePeak   = spikeTime(idxPeak   + preLength);

            spikeWidth = timePeak - timeTrough;

            medianForSplit = median(spikeWidth);
    %         topQuantileThr    = quantile(spikeWidth,0.75);
    %         bottomQuantileThr = quantile(spikeWidth,0.25);

            aboveMedBin    = spikeWidth > medianForSplit;
            belowMedBin    = spikeWidth < medianForSplit;

    %         aboveMedBin    = spikeWidth > topQuantileThr;
    %         belowMedBin    = spikeWidth < bottomQuantileThr;
    
            if isempty(pairNo)
                nexttile(2)
                histogram(spikeWidth)
                xline(medianForSplit,'r','median','LineWidth',1)
            end

        end

        if isempty(pairNo)
            title(['eSpike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ') metric: ' metric ', condition: ' state])
        end

        resAboveMed  = res1(aboveMedBin);
        resBelowMed  = res1(belowMedBin);
        
        medSplit.tarAboveMed{loopPairs} = resAboveMed;
        medSplit.tarBelowMed{loopPairs} = resBelowMed;

        waveAboveMed = waveformsMaxCell1(:,aboveMedBin);
        waveBelowMed = waveformsMaxCell1(:,belowMedBin);
        
        waveAboveMedRandIdx = randperm(size(waveAboveMed,2),Nwaveforms);
        waveBelowMedRandIdx = randperm(size(waveBelowMed,2),Nwaveforms);
        
        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRabove,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resAboveMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRbelow,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resBelowMed,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        for loopSurr = 1:500 
            [ccgRtemp, tR] = CCG([res1(surrIdx(loopSurr,:));res2],[ones(size(res1(surrIdx(loopSurr,:))));2*ones(size(res2))], ...
                    'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                    'norm', 'counts');

            ccgj(:,loopSurr) = ccgRtemp(:,1,2);
        end

        [aPoint, bPoint, aSim, bSim, ~, ~] = ... 
            computeBandsV2(ccgR(:,2,1),ccgj,tR);

        tRms = tR*1000;
        
        if isempty(pairNo)

            nexttile(4)
            patchline(spikeTime,waveAboveMed(:,waveAboveMedRandIdx(1)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(spikeTime,waveAboveMed(:,waveAboveMedRandIdx(k)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['above median (n = ' num2str(size(waveAboveMed,2)) ')'])
            
            nexttile(6)
            patchline(spikeTime,waveBelowMed(:,waveBelowMedRandIdx(1)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            hold on
            for k = 2:Nwaveforms
                patchline(spikeTime,waveBelowMed(:,waveBelowMedRandIdx(k)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
            end
            hold off
            set(gca,'YDir','reverse')
            title(['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
                        
            nexttile(8)
            plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',1)
            hold on
            plot(spikeTime,mean(waveAboveMed,2),'LineWidth',1);
            plot(spikeTime,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            set(gca,'YDir','reverse')
            legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
                   ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
                   ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
    
            nexttile(10)
    %         plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
            
            hold on
            plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'LineWidth',1)
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
            title(['normalized ACG: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
            legend('all','above median','below median')

        else

            nexttile(10)
            plot(spikeTime,spikeAvgMaxChanCell1,'k','LineWidth',1)
            hold on
            plot(spikeTime,mean(waveAboveMed,2),'LineWidth',1);
            plot(spikeTime,mean(waveBelowMed,2),'LineWidth',1);
            hold off
            xlim([-1 1])
            ylabel('[mV]')
            xlabel('[ms]')
            title(['avg. spike waveform: ' refStr])
            set(gca,'YDir','reverse')
            % legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
            %        ['above median (n = ' num2str(size(waveAboveMed,2)) ')'], ...
            %        ['below median (n = ' num2str(size(waveBelowMed,2)) ')'])
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')
    
            nexttile(12)
    %         plot(tR*1000,ccgR(:,1,2)/2,'LineWidth',1)
            % patch([tRms; flipud(tRms)],[aSim/sum(mean(ccgj,2));   flipud(bSim/sum(mean(ccgj,2)))],   [200 200 200]/255,'EdgeAlpha',0); % simultaneous bands
            patch([tRms; flipud(tRms)],[aPoint/sum(mean(ccgj,2)); flipud(bPoint/sum(mean(ccgj,2)))], [125 125 125]/255,'EdgeAlpha',0); % point-wise bands
            hold on
            plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'k','LineWidth',1)
    %         plot(tR*1000,ccgRabove(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveAboveMed,2))),'LineWidth',1)
    %         plot(tR*1000,ccgRbelow(:,1,2)*(0.5*(size(waveformsMaxCell1,2)/size(waveBelowMed,2))),'LineWidth',1)
            plot(tR*1000,ccgRabove(:,1,2)/sum(ccgRabove(:,1,2)),'LineWidth',1)
            plot(tR*1000,ccgRbelow(:,1,2)/sum(ccgRbelow(:,1,2)),'LineWidth',1)
            hold off
            xlim([-1 1])
            ylabel('Spike Probability')
            xlabel('[ms]')
            title({['CCG: ' pairStr],newline})
            % legend('all','above median','below median')
            box off
            set(gca,'FontSize',5)
            set(gca,'FontName','Arial')

        end
        
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
    
    % close 102
    % 
    % %% save median split subsets
    % save([datapath 'medSplitSets_' metric '_' state],'medSplit')
    
end

                         
                         
                         