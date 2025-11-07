function aggregateAnalysisLegacy

    clc; 

    %% hardcoding pairs of interest for waveform analysis

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

    pairsAna = [pairsAna; fliplr(pairsAna)];

    % constants
    session_name   = 'RoyMaze1';
    conn_type      = {'GapPairs','ExcPairs','InhPairs'};
    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    fpass          = 300;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    filterFlag     = false;
    Nwaveforms     = 10000;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

    preLength  = 36;
    postLength = 36;

    shankChanList ={1:8  , ...
                    9:16 , ...
                    17:24, ...
                    25:32, ...
                    33:40, ...
                    41:48, ...
                    49:56, ...
                    57:64};
    % probe info
    probe = buzsaki64probeLoc;

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    % pair data 
    datapath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';

    [jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);
    jit_var_out = addFlipped_Jit_var_out(jit_var_out);

    % get all pairs
    cell_pairs = {jit_var_out.cell_pair}';

    % convert to matrix
    cell_pairs_mat = [];
    for i = 1:size(cell_pairs,1)
        cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
    end

    %% group stat vars
    
    groupStat.mainRef             = [];
    groupStat.preRef              = [];
    groupStat.postRef             = [];
    
    groupStat.mainTar             = [];
    groupStat.preTar              = [];
    groupStat.postTar             = [];
    
    groupStat.mainRsqMax          = [];
    groupStat.preRsqMax           = [];
    groupStat.postRsqMax          = [];
    
    groupStat.mainMaxChan         = [];
    groupStat.preMaxChan          = [];
    groupStat.postMaxChan         = [];
    
    groupStat.mainFeatCCGall      = [];
    groupStat.preFeatCCGall       = [];
    groupStat.postFeatCCGall      = [];
    
    groupStat.mainFeatCCGnormAll  = [];
    groupStat.preFeatCCGnormAll   = [];
    groupStat.postFeatCCGnormAll  = [];
    
    groupStat.mainPkAgmonIdx      = [];
    groupStat.preFeatCCGagmonIdx  = [];
    groupStat.postFeatCCGagmonIdx = [];
    
    groupStat.mainFeatWAVEall     = [];
    groupStat.preFeatWAVEall      = [];
    groupStat.postFeatWAVEall     = [];
    
    groupStat.mainFeatCCGlatAll   = [];
    groupStat.preFeatCCGlatAll    = [];
    groupStat.postFeatCCGlatAll   = [];
    
    groupStat.mainFeatWAVElatAll  = [];
    groupStat.preFeatWAVElatAll   = [];
    groupStat.postFeatWAVElatAll  = [];
    
    groupStat.mainDistance        = [];
    groupStat.preDistance         = [];
    groupStat.postDistance        = [];

    %% get session maze time markers
    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    [data_dir, name, ~] = fileparts(spike_data_fullpath);
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');

    %% get channels from waveform files. Not computaitonal;y savvy but it's fastest solution
    extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
    RoySession1 = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);

    tWave = (-26:27)*(1/30);

    %%

    waveDataSpCounts = [];
    for i = 1:size(RoySession1.times,2)
        waveDataSpCounts(i) = size(RoySession1.times{i},1);
    end

    for j = 1:size(pairsAna,1)
        
        tic
        
        idx = j; % pair to run, right now we are working with pre-maze time period

        current_pair = pairsAna(idx,:);

        current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));

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

        chCell1 = RoySession1.maxWaveformCh(waveIdxCell1);
        chCell2 = RoySession1.maxWaveformCh(waveIdxCell2);

        cell1type = jit_var_out(current_pair_indices(2)).cell1type;
        cell2type = jit_var_out(current_pair_indices(2)).cell2type;

        cellID1 = current_pair(1);
        cellID2 = current_pair(2);

        % plot CCG 
        res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;

%         [chanDataAtCell1, onsetTime] = HiroLoadRawNSC(session_name,chCell1);
%         [chanDataAtCell2, onsetTime] = HiroLoadRawNSC(session_name,chCell2);
        
        [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,chCell1);
        [chanDataAtCell2, onsetTime] = HiroLoad300hz(session_name,chCell2);

        spikeTime = (1/30)*(-preLength+1:postLength);

        pairStr = [session_name ' - '  num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                                       num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];


        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);

        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(211,1),duration);

        res1sync = [];
        res2sync = [];

        for k = 1:211

            if (GSPExc(k) == 1) || (GSPInh(k) == 1)
                res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
                res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
            end

        end

        % remove sync spikes for waveforms
        res1nosync = setdiff(res1,res1sync);
        res2nosync = setdiff(res2,res2sync);

        if Nwaveforms < length(res1nosync)
            res1nosync = res1nosync(randperm(length(res1nosync),Nwaveforms));
        end
        if Nwaveforms < length(res2nosync)
            res2nosync = res2nosync(randperm(length(res2nosync),Nwaveforms));
        end

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;

        % randomly selected spikes
        spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;

%         [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
%         [spikeAvgMaxChanCell2, waveformsMaxCell2] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
% 
%         [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
%         [spikeAvgOpsChanCell2, waveformsOpsCell2] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        clear chanDataAtCell1
        clear chanDataAtCell2

        allChansAtCell1 = [shankChanList{n1shank}];
        allChansAtCell2 = [shankChanList{n2shank}];
        
        spikeAvgMaxShankCell1nosync = zeros(8,preLength+postLength);
        spikeAvgMaxShankCell2nosync = zeros(8,preLength+postLength);
        spikeAvgOpsShankCell1nosync = zeros(8,preLength+postLength);
        spikeAvgOpsShankCell2nosync = zeros(8,preLength+postLength);
        
        % home (hugs) shank
        for k = 1:8
            [chanDataAtCell1temp, onsetTime] = HiroLoad300hz(session_name,allChansAtCell1(k));
%             [chanDataAtCell1temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell1(k));
            spikeAvgMaxShankCell1nosync(k,:) = waveformAvg(chanDataAtCell1temp,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);

        end
        clear chanDataAtCell1temp

%         for k = 1:8
% 
%             [chanDataAtCell2temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell2(k));
%             spikeAvgMaxShankCell2nosync(k,:) = waveformAvg(chanDataAtCell2temp,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
% 
%         end
%         clear chanDataAtCell2temp
        
        % opposite shank
        for k = 1:8
            [chanDataAtCell2temp, onsetTime] = HiroLoad300hz(session_name,allChansAtCell2(k));
%             [chanDataAtCell2temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell2(k));
            spikeAvgOpsShankCell1nosync(k,:) = waveformAvg(chanDataAtCell2temp,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);

        end
        clear chanDataAtCell2temp

%         for k = 1:8
% 
%             [chanDataAtCell1temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell1(k));
%             spikeAvgOpsShankCell2nosync(k,:) = waveformAvg(chanDataAtCell1temp,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
% 
%         end
%         clear chanDataAtCell1temp

        ylims = get(gca,'ylim');

    %     if any(GSPExc)
    %         hold on;
    %         plot(tR(GSPExc == 1)*1000, 0.95*ylims(2),'k^');
    %     end
    % 
    %     if any(GSPInh)
    %         hold on;
    %         plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'kv');
    %     end


        %% pearson's R
        ccg = ccgR(:,1,2);

    %     ccgRange  = 106:136; % 1msec (0 to 1ms lag)
    %     waveRange = 36:66;   % 1msec (0 to 1ms lag)

        ccgRange  = 96:136; % 1.25msec (-0.25 to 1ms lag)
        waveRange = 26:66;  % 1.25msec (-0.25 to 1ms lag)

        RsqList = [];
        for k = 1:16
            if k <= 8
                wave4Ana = spikeAvgMaxShankCell1nosync(k,:)';
            elseif k > 8
                wave4Ana = spikeAvgOpsShankCell1nosync(k-8,:)';
            end
            
            R = corrcoef(ccg(ccgRange),wave4Ana(waveRange));
            
            RsqList(k) = R(1,2)^2;
            
        end 

        [RsqMax,maxMatchChIdx] = max(RsqList);
        
        % don't include for effect group
        if RsqMax < 0.5
            
            disp([num2str((j/size(pairsAna,1))*100) ' % done!'])
            toc
            
            continue 
        end
        
        allChans   = [allChansAtCell1 allChansAtCell2];
        maxMatchCh = allChans(maxMatchChIdx);
        
        % distance calc
        xSource   = probe.x(probe.chanNo == chCell1);
        ySource   = probe.y(probe.chanNo == chCell1);
        
        xReceiver = probe.x(probe.chanNo == maxMatchCh);
        yReceiver = probe.y(probe.chanNo == maxMatchCh);
        
        distance  = (((xSource - xReceiver)^2) + ((ySource - yReceiver)^2))^(1/2);
        
        % waveform for analysis based on max channel
        if maxMatchChIdx <= 8
            wave4Ana = spikeAvgMaxShankCell1nosync(maxMatchChIdx,:)';
        elseif maxMatchChIdx > 8
            wave4Ana = spikeAvgOpsShankCell1nosync(maxMatchChIdx-8,:)';
        end
        
        %% find peaks and troughs CCG
        peaksMask     = GSPExc(ccgRange);
        troughsMask   = GSPInh(ccgRange);
        tRanalysis    = tR(ccgRange)*1000;
        ccgAnalysis   = ccg(ccgRange);
        ccgjmAnalysis = ccgjm(ccgRange);
        JBSIEanalysis = JBSIE(ccgRange);
        JBSIIanalysis = JBSII(ccgRange);

        [pkCCG, pklocsCCG] = findpeaks(ccg(ccgRange));
        [trCCG, trlocsCCG] = findpeaks(-ccg(ccgRange));
        trCCG = -trCCG;

        % screen with mask
        pklocsCCGmask = zeros(length(peaksMask),1);
        pklocsCCGmask(pklocsCCG) = 1;

        trlocsCCGmask = zeros(length(troughsMask),1);
        trlocsCCGmask(trlocsCCG) = 1;

        pklocsCCGmask = peaksMask   & pklocsCCGmask;
        trlocsCCGmask = troughsMask & trlocsCCGmask;

        % identify multipeaks
        if (sum(pklocsCCGmask(10:12)) > 0) && (sum(pklocsCCGmask([1:9,13:41])) > 0)
            
            % zero lag peak
            pklocsCCGmask_0lag        = zeros(length(pklocsCCGmask),1);
            pklocsCCGmask_0lag(10:12) = pklocsCCGmask(10:12);
            
            [mainPk, ...
             mainTr, ...
             mainPkLat, ...
             mainTrLat, ...
             preFeatureCCGlat, ...
             preFeatureWAVElat, ...
             postFeatureCCGlat, ...
             postFeatureWAVElat, ...
             preFeatCCGidx, ...
             postFeatCCGidx, ...
             preFeatWAVEidx, ... 
             postFeatWAVEidx, ...
             preFeatCCG, ...
             postFeatCCG, ...
             preFeatWAVE, ... 
             postFeatWAVE, ...
             mainPkNorm, ...
             preFeatCCGnorm, ...
             postFeatCCGnorm, ...
             mainPkAgmonIdx, ...
             preFeatCCGagmonIdx, ...
             postFeatCCGagmonIdx] = featureLatencies(pklocsCCGmask_0lag, ...
                                                     trlocsCCGmask, ...
                                                     ccgAnalysis, ...
                                                     ccgjmAnalysis, ...
                                                     JBSIEanalysis, ...
                                                     JBSIIanalysis, ...
                                                     tRanalysis, ...
                                                     wave4Ana, ...
                                                     waveRange);
            
            %% pairing CCG an waveform features
            if ~isempty(mainPk) && ~isempty(mainTr)
                
                groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
                groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
                groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
                groupStat.mainPkAgmonIdx     = [groupStat.mainPkAgmonIdx mainPkAgmonIdx];
   
                groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
                groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
                
                groupStat.mainDistance = [groupStat.mainDistance distance];
                
                groupStat.mainRef = [groupStat.mainRef cellID1];
                groupStat.mainTar = [groupStat.mainTar cellID2];
                
                groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
                groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];

                if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
                    
                    groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
                    groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
                    groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
                    groupStat.preFeatCCGagmonIdx = [groupStat.preFeatCCGagmonIdx preFeatCCGagmonIdx];
                    
                    groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
                    groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
                     
                    groupStat.preDistance = [groupStat.preDistance distance];
                     
                    groupStat.preRef = [groupStat.preRef cellID1];
                    groupStat.preTar = [groupStat.preTar cellID2];

                    groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
                    groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];

                end

                if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
                    
                    groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
                    groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
                    groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
                    groupStat.postFeatCCGagmonIdx = [groupStat.postFeatCCGagmonIdx postFeatCCGagmonIdx];
                    
                    groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
                    groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
                    
                    groupStat.postDistance = [groupStat.postDistance distance];
                    
                    groupStat.postRef = [groupStat.postRef cellID1];
                    groupStat.postTar = [groupStat.postTar cellID2];

                    groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
                    groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];

                end
                
                yyaxis left
                plot(tR*1000,ccg,'b')
                hold on
                scatter(mainPkLat,mainPk,'b','LineWidth',2)
                scatter(preFeatureCCGlat,preFeatCCG,'b','LineWidth',2)
                scatter(postFeatureCCGlat,postFeatCCG,'b','LineWidth',2)
                hold off

                yyaxis right
                plot(spikeTime,wave4Ana,'r')
                hold on
                scatter(mainTrLat,mainTr,'r','LineWidth',2)
                scatter(preFeatureWAVElat,preFeatWAVE,'r','LineWidth',2)
                scatter(postFeatureWAVElat,postFeatWAVE,'r','LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlim([-3.5 3.5])

                titleChar = ['CCG and wave feature alighment - ' num2str(cellID1) cell1type ' v ' num2str(cellID2) cell2type' '- zero lag'];
                title(titleChar)

                disp([num2str((j/size(pairsAna,1))*100) ' % done!'])

                save_file = fullfile(datapath, titleChar);
                print(fig_use, save_file,'-dpng',resolution_use);
                
            end
            
            % lagged peaks
            pklocsCCGmask_non_0lag              = zeros(length(pklocsCCGmask),1);
            pklocsCCGmask_non_0lag([1:9,13:41]) = pklocsCCGmask([1:9,13:41]);
            
            [mainPk, ...
             mainTr, ...
             mainPkLat, ...
             mainTrLat, ...
             preFeatureCCGlat, ...
             preFeatureWAVElat, ...
             postFeatureCCGlat, ...
             postFeatureWAVElat, ...
             preFeatCCGidx, ...
             postFeatCCGidx, ...
             preFeatWAVEidx, ... 
             postFeatWAVEidx, ...
             preFeatCCG, ...
             postFeatCCG, ...
             preFeatWAVE, ... 
             postFeatWAVE, ...
             mainPkNorm, ...
             preFeatCCGnorm, ...
             postFeatCCGnorm, ...
             mainPkAgmonIdx, ...
             preFeatCCGagmonIdx, ...
             postFeatCCGagmonIdx] = featureLatencies(pklocsCCGmask_non_0lag, ...
                                                     trlocsCCGmask, ...
                                                     ccgAnalysis, ...
                                                     ccgjmAnalysis, ...
                                                     JBSIEanalysis, ...
                                                     JBSIIanalysis, ...
                                                     tRanalysis, ...
                                                     wave4Ana, ...
                                                     waveRange);
            
             %% pairing CCG an waveform features
            if ~isempty(mainPk) && ~isempty(mainTr)
                
                groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
                groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
                groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
                groupStat.mainPkAgmonIdx     = [groupStat.mainPkAgmonIdx mainPkAgmonIdx];
   
                groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
                groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
                
                groupStat.mainDistance = [groupStat.mainDistance distance];
                
                groupStat.mainRef = [groupStat.mainRef cellID1];
                groupStat.mainTar = [groupStat.mainTar cellID2];
                
                groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
                groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];

                if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
                    
                    groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
                    groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
                    groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
                    groupStat.preFeatCCGagmonIdx = [groupStat.preFeatCCGagmonIdx preFeatCCGagmonIdx];
                    
                    groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
                    groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
                     
                    groupStat.preDistance = [groupStat.preDistance distance];
                     
                    groupStat.preRef = [groupStat.preRef cellID1];
                    groupStat.preTar = [groupStat.preTar cellID2];

                    groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
                    groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];

                end

                if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
                    
                    groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
                    groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
                    groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
                    groupStat.postFeatCCGagmonIdx = [groupStat.postFeatCCGagmonIdx postFeatCCGagmonIdx];
                    
                    groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
                    groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
                    
                    groupStat.postDistance = [groupStat.postDistance distance];
                    
                    groupStat.postRef = [groupStat.postRef cellID1];
                    groupStat.postTar = [groupStat.postTar cellID2];

                    groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
                    groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];

                end
                
                yyaxis left
                plot(tR*1000,ccg,'b')
                hold on
                scatter(mainPkLat,mainPk,'b','LineWidth',2)
                scatter(preFeatureCCGlat,preFeatCCG,'b','LineWidth',2)
                scatter(postFeatureCCGlat,postFeatCCG,'b','LineWidth',2)
                hold off

                yyaxis right
                plot(spikeTime,wave4Ana,'r')
                hold on
                scatter(mainTrLat,mainTr,'r','LineWidth',2)
                scatter(preFeatureWAVElat,preFeatWAVE,'r','LineWidth',2)
                scatter(postFeatureWAVElat,postFeatWAVE,'r','LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlim([-3.5 3.5])

                titleChar = ['CCG and wave feature alighment - ' num2str(cellID1) cell1type ' v ' num2str(cellID2) cell2type ' - non-zero lag'];
                title(titleChar)

                disp([num2str((j/size(pairsAna,1))*100) ' % done!'])

                save_file = fullfile(datapath, titleChar);
                print(fig_use, save_file,'-dpng',resolution_use);
                
            end
            
        else 
            
            [mainPk, ...
             mainTr, ...
             mainPkLat, ...
             mainTrLat, ...
             preFeatureCCGlat, ...
             preFeatureWAVElat, ...
             postFeatureCCGlat, ...
             postFeatureWAVElat, ...
             preFeatCCGidx, ...
             postFeatCCGidx, ...
             preFeatWAVEidx, ... 
             postFeatWAVEidx, ...
             preFeatCCG, ...
             postFeatCCG, ...
             preFeatWAVE, ... 
             postFeatWAVE, ...
             mainPkNorm, ...
             preFeatCCGnorm, ...
             postFeatCCGnorm, ...
             mainPkAgmonIdx, ...
             preFeatCCGagmonIdx, ...
             postFeatCCGagmonIdx] = featureLatencies(pklocsCCGmask, ...
                                                     trlocsCCGmask, ...
                                                     ccgAnalysis, ...
                                                     ccgjmAnalysis, ...
                                                     JBSIEanalysis, ...
                                                     JBSIIanalysis, ...
                                                     tRanalysis, ...
                                                     wave4Ana, ...
                                                     waveRange);
            
             %% pairing CCG an waveform features
             if ~isempty(mainPk) && ~isempty(mainTr)
                
                groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
                groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
                groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
                groupStat.mainPkAgmonIdx     = [groupStat.mainPkAgmonIdx mainPkAgmonIdx];
   
                groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
                groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
                
                groupStat.mainDistance = [groupStat.mainDistance distance];
                
                groupStat.mainRef = [groupStat.mainRef cellID1];
                groupStat.mainTar = [groupStat.mainTar cellID2];
                
                groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
                groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];

                if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
                    
                    groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
                    groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
                    groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
                    groupStat.preFeatCCGagmonIdx = [groupStat.preFeatCCGagmonIdx preFeatCCGagmonIdx];
                    
                    groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
                    groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
                     
                    groupStat.preDistance = [groupStat.preDistance distance];
                     
                    groupStat.preRef = [groupStat.preRef cellID1];
                    groupStat.preTar = [groupStat.preTar cellID2];

                    groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
                    groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];

                end

                if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
                    
                    groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
                    groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
                    groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
                    groupStat.postFeatCCGagmonIdx = [groupStat.postFeatCCGagmonIdx postFeatCCGagmonIdx];
                    
                    groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
                    groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
                    
                    groupStat.postDistance = [groupStat.postDistance distance];
                    
                    groupStat.postRef = [groupStat.postRef cellID1];
                    groupStat.postTar = [groupStat.postTar cellID2];

                    groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
                    groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];

                end
                
                yyaxis left
                plot(tR*1000,ccg,'b')
                hold on
                scatter(mainPkLat,mainPk,'b','LineWidth',2)
                scatter(preFeatureCCGlat,preFeatCCG,'b','LineWidth',2)
                scatter(postFeatureCCGlat,postFeatCCG,'b','LineWidth',2)
                hold off

                yyaxis right
                plot(spikeTime,wave4Ana,'r')
                hold on
                scatter(mainTrLat,mainTr,'r','LineWidth',2)
                scatter(preFeatureWAVElat,preFeatWAVE,'r','LineWidth',2)
                scatter(postFeatureWAVElat,postFeatWAVE,'r','LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlim([-3.5 3.5])

                titleChar = ['CCG and wave feature alighment - ' num2str(cellID1) cell1type ' v ' num2str(cellID2) cell2type];
                title(titleChar)

                disp([num2str((j/size(pairsAna,1))*100) ' % done!'])

                save_file = fullfile(datapath, titleChar);
                print(fig_use, save_file,'-dpng',resolution_use);
                
            end
         
        end
        
        toc
        
        close all
        
        hcomb = figure(102);
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    end
    
    %% data struct for group stats
    save([datapath session_name 'GroupStats.mat'],'groupStat')
    
end

function [mainPk, ...
          mainTr, ...
          mainPkLat, ...
          mainTrLat, ...
          preFeatCCGlat, ...
          preFeatWAVElat, ...
          postFeatCCGlat, ...
          postFeatWAVElat, ...
          preFeatCCGidx, ...
          postFeatCCGidx, ...
          preFeatWAVEidx, ... 
          postFeatWAVEidx, ...
          preFeatCCG, ...
          postFeatCCG, ...
          preFeatWAVE, ... 
          postFeatWAVE, ...
          mainPkNorm, ...
          preFeatCCGnorm, ...
          postFeatCCGnorm, ...
          mainPkAgmonIdx, ...
          preFeatCCGagmonIdx, ...
          postFeatCCGagmonIdx] = featureLatencies(pklocsCCGmask, ...
                                              trlocsCCGmask, ...
                                              ccgAnalysis, ...
                                              ccgjmAnalysis, ...
                                              JBSIEanalysis, ...
                                              JBSIIanalysis, ...
                                              tRanalysis, ...
                                              wave4Ana, ...
                                              waveRange)
                                          
    % normalize ccg 
    ccgNorm = (ccgAnalysis-ccgjmAnalysis)/sum(ccgAnalysis);
                                          
    wave4Ana = wave4Ana(waveRange);
                                          
    screenedPkIdx = find(pklocsCCGmask);
    screenedTrIdx = find(trlocsCCGmask);

    [mainPk,idx]   = max(ccgAnalysis(screenedPkIdx));
    idxMainPk      = screenedPkIdx(idx);
    mainPkLat      = tRanalysis(screenedPkIdx(idx));
    mainPkNorm     = ccgNorm(screenedPkIdx(idx));  % normalized 
    mainPkAgmonIdx = JBSIEanalysis(screenedPkIdx(idx));

    preFeatCCGidx  = find(idxMainPk > screenedTrIdx,1,'last'); 
    postFeatCCGidx = find(idxMainPk < screenedTrIdx,1);
    
    % pre/post troughs    
    if ~isempty(preFeatCCGidx)
        preFeatCCG         = ccgAnalysis(screenedTrIdx(preFeatCCGidx));
        preFeatCCGnorm     = ccgNorm(screenedTrIdx(preFeatCCGidx)); % normalized 
        preFeatCCGagmonIdx = JBSIEanalysis(screenedTrIdx(preFeatCCGidx));
        preFeatCCGlat      = tRanalysis(screenedTrIdx(preFeatCCGidx));
    else
        preFeatCCG         = [];
        preFeatCCGnorm     = [];
        preFeatCCGagmonIdx = [];
        preFeatCCGlat      = [];
    end

    if ~isempty(postFeatCCGidx)
        postFeatCCG         = ccgAnalysis(screenedTrIdx(postFeatCCGidx));
        postFeatCCGnorm     = ccgNorm(screenedTrIdx(postFeatCCGidx)); % normalized 
        postFeatCCGagmonIdx = JBSIEanalysis(screenedTrIdx(postFeatCCGidx));
        postFeatCCGlat      = tRanalysis(screenedTrIdx(postFeatCCGidx));
    else
        postFeatCCG         = [];
        postFeatCCGnorm     = [];
        postFeatCCGagmonIdx = [];
        postFeatCCGlat      = [];
    end

    %% sign flipped on waveform because it's an eSpike

    % troughs
    [trWAVE, trlocsWAVE] = findpeaks(-wave4Ana);
    trWAVE = -trWAVE;

    [mainTr,idx] = min(trWAVE);
    idxMainTr    = trlocsWAVE(idx);
    mainTrLat    = tRanalysis(idxMainTr);

    % peaks
    [pkWAVE, pklocsWAVE] = findpeaks(wave4Ana);
    preFeatWAVEidx  = find(idxMainTr > pklocsWAVE,1,'last'); 
    postFeatWAVEidx = find(idxMainTr < pklocsWAVE,1);
      
    % pre/post peaks
    if ~isempty(preFeatWAVEidx)
         preFeatWAVE       = wave4Ana(pklocsWAVE(preFeatWAVEidx));
         preFeatWAVElat    = tRanalysis(pklocsWAVE(preFeatWAVEidx));
    else
         preFeatWAVE       = [];
         preFeatWAVElat    = [];
    end

    if ~isempty(postFeatWAVEidx)
        postFeatWAVE       = wave4Ana(pklocsWAVE(postFeatWAVEidx));
        postFeatWAVElat    = tRanalysis(pklocsWAVE(postFeatWAVEidx));
    else
        postFeatWAVE       = [];
        postFeatWAVElat    = [];
    end
    
    %% if CCG peak and wave trough are more than 0.2msec apart set
    if abs(mainTrLat-mainPkLat) > 0.2
        mainTr              = [];
        mainPk              = [];
        mainPkNorm          = [];
        mainPkAgmonIdx      = [];
        
        mainTrLat           = [];
        mainPkLat           = [];
        
        preFeatCCG          = [];
        preFeatCCGlat       = [];
        preFeatCCGnorm      = [];
        preFeatCCGagmonIdx  = [];
        
        postFeatCCG         = [];
        postFeatCCGlat      = [];
        postFeatCCGnorm     = [];
        postFeatCCGagmonIdx = [];
        
        preFeatWAVE         = [];
        preFeatWAVElat      = [];
        
        postFeatWAVE        = [];
        postFeatWAVElat     = [];
    end

end