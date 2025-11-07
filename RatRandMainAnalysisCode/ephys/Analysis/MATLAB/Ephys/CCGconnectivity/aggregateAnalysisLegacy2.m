function aggregateAnalysis(datasetID,animalIdx,anaDirection)
    
    clc; 
        
    %% hardcoding pairs of interest for waveform analysis
    
    fig_use        = 102;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    % pair data 
    if strcmp(datasetID,'RatRoy')
        datapath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
        load([datapath 'groupStatsRatRoy'])
    elseif strcmp(datasetID,'RatAGday1')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/';
        load([datapath 'groupStatsRatAGday1'])
    elseif strcmp(datasetID,'RatAGday2')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/';
        load([datapath 'groupStatsRatAGday2'])
    elseif strcmp(datasetID,'RatN')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/';
        load([datapath 'groupStatsRatN'])
    elseif strcmp(datasetID,'RatS')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/';
        load([datapath 'groupStatsRatS'])
    elseif strcmp(datasetID,'RatU')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/';
        load([datapath 'groupStatsRatU'])
    elseif strcmp(datasetID,'Steinmetz')
        datapath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';
        load([datapath 'groupStatsSteinmetzDataset'])
    elseif strcmp(datasetID,'AllenInstitute')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
        animalsList = [715093703
                       719161530
                       721123822
                       732592105
                       737581020
                       739448407
                       742951821
                       743475441
                       744228101
                       746083955
                       750332458
                       750749662
                       751348571
                       754312389
                       754829445
                       755434585
                       756029989
                       757216464
                       757970808
                       758798717
                       759883607
                       760345702
                       760693773
                       761418226
                       762120172
                       762602078
                       763673393
                       766640955
                       767871931
                       768515987
                       771160300
                       771990200
                       773418906
                       774875821
                       778240327
                       778998620
                       779839471
                       781842082
                       786091066
                       787025148
                       789848216
                       791319847
                       793224716
                       794812542
                       797828357
                       798911424
                       799864342
                       816200189
                       819186360
                       819701982
                       821695405
                       829720705
                       831882777
                       835479236
                       839068429
                       839557629
                       840012044
                       847657808];
        load([datapath 'groupStatsMouse' num2str(animalsList(animalIdx))])
        pairGroupStatTable = pairGroupStatsTable;
        clear pairGroupStatsTable
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

    %%

    for loopPairs = 1:size(pairGroupStatTable,1)
        
        %% pearson's R
        if strcmp(anaDirection,'forward')
            ccg = pairGroupStatTable.pairRawCCG{loopPairs,1};

            refID = pairGroupStatTable.refNeuronID(loopPairs);
            tarID = pairGroupStatTable.tarNeuronID(loopPairs);
            
            if strcmp(datasetID,'Steinmetz') 
                refWaveAtRefShank = pairGroupStatTable.refWaveforms{loopPairs,1}';
                refWaveAtTarShank = pairGroupStatTable.refWaveforms{loopPairs,1}';
                
                channelSetRef = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},2);
                channelSetTar = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},2);
            elseif strcmp(datasetID,'AllenInstitute')
                refWaveAtRefShank = pairGroupStatTable.refWaveforms{loopPairs,1};
                refWaveAtTarShank = pairGroupStatTable.refWaveforms{loopPairs,1};
                
                channelSetRef = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},1);
                channelSetTar = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},1);
            else
                refWaveAtRefShank = pairGroupStatTable.refWaveAtRefShank{loopPairs,1};
                refWaveAtTarShank = pairGroupStatTable.refWaveAtTarShank{loopPairs,1};

                channelSetRef = pairGroupStatTable.channelSetRef{loopPairs,1};
                channelSetTar = pairGroupStatTable.channelSetTar{loopPairs,1};
            end
        elseif strcmp(anaDirection,'reverse')
            ccg = flipud(pairGroupStatTable.pairRawCCG{loopPairs,1});
            
            tarID = pairGroupStatTable.refNeuronID(loopPairs);
            refID = pairGroupStatTable.tarNeuronID(loopPairs);
            
            if strcmp(datasetID,'Steinmetz') 
                refWaveAtTarShank = pairGroupStatTable.tarWaveforms{loopPairs,1}';
                refWaveAtRefShank = pairGroupStatTable.tarWaveforms{loopPairs,1}';
                
                channelSetTar = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},2);
                channelSetRef = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},2);
            elseif strcmp(datasetID,'AllenInstitute')
                refWaveAtTarShank = pairGroupStatTable.tarWaveforms{loopPairs,1};
                refWaveAtRefShank = pairGroupStatTable.tarWaveforms{loopPairs,1};
                
                channelSetTar = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},1);
                channelSetRef = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},1);    
            else
                refWaveAtTarShank = pairGroupStatTable.tarWaveAtRefShank{loopPairs,1};
                refWaveAtRefShank = pairGroupStatTable.tarWaveAtTarShank{loopPairs,1};

                channelSetTar = pairGroupStatTable.channelSetRef{loopPairs,1};
                channelSetRef = pairGroupStatTable.channelSetTar{loopPairs,1};
            end
        end

    %     ccgRange  = 106:136; % 1msec (0 to 1ms lag)
    %     waveRange = 36:66;   % 1msec (0 to 1ms lag)

        
        if strcmp(datasetID,'Steinmetz')
            ccgRange  = (96:136)'; % 1.25msec (-0.25 to 1ms lag)
            waveRange = (31:71)';  % 1.25msec (-0.25 to 1ms lag)
        elseif strcmp(datasetID,'AllenInstitute')
%             ccgRange  = (96:136)'; % 1.25msec (-0.25 to 1ms lag)
%             waveRange = (10:50)';  % 1.25msec (-0.25 to 1ms lag)
            ccgRange  = (87:168)'; % 1.25msec (-0.25 to 1ms lag)
            waveRange = (1:82)';  % 1.25msec (-0.25 to 1ms lag)
        else
            ccgRange  = (96:136)'; % 1.25msec (-0.25 to 1ms lag)
            waveRange = (26:66)';  % 1.25msec (-0.33 to 1ms lag)
        end

        RsqList = [];
        for k = 1:(size(refWaveAtRefShank,1) + size(refWaveAtTarShank,1))
            if k <= size(refWaveAtRefShank,1)
                wave4Ana = refWaveAtRefShank(k,:)';
            elseif k > size(refWaveAtRefShank,1)
                wave4Ana = refWaveAtTarShank(k-size(refWaveAtRefShank,1),:)';
            end
            
            R = corrcoef(ccg(ccgRange),wave4Ana(waveRange));
            
            if R(1,2) >= 0
                RsqList(k) = 0;
            elseif R(1,2) < 0
                RsqList(k) = R(1,2)^2;
            end
            
        end 

        [RsqMax,maxMatchChIdx] = max(RsqList);
        
        % don't include for effect group
        if RsqMax < 0.5
            
%             disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
            disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
            
            continue 
        end
        
        allChans   = [channelSetRef channelSetTar];
        maxMatchCh = allChans(maxMatchChIdx);
        
        % distance calc
%         xSource   = probe.x(probe.chanNo == chCell1);
%         ySource   = probe.y(probe.chanNo == chCell1);
%         
%         xReceiver = probe.x(probe.chanNo == maxMatchCh);
%         yReceiver = probe.y(probe.chanNo == maxMatchCh);
%         
%         distance  = (((xSource - xReceiver)^2) + ((ySource - yReceiver)^2))^(1/2);

        distance = pairGroupStatTable.pairDistance(loopPairs);
        
        % waveform for analysis based on max channel
        if maxMatchChIdx <= size(refWaveAtRefShank,1)
            wave4Ana = refWaveAtRefShank(maxMatchChIdx,:)';
        elseif maxMatchChIdx > size(refWaveAtRefShank,1)
            wave4Ana = refWaveAtTarShank(maxMatchChIdx-size(refWaveAtRefShank,1),:)';
        end
        
        %% find peaks and troughs CCG
        if strcmp(anaDirection,'forward')
            peaksMask     = pairGroupStatTable.GSPExc{loopPairs,1};         peaksMask     = peaksMask(ccgRange);
            troughsMask   = pairGroupStatTable.GSPInh{loopPairs,1};         troughsMask   = troughsMask(ccgRange);
            ccgjmAnalysis = pairGroupStatTable.jitterMean{loopPairs,1};     ccgjmAnalysis = ccgjmAnalysis(ccgRange);

            refCellType   = pairGroupStatTable.refCellExplorerType{loopPairs,1};
            tarCellType   = pairGroupStatTable.tarCellExplorerType{loopPairs,1};
        elseif strcmp(anaDirection,'reverse')
            peaksMask     = flipud(pairGroupStatTable.GSPExc{loopPairs,1});         peaksMask     = peaksMask(ccgRange);
            troughsMask   = flipud(pairGroupStatTable.GSPInh{loopPairs,1});         troughsMask   = troughsMask(ccgRange);
            ccgjmAnalysis = flipud(pairGroupStatTable.jitterMean{loopPairs,1});     ccgjmAnalysis = ccgjmAnalysis(ccgRange);
            
            tarCellType   = pairGroupStatTable.refCellExplorerType{loopPairs,1};
            refCellType   = pairGroupStatTable.tarCellExplorerType{loopPairs,1};
        end
        
        tR            = pairGroupStatTable.CCGbinLagTimes{loopPairs,1}*1000;
        tRanalysis    = pairGroupStatTable.CCGbinLagTimes{loopPairs,1}(ccgRange);
        ccgAnalysis   = ccg(ccgRange);
        tWave         = pairGroupStatTable.tWave{loopPairs,1}';
        
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
        
        % reject pairs with no significant bins in the analysis range
        if ~((sum(pklocsCCGmask(10:12)) > 0) && (sum(pklocsCCGmask([1:9,13:41])) > 0))
            
%                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
            disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])

            continue
        end
        
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
             postFeatCCGnorm] = featureLatencies(pklocsCCGmask_0lag, ...
                                                     trlocsCCGmask, ...
                                                     ccgAnalysis, ...
                                                     ccgjmAnalysis, ...
                                                     tRanalysis, ...
                                                     wave4Ana, ...
                                                     waveRange);
            
            %% pairing CCG an waveform features
            if ~isempty(mainPk) && ~isempty(mainTr)
                
                groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
                groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
                groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
   
                groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
                groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
                
                groupStat.mainDistance = [groupStat.mainDistance distance];
                
                groupStat.mainRef = [groupStat.mainRef refID];
                groupStat.mainTar = [groupStat.mainTar tarID];
                
                groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
                groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];

                if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
                    
                    groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
                    groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
                    groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
                    
                    groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
                    groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
                     
                    groupStat.preDistance = [groupStat.preDistance distance];
                     
                    groupStat.preRef = [groupStat.preRef refID];
                    groupStat.preTar = [groupStat.preTar tarID];

                    groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
                    groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];

                end

                if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
                    
                    groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
                    groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
                    groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
                    
                    groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
                    groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
                    
                    groupStat.postDistance = [groupStat.postDistance distance];
                    
                    groupStat.postRef = [groupStat.postRef refID];
                    groupStat.postTar = [groupStat.postTar tarID];

                    groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
                    groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];

                end
                
                yyaxis left
                plot(tR,ccg,'b')
                hold on
                scatter(mainPkLat*1000,mainPk,'b','LineWidth',2)
                scatter(preFeatureCCGlat*1000,preFeatCCG,'b','LineWidth',2)
                scatter(postFeatureCCGlat*1000,postFeatCCG,'b','LineWidth',2)
                hold off

                yyaxis right
                plot(tWave,wave4Ana,'r')
                hold on
                scatter(mainTrLat*1000,mainTr,'r','LineWidth',2)
                scatter(preFeatureWAVElat*1000,preFeatWAVE,'r','LineWidth',2)
                scatter(postFeatureWAVElat*1000,postFeatWAVE,'r','LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlim([-3.5 3.5])

                titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
                title(titleChar)

%                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
                disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])

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
             postFeatCCGnorm] = featureLatencies(pklocsCCGmask_non_0lag, ...
                                                 trlocsCCGmask, ...
                                                 ccgAnalysis, ...
                                                 ccgjmAnalysis, ...
                                                 tRanalysis, ...
                                                 wave4Ana, ...
                                                 waveRange);
            
             %% pairing CCG an waveform features
            if ~isempty(mainPk) && ~isempty(mainTr)
                
                groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
                groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
                groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
   
                groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
                groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
                
                groupStat.mainDistance = [groupStat.mainDistance distance];
                
                groupStat.mainRef = [groupStat.mainRef refID];
                groupStat.mainTar = [groupStat.mainTar tarID];
                
                groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
                groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];

                if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
                    
                    groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
                    groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
                    groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
                    
                    groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
                    groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
                     
                    groupStat.preDistance = [groupStat.preDistance distance];
                     
                    groupStat.preRef = [groupStat.preRef refID];
                    groupStat.preTar = [groupStat.preTar tarID];

                    groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
                    groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];

                end

                if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
                    
                    groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
                    groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
                    groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
                    
                    groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
                    groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
                    
                    groupStat.postDistance = [groupStat.postDistance distance];
                    
                    groupStat.postRef = [groupStat.postRef refID];
                    groupStat.postTar = [groupStat.postTar tarID];

                    groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
                    groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];

                end
                
                yyaxis left
                plot(tR,ccg,'b')
                hold on
                scatter(mainPkLat*1000,mainPk,'b','LineWidth',2)
                scatter(preFeatureCCGlat*1000,preFeatCCG,'b','LineWidth',2)
                scatter(postFeatureCCGlat*1000,postFeatCCG,'b','LineWidth',2)
                hold off

                yyaxis right
                plot(tWave,wave4Ana,'r')
                hold on
                scatter(mainTrLat*1000,mainTr,'r','LineWidth',2)
                scatter(preFeatureWAVElat*1000,preFeatWAVE,'r','LineWidth',2)
                scatter(postFeatureWAVElat*1000,postFeatWAVE,'r','LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlim([-3.5 3.5])

                titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
                title(titleChar)

%                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
                disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])

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
             postFeatCCGnorm] = featureLatencies(pklocsCCGmask, ...
                                                     trlocsCCGmask, ...
                                                     ccgAnalysis, ...
                                                     ccgjmAnalysis, ...
                                                     tRanalysis, ...
                                                     wave4Ana, ...
                                                     waveRange);
            
             %% pairing CCG an waveform features
             if ~isempty(mainPk) && ~isempty(mainTr)
                
                groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
                groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
                groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
   
                groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
                groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
                
                groupStat.mainDistance = [groupStat.mainDistance distance];
                
                groupStat.mainRef = [groupStat.mainRef refID];
                groupStat.mainTar = [groupStat.mainTar tarID];
                
                groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
                groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];

                if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
                    
                    groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
                    groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
                    groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
                    
                    groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
                    groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
                     
                    groupStat.preDistance = [groupStat.preDistance distance];
                     
                    groupStat.preRef = [groupStat.preRef refID];
                    groupStat.preTar = [groupStat.preTar tarID];

                    groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
                    groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];

                end

                if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
                    
                    groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
                    groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
                    groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
                    
                    groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
                    groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
                    
                    groupStat.postDistance = [groupStat.postDistance distance];
                    
                    groupStat.postRef = [groupStat.postRef refID];
                    groupStat.postTar = [groupStat.postTar tarID];

                    groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
                    groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];

                end
                
                yyaxis left
                plot(tR,ccg,'b')
                hold on
                scatter(mainPkLat*1000,mainPk,'b','LineWidth',2)
                scatter(preFeatureCCGlat*1000,preFeatCCG,'b','LineWidth',2)
                scatter(postFeatureCCGlat*1000,postFeatCCG,'b','LineWidth',2)
                hold off

                yyaxis right
                plot(tWave,wave4Ana,'r')
                hold on
                scatter(mainTrLat*1000,mainTr,'r','LineWidth',2)
                scatter(preFeatureWAVElat*1000,preFeatWAVE,'r','LineWidth',2)
                scatter(postFeatureWAVElat*1000,postFeatWAVE,'r','LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlim([-3.5 3.5])

                titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
                title(titleChar)

%                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
                disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
                
                save_file = fullfile(datapath, titleChar);
                print(fig_use, save_file,'-dpng',resolution_use);
                
            end
         
        end
                
        close all
        
        hcomb = figure(102);
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    end
    
    %% data struct for group stats
    save([datapath datasetID anaDirection 'ProcessedGroupStats.mat'],'groupStat')
    
    close all
    
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
          postFeatCCGnorm] = featureLatencies(pklocsCCGmask, ...
                                              trlocsCCGmask, ...
                                              ccgAnalysis, ...
                                              ccgjmAnalysis, ...
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

    preFeatCCGidx  = find(idxMainPk > screenedTrIdx,1,'last'); 
    postFeatCCGidx = find(idxMainPk < screenedTrIdx,1);
    
    % pre/post troughs    
    if ~isempty(preFeatCCGidx)
        preFeatCCG         = ccgAnalysis(screenedTrIdx(preFeatCCGidx));
        preFeatCCGnorm     = ccgNorm(screenedTrIdx(preFeatCCGidx)); % normalized 
        preFeatCCGlat      = tRanalysis(screenedTrIdx(preFeatCCGidx));
    else
        preFeatCCG         = [];
        preFeatCCGnorm     = [];
        preFeatCCGlat      = [];
    end

    if ~isempty(postFeatCCGidx)
        postFeatCCG         = ccgAnalysis(screenedTrIdx(postFeatCCGidx));
        postFeatCCGnorm     = ccgNorm(screenedTrIdx(postFeatCCGidx)); % normalized 
        postFeatCCGlat      = tRanalysis(screenedTrIdx(postFeatCCGidx));
    else
        postFeatCCG         = [];
        postFeatCCGnorm     = [];
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
        
        mainTrLat           = [];
        mainPkLat           = [];
        
        preFeatCCG          = [];
        preFeatCCGlat       = [];
        preFeatCCGnorm      = [];
        
        postFeatCCG         = [];
        postFeatCCGlat      = [];
        postFeatCCGnorm     = [];
        
        preFeatWAVE         = [];
        preFeatWAVElat      = [];
        
        postFeatWAVE        = [];
        postFeatWAVElat     = [];
    end

end