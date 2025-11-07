function aggregateAnalysisV2(datasetID,animalIdx,anaDirection,GJanaFlag)
    
    clc; 
    
    warning('off','signal:findpeaks:largeMinPeakHeight')
    
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
    
    fs = 3e4;
    
    % pair data 
    if strcmp(datasetID,'RatRoy')
        datapath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
        if GJanaFlag
            load([datapath 'groupStatsRatRoyPutativeGJ'])
        else
            load([datapath 'groupStatsRatRoy'])
        end
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        noPairCompTable = countPairsHiro;
    elseif strcmp(datasetID,'RatAGday1')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/';
        if GJanaFlag
            load([datapath 'groupStatsRatAGday1putativeGJ'])
        else
            load([datapath 'groupStatsRatAGday1'])
        end 
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        noPairCompTable = countPairsUtku;
    elseif strcmp(datasetID,'RatAGday2')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/';
        if GJanaFlag
            load([datapath 'groupStatsRatAGday2putativeGJ'])
        else
            load([datapath 'groupStatsRatAGday2'])
        end
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        noPairCompTable = countPairsUtku;
    elseif strcmp(datasetID,'RatN')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/';
        if GJanaFlag
            load([datapath 'groupStatsRatNputativeGJ'])
        else
            load([datapath 'groupStatsRatN'])
        end
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        noPairCompTable = countPairsBapun;
    elseif strcmp(datasetID,'RatS')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/';
        if GJanaFlag
            load([datapath 'groupStatsRatSputativeGJ'])
        else
            load([datapath 'groupStatsRatS'])
        end
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        noPairCompTable = countPairsBapun;
    elseif strcmp(datasetID,'RatU')
        datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/';
        if GJanaFlag
            load([datapath 'groupStatsRatUputativeGJ'])
        else
            load([datapath 'groupStatsRatU'])
        end
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        noPairCompTable = countPairsBapun;
    elseif strcmp(datasetID,'Steinmetz')
        datapath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';
        load([datapath 'groupStatsSteinmetzDataset'])
        
        % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
    elseif strcmp(datasetID,'AllenInstitute')
%         datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
        datapath = '/media/nasko/WD_BLACK31/BOTtemp/';
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
        
        display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(animalIdx))])
        
        if GJanaFlag
            load([datapath 'groupStatsMousePutativeGJ' num2str(animalsList(animalIdx))])
            pairGroupStatTable = pairGroupStatsTable;
            clear pairGroupStatsTable
        else
            load([datapath 'groupStatsMouse' num2str(animalsList(animalIdx))])
            pairGroupStatTable = pairGroupStatsTable;
            clear pairGroupStatsTable
        end
        
         % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
        
        pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable);
        
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/noPairCompTable.mat');
    end

    %%

    for loopPairs = 1:size(pairGroupStatTable,1)
       
        display(['pair: ' num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1))])
        if strcmp(anaDirection,'forward')
            ccg = pairGroupStatTable.pairRawCCG{loopPairs,1};

            refID = pairGroupStatTable.refNeuronID(loopPairs);
            tarID = pairGroupStatTable.tarNeuronID(loopPairs);
            
            refCellType = pairGroupStatTable.refCellExplorerType{loopPairs,1};
            tarCellType = pairGroupStatTable.tarCellExplorerType{loopPairs,1};
            
            acg_narrow = CCG(pairGroupStatTable.refSpikeTimes{loopPairs,1},ones(size(pairGroupStatTable.refSpikeTimes{loopPairs,1})),'binSize',0.0005,'duration',0.100,'norm','rate','Fs',1/fs);
            fit_params_out = fit_ACG(acg_narrow,false);
            pairGroupStatTable.refACGtauRise(loopPairs) = fit_params_out.acg_tau_rise;
            
            acg_narrow = CCG(pairGroupStatTable.tarSpikeTimes{loopPairs,1},ones(size(pairGroupStatTable.tarSpikeTimes{loopPairs,1})),'binSize',0.0005,'duration',0.100,'norm','rate','Fs',1/fs);
            fit_params_out = fit_ACG(acg_narrow,false);
            pairGroupStatTable.tarACGtauRise(loopPairs) = fit_params_out.acg_tau_rise;
            
            %%
            
%             if strcmp(datasetID,'Steinmetz') 
%                 refWaveAtRefShank = pairGroupStatTable.refWaveforms{loopPairs,1}';
%                 refWaveAtTarShank = pairGroupStatTable.refWaveforms{loopPairs,1}';
%                 
%                 channelSetRef = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},2);
%                 channelSetTar = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},2);
%             elseif strcmp(datasetID,'AllenInstitute')
%                 refWaveAtRefShank = pairGroupStatTable.refWaveforms{loopPairs,1};
%                 refWaveAtTarShank = pairGroupStatTable.refWaveforms{loopPairs,1};
%                 
%                 channelSetRef = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},1);
%                 channelSetTar = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},1);
%             else
%                 refWaveAtRefShank = pairGroupStatTable.refWaveAtRefShank{loopPairs,1};
%                 refWaveAtTarShank = pairGroupStatTable.refWaveAtTarShank{loopPairs,1};
% 
%                 channelSetRef = pairGroupStatTable.channelSetRef{loopPairs,1};
%                 channelSetTar = pairGroupStatTable.channelSetTar{loopPairs,1};
%             end
%         elseif strcmp(anaDirection,'reverse')
%             ccg = flipud(pairGroupStatTable.pairRawCCG{loopPairs,1});
%             
%             tarID = pairGroupStatTable.refNeuronID(loopPairs);
%             refID = pairGroupStatTable.tarNeuronID(loopPairs);
%             
%             if strcmp(datasetID,'Steinmetz') 
%                 refWaveAtTarShank = pairGroupStatTable.tarWaveforms{loopPairs,1}';
%                 refWaveAtRefShank = pairGroupStatTable.tarWaveforms{loopPairs,1}';
%                 
%                 channelSetTar = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},2);
%                 channelSetRef = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},2);
%             elseif strcmp(datasetID,'AllenInstitute')
%                 refWaveAtTarShank = pairGroupStatTable.tarWaveforms{loopPairs,1};
%                 refWaveAtRefShank = pairGroupStatTable.tarWaveforms{loopPairs,1};
%                 
%                 channelSetTar = 1:size(pairGroupStatTable.refWaveforms{loopPairs,1},1);
%                 channelSetRef = 1:size(pairGroupStatTable.tarWaveforms{loopPairs,1},1);    
%             else
%                 refWaveAtTarShank = pairGroupStatTable.tarWaveAtRefShank{loopPairs,1};
%                 refWaveAtRefShank = pairGroupStatTable.tarWaveAtTarShank{loopPairs,1};
% 
%                 channelSetTar = pairGroupStatTable.channelSetRef{loopPairs,1};
%                 channelSetRef = pairGroupStatTable.channelSetTar{loopPairs,1};
%             end
        end
        
    % plot waves at other shank
    
%         tiledlayout(2,1)
%         
%         refID       = pairGroupStatTable.refNeuronID(loopPairs,1);
%         tarID       = pairGroupStatTable.tarNeuronID(loopPairs,1);
%         refCellType = pairGroupStatTable.refCellExplorerType{loopPairs,1};
%         tarCellType = pairGroupStatTable.tarCellExplorerType{loopPairs,1};
%         
%         tarChan     = pairGroupStatTable.tarChannel(loopPairs,1);
%         
%         tR     = pairGroupStatTable.CCGbinLagTimes{1,1}*1000;
%         ccg    = pairGroupStatTable.pairRawCCG{loopPairs,1};
%         
%         GSPExc = pairGroupStatTable.GSPExc{loopPairs,1};
%         GSPInh = pairGroupStatTable.GSPInh{loopPairs,1};
%         
%         tWave  = pairGroupStatTable.tWave{1,1};
%         refWaveAtTarShank = pairGroupStatTable.refWaveAtTarShank{loopPairs,1};
%         channelSetTar = pairGroupStatTable.channelSetTar{loopPairs,1};
%         
%         nexttile
%         plot(tR,ccg,'LineWidth',3)
%         hold on
%         xlim([-1 1])
%         ylims = get(gca,'ylim');
%         title([num2str(refID) refCellType ' v ' num2str(tarID) tarCellType ' target chan: ' num2str(tarChan)]);
% 
%         if any(GSPExc)
%             hold on;
%             plot(tR(GSPExc == 1), 0.95*ylims(2), 'r^');
%         end
%         if any(GSPInh)
%             hold on;
%             plot(tR(GSPInh == 1), 0.95*ylims(2),'bv');
%         end
%         hold off
%         
%         nexttile
%         hold on
%         plot(tWave,refWaveAtTarShank,'LineWidth',3)
%         hold off
%         set(gca, 'YDir','reverse')
%         legend(num2str(channelSetTar'))
%          
%         xlim([-1 1])
        
%         titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
        
    %     ccgRange  = 106:136; % 1msec (0 to 1ms lag)
    %     waveRange = 36:66;   % 1msec (0 to 1ms lag)
        
        

        %% measure spike widths (trough-to-peak)
        if ~strcmp(datasetID,'AllenInstitute') || (~strcmp(datasetID,'RatAGday1') || ~strcmp(datasetID,'RatAGday2'))
            
            preLength   = find(pairGroupStatTable.tWave{1,1} == 0); % where the trough is
            metricParam = 15;  

            spikeTime = pairGroupStatTable.tWave{1,:};

            %
            if strcmp(datasetID,'Steinmetz')
                waveTemp  = pairGroupStatTable.refWaveforms{loopPairs,:}(:,pairGroupStatTable.refChannel(loopPairs));
            else
                waveTemp  = pairGroupStatTable.refWaveforms{loopPairs,:}(pairGroupStatTable.refChannel(loopPairs),:);
            end

            [~,idxTrough]     = min(waveTemp(preLength-metricParam:preLength+metricParam));
            [~,idxPeak]       = max(waveTemp(preLength:preLength+metricParam));

            timeTrough = spikeTime(idxTrough + preLength-metricParam);
            timePeak   = spikeTime(idxPeak   + preLength);

            pairGroupStatTable.refTroughToPeak(loopPairs) = timePeak - timeTrough;

            if strcmp(datasetID,'Steinmetz')
                waveTemp  = pairGroupStatTable.tarWaveforms{loopPairs,:}(:,pairGroupStatTable.tarChannel(loopPairs));
            else
                waveTemp  = pairGroupStatTable.tarWaveforms{loopPairs,:}(pairGroupStatTable.tarChannel(loopPairs),:);
            end

            [~,idxTrough]     = min(waveTemp(preLength-metricParam:preLength+metricParam));
            [~,idxPeak]       = max(waveTemp(preLength:preLength+metricParam));

            timeTrough = spikeTime(idxTrough + preLength-metricParam);
            timePeak   = spikeTime(idxPeak   + preLength);

            pairGroupStatTable.tarTroughToPeak(loopPairs) = timePeak - timeTrough;
            
        end
        
        if strcmp(datasetID,'RatAGday1') 
            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/AG_2019-12-23_NSD.cell_metrics.cellinfo.mat')
            
            pairGroupStatTable.refTroughToPeak(loopPairs) = cell_metrics.troughToPeak(refID);
            pairGroupStatTable.tarTroughToPeak(loopPairs) = cell_metrics.troughToPeak(tarID);
        elseif strcmp(datasetID,'RatAGday2')
            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/AG_2019-12-27_NSD.cell_metrics.cellinfo.mat')
            
            pairGroupStatTable.refTroughToPeak(loopPairs) = cell_metrics.troughToPeak(refID);
            pairGroupStatTable.tarTroughToPeak(loopPairs) = cell_metrics.troughToPeak(tarID);
        end
        %% find CCG peaks and troughs (not putative GJ)
        
        if ~GJanaFlag
            jitterMean   = pairGroupStatTable.jitterMean{loopPairs,1};
            ccgCorr      = ccg - jitterMean;
            thresh       = max(abs(ccgCorr))/4;

            patchLen = 5;
            patch    = ones(patchLen,1); 
            patchAdj = floor(length(patch)/2);

    %         if (median(ccg) >= 0) || (median(ccg) <= 1) % dealing with false positives
                [~,idx]  = findpeaks(smoothdata(ccg,"gauss",5), ...
                                    'MinPeakProminence',thresh, ...
                                    'MinPeakWidth',2, ...
                                    'MinPeakDistance',9, ... % minimum 300 miscrosec
                                    'Annotate','extents');
    %         else
    %             [~,idx]  = findpeaks(ccgCorr,'MinPeakProminence',thresh,'MinPeakWidth',2);
    %         end
            filtIdx      = zeros(211,1);
            filtIdx(idx) = 1;
            filtIdx      = filtIdx & pairGroupStatTable.GSPExc{loopPairs,1};
            latPeakCCGactual = pairGroupStatTable.CCGbinLagTimes{1,1}(filtIdx)*1000;
            peakCCGactual    = ccg(filtIdx);
            filtIdx      = conv(filtIdx,patch);
            filtIdx      = find(filtIdx) - patchAdj;
            peakCCG      = ccg(filtIdx);
            idxPeakCCG   = filtIdx;
            latPeakCCG   = pairGroupStatTable.CCGbinLagTimes{1,1}(filtIdx)*1000;

    %         if (median(ccg) >= 0) || (median(ccg) <= 1) % dealing with false positives
                [~,idx]  = findpeaks(-smoothdata(ccg,"gauss",5), ... 
                                    'MinPeakProminence',thresh, ... 
                                    'MinPeakWidth',2, ... 
                                    'MinPeakDistance',9, ... % minimum 300 miscrosec
                                    'Annotate','extents');
    %         else
    %             [~,idx]  = findpeaks(-ccgCorr,'MinPeakProminence',thresh,'MinPeakHeigh',2,'Annotate','extents');
    %         end
            filtIdx      = zeros(211,1);
            filtIdx(idx) = 1;
            filtIdx      = filtIdx & pairGroupStatTable.GSPInh{loopPairs,1};
            latTroughCCGactual = pairGroupStatTable.CCGbinLagTimes{1,1}(filtIdx)*1000;
            troughCCGactual    = ccg(filtIdx);
            filtIdx      = conv(filtIdx,patch);
            filtIdx      = find(filtIdx) - patchAdj;
            troughCCG    = ccg(filtIdx);
            idxTroughCCG = filtIdx;
            latTroughCCG = pairGroupStatTable.CCGbinLagTimes{1,1}(filtIdx)*1000;

            featuresCCGidx  = [idxPeakCCG; idxTroughCCG];
            featuresCCGval  = [peakCCG;    troughCCG];
            featuresCCGlat  = [latPeakCCG; latTroughCCG];
            featuresCCGactualLat   = [latPeakCCGactual; latTroughCCGactual];
    %         if length(featuresCCGactualLat) < 2 % require at least two significant CCG bins after findpeaks screening (this is conservative!)
    %            continue 
    %         end
            featuresCCGactualVal   = [peakCCGactual;    troughCCGactual];
            featuresCCGactualGroup = [];
            for k = 1:length(featuresCCGactualLat)
                featuresCCGactualGroup = [featuresCCGactualGroup; k*ones(patchLen,1)];
            end
            featuresCCGactualType = [repmat({'peak'}, 1, length(peakCCGactual)), repmat({'trough'}, 1, length(troughCCGactual))];
            featuresCCGtype       = [repmat({'peak'}, 1, length(idxPeakCCG)),    repmat({'trough'}, 1, length(idxTroughCCG))]; 

    %         latPeakCCGactual((latPeakCCGactual < 0.1) | (latPeakCCGactual > 1)) = [];
    %         pairGroupStatTable.featuresCCGactualLat{loopPairs,1} = latPeakCCGactual;
        end
        %% number of pairs
        if strcmp(datasetID,'RatRoy')
            pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons;
            pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP;
            pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII;
            pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI;
        elseif strcmp(datasetID,'RatAGday1')
            if strcmp(pairGroupStatTable.refBrainRegion(loopPairs,:),'CA1')
               pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(1,1);
               pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP(1,1);
               pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(1,1);
               pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(1,1);
            elseif strcmp(pairGroupStatTable.refBrainRegion(loopPairs,:),'CA3')
               pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(4,1);
               pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP(4,1);
               pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(4,1);
               pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(4,1);
            end
        elseif strcmp(datasetID,'RatAGday2')
            pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(5,1);
            pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP(5,1);
            pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(5,1);
            pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(5,1);
        elseif strcmp(datasetID,'RatN')
            pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(1,1);
            pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP(1,1);
            pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(1,1);
            pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(1,1);
        elseif strcmp(datasetID,'RatS')
            pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(2,1) + noPairCompTable.nPairComparisons(3,1);
            pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP(2,1)      + noPairCompTable.nPairCompPP(3,1);
            pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(2,1)      + noPairCompTable.nPairCompII(3,1);
            pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(2,1)  + noPairCompTable.nPairCompIPorPI(3,1);
        elseif strcmp(datasetID,'RatU')
            pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(4,1) + noPairCompTable.nPairComparisons(5,1);
             pairGroupStatTable.nPairCompPP(loopPairs,:)     = noPairCompTable.nPairCompPP(4,1)      + noPairCompTable.nPairCompPP(5,1);
            pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(4,1)      + noPairCompTable.nPairCompII(5,1);
            pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(4,1)  + noPairCompTable.nPairCompIPorPI(5,1);
        elseif strcmp(datasetID,'AllenInstitute')
            idxTemp = find( ...
                           ((noPairCompTable.animalID == pairGroupStatTable.animalID(loopPairs,1)) & ...
                           strcmp(string(noPairCompTable.probeName),pairGroupStatTable.probeID{loopPairs,1})) & ...
                           strcmp(string(noPairCompTable.brainRegion),pairGroupStatTable.brainRegion{loopPairs,1}));
            pairGroupStatTable.nPairComparisons(loopPairs,:) = noPairCompTable.nPairComparisons(idxTemp,1);
            pairGroupStatTable.nPairCompPP(loopPairs,:)      = noPairCompTable.nPairCompPP(idxTemp,1);
            pairGroupStatTable.nPairCompII(loopPairs,:)      = noPairCompTable.nPairCompII(idxTemp,1);
            pairGroupStatTable.nPairCompIPorPI(loopPairs,:)  = noPairCompTable.nPairCompIPorPI(idxTemp,1);
        end
        
        
        %% exquisite detection (not putative GJ)
        if ~GJanaFlag
            if (sum(intersect(85:127,find(pairGroupStatTable.GSPExc{loopPairs,1})) >= 2)  && ...
                sum(intersect(85:127,find(pairGroupStatTable.GSPInh{loopPairs,1})) >= 1)) && ...
               (median(pairGroupStatTable.pairRawCCG{loopPairs,1}(~pairGroupStatTable.GSPExc{loopPairs,1} | ...
                                                                  ~pairGroupStatTable.GSPInh{loopPairs,1})) >= 5)
                pairGroupStatTable.flagExq{loopPairs,1} = true;
            else
                pairGroupStatTable.flagExq{loopPairs,1} = false;
            end

            if (sum(intersect(85:127,find(pairGroupStatTable.GSPExc{loopPairs,1})) >= 2)  && ...
                sum(intersect(85:127,find(pairGroupStatTable.GSPInh{loopPairs,1})) >= 1)) && ...
               (median(pairGroupStatTable.pairRawCCG{loopPairs,1}(~pairGroupStatTable.GSPExc{loopPairs,1} | ...
                                                                  ~pairGroupStatTable.GSPInh{loopPairs,1})) <= 5)
                pairGroupStatTable.flagExqNoBaseline{loopPairs,1} = true;
            else
                pairGroupStatTable.flagExqNoBaseline{loopPairs,1} = false;
            end
        end
%         %% excitatory monosynapse detection
%         if sum(intersect(127:211,idxPeakCCG)) >= 3
%             pairGroupStatTable.flagExcMonSyp{loopPairs,1} = true; 
% %             plot(pairGroupStatTable.CCGbinLagTimes{1,1}*1000,ccg)
%         else
%             pairGroupStatTable.flagExcMonSyp{loopPairs,1} = false; 
%         end
%         
%         %% inhibitory monosynapse detection
%         if sum(intersect(127:211,idxTroughCCG)) >= 3
%             pairGroupStatTable.flagInhMonSyp{loopPairs,1} = true;
% %             plot(pairGroupStatTable.CCGbinLagTimes{1,1}*1000,ccg)
%         else
%             pairGroupStatTable.flagInhMonSyp{loopPairs,1} = false;
%         end
        
        if GJanaFlag
            
            ref      = pairGroupStatTable.refSpikeTimes{loopPairs};
            tar      = pairGroupStatTable.tarSpikeTimes{loopPairs};
            binSize  = 1/fs;
            duration = 0.007;
            
            % make sure it is a column vector
            ref = ref(:);
            tar = tar(:);
            
            [ccgR, tR] = CCG([ref;tar],[ones(size(ref));2*ones(size(tar))], ...
                'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                'norm', 'counts');
            
            ccgHD = ccgR(:,1,2);
            tHD   = tR*1000;
            
            pairGroupStatTable.ccgHD{loopPairs} = ccgHD;
            pairGroupStatTable.tHD{loopPairs}   = tHD;
            
            %% remove symmetry rule Nov/Dec 2024
%             symmetryCorr = corrcoef(ccgHD(106:136),flipud(ccgHD(76:106)));
% 
%             if (symmetryCorr(2,1) > 0.75) && median(ccgHD) >= 5
%                 pairGroupStatTable.flagGJ{loopPairs,1} = true;
% %                 GJpairList = [GJpairList loopPairs];
% %                 plot(pairGroupStatTable.CCGbinLagTimes{1,1}*1000,ccg)
%             else
%                 pairGroupStatTable.flagGJ{loopPairs,1} = false;
%             end  

            %% require only central millisecond bin rule from now on
            
            if (sum(pairGroupStatTable.GSPExc{loopPairs} == [0,0,0,0,1,0,0,0,0]') == 9) && median(ccgHD) >= 5
                pairGroupStatTable.flagGJ{loopPairs,1} = true;
            else
                pairGroupStatTable.flagGJ{loopPairs,1} = false;
            end  
        else %% gap junction detection for Exq screen (old!! and not putative GJ)
            if ~isempty(intersect(104:108,idxPeakCCG)) %&& ((strcmp(refCellType(1),'i') && strcmp(tarCellType(1),'i')) || (strcmp(refCellType(1),'p') && strcmp(tarCellType(1),'p')))
    %             % Calculate the sortedness index
    %             differences = diff(-pairGroupStatTable.jitterMean{loopPairs,1}(106:136));
    %             sortednessIndexRHS = sum(differences > 0) / length(differences);
    % 
    %             differences = diff(-flipud(pairGroupStatTable.jitterMean{loopPairs,1}(76:106)));
    %             sortednessIndexLHS = sum(differences > 0) / length(differences);
                
                % symmetry measure
                symmetryCorr = corrcoef(ccg(106:136),flipud(ccg(76:106)));

                if (symmetryCorr(2,1) > 0.9) %&& ((sortednessIndexRHS > 0.75) && (sortednessIndexLHS > 0.75))
                    pairGroupStatTable.flagGJ{loopPairs,1} = true;
    %                 GJpairList = [GJpairList loopPairs];
    %                 plot(pairGroupStatTable.CCGbinLagTimes{1,1}*1000,ccg)
                else
                    pairGroupStatTable.flagGJ{loopPairs,1} = false;
                end
            else
                pairGroupStatTable.flagGJ{loopPairs,1} = false;    
            end
        end
         
%         %% find aligned channel
%         pairGroupStatTable.flagAligned(loopPairs,1) = false; 
%         
%         for loopChannels = 1:size(refWaveAtTarShank,1)
%             
%             if isempty([latPeakWaveList{loopChannels} latTroughWaveList{loopChannels}])
%                 continue
%             elseif isempty(latPeakWaveList{loopChannels}) && ~isempty(latTroughWaveList{loopChannels})
%                 
%                 featuresWaveIdx = idxTroughWaveList{loopChannels};
%                 featuresWaveLat = latTroughWaveList{loopChannels};
%                 featuresWaveVal = refWaveAtTarShank(loopChannels,idxTroughWaveList{loopChannels});
%                 
%             elseif ~isempty(latPeakWaveList{loopChannels}) && isempty(latTroughWaveList{loopChannels})
%                
%                 featuresWaveIdx = idxPeakWaveList{loopChannels};
%                 featuresWaveLat = latPeakWaveList{loopChannels};
%                 featuresWaveVal = refWaveAtTarShank(loopChannels,idxPeakWaveList{loopChannels});
%                 
%             else
%                 
%                 featuresWaveIdx = [idxPeakWaveList{loopChannels} idxTroughWaveList{loopChannels}];
%                 featuresWaveLat = [latPeakWaveList{loopChannels} latTroughWaveList{loopChannels}];
%                 featuresWaveVal = refWaveAtTarShank(loopChannels,[idxPeakWaveList{loopChannels} idxTroughWaveList{loopChannels}]);
%                 
%             end
%             
%             % alignment logic 
%             
%             if ~isempty(intersect(featuresCCGlat,featuresWaveLat))
% %                 intersect(featuresCCGlat,featuresWaveLat) 
%                 
%                 tiledlayout(2,1)
%                 
% %                 nexttile
% %                 plot(pairGroupStatTable.CCGbinLagTimes{1,1}*1000,ccg)
% %                 hold on 
% %                 plot([latPeakCCG; latTroughCCG],[peakCCG; troughCCG],'o')
% %                 hold off
% %                 xlim([-3.5 3.5])
% %                 
% %                 nexttile
% %                 plot(pairGroupStatTable.tWave{1,1},refWaveAtTarShank(loopChannels,:))
% %                 hold on
% %                 plot(featuresWaveLat,featuresWaveVal,'o')
% %                 hold off
% %                 xlim([-3.5 3.5])
%                 
%             end
%             
%             chansAligned       = [chansAligned   loopChannels];
%             numFeatAligned     = [numFeatAligned length(intersect(featuresCCGlat,featuresWaveLat))];
%             waveAmpAligned     = [waveAmpAligned max(abs(featuresWaveVal))];
% 
%         end
%         
%         % no alignment so skip to next pair (at two features required)
%         if isempty(numFeatAligned)
%             continue
%         elseif max(numFeatAligned) < 1
%             continue
%         end
%         
%         chansAlignedFilt   = chansAligned(numFeatAligned   == max(numFeatAligned));
%         waveAmpAlignedFilt = waveAmpAligned(numFeatAligned == max(numFeatAligned));
% 
%         chansAlignedFinal  = chansAlignedFilt(waveAmpAlignedFilt == max(waveAmpAlignedFilt));
%         
%         find wave peaks and troughs in each channel
%         
%         chansAligned   = [];
%         numFeatAligned = [];
%         waveAmpAligned = [];
%         
%         if strcmp(datasetID,'AllenInstitute') || strcmp(datasetID,'Steinmetz')
%             thresh = max(abs(refWaveAtTarShank(pairGroupStatTable.refChannel(loopPairs,:),:)))/100;
%         else
%             thresh = max(abs(pairGroupStatTable.refWaveAtRefChan{loopPairs,1}))/100;
%         end
%         
%         for loopChannels = 1:size(refWaveAtTarShank,1)
%             
%             [~,idx]    = findpeaks(refWaveAtTarShank(loopChannels,:), ...
%                                                     'MinPeakDistance',12, ... % minimum 400 miscrosec
%                                                     'MinPeakHeight',thresh);
%             
%             % analyze waveform alignments between -0.25 and 1ms                                     
%             idx((pairGroupStatTable.tWave{1,1}(idx) <= -0.25) | (pairGroupStatTable.tWave{1,1}(idx) >= 1)) = [];
%             
%             % SNR filter 
%             idx(abs((refWaveAtTarShank(loopChannels,idx) - mean(refWaveAtTarShank(loopChannels,:)))/std(refWaveAtTarShank(loopChannels,:))) < 1) = [];
%             % remove noise spikes in rat N's distal STA
%             if strcmp(datasetID,'RatN')
%                 absDiff = abs(diff(refWaveAtTarShank(loopChannels,:)));
%                 idx(find(((absDiff(idx-1) - mean(absDiff))/std(absDiff)) > 2)) = [];
%             end
%             
%             if ~isempty(idx)
%                 latPeakWaveList{loopChannels} = pairGroupStatTable.tWave{1,1}(idx);
%                 idxPeakWaveList{loopChannels} = idx;
%             else
%                 latPeakWaveList{loopChannels} = {};
%                 idxPeakWaveList{loopChannels} = {};
%             end
%             
%             [~,idx]    = findpeaks(-refWaveAtTarShank(loopChannels,:), ...
%                                                      'MinPeakDistance',12, ... % minimum 400 miscrosec
%                                                      'MinPeakHeight',thresh);
%                                                  
%             % analyze waveform alignments between -0.25 and 1ms                                     
%             idx((pairGroupStatTable.tWave{1,1}(idx) <= -0.25) | (pairGroupStatTable.tWave{1,1}(idx) >= 1)) = [];
%             
%             % SNR filter 
%             idx(abs((refWaveAtTarShank(loopChannels,idx) - mean(refWaveAtTarShank(loopChannels,:)))/std(refWaveAtTarShank(loopChannels,:))) < 1) = [];
%             % remove noise spikes in rat N's distal STA
%             if strcmp(datasetID,'RatN')
%                 absDiff = abs(diff(refWaveAtTarShank(loopChannels,:)));
%                 idx(find(((absDiff(idx-1) - mean(absDiff))/std(absDiff)) > 2)) = [];
%             end
%             
%             if ~isempty(idx)
%                 latTroughWaveList{loopChannels} = pairGroupStatTable.tWave{1,1}(idx);
%                 idxTroughWaveList{loopChannels} = idx;
%             else
%                 latTroughWaveList{loopChannels} = {};
%                 idxTroughWaveList{loopChannels} = {};
%             end
%             
%         end
%         

%         
%         %% find aligned features
%                 
%         if ~isempty(latPeakWaveList{chansAlignedFinal}) && isempty(latTroughWaveList{chansAlignedFinal})
%             featuresWaveLatActual = latPeakWaveList{chansAlignedFinal};
%         elseif isempty(latPeakWaveList{chansAlignedFinal}) && ~isempty(latTroughWaveList{chansAlignedFinal})
%             featuresWaveLatActual = latTroughWaveList{chansAlignedFinal};
%         elseif ~isempty(latPeakWaveList{chansAlignedFinal}) && ~isempty(latTroughWaveList{chansAlignedFinal})
%             featuresWaveLatActual = [latPeakWaveList{chansAlignedFinal} latTroughWaveList{chansAlignedFinal}];
%         end
%         if ~isempty(idxPeakWaveList{chansAlignedFinal}) && isempty(idxTroughWaveList{chansAlignedFinal})
%             featuresWaveVal  = refWaveAtTarShank(chansAlignedFinal,idxPeakWaveList{chansAlignedFinal});
%         elseif isempty(idxPeakWaveList{chansAlignedFinal}) && ~isempty(idxTroughWaveList{chansAlignedFinal})
%             featuresWaveVal  = refWaveAtTarShank(chansAlignedFinal,idxTroughWaveList{chansAlignedFinal});
%         elseif ~isempty(idxPeakWaveList{chansAlignedFinal}) && ~isempty(idxTroughWaveList{chansAlignedFinal})
%             featuresWaveVal  = [refWaveAtTarShank(chansAlignedFinal,idxPeakWaveList{chansAlignedFinal}) ...
%                                 refWaveAtTarShank(chansAlignedFinal,idxTroughWaveList{chansAlignedFinal})];
%         end
%         featuresWaveType = [repmat({'peak'}, 1, length(latPeakWaveList{chansAlignedFinal})), repmat({'trough'}, 1, length(latTroughWaveList{chansAlignedFinal}))];
%         
%         % ALIGNMENT DETECTION
%         [~,CCGalignmentIdx,WaveAlignmentIdx] = intersect(featuresCCGlat,featuresWaveLatActual);
%         
%         if (length(CCGalignmentIdx) < 2) || (length(WaveAlignmentIdx) < 2) % (length(CCGalignmentIdx) < 1) || (length(WaveAlignmentIdx) < 1)
%             continue
%         elseif  (sum(strcmp({featuresCCGtype{CCGalignmentIdx}},{featuresWaveType{WaveAlignmentIdx}})) == length(CCGalignmentIdx)) || ...
%                 ((length(CCGalignmentIdx) - sum(strcmp({featuresCCGtype{CCGalignmentIdx}},{featuresWaveType{WaveAlignmentIdx}}))) < 2)
%             continue
%         end
%         
% %         featuresCCGactualLat(featuresCCGactualGroup(CCGalignmentIdx));
% %         featuresCCGactualVal(featuresCCGactualGroup(CCGalignmentIdx));
% %         featuresCCGtype(CCGalignmentIdx);
% %         
% %         featuresWaveLatActual(WaveAlignmentIdx)
% %         featuresWaveVal(WaveAlignmentIdx)
% %         featuresWaveType(WaveAlignmentIdx)
%         
%         % at this point alignment is successful 
%         pairGroupStatTable.flagAligned(loopPairs,1) = true; 
%         
%         % alignment stats!!
%         pairGroupStatTable.channelOfAlignment(loopPairs,1)     = chansAlignedFinal;
%         pairGroupStatTable.CCGAlignedLatencies{loopPairs,1}    = featuresCCGactualLat( featuresCCGactualGroup(CCGalignmentIdx));
%         pairGroupStatTable.CCGAlignedAmplitude{loopPairs,1}    = featuresCCGactualVal( featuresCCGactualGroup(CCGalignmentIdx));
%         pairGroupStatTable.CCGAlignedFeatureType{loopPairs,1}  = featuresCCGactualType(featuresCCGactualGroup(CCGalignmentIdx));
%         pairGroupStatTable.WaveAlignedLatencies{loopPairs,1}   = featuresWaveLatActual(WaveAlignmentIdx);
%         pairGroupStatTable.WaveAlignedAmplitude{loopPairs,1}   = featuresWaveVal(      WaveAlignmentIdx);
%         pairGroupStatTable.WaveAlignedFeatureType{loopPairs,1} = featuresWaveType(     WaveAlignmentIdx);
%         
%         tiledlayout(2,1)
%                 
%         nexttile
%         plot(pairGroupStatTable.CCGbinLagTimes{1,1}*1000,ccg)
%         hold on 
%         plot(featuresCCGactualLat(featuresCCGactualGroup(CCGalignmentIdx)), ...
%              featuresCCGactualVal(featuresCCGactualGroup(CCGalignmentIdx)),'o')
%         text(featuresCCGactualLat(featuresCCGactualGroup(CCGalignmentIdx)), ...
%              featuresCCGactualVal(featuresCCGactualGroup(CCGalignmentIdx)), ...
%              featuresCCGactualType(featuresCCGactualGroup(CCGalignmentIdx)))
%         hold off
%         xlim([-3.5 3.5])
% %         titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
%         title(num2str(pairGroupStatTable.pairDistance(loopPairs,1)))
% %         title(titleChar)
% 
%         nexttile
%         plot(pairGroupStatTable.tWave{1,1},refWaveAtTarShank(chansAlignedFinal,:))
%         hold on
%         plot(featuresWaveLatActual(WaveAlignmentIdx), ...
%              featuresWaveVal(      WaveAlignmentIdx),'o')
%         text(featuresWaveLatActual(WaveAlignmentIdx), ...
%              featuresWaveVal(      WaveAlignmentIdx), ...
%              featuresWaveType(     WaveAlignmentIdx))
%         hold off
%         xlim([-3.5 3.5])
%         
%         drawnow
        
%         save_file = fullfile(datapath, titleChar);
%         print(fig_use, save_file,'-dpng',resolution_use);

        
        %% test opposite max chan
        
        
        
%         idxPeakWaveList{chansAlignedFinal} 
%         idxTroughWaveList{chansAlignedFinal}
       
%         if strcmp(datasetID,'Steinmetz')
%             ccgRange  = (96:136)'; % 1.25msec (-0.25 to 1ms lag)
%             waveRange = (31:71)';  % 1.25msec (-0.25 to 1ms lag)
%         elseif strcmp(datasetID,'AllenInstitute')
% %             ccgRange  = (96:136)'; % 1.25msec (-0.25 to 1ms lag)
% %             waveRange = (10:50)';  % 1.25msec (-0.25 to 1ms lag)
%             ccgRange  = (87:168)'; % 1.25msec (-0.25 to 1ms lag)
%             waveRange = (1:82)';  % 1.25msec (-0.25 to 1ms lag)
%         else
%             ccgRange  = (96:136)'; % 1.25msec (-0.25 to 1ms lag)
%             waveRange = (26:66)';  % 1.25msec (-0.33 to 1ms lag)
%         end
% 
%         RsqList = [];
%         for k = 1:(size(refWaveAtRefShank,1) + size(refWaveAtTarShank,1))
%             if k <= size(refWaveAtRefShank,1)
%                 wave4Ana = refWaveAtRefShank(k,:)';
%             elseif k > size(refWaveAtRefShank,1)
%                 wave4Ana = refWaveAtTarShank(k-size(refWaveAtRefShank,1),:)';
%             end
%             
%             R = corrcoef(ccg(ccgRange),wave4Ana(waveRange));
%             
%             if R(1,2) >= 0
%                 RsqList(k) = 0;
%             elseif R(1,2) < 0
%                 RsqList(k) = R(1,2)^2;
%             end
%             
%         end 
% 
%         [RsqMax,maxMatchChIdx] = max(RsqList);
%         
%         % don't include for effect group
%         if RsqMax < 0.5
%             
% %             disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
%             disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
%             
%             continue 
%         end
%         
%         allChans   = [channelSetRef channelSetTar];
%         maxMatchCh = allChans(maxMatchChIdx);
%         
%         % distance calc
% %         xSource   = probe.x(probe.chanNo == chCell1);
% %         ySource   = probe.y(probe.chanNo == chCell1);
% %         
% %         xReceiver = probe.x(probe.chanNo == maxMatchCh);
% %         yReceiver = probe.y(probe.chanNo == maxMatchCh);
% %         
% %         distance  = (((xSource - xReceiver)^2) + ((ySource - yReceiver)^2))^(1/2);
% 
%         distance = pairGroupStatTable.pairDistance(loopPairs);
%         
%         % waveform for analysis based on max channel
%         if maxMatchChIdx <= size(refWaveAtRefShank,1)
%             wave4Ana = refWaveAtRefShank(maxMatchChIdx,:)';
%         elseif maxMatchChIdx > size(refWaveAtRefShank,1)
%             wave4Ana = refWaveAtTarShank(maxMatchChIdx-size(refWaveAtRefShank,1),:)';
%         end
%         
%         %% find peaks and troughs CCG
%         if strcmp(anaDirection,'forward')
%             peaksMask     = pairGroupStatTable.GSPExc{loopPairs,1};         peaksMask     = peaksMask(ccgRange);
%             troughsMask   = pairGroupStatTable.GSPInh{loopPairs,1};         troughsMask   = troughsMask(ccgRange);
%             ccgjmAnalysis = pairGroupStatTable.jitterMean{loopPairs,1};     ccgjmAnalysis = ccgjmAnalysis(ccgRange);
% 
%             refCellType   = pairGroupStatTable.refCellExplorerType{loopPairs,1};
%             tarCellType   = pairGroupStatTable.tarCellExplorerType{loopPairs,1};
%         elseif strcmp(anaDirection,'reverse')
%             peaksMask     = flipud(pairGroupStatTable.GSPExc{loopPairs,1});         peaksMask     = peaksMask(ccgRange);
%             troughsMask   = flipud(pairGroupStatTable.GSPInh{loopPairs,1});         troughsMask   = troughsMask(ccgRange);
%             ccgjmAnalysis = flipud(pairGroupStatTable.jitterMean{loopPairs,1});     ccgjmAnalysis = ccgjmAnalysis(ccgRange);
%             
%             tarCellType   = pairGroupStatTable.refCellExplorerType{loopPairs,1};
%             refCellType   = pairGroupStatTable.tarCellExplorerType{loopPairs,1};
%         end
%         
%         tR            = pairGroupStatTable.CCGbinLagTimes{loopPairs,1}*1000;
%         tRanalysis    = pairGroupStatTable.CCGbinLagTimes{loopPairs,1}(ccgRange);
%         ccgAnalysis   = ccg(ccgRange);
%         tWave         = pairGroupStatTable.tWave{loopPairs,1}';
%         
%         [pkCCG, pklocsCCG] = findpeaks(ccg(ccgRange));
%         [trCCG, trlocsCCG] = findpeaks(-ccg(ccgRange));
%         trCCG = -trCCG;
% 
%         % screen with mask
%         pklocsCCGmask = zeros(length(peaksMask),1);
%         pklocsCCGmask(pklocsCCG) = 1;
% 
%         trlocsCCGmask = zeros(length(troughsMask),1);
%         trlocsCCGmask(trlocsCCG) = 1;
% 
%         pklocsCCGmask = peaksMask   & pklocsCCGmask;
%         trlocsCCGmask = troughsMask & trlocsCCGmask;
%         
%         % reject pairs with no significant bins in the analysis range
%         if ~((sum(pklocsCCGmask(10:12)) > 0) && (sum(pklocsCCGmask([1:9,13:41])) > 0))
%             
% %                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
%             disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
% 
%             continue
%         end
%         
%         % identify multipeaks
%         if (sum(pklocsCCGmask(10:12)) > 0) && (sum(pklocsCCGmask([1:9,13:41])) > 0)
%             
%             % zero lag peak
%             pklocsCCGmask_0lag        = zeros(length(pklocsCCGmask),1);
%             pklocsCCGmask_0lag(10:12) = pklocsCCGmask(10:12);
%             
%             [mainPk, ...
%              mainTr, ...
%              mainPkLat, ...
%              mainTrLat, ...
%              preFeatureCCGlat, ...
%              preFeatureWAVElat, ...
%              postFeatureCCGlat, ...
%              postFeatureWAVElat, ...
%              preFeatCCGidx, ...
%              postFeatCCGidx, ...
%              preFeatWAVEidx, ... 
%              postFeatWAVEidx, ...
%              preFeatCCG, ...
%              postFeatCCG, ...
%              preFeatWAVE, ... 
%              postFeatWAVE, ...
%              mainPkNorm, ...
%              preFeatCCGnorm, ...
%              postFeatCCGnorm] = featureLatencies(pklocsCCGmask_0lag, ...
%                                                      trlocsCCGmask, ...
%                                                      ccgAnalysis, ...
%                                                      ccgjmAnalysis, ...
%                                                      tRanalysis, ...
%                                                      wave4Ana, ...
%                                                      waveRange);
%             
%             %% pairing CCG an waveform features
%             if ~isempty(mainPk) && ~isempty(mainTr)
%                 
%                 groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
%                 groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
%                 groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
%    
%                 groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
%                 groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
%                 
%                 groupStat.mainDistance = [groupStat.mainDistance distance];
%                 
%                 groupStat.mainRef = [groupStat.mainRef refID];
%                 groupStat.mainTar = [groupStat.mainTar tarID];
%                 
%                 groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
%                 groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];
% 
%                 if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
%                     
%                     groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
%                     groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
%                     groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
%                     
%                     groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
%                     groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
%                      
%                     groupStat.preDistance = [groupStat.preDistance distance];
%                      
%                     groupStat.preRef = [groupStat.preRef refID];
%                     groupStat.preTar = [groupStat.preTar tarID];
% 
%                     groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
%                     groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];
% 
%                 end
% 
%                 if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
%                     
%                     groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
%                     groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
%                     groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
%                     
%                     groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
%                     groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
%                     
%                     groupStat.postDistance = [groupStat.postDistance distance];
%                     
%                     groupStat.postRef = [groupStat.postRef refID];
%                     groupStat.postTar = [groupStat.postTar tarID];
% 
%                     groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
%                     groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];
% 
%                 end
%                 
%                 yyaxis left
%                 plot(tR,ccg,'b')
%                 hold on
%                 scatter(mainPkLat*1000,mainPk,'b','LineWidth',2)
%                 scatter(preFeatureCCGlat*1000,preFeatCCG,'b','LineWidth',2)
%                 scatter(postFeatureCCGlat*1000,postFeatCCG,'b','LineWidth',2)
%                 hold off
% 
%                 yyaxis right
%                 plot(tWave,wave4Ana,'r')
%                 hold on
%                 scatter(mainTrLat*1000,mainTr,'r','LineWidth',2)
%                 scatter(preFeatureWAVElat*1000,preFeatWAVE,'r','LineWidth',2)
%                 scatter(postFeatureWAVElat*1000,postFeatWAVE,'r','LineWidth',2)
%                 hold off
%                 set(gca, 'YDir','reverse')
%                 xlim([-3.5 3.5])
% 
%                 titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
%                 title(titleChar)
% 
% %                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
%                 disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
% 
%                 save_file = fullfile(datapath, titleChar);
%                 print(fig_use, save_file,'-dpng',resolution_use);
%                 
%             end
%             
%             % lagged peaks
%             pklocsCCGmask_non_0lag              = zeros(length(pklocsCCGmask),1);
%             pklocsCCGmask_non_0lag([1:9,13:41]) = pklocsCCGmask([1:9,13:41]);
%             
%             [mainPk, ...
%              mainTr, ...
%              mainPkLat, ...
%              mainTrLat, ...
%              preFeatureCCGlat, ...
%              preFeatureWAVElat, ...
%              postFeatureCCGlat, ...
%              postFeatureWAVElat, ...
%              preFeatCCGidx, ...
%              postFeatCCGidx, ...
%              preFeatWAVEidx, ... 
%              postFeatWAVEidx, ...
%              preFeatCCG, ...
%              postFeatCCG, ...
%              preFeatWAVE, ... 
%              postFeatWAVE, ...
%              mainPkNorm, ...
%              preFeatCCGnorm, ...
%              postFeatCCGnorm] = featureLatencies(pklocsCCGmask_non_0lag, ...
%                                                  trlocsCCGmask, ...
%                                                  ccgAnalysis, ...
%                                                  ccgjmAnalysis, ...
%                                                  tRanalysis, ...
%                                                  wave4Ana, ...
%                                                  waveRange);
%             
%              %% pairing CCG an waveform features
%             if ~isempty(mainPk) && ~isempty(mainTr)
%                 
%                 groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
%                 groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
%                 groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
%    
%                 groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
%                 groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
%                 
%                 groupStat.mainDistance = [groupStat.mainDistance distance];
%                 
%                 groupStat.mainRef = [groupStat.mainRef refID];
%                 groupStat.mainTar = [groupStat.mainTar tarID];
%                 
%                 groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
%                 groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];
% 
%                 if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
%                     
%                     groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
%                     groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
%                     groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
%                     
%                     groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
%                     groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
%                      
%                     groupStat.preDistance = [groupStat.preDistance distance];
%                      
%                     groupStat.preRef = [groupStat.preRef refID];
%                     groupStat.preTar = [groupStat.preTar tarID];
% 
%                     groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
%                     groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];
% 
%                 end
% 
%                 if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
%                     
%                     groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
%                     groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
%                     groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
%                     
%                     groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
%                     groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
%                     
%                     groupStat.postDistance = [groupStat.postDistance distance];
%                     
%                     groupStat.postRef = [groupStat.postRef refID];
%                     groupStat.postTar = [groupStat.postTar tarID];
% 
%                     groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
%                     groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];
% 
%                 end
%                 
%                 yyaxis left
%                 plot(tR,ccg,'b')
%                 hold on
%                 scatter(mainPkLat*1000,mainPk,'b','LineWidth',2)
%                 scatter(preFeatureCCGlat*1000,preFeatCCG,'b','LineWidth',2)
%                 scatter(postFeatureCCGlat*1000,postFeatCCG,'b','LineWidth',2)
%                 hold off
% 
%                 yyaxis right
%                 plot(tWave,wave4Ana,'r')
%                 hold on
%                 scatter(mainTrLat*1000,mainTr,'r','LineWidth',2)
%                 scatter(preFeatureWAVElat*1000,preFeatWAVE,'r','LineWidth',2)
%                 scatter(postFeatureWAVElat*1000,postFeatWAVE,'r','LineWidth',2)
%                 hold off
%                 set(gca, 'YDir','reverse')
%                 xlim([-3.5 3.5])
% 
%                 titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
%                 title(titleChar)
% 
% %                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
%                 disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
% 
%                 save_file = fullfile(datapath, titleChar);
%                 print(fig_use, save_file,'-dpng',resolution_use);
%                 
%             end
%             
%         else 
%             
%             [mainPk, ...
%              mainTr, ...
%              mainPkLat, ...
%              mainTrLat, ...
%              preFeatureCCGlat, ...
%              preFeatureWAVElat, ...
%              postFeatureCCGlat, ...
%              postFeatureWAVElat, ...
%              preFeatCCGidx, ...
%              postFeatCCGidx, ...
%              preFeatWAVEidx, ... 
%              postFeatWAVEidx, ...
%              preFeatCCG, ...
%              postFeatCCG, ...
%              preFeatWAVE, ... 
%              postFeatWAVE, ...
%              mainPkNorm, ...
%              preFeatCCGnorm, ...
%              postFeatCCGnorm] = featureLatencies(pklocsCCGmask, ...
%                                                      trlocsCCGmask, ...
%                                                      ccgAnalysis, ...
%                                                      ccgjmAnalysis, ...
%                                                      tRanalysis, ...
%                                                      wave4Ana, ...
%                                                      waveRange);
%             
%              %% pairing CCG an waveform features
%              if ~isempty(mainPk) && ~isempty(mainTr)
%                 
%                 groupStat.mainFeatCCGall  = [groupStat.mainFeatCCGall  mainPk];
%                 groupStat.mainFeatWAVEall = [groupStat.mainFeatWAVEall mainTr];
%                 groupStat.mainFeatCCGnormAll = [groupStat.mainFeatCCGnormAll mainPkNorm];
%    
%                 groupStat.mainFeatCCGlatAll  = [groupStat.mainFeatCCGlatAll  mainPkLat];
%                 groupStat.mainFeatWAVElatAll = [groupStat.mainFeatWAVElatAll mainTrLat];
%                 
%                 groupStat.mainDistance = [groupStat.mainDistance distance];
%                 
%                 groupStat.mainRef = [groupStat.mainRef refID];
%                 groupStat.mainTar = [groupStat.mainTar tarID];
%                 
%                 groupStat.mainRsqMax  = [groupStat.mainRsqMax  RsqMax];
%                 groupStat.mainMaxChan = [groupStat.mainMaxChan maxMatchCh];
% 
%                 if ~isempty(preFeatCCGidx) && ~isempty(preFeatWAVEidx)
%                     
%                     groupStat.preFeatCCGall   = [groupStat.preFeatCCGall  preFeatCCG];
%                     groupStat.preFeatWAVEall  = [groupStat.preFeatWAVEall preFeatWAVE];
%                     groupStat.preFeatCCGnormAll  = [groupStat.preFeatCCGnormAll preFeatCCGnorm];
%                     
%                     groupStat.preFeatCCGlatAll   = [groupStat.preFeatCCGlatAll  preFeatureCCGlat];
%                     groupStat.preFeatWAVElatAll  = [groupStat.preFeatWAVElatAll preFeatureWAVElat];
%                      
%                     groupStat.preDistance = [groupStat.preDistance distance];
%                      
%                     groupStat.preRef = [groupStat.preRef refID];
%                     groupStat.preTar = [groupStat.preTar tarID];
% 
%                     groupStat.preRsqMax  = [groupStat.preRsqMax  RsqMax];
%                     groupStat.preMaxChan = [groupStat.preMaxChan maxMatchCh];
% 
%                 end
% 
%                 if ~isempty(postFeatCCGidx) && ~isempty(postFeatWAVEidx)
%                     
%                     groupStat.postFeatCCGall  = [groupStat.postFeatCCGall  postFeatCCG];
%                     groupStat.postFeatWAVEall = [groupStat.postFeatWAVEall postFeatWAVE];
%                     groupStat.postFeatCCGnormAll  = [groupStat.postFeatCCGnormAll postFeatCCGnorm];
%                     
%                     groupStat.postFeatCCGlatAll  = [groupStat.postFeatCCGlatAll  postFeatureCCGlat];
%                     groupStat.postFeatWAVElatAll = [groupStat.postFeatWAVElatAll postFeatureWAVElat];
%                     
%                     groupStat.postDistance = [groupStat.postDistance distance];
%                     
%                     groupStat.postRef = [groupStat.postRef refID];
%                     groupStat.postTar = [groupStat.postTar tarID];
% 
%                     groupStat.postRsqMax  = [groupStat.postRsqMax RsqMax];
%                     groupStat.postMaxChan = [groupStat.postMaxChan maxMatchCh];
% 
%                 end
%                 
%                 yyaxis left
%                 plot(tR,ccg,'b')
%                 hold on
%                 scatter(mainPkLat*1000,mainPk,'b','LineWidth',2)
%                 scatter(preFeatureCCGlat*1000,preFeatCCG,'b','LineWidth',2)
%                 scatter(postFeatureCCGlat*1000,postFeatCCG,'b','LineWidth',2)
%                 hold off
% 
%                 yyaxis right
%                 plot(tWave,wave4Ana,'r')
%                 hold on
%                 scatter(mainTrLat*1000,mainTr,'r','LineWidth',2)
%                 scatter(preFeatureWAVElat*1000,preFeatWAVE,'r','LineWidth',2)
%                 scatter(postFeatureWAVElat*1000,postFeatWAVE,'r','LineWidth',2)
%                 hold off
%                 set(gca, 'YDir','reverse')
%                 xlim([-3.5 3.5])
% 
%                 titleChar = ['CCG and wave feature alighment - ' num2str(refID) refCellType ' v ' num2str(tarID) tarCellType];
%                 title(titleChar)
% 
% %                 disp([num2str((loopPairs/size(pairGroupStatTable,1))*100) ' % done!'])
%                 disp([num2str(loopPairs) '/' num2str(size(pairGroupStatTable,1)) ' done!'])
%                 
%                 save_file = fullfile(datapath, titleChar);
%                 print(fig_use, save_file,'-dpng',resolution_use);
%                 
%             end
%          
%         end
%                 
%         close all
%         
%         hcomb = figure(102);
%         arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% 
    end
    
    %% data struct for group stats
    if GJanaFlag
        if strcmp(datasetID,'AllenInstitute')
            save([datapath datasetID num2str(animalsList(animalIdx)) 'putativeGJprocessedGroupStats.mat'],'pairGroupStatTable')
        else
            save([datapath datasetID 'putativeGJprocessedGroupStats.mat'],'pairGroupStatTable')
        end
    else
        if strcmp(datasetID,'AllenInstitute')
            save([datapath datasetID num2str(animalsList(animalIdx)) 'processedGroupStats.mat'],'pairGroupStatTable')
        else
            save([datapath datasetID 'processedGroupStats.mat'],'pairGroupStatTable')
        end
    end
        
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

function pairGroupStatTable = addReversePairs(datasetID,pairGroupStatTable)
    
    reversePairGroupStatTable = table;

    if strcmp(datasetID,'RatRoy')

        reversePairGroupStatTable.brainRegion             = pairGroupStatTable.brainRegion;
        reversePairGroupStatTable.refNeuronID             = pairGroupStatTable.tarNeuronID;
        reversePairGroupStatTable.tarNeuronID             = pairGroupStatTable.refNeuronID;
        reversePairGroupStatTable.refSpikeTimes           = pairGroupStatTable.tarSpikeTimes;
        reversePairGroupStatTable.tarSpikeTimes           = pairGroupStatTable.refSpikeTimes;
        reversePairGroupStatTable.refChannel              = pairGroupStatTable.tarChannel;
        reversePairGroupStatTable.tarChannel              = pairGroupStatTable.refChannel;
        reversePairGroupStatTable.refShank                = pairGroupStatTable.tarShank;
        reversePairGroupStatTable.tarShank                = pairGroupStatTable.refShank;
        reversePairGroupStatTable.channelSetRef           = pairGroupStatTable.channelSetTar;
        reversePairGroupStatTable.channelSetTar           = pairGroupStatTable.channelSetRef;
        reversePairGroupStatTable.refCellExplorerType     = pairGroupStatTable.tarCellExplorerType;
        reversePairGroupStatTable.tarCellExplorerType     = pairGroupStatTable.refCellExplorerType;
        reversePairGroupStatTable.refWaveforms            = pairGroupStatTable.tarWaveforms;
        reversePairGroupStatTable.tarWaveforms            = pairGroupStatTable.refWaveforms;
        reversePairGroupStatTable.tWave                   = pairGroupStatTable.tWave;
        reversePairGroupStatTable.pairDistance            = pairGroupStatTable.pairDistance;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.pairRawCCG{i,1}     = flipud(pairGroupStatTable.pairRawCCG{i,1});
        end
        reversePairGroupStatTable.refRawACG               = pairGroupStatTable.tarRawACG;
        reversePairGroupStatTable.tarRawACG               = pairGroupStatTable.refRawACG;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.jitterMean{i,1}     = flipud(pairGroupStatTable.jitterMean{i,1});
            reversePairGroupStatTable.CCGbinLagTimes{i,1} = flipud(pairGroupStatTable.CCGbinLagTimes{i,1});
            reversePairGroupStatTable.GSPExc{i,1}         = flipud(pairGroupStatTable.GSPExc{i,1});
            reversePairGroupStatTable.GSPInh{i,1}         = flipud(pairGroupStatTable.GSPInh{i,1});
            reversePairGroupStatTable.pvalE{i,1}          = flipud(pairGroupStatTable.pvalE{i,1});
            reversePairGroupStatTable.pvalI{i,1}          = flipud(pairGroupStatTable.pvalI{i,1});
            reversePairGroupStatTable.LSPExc{i,1}         = flipud(pairGroupStatTable.LSPExc{i,1});
            reversePairGroupStatTable.LSPInh{i,1}         = flipud(pairGroupStatTable.LSPInh{i,1});
        end
        reversePairGroupStatTable.refSnippetsAtRefChan    = pairGroupStatTable.tarSnippetsAtTarChan;
        reversePairGroupStatTable.refSnippetsAtTarChan    = pairGroupStatTable.tarSnippetsAtRefChan;
        reversePairGroupStatTable.tarSnippetsAtTarChan    = pairGroupStatTable.refSnippetsAtRefChan;
        reversePairGroupStatTable.tarSnippetsAtRefChan    = pairGroupStatTable.refSnippetsAtTarChan;
        reversePairGroupStatTable.refWaveAtRefChan        = pairGroupStatTable.tarWaveAtTarChan;
        reversePairGroupStatTable.refWaveAtTarChan        = pairGroupStatTable.tarWaveAtRefChan;
        reversePairGroupStatTable.tarWaveAtTarChan        = pairGroupStatTable.refWaveAtRefChan;
        reversePairGroupStatTable.tarWaveAtRefChan        = pairGroupStatTable.refWaveAtTarChan;
        reversePairGroupStatTable.refWaveAtRefShank       = pairGroupStatTable.tarWaveAtTarShank;
        reversePairGroupStatTable.refWaveAtTarShank       = pairGroupStatTable.tarWaveAtRefShank;
        reversePairGroupStatTable.tarWaveAtTarShank       = pairGroupStatTable.refWaveAtRefShank;
        reversePairGroupStatTable.tarWaveAtRefShank       = pairGroupStatTable.refWaveAtTarShank;

        pairGroupStatTable = [pairGroupStatTable; reversePairGroupStatTable];
        
    elseif strcmp(datasetID,'RatAGday1') || strcmp(datasetID,'RatAGday2')
        
        reversePairGroupStatTable.refBrainRegion          = pairGroupStatTable.tarBrainRegion;
        reversePairGroupStatTable.tarBrainRegion          = pairGroupStatTable.refBrainRegion;
        reversePairGroupStatTable.refClusterID            = pairGroupStatTable.tarClusterID;
        reversePairGroupStatTable.tarClusterID            = pairGroupStatTable.refClusterID;
        reversePairGroupStatTable.refNeuronID             = pairGroupStatTable.tarNeuronID;
        reversePairGroupStatTable.tarNeuronID             = pairGroupStatTable.refNeuronID;
        reversePairGroupStatTable.refSpikeTimes           = pairGroupStatTable.tarSpikeTimes;
        reversePairGroupStatTable.tarSpikeTimes           = pairGroupStatTable.refSpikeTimes;
        reversePairGroupStatTable.refChannel              = pairGroupStatTable.tarChannel;
        reversePairGroupStatTable.tarChannel              = pairGroupStatTable.refChannel;
        reversePairGroupStatTable.refShank                = pairGroupStatTable.tarShank;
        reversePairGroupStatTable.tarShank                = pairGroupStatTable.refShank;
        reversePairGroupStatTable.channelSetRef           = pairGroupStatTable.channelSetTar;
        reversePairGroupStatTable.channelSetTar           = pairGroupStatTable.channelSetRef;
        reversePairGroupStatTable.refProbe                = pairGroupStatTable.tarProbe;
        reversePairGroupStatTable.tarProbe                = pairGroupStatTable.refProbe;
        reversePairGroupStatTable.refCellExplorerType     = pairGroupStatTable.tarCellExplorerType;
        reversePairGroupStatTable.tarCellExplorerType     = pairGroupStatTable.refCellExplorerType;
        reversePairGroupStatTable.refWaveforms            = pairGroupStatTable.tarWaveforms;
        reversePairGroupStatTable.tarWaveforms            = pairGroupStatTable.refWaveforms;
        reversePairGroupStatTable.tWave                   = pairGroupStatTable.tWave;
        reversePairGroupStatTable.pairDistance            = pairGroupStatTable.pairDistance;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.pairRawCCG{i,1}     = flipud(pairGroupStatTable.pairRawCCG{i,1});
        end
        reversePairGroupStatTable.refRawACG               = pairGroupStatTable.tarRawACG;
        reversePairGroupStatTable.tarRawACG               = pairGroupStatTable.refRawACG;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.jitterMean{i,1}     = flipud(pairGroupStatTable.jitterMean{i,1});
            reversePairGroupStatTable.CCGbinLagTimes{i,1} = flipud(pairGroupStatTable.CCGbinLagTimes{i,1});
            reversePairGroupStatTable.GSPExc{i,1}         = flipud(pairGroupStatTable.GSPExc{i,1});
            reversePairGroupStatTable.GSPInh{i,1}         = flipud(pairGroupStatTable.GSPInh{i,1});
            reversePairGroupStatTable.pvalE{i,1}          = flipud(pairGroupStatTable.pvalE{i,1});
            reversePairGroupStatTable.pvalI{i,1}          = flipud(pairGroupStatTable.pvalI{i,1});
            reversePairGroupStatTable.LSPExc{i,1}         = flipud(pairGroupStatTable.LSPExc{i,1});
            reversePairGroupStatTable.LSPInh{i,1}         = flipud(pairGroupStatTable.LSPInh{i,1});
        end
        reversePairGroupStatTable.refSnippetsAtRefChan    = pairGroupStatTable.tarSnippetsAtTarChan;
        reversePairGroupStatTable.refSnippetsAtTarChan    = pairGroupStatTable.tarSnippetsAtRefChan;
        reversePairGroupStatTable.tarSnippetsAtTarChan    = pairGroupStatTable.refSnippetsAtRefChan;
        reversePairGroupStatTable.tarSnippetsAtRefChan    = pairGroupStatTable.refSnippetsAtTarChan;
        reversePairGroupStatTable.refWaveAtRefChan        = pairGroupStatTable.tarWaveAtTarChan;
        reversePairGroupStatTable.refWaveAtTarChan        = pairGroupStatTable.tarWaveAtRefChan;
        reversePairGroupStatTable.tarWaveAtTarChan        = pairGroupStatTable.refWaveAtRefChan;
        reversePairGroupStatTable.tarWaveAtRefChan        = pairGroupStatTable.refWaveAtTarChan;
        reversePairGroupStatTable.refWaveAtRefShank       = pairGroupStatTable.tarWaveAtTarShank;
        reversePairGroupStatTable.refWaveAtTarShank       = pairGroupStatTable.tarWaveAtRefShank;
        reversePairGroupStatTable.tarWaveAtTarShank       = pairGroupStatTable.refWaveAtRefShank;
        reversePairGroupStatTable.tarWaveAtRefShank       = pairGroupStatTable.refWaveAtTarShank;

        pairGroupStatTable = [pairGroupStatTable; reversePairGroupStatTable];
        
    elseif strcmp(datasetID,'RatN') || (strcmp(datasetID,'RatS') || strcmp(datasetID,'RatU'))
        
        reversePairGroupStatTable.brainRegion             = pairGroupStatTable.brainRegion;
        reversePairGroupStatTable.refNeuronID             = pairGroupStatTable.tarNeuronID;
        reversePairGroupStatTable.tarNeuronID             = pairGroupStatTable.refNeuronID;
        reversePairGroupStatTable.refSpikeTimes           = pairGroupStatTable.tarSpikeTimes;
        reversePairGroupStatTable.tarSpikeTimes           = pairGroupStatTable.refSpikeTimes;
        reversePairGroupStatTable.refChannel              = pairGroupStatTable.tarChannel;
        reversePairGroupStatTable.tarChannel              = pairGroupStatTable.refChannel;
        reversePairGroupStatTable.refShank                = pairGroupStatTable.tarShank;
        reversePairGroupStatTable.tarShank                = pairGroupStatTable.refShank;
        reversePairGroupStatTable.channelSetRef           = pairGroupStatTable.channelSetTar;
        reversePairGroupStatTable.channelSetTar           = pairGroupStatTable.channelSetRef;
        reversePairGroupStatTable.refProbe                = pairGroupStatTable.tarProbe;
        reversePairGroupStatTable.tarProbe                = pairGroupStatTable.refProbe;
        reversePairGroupStatTable.refCellExplorerType     = pairGroupStatTable.tarCellExplorerType;
        reversePairGroupStatTable.tarCellExplorerType     = pairGroupStatTable.refCellExplorerType;
        reversePairGroupStatTable.refWaveforms            = pairGroupStatTable.tarWaveforms;
        reversePairGroupStatTable.tarWaveforms            = pairGroupStatTable.refWaveforms;
        reversePairGroupStatTable.tWave                   = pairGroupStatTable.tWave;
        reversePairGroupStatTable.pairDistance            = pairGroupStatTable.pairDistance;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.pairRawCCG{i,1}     = flipud(pairGroupStatTable.pairRawCCG{i,1});
        end
        reversePairGroupStatTable.refRawACG               = pairGroupStatTable.tarRawACG;
        reversePairGroupStatTable.tarRawACG               = pairGroupStatTable.refRawACG;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.jitterMean{i,1}     = flipud(pairGroupStatTable.jitterMean{i,1});
            reversePairGroupStatTable.CCGbinLagTimes{i,1} = flipud(pairGroupStatTable.CCGbinLagTimes{i,1});
            reversePairGroupStatTable.GSPExc{i,1}         = flipud(pairGroupStatTable.GSPExc{i,1});
            reversePairGroupStatTable.GSPInh{i,1}         = flipud(pairGroupStatTable.GSPInh{i,1});
            reversePairGroupStatTable.pvalE{i,1}          = flipud(pairGroupStatTable.pvalE{i,1});
            reversePairGroupStatTable.pvalI{i,1}          = flipud(pairGroupStatTable.pvalI{i,1});
            reversePairGroupStatTable.LSPExc{i,1}         = flipud(pairGroupStatTable.LSPExc{i,1});
            reversePairGroupStatTable.LSPInh{i,1}         = flipud(pairGroupStatTable.LSPInh{i,1});
        end
        reversePairGroupStatTable.refSnippetsAtRefChan    = pairGroupStatTable.tarSnippetsAtTarChan;
        reversePairGroupStatTable.refSnippetsAtTarChan    = pairGroupStatTable.tarSnippetsAtRefChan;
        reversePairGroupStatTable.tarSnippetsAtTarChan    = pairGroupStatTable.refSnippetsAtRefChan;
        reversePairGroupStatTable.tarSnippetsAtRefChan    = pairGroupStatTable.refSnippetsAtTarChan;
        reversePairGroupStatTable.refWaveAtRefChan        = pairGroupStatTable.tarWaveAtTarChan;
        reversePairGroupStatTable.refWaveAtTarChan        = pairGroupStatTable.tarWaveAtRefChan;
        reversePairGroupStatTable.tarWaveAtTarChan        = pairGroupStatTable.refWaveAtRefChan;
        reversePairGroupStatTable.tarWaveAtRefChan        = pairGroupStatTable.refWaveAtTarChan;
        reversePairGroupStatTable.refWaveAtRefShank       = pairGroupStatTable.tarWaveAtTarShank;
        reversePairGroupStatTable.refWaveAtTarShank       = pairGroupStatTable.tarWaveAtRefShank;
        reversePairGroupStatTable.tarWaveAtTarShank       = pairGroupStatTable.refWaveAtRefShank;
        reversePairGroupStatTable.tarWaveAtRefShank       = pairGroupStatTable.refWaveAtTarShank;

        pairGroupStatTable = [pairGroupStatTable; reversePairGroupStatTable];
        
    elseif strcmp(datasetID,'Steinmetz')
        
        reversePairGroupStatTable.probeID                 = pairGroupStatTable.probeID;
        reversePairGroupStatTable.brainRegion             = pairGroupStatTable.brainRegion;
        reversePairGroupStatTable.refNeuronID             = pairGroupStatTable.tarNeuronID;
        reversePairGroupStatTable.tarNeuronID             = pairGroupStatTable.refNeuronID;
        reversePairGroupStatTable.refSpikeTimes           = pairGroupStatTable.tarSpikeTimes;
        reversePairGroupStatTable.tarSpikeTimes           = pairGroupStatTable.refSpikeTimes;
        reversePairGroupStatTable.refChannel              = pairGroupStatTable.tarChannel;
        reversePairGroupStatTable.tarChannel              = pairGroupStatTable.refChannel;
        reversePairGroupStatTable.refCellExplorerType     = pairGroupStatTable.tarCellExplorerType;
        reversePairGroupStatTable.tarCellExplorerType     = pairGroupStatTable.refCellExplorerType;
        reversePairGroupStatTable.refWaveforms            = pairGroupStatTable.tarWaveforms;
        reversePairGroupStatTable.tarWaveforms            = pairGroupStatTable.refWaveforms;
        reversePairGroupStatTable.tWave                   = pairGroupStatTable.tWave;
        reversePairGroupStatTable.pairDistance            = pairGroupStatTable.pairDistance;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.pairRawCCG{i,1}     = flipud(pairGroupStatTable.pairRawCCG{i,1});
        end
        reversePairGroupStatTable.refRawACG               = pairGroupStatTable.tarRawACG;
        reversePairGroupStatTable.tarRawACG               = pairGroupStatTable.refRawACG;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.jitterMean{i,1}     = flipud(pairGroupStatTable.jitterMean{i,1});
            reversePairGroupStatTable.CCGbinLagTimes{i,1} = flipud(pairGroupStatTable.CCGbinLagTimes{i,1});
            reversePairGroupStatTable.GSPExc{i,1}         = flipud(pairGroupStatTable.GSPExc{i,1});
            reversePairGroupStatTable.GSPInh{i,1}         = flipud(pairGroupStatTable.GSPInh{i,1});
            reversePairGroupStatTable.pvalE{i,1}          = flipud(pairGroupStatTable.pvalE{i,1});
            reversePairGroupStatTable.pvalI{i,1}          = flipud(pairGroupStatTable.pvalI{i,1});
            reversePairGroupStatTable.LSPExc{i,1}         = flipud(pairGroupStatTable.LSPExc{i,1});
            reversePairGroupStatTable.LSPInh{i,1}         = flipud(pairGroupStatTable.LSPInh{i,1});
        end
        
        pairGroupStatTable = [pairGroupStatTable; reversePairGroupStatTable];
        
    elseif strcmp(datasetID,'AllenInstitute')
        
        reversePairGroupStatTable.animalID                = pairGroupStatTable.animalID;
        reversePairGroupStatTable.probeID                 = pairGroupStatTable.probeID;
        reversePairGroupStatTable.brainRegion             = pairGroupStatTable.brainRegion;
        reversePairGroupStatTable.refNeuronID             = pairGroupStatTable.tarNeuronID;
        reversePairGroupStatTable.tarNeuronID             = pairGroupStatTable.refNeuronID;
        reversePairGroupStatTable.refSpikeTimes           = pairGroupStatTable.tarSpikeTimes;
        reversePairGroupStatTable.tarSpikeTimes           = pairGroupStatTable.refSpikeTimes;
        reversePairGroupStatTable.refChannel              = pairGroupStatTable.tarChannel;
        reversePairGroupStatTable.tarChannel              = pairGroupStatTable.refChannel;
        reversePairGroupStatTable.refCellExplorerType     = pairGroupStatTable.tarCellExplorerType;
        reversePairGroupStatTable.tarCellExplorerType     = pairGroupStatTable.refCellExplorerType;
        reversePairGroupStatTable.refOptoType2X           = pairGroupStatTable.tarOptoType2X;
        reversePairGroupStatTable.tarOptoType2X           = pairGroupStatTable.refOptoType2X;
        reversePairGroupStatTable.refOptoType5X           = pairGroupStatTable.tarOptoType5X;
        reversePairGroupStatTable.tarOptoType5X           = pairGroupStatTable.refOptoType5X;
        reversePairGroupStatTable.refWaveforms            = pairGroupStatTable.tarWaveforms;
        reversePairGroupStatTable.tarWaveforms            = pairGroupStatTable.refWaveforms;
        reversePairGroupStatTable.tWave                   = pairGroupStatTable.tWave;
        reversePairGroupStatTable.pairDistance            = pairGroupStatTable.pairDistance;
        reversePairGroupStatTable.refSNR                  = pairGroupStatTable.tarSNR;
        reversePairGroupStatTable.tarSNR                  = pairGroupStatTable.refSNR;
        reversePairGroupStatTable.refFiringRate           = pairGroupStatTable.tarFiringRate;
        reversePairGroupStatTable.tarFiringRate           = pairGroupStatTable.refFiringRate;
        reversePairGroupStatTable.refTroughToPeakLength   = pairGroupStatTable.tarTroughToPeakLength;
        reversePairGroupStatTable.tarTroughToPeakLength   = pairGroupStatTable.refTroughToPeakLength;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.pairRawCCG{i,1}     = flipud(pairGroupStatTable.pairRawCCG{i,1});
        end
        reversePairGroupStatTable.refRawACG               = pairGroupStatTable.tarRawACG;
        reversePairGroupStatTable.tarRawACG               = pairGroupStatTable.refRawACG;
        for i = 1:size(pairGroupStatTable,1)
            reversePairGroupStatTable.jitterMean{i,1}     = flipud(pairGroupStatTable.jitterMean{i,1});
            reversePairGroupStatTable.CCGbinLagTimes{i,1} = flipud(pairGroupStatTable.CCGbinLagTimes{i,1});
            reversePairGroupStatTable.GSPExc{i,1}         = flipud(pairGroupStatTable.GSPExc{i,1});
            reversePairGroupStatTable.GSPInh{i,1}         = flipud(pairGroupStatTable.GSPInh{i,1});
            reversePairGroupStatTable.pvalE{i,1}          = flipud(pairGroupStatTable.pvalE{i,1});
            reversePairGroupStatTable.pvalI{i,1}          = flipud(pairGroupStatTable.pvalI{i,1});
            reversePairGroupStatTable.LSPExc{i,1}         = flipud(pairGroupStatTable.LSPExc{i,1});
            reversePairGroupStatTable.LSPInh{i,1}         = flipud(pairGroupStatTable.LSPInh{i,1});
        end
        
        pairGroupStatTable = [pairGroupStatTable; reversePairGroupStatTable];
        
    end

end