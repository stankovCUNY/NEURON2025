function UtkuScreenConvV2(dayNo,nullFlag)

    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    if strcmp(dayNo,'1')
        data_dir = [UtkuPath 'AG_2019-12-23_NSD' '/'];
    elseif strcmp(dayNo,'2')
        data_dir = [UtkuPath 'AG_2019-12-27_NSD' '/'];
    elseif strcmp(dayNo,'B1')
        data_dir = [UtkuPath 'B1_2022-06-24_NSD_CA1_24hrs/2022-06-24_NSD_CA1_24hrs/2022-06-24_NSD_CA1_24hrscrs-merged_cleaned.GUI' '/'];
    end    

    spikeTimes = []; % spiket for Eran
    spikeInd   = []; % spikeind for Eran
    
    SampleRate = 30000;                             % SampleRate for Eran
    alpha      = 0.05;                              % alpah for Eran
    wintype    = 'gauss';                           % wintype for Eran

    if strcmp(dayNo,'1') || strcmp(dayNo,'2') % healthy rat AG for days 1 and 2
        
        UtkuData = loadUtkuV2(dayNo);
        
        for i = 1:length(UtkuData.cell_metrics.UID) 
            spikeTimes = [UtkuData.cell_metrics.spikes.times{1,i}; spikeTimes];
            spikeInd   = [ones(length(UtkuData.cell_metrics.spikes.times{1,i}),1)*UtkuData.cell_metrics.UID(i) ; spikeInd];
        end

        cellID     = UtkuData.cell_metrics.UID;         % Cells for Eran
        shankID    = UtkuData.cell_metrics.shankID;     % shank for Eran
        
    elseif strcmp(dayNo,'B1') % aged rat B1
        UtkuScreenConvV2('2',true)
        spikeTimes = double(readNPY([data_dir 'spike_times.npy']))/SampleRate;
        spikeInd   = double(readNPY([data_dir 'spike_clusters.npy']));
        
        cellID     = unique(spikeInd)';
        
        clusterinfo   = ratB1clusterInfo;
        channelShanks = double(readNPY([data_dir 'channel_shanks.npy'])) + 1;
        
        shankID = [];
        
        for i = 1:length(cellID)
                        
            tempIdx = find(cellID(i) == clusterinfo.cluster_id);
            
            unitQuality{i} = string(clusterinfo.group(tempIdx));
            channel        = clusterinfo.ch(tempIdx) + 1;
            shankID(i)     = channelShanks(channel);
            
        end
        
        % track "mua" and "noise" cells
        badUnitQualityBin = zeros(length(unitQuality),1);
        for i = 1:length(unitQuality) 
        
            if strcmp(unitQuality{i},"mua") || strcmp(unitQuality{i},"noise")
                badUnitQualityBin(i) = 1; 
            end
            
        end
        
        % remove "mua" and "noise" cells
        cellID  = cellID(~badUnitQualityBin);
        shankID = shankID(~badUnitQualityBin);
        
        % track "mua" and "noise" spike times
        goodCellSpikeIdx = zeros(length(spikeTimes),1);
        for i = 1:length(cellID)
            goodCellSpikeIdx = goodCellSpikeIdx | (cellID(i) == spikeInd);
        end
        
        % remove "mua" and "noise" spike times
        spikeTimes = spikeTimes(goodCellSpikeIdx);
        spikeInd   = spikeInd(goodCellSpikeIdx);
        
    end
    
    if nullFlag
        one_ms = 0.001;
        delta = 5*one_ms;
        spikeTimes = floor(spikeTimes/delta)*delta + rand(size(spikeTimes))*delta;
    end

    % spikeind for Eran
    % cluID = zeros(size(UtkuData.spindices,1),1);
    % for i = 1:size(UtkuData.spindices,1)
    %     
    %     tempCellID = UtkuData.spindices(i,2);
    %     cluID(i)   = UtkuData.cluID(tempCellID);
    %     
    % end

    %% jscale 1
    jscale = 1; % jscale for Eran

    [ExcPairs, ...
     InhPairs, ...
     GapPairs, ...
     Rzero] = EranConv_group_forUtku(spikeTimes, ...
                                     spikeInd, ...
                                     cellID, ...
                                     SampleRate, ...
                                     jscale, ...
                                     alpha,...
                                     shankID, ...
                                     wintype);

    pairs.ExcPairs   = ExcPairs;
    pairs.InhPairs   = InhPairs;
    pairs.GapPairs   = GapPairs;
    pairs.RZero      = Rzero;
    pairs.jscale     = jscale;
    pairs.nullFlag   = nullFlag;
    pairs.spikeTimes = spikeTimes;
    pairs.spikeInd   = spikeInd;
    
    if nullFlag
        save(fullfile(data_dir,['UtkuDay' dayNo 'null_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
    else
        save(fullfile(data_dir,['UtkuDay' dayNo '_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
    end

    %% jscale 5
    jscale = 5;  % jscale for Eran

    [ExcPairs, ...
     InhPairs, ...
     GapPairs, ...
     Rzero] = EranConv_group_forUtku(spikeTimes, ...
                                     spikeInd, ...
                                     cellID, ...
                                     SampleRate, ...
                                     jscale, ...
                                     alpha,...
                                     shankID, ...
                                     wintype);

    pairs.ExcPairs   = ExcPairs;
    pairs.InhPairs   = InhPairs;
    pairs.GapPairs   = GapPairs;
    pairs.RZero      = Rzero;
    pairs.jscale     = jscale;
    pairs.nullFlag   = nullFlag;
    pairs.spikeTimes = spikeTimes;
    pairs.spikeInd   = spikeInd;
    
    if nullFlag
        save(fullfile(data_dir,['UtkuDay' dayNo 'null_jscale' num2str(jscale) '_alpha' ...
                num2str(round(alpha*100)) '_pairs']), ...
                'pairs', 'jscale', 'alpha')    
    else
        save(fullfile(data_dir,['UtkuDay' dayNo '_jscale' num2str(jscale) '_alpha' ...
                num2str(round(alpha*100)) '_pairs']), ...
                'pairs', 'jscale', 'alpha')    
    end

    %%

%     dataJscale1 = load(['UtkuDay' dayNo '_jscale1_alpha5_pairs.mat']);
%     dataJscale5 = load(['UtkuDay' dayNo '_jscale5_alpha5_pairs.mat']);
% 
%     ExcPairs = [dataJscale1.pairs.ExcPairs(:,1:2); dataJscale5.pairs.ExcPairs(:,1:2)];
%     InhPairs = [dataJscale1.pairs.InhPairs(:,1:2); dataJscale5.pairs.InhPairs(:,1:2)];
%     GapPairs = [dataJscale1.pairs.GapPairs(:,1:2); dataJscale5.pairs.GapPairs(:,1:2)];
% 
%     [~,idxUnique,~] = unique(ExcPairs,'rows');
%     ExcPairs(idxUnique,:) = [];
% 
%     [~,idxUnique,~] = unique(InhPairs,'rows');
%     InhPairs(idxUnique,:) = [];
% 
%     [~,idxUnique,~] = unique(GapPairs,'rows');
%     GapPairs(idxUnique,:) = [];
% 
%     allPairs = [ExcPairs; InhPairs; GapPairs];
%     [~,idxUnique,~] = unique(allPairs,'rows');
%     allPairs = allPairs(idxUnique,:);
end

