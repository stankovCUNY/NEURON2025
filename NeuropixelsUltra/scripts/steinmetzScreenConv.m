function steinmetzScreenConv

    dataPath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';
    
    brainRegion = steinmetzBrainArea;
    brainRegionList = unique(brainRegion);

    for region = 1:length(brainRegionList)
        
        display([num2str((region/length(brainRegionList))*100) '% done'])
        
        if strcmp(brainRegionList{region},'MOs2/3')
            brainRegionSaveStr = 'MOs2and3';
        elseif strcmp(brainRegionList{region},'ORBl2/3')
            brainRegionSaveStr = 'ORBl2and3';
        elseif strcmp(brainRegionList{region},'VISam2/3')
            brainRegionSaveStr = 'VISam2and3';
        elseif strcmp(brainRegionList{region},'VISpm2/3')
            brainRegionSaveStr = 'VISpm2and3';
        else
            brainRegionSaveStr = brainRegionList{region};
        end
        
        load([dataPath 'neuronsSteinmetzDataArea' brainRegionSaveStr '.mat'])
    
        spikeTimes = []; % spiket for Eran
        spikeInd   = []; % spikeind for Eran

        for i = 1:length(neurons.unitID) 
            spikeTimes = [neurons.spikeTimes{i}'; spikeTimes];
            spikeInd   = [ones(length(neurons.spikeTimes{1,i}),1)*(double(neurons.unitID(i))); spikeInd];
        end

        cellID     = double(neurons.unitID);             % Cells for Eran
        probeID    = double(cell2num(neurons.probeID));	 % shank for Eran, set as 1 bc of neuropixel
        SampleRate = 30000;                              % SampleRate for Eran
        alpha      = 0.05;                               % alpah for Eran
        wintype    = 'gauss';                            % wintype for Eran

        %% jscale 1
        jscale = 1; % jscale for Eran

        [ExcPairs, ...
         InhPairs, ...
         GapPairs, ...
         Rzero] = EranConv_group_forUtku(spikeTimes', ...
                                         spikeInd', ...
                                         cellID, ...
                                         SampleRate, ...
                                         jscale, ...
                                         alpha,...
                                         probeID, ...
                                         wintype);

        pairs.ExcPairs = ExcPairs;
        pairs.InhPairs = InhPairs;
        pairs.GapPairs = GapPairs;
        pairs.RZero    = Rzero;
        pairs.jscale   = jscale;

        save(fullfile(dataPath,['SteinmetzMouseRegion' brainRegionSaveStr '_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')

        %% jscale 5
        jscale = 5;  % jscale for Eran

        [ExcPairs, ...
         InhPairs, ...
         GapPairs, ...
         Rzero] = EranConv_group_forUtku(spikeTimes', ...
                                         spikeInd', ...
                                         cellID, ...
                                         SampleRate, ...
                                         jscale, ...
                                         alpha,...
                                         probeID, ...
                                         wintype);

        pairs.ExcPairs = ExcPairs;
        pairs.InhPairs = InhPairs;
        pairs.GapPairs = GapPairs;
        pairs.RZero    = Rzero;
        pairs.jscale   = jscale;

        save(fullfile(dataPath,['SteinmetzMouseRegion' brainRegionSaveStr '_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
            
    end
 end