function BapunRunCCGputativeGJ(RatID,nullFlag)
    
    flagAdjusted = false;

    pairType      = 'all';
    [filteredPairs,nullSpikeTimes] = BapunScreenJitter(RatID,pairType,nullFlag);
    
    pairGroupStatTable = table;
    
    if strcmp(RatID,'N')
        
        neurons     = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        probegroup  = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
        snippetsDir = '/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/snippetsPerChanPerUnit/';
        
        % unit global snapshots
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/unitsGlobalSnapshots.mat')
        
        % load adjusted spike times
        if flagAdjusted
            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatN_adjustedSpikeTimes.mat');
        end
        
        if nullFlag
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/Bapuns_RatNnull_jscale1_alpha5_pairs.mat');
        else
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/Bapuns_RatN_jscale1_alpha5_pairs.mat');
        end
        
        chanID_xml = [13  2   14  1   15  0   12  3   11  4   10  5   9   6   8   7   ...
                      29  18  30  17  31  16  28  19  27  20  26  21  25  22  24  23  ...
                      45  34  46  33  47  32  44  35  43  36  42  37  41  38  40  39  ...
                      61  50  62  49  63  48  60  51  59  52  58  53  57  54  56  55  ...
                      77  66  78  65  79  64  76  67  75  68  74  69  73  70  72  71  ...
                      93  82  94  81  95  80  92  83  91  84  90  85  89  86  88  87  ...
                      109 98  110 97  111 96  108 99  107 100 106 101 105 102 104 103 ...
                      125 114 126 113 127 112 124 115 123 116 122 117 121 118 120 119
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
        
        neurons     = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        probegroup  = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');
        snippetsDir = '/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/snippetsPerChanPerUnit/';
        
        % unit global snapshots
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/unitsGlobalSnapshots.mat')
        
        % load adjusted spike times
        if flagAdjusted
            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatS_adjustedSpikeTimes.mat');
        end
        
        if nullFlag
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/Bapuns_RatSnull_jscale1_alpha5_pairs.mat');
        else
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/Bapuns_RatS_jscale1_alpha5_pairs.mat');
        end
        
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
        
        neurons     = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
        probegroup  = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');
        snippetsDir = '/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/snippetsPerChanPerUnit/';
        
        % unit global snapshots
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/unitsGlobalSnapshots.mat')
        
        % load adjusted spike times
        if flagAdjusted
            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatU_adjustedSpikeTimes.mat');
        end
        
        if nullFlag
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/Bapuns_RatUnull_jscale1_alpha5_pairs.mat');
        else
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/Bapuns_RatU_jscale1_alpha5_pairs.mat');
        end
        
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
          
    if ~isempty(dataJscale1.pairs.ExcPairs) 
        ExcPairs = dataJscale1.pairs.ExcPairs(:,1:2);
    else
        ExcPairs = [];
    end 

    if ~isempty(dataJscale1.pairs.InhPairs) 
        InhPairs = dataJscale1.pairs.InhPairs(:,1:2);
    else
        InhPairs = [];
    end

    if ~isempty(dataJscale1.pairs.GapPairs) 
        GapPairs = dataJscale1.pairs.GapPairs(:,1:2);
    else
        GapPairs = [];
    end

    filteredPairs = [ExcPairs; InhPairs; GapPairs];

    [~,idxUnique,~] = unique(filteredPairs,'rows');
    filteredPairs = filteredPairs(idxUnique,:);
    
    matWaveIdx = find(cell2mat(probegroup.connected));
    
    %%

    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    lagLimit       = duration/2;
    fs             = 30000;
    fpass          = 300;
    binSize        = 1/1000;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    Nwaveforms     = 100;
    filterFlag     = true;
    
    tWave          = (-35:36)*(1/30);
    
    figPath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    dataPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    for loopPairs = 1:size(filteredPairs,1)
        
        display([num2str(100*(loopPairs/size(filteredPairs,1))) '% done'])
        
        cellID1 = filteredPairs(loopPairs,1);
        cellID2 = filteredPairs(loopPairs,2);
        
        cell1type = neurons.neuron_type(filteredPairs(loopPairs,1),:);
        cell2type = neurons.neuron_type(filteredPairs(loopPairs,2),:);
        
        % reject MUAs at this stage
        if strcmp(cell1type,'mua  ') || strcmp(cell2type,'mua  ')
            continue
        end
        
        if flagAdjusted            
            res1 = adjustedSpikeTimes.spikeTimes{find(adjustedSpikeTimes.cellID == filteredPairs(loopPairs,1)),1};
            res2 = adjustedSpikeTimes.spikeTimes{find(adjustedSpikeTimes.cellID == filteredPairs(loopPairs,2)),1};
        elseif nullFlag
            res1    = nullSpikeTimes.spikeTimes((nullSpikeTimes.spikeInd) == filteredPairs(loopPairs,1));
            res2    = nullSpikeTimes.spikeTimes((nullSpikeTimes.spikeInd) == filteredPairs(loopPairs,2));
        else
            res1    = neurons.spiketrains{1,filteredPairs(loopPairs,1)};
            res2    = neurons.spiketrains{1,filteredPairs(loopPairs,2)};
        end
        
        if strcmp(RatID,'N') % sorted by shank, need to find max channel on my own bc channel markers are unreliable
            shank1  = double(neurons.shank_ids(filteredPairs(loopPairs,1))) + 1;
            shank2  = double(neurons.shank_ids(filteredPairs(loopPairs,2))) + 1;
            
            [~,ch1] = max(max(abs(unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,1),1}(channelSet{shank1},:)),[],2));
            [~,ch2] = max(max(abs(unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,2),1}(channelSet{shank2},:)),[],2));
            
            ch1 = ch1 + 16*(shank1-1);
            ch2 = ch2 + 16*(shank2-1);
        elseif strcmp(RatID,'S') % need to find max channel on my own bc channel markers are unreliable
            [~,ch1] = max(max(abs(unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,1),1}),[],2));
            [~,ch2] = max(max(abs(unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,2),1}),[],2));
            shank1  = probegroup.shank_id(ch1 == (cell2num(probegroup.channel_id) + 1)) + 1;
            shank2  = probegroup.shank_id(ch2 == (cell2num(probegroup.channel_id) + 1)) + 1;
        elseif strcmp(RatID,'U') % channel markers are reliable
            ch1     = double(neurons.peak_channels(filteredPairs(loopPairs,1))) + 1;
            ch2     = double(neurons.peak_channels(filteredPairs(loopPairs,2))) + 1;
            shank1  = double(neurons.shank_ids(filteredPairs(loopPairs,1))) + 1;
            shank2  = double(neurons.shank_ids(filteredPairs(loopPairs,2))) + 1;
        end
        
        probe1  = probegroup.probe_id(ch1 == (cell2num(probegroup.channel_id) + 1)) + 1;
        probe2  = probegroup.probe_id(ch2 == (cell2num(probegroup.channel_id) + 1)) + 1;
        
        y1 = double(probegroup.y(ch1 == (cell2num(probegroup.channel_id) + 1)));
        y2 = double(probegroup.y(ch2 == (cell2num(probegroup.channel_id) + 1)));

        x1 = double(probegroup.x(ch1 == (cell2num(probegroup.channel_id) + 1)));
        x2 = double(probegroup.x(ch2 == (cell2num(probegroup.channel_id) + 1)));

        d = sqrt((x2-x1)^2 + (y2-y1)^2);
        
        if d <= 55
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
        
        % retrieve waveform snippets
        if flagAdjusted
            snippetsAtRefChan = load([snippetsDir 'ch' num2str(ch1) 'adjustedSnippets.mat']);
            snippetsAtTarChan = load([snippetsDir 'ch' num2str(ch2) 'adjustedSnippets.mat']);
        else
            snippetsAtRefChan = load([snippetsDir 'ch' num2str(ch1) 'snippets.mat']);
            snippetsAtTarChan = load([snippetsDir 'ch' num2str(ch2) 'snippets.mat']);
        end
        
        pairStr = ['Rat' RatID ' (GJ screen) - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' ...
                                                   num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')' ...
                                                   ' distance: ' num2str(d) ' microns'];
        

        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
        
        if ~(GSPExc(round(length(tR)/2)) == 1) 

           close all

           hcomb = figure(102);
           arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

           continue 
        end
                        
        % removing the synchronous spikes                
        [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(res1',res2',binSize,lagLimit);
                        
        res1sync = [];
        res2sync = [];
        
        for k = 1:length(tR)
            if (GSPExc(k) == 1) || (GSPInh(k) == 1)
                if isempty(spikeTimesInBin{k})
                    continue
                else
                    res1sync = [res1sync; spikeTimesInBin{k}(:,1)];
                    res2sync = [res2sync; spikeTimesInBin{k}(:,2)];
                end
            end
        end
        
        % not all spike times are in the data times series 
        refSpikeTimesSnippets = snippetsAtRefChan.snippetsData.spikeTimes{snippetsAtRefChan.snippetsData.neuronID == cellID1,1};
        refSpikeTimesSnippets = refSpikeTimesSnippets(1:size(snippetsAtRefChan.snippetsData.waveforms{snippetsAtRefChan.snippetsData.neuronID == cellID1,1},2));
        
        tarSpikeTimesSnippets = snippetsAtRefChan.snippetsData.spikeTimes{snippetsAtRefChan.snippetsData.neuronID == cellID2,1};
        tarSpikeTimesSnippets = tarSpikeTimesSnippets(1:size(snippetsAtRefChan.snippetsData.waveforms{snippetsAtRefChan.snippetsData.neuronID == cellID2,1},2));
        
        % remove sync spikes for waveforms
        [~,filterIdxRef] = setdiff(refSpikeTimesSnippets,res1sync); filterIdxRef = sort(filterIdxRef);
        [~,filterIdxTar] = setdiff(tarSpikeTimesSnippets,res2sync); filterIdxTar = sort(filterIdxTar);
        
        % extract specific neuron and convert from nanovolts to millivolts
        refSnippetsAtRefChan = snippetsAtRefChan.snippetsData.waveforms{snippetsAtRefChan.snippetsData.neuronID == cellID1,1}/1e3; 
        refSnippetsAtTarChan = snippetsAtTarChan.snippetsData.waveforms{snippetsAtTarChan.snippetsData.neuronID == cellID1,1}/1e3; 
        tarSnippetsAtTarChan = snippetsAtTarChan.snippetsData.waveforms{snippetsAtRefChan.snippetsData.neuronID == cellID2,1}/1e3; 
        tarSnippetsAtRefChan = snippetsAtRefChan.snippetsData.waveforms{snippetsAtTarChan.snippetsData.neuronID == cellID2,1}/1e3; 
        
        % baseline snippets
        refSnippetsAtRefChan = refSnippetsAtRefChan - refSnippetsAtRefChan(1,:);
        refSnippetsAtTarChan = refSnippetsAtTarChan - refSnippetsAtTarChan(1,:);
        tarSnippetsAtTarChan = tarSnippetsAtTarChan - tarSnippetsAtTarChan(1,:);
        tarSnippetsAtRefChan = tarSnippetsAtRefChan - tarSnippetsAtRefChan(1,:);
        
        % filter out synchronous waveforms 
        refSnippetsAtRefChan = refSnippetsAtRefChan(:,filterIdxRef);
        refSnippetsAtTarChan = refSnippetsAtTarChan(:,filterIdxRef);
        tarSnippetsAtTarChan = tarSnippetsAtTarChan(:,filterIdxTar);
        tarSnippetsAtRefChan = tarSnippetsAtRefChan(:,filterIdxTar);
        
        % mean specific waveforms snippets
        refWaveAtRefChan     = mean(refSnippetsAtRefChan,2);
        refWaveAtTarChan     = mean(refSnippetsAtTarChan,2);
        tarWaveAtTarChan     = mean(tarSnippetsAtTarChan,2);
        tarWaveAtRefChan     = mean(tarSnippetsAtRefChan,2);
        
        channelSetRef = channelSet{shank1};
        channelSetTar = channelSet{shank2};
        
        refWaveAtRefShank = zeros(length(channelSetRef),72);
        refWaveAtTarShank = zeros(length(channelSetRef),72);
        tarWaveAtTarShank = zeros(length(channelSetRef),72);
        tarWaveAtRefShank = zeros(length(channelSetRef),72);

        % ref shank snippets
        parfor loopChans = 1:length(channelSetRef)
            snippetsTemp = load([snippetsDir 'ch' num2str(channelSetRef(loopChans)) 'snippets.mat']);
            
            snippetsRefTemp = snippetsTemp.snippetsData.waveforms{snippetsTemp.snippetsData.neuronID == cellID1,1}/1e3;
            snippetsRefTemp = snippetsRefTemp(:,filterIdxRef);
            snippetsRefTemp = mean(snippetsRefTemp - snippetsRefTemp(1,:),2);
            refWaveAtRefShank(loopChans,:) = snippetsRefTemp;
            
            snippetsTarTemp = snippetsTemp.snippetsData.waveforms{snippetsTemp.snippetsData.neuronID == cellID2,1}/1e3;
            snippetsTarTemp = snippetsTarTemp(:,filterIdxTar);
            snippetsTarTemp = mean(snippetsTarTemp - snippetsTarTemp(1,:),2);
            tarWaveAtRefShank(loopChans,:) = snippetsTarTemp;
        end
        
        % tar shank snippets
        parfor loopChans = 1:length(channelSetTar)
            snippetsTemp = load([snippetsDir 'ch' num2str(channelSetTar(loopChans)) 'snippets.mat']);
            
            snippetsRefTemp = snippetsTemp.snippetsData.waveforms{snippetsTemp.snippetsData.neuronID == cellID1,1}/1e3;
            snippetsRefTemp = snippetsRefTemp(:,filterIdxRef);
            snippetsRefTemp = mean(snippetsRefTemp - snippetsRefTemp(1,:),2);
            refWaveAtTarShank(loopChans,:) = snippetsRefTemp;
            
            snippetsTarTemp = snippetsTemp.snippetsData.waveforms{snippetsTemp.snippetsData.neuronID == cellID2,1}/1e3;
            snippetsTarTemp = snippetsTarTemp(:,filterIdxTar);
            snippetsTarTemp = mean(snippetsTarTemp - snippetsTarTemp(1,:),2);
            tarWaveAtTarShank(loopChans,:) = snippetsTarTemp;
        end
                        
        if plotFlag
            
            title(pairStr)

            ylims = get(gca,'ylim');

            if any(GSPExc)
                hold on;
                plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^');
            end
            if any(GSPInh)
                hold on;
                plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv');
            end
            xlim([-3.5 3.5])
            
            [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
                
            title(pairStr)

            ylims = get(gca,'ylim');

            if any(GSPExc)
                hold on;
                plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^');
            end
            if any(GSPInh)
                hold on;
                plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv');
            end
            xlim([-3.5 3.5])
            
%             nexttile
%             plot(tWave,refWaveAtRefChan,'LineWidth',2,'Color','#0072BD')
%             hold on
%             plot(tWave,refWaveAtTarChan,'--','LineWidth',2,'Color','#0072BD')
%             hold off
%             legend('Ref e-spike at ref max chan', ...
%                    'Ref e-spike at tar max chan', ...
%                    'Location','Northwest')
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' average'])
%             
%             nexttile
%             plot(tWave,tarWaveAtTarChan,'LineWidth',2,'Color','#D95319')
%             hold on
%             plot(tWave,tarWaveAtRefChan,'--','LineWidth',2,'Color','#D95319')
%             hold off
%             legend('Tar e-spike at tar max chan', ...
%                    'Tar e-spike at ref max chan', ...
%                    'Location','Northwest')
%             xlabel('Time revere order [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'XDir','reverse')
%             set(gca, 'YDir','reverse')
%             title(['Target e-spike: ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')' ' average'])
%             
%             nexttile
%             patchline(tWave,refSnippetsAtRefChan(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(tWave,refSnippetsAtRefChan(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
%             end
%             hold off
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' at ref max channel (100 random spikes)'])
%             
%             nexttile
%             patchline(tWave,tarSnippetsAtTarChan(:,1),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(tWave,tarSnippetsAtTarChan(:,k),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
%             end
%             hold off
%             xlabel('Time revere order [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'XDir','reverse')
%             set(gca, 'YDir','reverse')
%             title(['Target e-spike: ' num2str(cellID2) cell2type ' at tar max channel (100 random spikes)'])
%             
%             nexttile
%             patchline(tWave,refSnippetsAtTarChan(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(tWave,refSnippetsAtTarChan(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             end
%             hold off
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' at tar max channel (100 random spikes)'])
%             
%             nexttile
%             patchline(tWave,tarSnippetsAtRefChan(:,1),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(tWave,tarSnippetsAtRefChan(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             end
%             hold off
%             xlabel('Time revere order [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'XDir','reverse')
%             set(gca, 'YDir','reverse')
%             title(['Target e-spike: ' num2str(cellID2) cell2type ' at ref max channel (100 random spikes)'])
%             
%             nexttile
%             plot(tWave,refWaveAtTarShank)
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
% %             legend(cellstr(num2str(cell2num(probegroup.channel_id(find((probegroup.shank_id + 1) == shank1)))' + 1, 'ch%-d')),'Location','West')
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' on all ' num2str(cellID2) cell2type ' shank ' num2str(shank2) ' channels'])
% 
%             nexttile
%             plot(tWave,tarWaveAtRefShank)
%             xlabel('Time revere order [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'XDir','reverse')
%             set(gca, 'YDir','reverse')
% %             legend(cellstr(num2str(cell2num(probegroup.channel_id(find((probegroup.shank_id + 1) == shank2)))' + 1, 'ch%-d')),'Location','West')
%             title(['Target e-spike: ' num2str(cellID2) cell2type ' on all ' num2str(cellID1) cell1type ' shank ' num2str(shank1) ' channels'])
%             
%             save_file = fullfile(figPath, pairStr);
%             print(fig_use, save_file,'-djpeg',resolution_use);
            
            close all

            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
            
        end
        
        %%
                
        pairGroupStatsTempTable = table;

        pairGroupStatsTempTable.brainRegion   = "CA1";

        pairGroupStatsTempTable.refNeuronID   = cellID1;
        pairGroupStatsTempTable.tarNeuronID   = cellID2;

        pairGroupStatsTempTable.refSpikeTimes = {res1};
        pairGroupStatsTempTable.tarSpikeTimes = {res2};

        pairGroupStatsTempTable.refChannel    = ch1;
        pairGroupStatsTempTable.tarChannel    = ch2;

        pairGroupStatsTempTable.refShank      = shank1;
        pairGroupStatsTempTable.tarShank      = shank2;

        pairGroupStatsTempTable.channelSetRef = {channelSetRef};
        pairGroupStatsTempTable.channelSetTar = {channelSetTar};

        pairGroupStatsTempTable.refProbe      = probe1;
        pairGroupStatsTempTable.tarProbe      = probe2;

        pairGroupStatsTempTable.refCellExplorerType = {cell1type};
        pairGroupStatsTempTable.tarCellExplorerType = {cell2type};

        pairGroupStatsTempTable.refWaveforms = {unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,1),1}};
        pairGroupStatsTempTable.tarWaveforms = {unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,2),1}};
        pairGroupStatsTempTable.tWave        = {tWave};

        pairGroupStatsTempTable.pairDistance = d;

        pairGroupStatsTempTable.pairRawCCG     = {ccgR(:,1,2)};
        pairGroupStatsTempTable.refRawACG      = {ccgR(:,1,1)};
        pairGroupStatsTempTable.tarRawACG      = {ccgR(:,2,2)};
        pairGroupStatsTempTable.jitterMean     = {ccgjm};
        pairGroupStatsTempTable.CCGbinLagTimes = {tR};
        pairGroupStatsTempTable.GSPExc         = {GSPExc};
        pairGroupStatsTempTable.GSPInh         = {GSPInh};
        pairGroupStatsTempTable.pvalE          = {pvalE};
        pairGroupStatsTempTable.pvalI          = {pvalI};
        pairGroupStatsTempTable.LSPExc         = {LSPExc};
        pairGroupStatsTempTable.LSPInh         = {LSPInh};

        pairGroupStatsTempTable.refSnippetsAtRefChan = {refSnippetsAtRefChan};
        pairGroupStatsTempTable.refSnippetsAtTarChan = {refSnippetsAtTarChan};
        pairGroupStatsTempTable.tarSnippetsAtTarChan = {tarSnippetsAtTarChan};
        pairGroupStatsTempTable.tarSnippetsAtRefChan = {tarSnippetsAtRefChan};

        pairGroupStatsTempTable.refWaveAtRefChan     = {refWaveAtRefChan};
        pairGroupStatsTempTable.refWaveAtTarChan     = {refWaveAtTarChan};
        pairGroupStatsTempTable.tarWaveAtTarChan     = {tarWaveAtTarChan};
        pairGroupStatsTempTable.tarWaveAtRefChan     = {tarWaveAtRefChan};

        pairGroupStatsTempTable.refWaveAtRefShank    = {refWaveAtRefShank};
        pairGroupStatsTempTable.refWaveAtTarShank    = {refWaveAtTarShank};
        pairGroupStatsTempTable.tarWaveAtTarShank    = {tarWaveAtTarShank};
        pairGroupStatsTempTable.tarWaveAtRefShank    = {tarWaveAtRefShank};

        pairGroupStatTable = [pairGroupStatTable; pairGroupStatsTempTable];
        
    end
    
    if nullFlag
        save([dataPath 'Rat' RatID '/' 'nullGroupStatsRat'       RatID 'putativeGJ.mat'],'pairGroupStatTable')
    else
        save([dataPath 'Rat' RatID '/' 'groupStatsRatPutativeGJ' RatID 'putativeGJ.mat'],'pairGroupStatTable')
    end
    
    close all
    
end