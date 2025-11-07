function UtkuRunCCGputativeGJ(dayNo,nullFlag)

    flagAdjusted = false;

    UtkuData      = loadUtkuV2(dayNo); % NOTE: spike time data is in minutes!!!
    pairType      = 'exquisite';
    
    [filteredPairs,nullSpikeTimes] = UtkuScreenJitterV2(dayNo,pairType,nullFlag);
    
    pairGroupStatTable = table;
    
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
    filterFlag     = false;
        
    channelSet    = {[1:32],[33:64],[65:96],[97,128],[129:160],[161:192]};
    channelShanks = [ones(32,1); 2*ones(32,1); 3*ones(32,1); 4*ones(32,1); 5*ones(32,1); 6*ones(32,1);];
    
    preLength  = 36;
    postLength = 36;
    
    % organizing and plotting variabls
    waveTimeShort   = UtkuData.cell_metrics.waveforms.time{1,1};     % 48 time points
    waveTimeLong    = UtkuData.cell_metrics.waveforms.time_all{1,1}; % 72 time points
    electrodeGroups = UtkuData.cell_metrics.general.electrodeGroups; 

    figpath     = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    dataPath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    day = str2double(dayNo);
    if day == 1
        chanDataDir = ['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/' 'AG_2019-12_23_NSD' '/'];
        chanCoords  = readtable('AG_2019-12-23_NSD.Probe.xlsx'); 
        snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        
        % unit global snapshots
        load([dataPath 'AG_2019-12-23_NSD' '/' 'unitsGlobalSnapshots.mat'])
        
        if nullFlag
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/UtkuDay1null_jscale1_alpha5_pairs.mat');
        else
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/AG1_jscale1_alpha5_pairs.mat');
        end
        
    elseif day == 2
        chanDataDir = ['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/' 'AG_2019-12_27_NSD' '/'];
        chanCoords  = readtable('AG_2019-12-27_NSD.Probe.xlsx'); 
        snippetsDir = [chanDataDir 'snippetsPerChanPerUnit/'];
        
        % unit global snapshots
        load([dataPath 'AG_2019-12-27_NSD' '/' 'unitsGlobalSnapshots.mat'])
        
        if nullFlag
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/UtkuDay2null_jscale1_alpha5_pairs.mat');
        else
            dataJscale1 = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/AG2_jscale1_alpha5_pairs.mat');
        end
        
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
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    for loopPairs = 1:size(filteredPairs,1)
        
        display([num2str(100*(loopPairs/size(filteredPairs,1))) '% done'])
        
        if nullFlag
            res1    = nullSpikeTimes.spikeTimes((nullSpikeTimes.spikeInd) == filteredPairs(loopPairs,1));
            res2    = nullSpikeTimes.spikeTimes((nullSpikeTimes.spikeInd) == filteredPairs(loopPairs,2));
        else
            res1    = UtkuData.cell_metrics.spikes.times{1,filteredPairs(loopPairs,1)};
            res2    = UtkuData.cell_metrics.spikes.times{1,filteredPairs(loopPairs,2)};
        end
        
        cellID1 = UtkuData.cell_metrics.UID(filteredPairs(loopPairs,1));
        cellID2 = UtkuData.cell_metrics.UID(filteredPairs(loopPairs,2));

        ch1     = UtkuData.cell_metrics.maxWaveformCh1(filteredPairs(loopPairs,1));
        ch2     = UtkuData.cell_metrics.maxWaveformCh1(filteredPairs(loopPairs,2));

        clu1    = UtkuData.cell_metrics.cluID(filteredPairs(loopPairs,1));
        clu2    = UtkuData.cell_metrics.cluID(filteredPairs(loopPairs,2));
        
        shank1  = UtkuData.cell_metrics.shankID(filteredPairs(loopPairs,1));
        shank2  = UtkuData.cell_metrics.shankID(filteredPairs(loopPairs,2));
        
        if shank1 <= 4; probe1 = 1; elseif shank1 > 4; probe1 = 2; end
        if shank2 <= 4; probe2 = 1; elseif shank2 > 4; probe2 = 2; end
        
        cell1type = UtkuData.cell_metrics.putativeCellType{1,filteredPairs(loopPairs,1)};
        cell2type = UtkuData.cell_metrics.putativeCellType{1,filteredPairs(loopPairs,2)};
        
        region1 = UtkuData.cell_metrics.brainRegion{1,filteredPairs(loopPairs,1)}; 
        region2 = UtkuData.cell_metrics.brainRegion{1,filteredPairs(loopPairs,2)};
        
        z1 = double(chanCoords.Z(ch1 == (chanCoords.ChannelNumberComingOutFromProbe + 1)));
        z2 = double(chanCoords.Z(ch2 == (chanCoords.ChannelNumberComingOutFromProbe + 1)));

        x1 = double(chanCoords.X(ch1 == (chanCoords.ChannelNumberComingOutFromProbe + 1)));
        x2 = double(chanCoords.X(ch2 == (chanCoords.ChannelNumberComingOutFromProbe + 1)));

        d = sqrt((x2-x1)^2 + (z2-z1)^2);
        
        if d <= 55
            continue
        end
        
        cell1waveform = UtkuData.cell_metrics.waveforms.filt_all{1,filteredPairs(loopPairs,1)};
        cell2waveform = UtkuData.cell_metrics.waveforms.filt_all{1,filteredPairs(loopPairs,2)};
        
        waveform1Idx = find(electrodeGroups{1,shank1} == ch1);
        waveform2Idx = find(electrodeGroups{1,shank2} == ch2);
        
%         spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;
        
        if strcmp(cell1type,'Pyramidal Cell')
            cell1type = 'p';
        elseif strcmp(cell1type,'Narrow Interneuron')
            cell1type = 'i-narrow';
        elseif strcmp(cell1type,'Wide Interneuron')
            cell1type = 'i-wide';
        end
        
        if strcmp(cell2type,'Pyramidal Cell')
            cell2type = 'p';
        elseif strcmp(cell2type,'Narrow Interneuron')
            cell2type = 'i-narrow';
        elseif strcmp(cell2type,'Wide Interneuron')
            cell2type = 'i-wide';
        end
        
        % retrieve waveform snippets
        if flagAdjusted
            snippetsAtRefChan = load([snippetsDir 'shank' num2str(channelShanks(ch1)) '/ch' num2str(ch1) 'adjustedSnippets.mat']);
            snippetsAtTarChan = load([snippetsDir 'shank' num2str(channelShanks(ch2)) '/ch' num2str(ch2) 'adjustedSnippets.mat']);
        else           
            snippetsAtRefChan = load([snippetsDir 'shank' num2str(channelShanks(ch1)) '/ch' num2str(ch1) 'snippets.mat']);
            snippetsAtTarChan = load([snippetsDir 'shank' num2str(channelShanks(ch2)) '/ch' num2str(ch2) 'snippets.mat']);
        end
        
%         tiledlayout(5,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
        
        pairStr = ['AG' dayNo ' (GJ screen) - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', prb:' num2str(probe1) ', clu: ' num2str(clu1) ', reg: ' region1 ')' ' v ' ...
                                    num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', prb:' num2str(probe2) ', clu: ' num2str(clu2) ', reg: ' region2 ')'...
                                    ' dist: ' num2str(d) ' microns'];
%         nexttile([1 ,2])
        nexttile
        
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
        
        % removing the synchronous spikes                
        [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(res1,res2,binSize,lagLimit);
                        
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
            
            snippetsTemp = load([snippetsDir 'shank' num2str(channelShanks(channelSetRef(loopChans))) '/ch' num2str(channelSetRef(loopChans)) 'snippets.mat']);
            
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
            snippetsTemp = load([snippetsDir 'shank' num2str(channelShanks(channelSetTar(loopChans))) '/ch' num2str(channelSetTar(loopChans)) 'snippets.mat']);
            
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

%             pairStr = [num2str(cellID1) cell1type ' v ' ...
%                        num2str(cellID2) cell2type];

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
            xlim([-(duration/2)*1000 (duration/2)*1000])
            
%             nexttile
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
            xlim([-(duration/2)*1000 (duration/2)*1000])
            
%             nexttile
%             plot(waveTimeLong,refWaveAtRefChan,'LineWidth',2,'Color','#0072BD')
%             hold on
%             plot(waveTimeLong,refWaveAtTarChan,'--','LineWidth',2,'Color','#0072BD')
%             hold off
%             legend('Ref e-spike at ref max chan', ...
%                    'Ref e-spike at tar max chan', ...
%                    'Location','Northwest')
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', prb:' num2str(probe1) ', clu: ' num2str(clu1) ', reg: ' region1 ')' ' average'])
%             
%             nexttile
%             plot(waveTimeLong,tarWaveAtTarChan,'LineWidth',2,'Color','#D95319')
%             hold on
%             plot(waveTimeLong,tarWaveAtRefChan,'--','LineWidth',2,'Color','#D95319')
%             hold off
%             legend('Tar e-spike at tar max chan', ...
%                    'Tar e-spike at ref max chan', ...
%                    'Location','Northwest')
%             xlabel('Time revere order [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'XDir','reverse')
%             set(gca, 'YDir','reverse')
%             title(['Target e-spike: ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', prb:' num2str(probe2) ', clu: ' num2str(clu2) ', reg: ' region2 ')' ' average'])
%             
%             nexttile
%             patchline(waveTimeLong,refSnippetsAtRefChan(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(waveTimeLong,refSnippetsAtRefChan(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
%             end
%             hold off
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' at ref max channel (100 random spikes)'])
%             
%             nexttile
%             patchline(waveTimeLong,tarSnippetsAtTarChan(:,1),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(waveTimeLong,tarSnippetsAtTarChan(:,k),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
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
%             patchline(waveTimeLong,refSnippetsAtTarChan(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(waveTimeLong,refSnippetsAtTarChan(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             end
%             hold off
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' at tar max channel (100 random spikes)'])
%             
%             nexttile
%             patchline(waveTimeLong,tarSnippetsAtRefChan(:,1),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
%             hold on
%             for k = 2:Nwaveforms
%                 patchline(waveTimeLong,tarSnippetsAtRefChan(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
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
%             plot(waveTimeLong,refWaveAtTarShank)
%             xlabel('Time [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
% %             legend(cellstr(num2str(cell2num(probegroup.channel_id(find((probegroup.shank_id + 1) == shank1)))' + 1, 'ch%-d')),'Location','West')
%             set(gca, 'YDir','reverse')
%             title(['Reference e-spike: ' num2str(cellID1) cell1type ' on all ' num2str(cellID2) cell2type ' shank ' num2str(shank2) ' channels'])
% 
%             nexttile
%             plot(waveTimeLong,tarWaveAtRefShank)
%             xlabel('Time revere order [ms]')
%             ylabel('[mV]')
%             xlim([-3.5 3.5])
%             set(gca, 'XDir','reverse')
%             set(gca, 'YDir','reverse')
% %             legend(cellstr(num2str(cell2num(probegroup.channel_id(find((probegroup.shank_id + 1) == shank2)))' + 1, 'ch%-d')),'Location','West')
%             title(['Target e-spike: ' num2str(cellID2) cell2type ' on all ' num2str(cellID1) cell1type ' shank ' num2str(shank1) ' channels'])
            
            save_file = fullfile(figpath, ['waveform ana - ' pairStr]);
            print(fig_use, save_file,'-djpeg',resolution_use);

            close all

            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
            
        end
        
        %%
                
        pairGroupStatsTempTable = table;

        pairGroupStatsTempTable.refBrainRegion = region1;
        pairGroupStatsTempTable.tarBrainRegion = region2;
        
        pairGroupStatsTempTable.refClusterID   = clu1;
        pairGroupStatsTempTable.tarClusterID   = clu2;

        pairGroupStatsTempTable.refNeuronID    = cellID1;
        pairGroupStatsTempTable.tarNeuronID    = cellID2;

        pairGroupStatsTempTable.refSpikeTimes  = {res1};
        pairGroupStatsTempTable.tarSpikeTimes  = {res2};

        pairGroupStatsTempTable.refChannel     = ch1;
        pairGroupStatsTempTable.tarChannel     = ch2;

        pairGroupStatsTempTable.refShank       = shank1;
        pairGroupStatsTempTable.tarShank       = shank2;

        pairGroupStatsTempTable.channelSetRef  = {channelSetRef};
        pairGroupStatsTempTable.channelSetTar  = {channelSetTar};

        pairGroupStatsTempTable.refProbe       = probe1;
        pairGroupStatsTempTable.tarProbe       = probe2;

        pairGroupStatsTempTable.refCellExplorerType = {cell1type};
        pairGroupStatsTempTable.tarCellExplorerType = {cell2type};

        pairGroupStatsTempTable.refWaveforms = {unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,1),1}};
        pairGroupStatsTempTable.tarWaveforms = {unitsGlobalSnapshots.waveforms{filteredPairs(loopPairs,2),1}};
        pairGroupStatsTempTable.tWave        = {waveTimeLong};

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
        if day == 1
            save([dataPath 'AG_2019-12-23_NSD' '/' 'nullGroupStatsRatAGday' dayNo 'putativeGJ.mat'],'pairGroupStatTable')
        elseif day == 2
            save([dataPath 'AG_2019-12-27_NSD' '/' 'nullGroupStatsRatAGday' dayNo 'putativeGJ.mat'],'pairGroupStatTable')
        end 
    else
        if day == 1
            save([dataPath 'AG_2019-12-23_NSD' '/' 'groupStatsRatAGday' dayNo 'putativeGJ.mat'],'pairGroupStatTable')
        elseif day == 2
            save([dataPath 'AG_2019-12-27_NSD' '/' 'groupStatsRatAGday' dayNo 'putativeGJ.mat'],'pairGroupStatTable')
        end 
    end
    
    
    
end