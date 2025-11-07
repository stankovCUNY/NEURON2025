function HiroRunCCGputativeGJ(nullFlag)

    pairGroupStatTable = table;
    
    %% hardcoding pairs of interest for waveform analysis

    % constants
    session_name   = 'RoyMaze1';
    % session_name   = 'KevinMaze1';
    conn_type      = {'GapPairs','ExcPairs','InhPairs'};
    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    lagLimit       = duration/2;
    fpass          = 300;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    filterFlag     = false;
    Nwaveforms     = 100;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using
    
    tWave          = (-35:36)*(1/30);
    
    figPath  = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/figures/';
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';

    preLength  = 36;
    postLength = 36;

    if strcmp(session_name,'RoyMaze1')

        fs        = 30000;
        binSize   = 1/1000;
        spikeTime = (1/30)*(-preLength+1:postLength);
                
        snippetsDir = '/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/snippetsPerChanPerUnit/';

    elseif strcmp(session_name,'KevinMaze1')

        fs        = 32000;
        binSize  = 1/fs;
        spikeTime = (1/32)*(-preLength+1:postLength);

    end

    shankChanList ={1:8  , ...
                    9:16 , ...
                    17:24, ...
                    25:32, ...
                    33:40, ...
                    41:48, ...
                    49:56, ...
                    57:64};

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    pos = [70 230 1860 2660]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';

    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    % pair data 
    datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";
    
    % jitter data adjust for null (double jitter)
    if nullFlag
        dataJscale1 = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoyMaze1_jscale1_alpha5_pairs_NULL.mat');
        
        session_name = [session_name 'NULL'];
        [jit_var_out,~] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);
        
        if ~isempty(dataJscale1.pairs(1).ExcPairs) 
            ExcPairs = dataJscale1.pairs(1).ExcPairs(:,1:2);
        else
            ExcPairs = [];
        end 

        if ~isempty(dataJscale1.pairs(1).InhPairs) 
            InhPairs = dataJscale1.pairs(1).InhPairs(:,1:2);
        else
            InhPairs = [];
        end

        if ~isempty(dataJscale1.pairs(1).GapPairs) 
            GapPairs = dataJscale1.pairs(1).GapPairs(:,1:2);
        else
            GapPairs = [];
        end
        
    else
        ExcPairs = [];
        InhPairs = []; 
        GapPairs = [];
        
        ExcPairsTemp = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoyMaze1_ExcPairs_jscale1_alpha5_jitter_stats.mat');
        InhPairsTemp = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoyMaze1_InhPairs_jscale1_alpha5_jitter_stats.mat');
        GapPairsTemp = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoyMaze1_GapPairs_jscale1_alpha5_jitter_stats.mat');
        
        jit_var_out = [ExcPairsTemp.jit_var_out InhPairsTemp.jit_var_out GapPairsTemp.jit_var_out];
        
        for loopRows = 1:size(ExcPairsTemp.jit_var_out,2)
            ExcPairs = [ExcPairs; ExcPairsTemp.jit_var_out(1,loopRows).cell_pair];
        end 
        for loopRows = 1:size(InhPairsTemp.jit_var_out,2)
            InhPairs = [InhPairs; InhPairsTemp.jit_var_out(1,loopRows).cell_pair];
        end 
        for loopRows = 1:size(GapPairsTemp.jit_var_out,2)
            GapPairs = [GapPairs; GapPairsTemp.jit_var_out(1,loopRows).cell_pair];
        end 
    end

    pairsAna = [ExcPairs; InhPairs; GapPairs];

    [~,idxUnique,~] = unique(pairsAna,'rows');
    pairsAna = pairsAna(idxUnique,:);
    
    jit_var_out = jit_var_out(idxUnique);

    %% get session maze time markers
    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    [data_dir, name, ~] = fileparts(spike_data_fullpath);
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-spikes.mat'));

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')
    
    chanCoords = buzsaki64probeLoc;
    
    %%
    
%     if nullFlag
%         pairsAna = cell_pairs_mat;
%     end
    
    cell_pairs_mat = pairsAna;

    for loopPairs = 1:size(pairsAna,1)
        
        display([num2str(100*(loopPairs/size(pairsAna,1))) '% done'])

        idx = loopPairs; % pair to run, right now we are working with pre-maze time period

        current_pair = pairsAna(idx,:);

        current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));

        shank1 = jit_var_out(current_pair_indices).cell1shank;
        shank2 = jit_var_out(current_pair_indices).cell2shank;

        cellID1 = current_pair(1);
        cellID2 = current_pair(2);
        
        refIdx = find(unitsGlobalSnapshots.neuronID == cellID1);
        tarIdx = find(unitsGlobalSnapshots.neuronID == cellID2);

        refWave = unitsGlobalSnapshots.waveforms{refIdx};
        tarWave = unitsGlobalSnapshots.waveforms{tarIdx};

        [~,ch1] = min(min(refWave'));
        [~,ch2] = min(min(tarWave'));

        cell1type = jit_var_out(current_pair_indices).cell1type;
        cell2type = jit_var_out(current_pair_indices).cell2type;

        res1 = jit_var_out(current_pair_indices).cell1_spike_times;
        res2 = jit_var_out(current_pair_indices).cell2_spike_times;
        
        y1 = chanCoords.y(ch1);
        y2 = chanCoords.y(ch2);

        x1 = chanCoords.x(ch1);
        x2 = chanCoords.x(ch2);

        d = sqrt((x2-x1)^2 + (y2-y1)^2);
        
        if d <= 55
            continue
        end
        
        snippetsAtRefChan = load([snippetsDir 'ch' num2str(ch1) 'snippets.mat']);
        snippetsAtTarChan = load([snippetsDir 'ch' num2str(ch2) 'snippets.mat']);

        %%
%         tiledlayout(5,2, 'Padding', 'none', 'TileSpacing', 'compact')
        
        pairStr = [session_name ' (GJ screen) - ' num2str(cellID1) cell1type ' (ch ' num2str(ch1) ', sh ' num2str(shank1) ')' ' v ' ...
                                                  num2str(cellID2) cell2type ' (ch ' num2str(ch2) ', sh ' num2str(shank2) ')' ...
                                                ' distance ' num2str(d) ' microns'];

        nexttile

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
                            
        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(length(GSPExc),1),duration);

        % tracking spike in bins
        [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(res1',res2',binSize,lagLimit);
        
        if ~(GSPExc(round(length(tR)/2)) == 1) 

           close all

           hcomb = figure(102);
           arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

           continue 
        end
        
        % removing the synchronous spikes
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

        % remove sync spikes for waveforms
        res1nosync = setdiff(res1,res1sync);
        res2nosync = setdiff(res2,res2sync);
        
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
        
        channelSetRef = shankChanList{shank1};
        channelSetTar = shankChanList{shank2};
        
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
            
            save_file = fullfile(dataPath, pairStr);
            print(fig_use, save_file,'-djpeg',resolution_use);
            
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

        pairGroupStatsTempTable.refCellExplorerType = {cell1type};
        pairGroupStatsTempTable.tarCellExplorerType = {cell2type};

        pairGroupStatsTempTable.refWaveforms = {unitsGlobalSnapshots.waveforms{unitsGlobalSnapshots.neuronID == pairsAna(loopPairs,1),1}};
        pairGroupStatsTempTable.tarWaveforms = {unitsGlobalSnapshots.waveforms{unitsGlobalSnapshots.neuronID == pairsAna(loopPairs,2),1}};
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
        save([dataPath 'groupStatsRatRoyNULLputativeGJ.mat'],'pairGroupStatTable')
    else
        save([dataPath 'groupStatsRatRoyPutativeGJ.mat'],'pairGroupStatTable')
    end
       
    close all
    
end


                         
                         
                         