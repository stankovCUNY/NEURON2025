function BOTrunCCGputativeGJ(nullFlag)
     
     % jitter constants
     jscale         = 1;
     alpha_name     = 5;
     duration       = 0.007;
     fs             = 30000;
     binSize        = 1/1000;
     fig_use        = 102;
     njitter        = 500;
     alpha          = 0.05;
     for_grant      = false;
     plotFlag       = true;
     plotOptoCCG    = false;
     plotUnionCCG   = false;
     resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
     tWave          = (-19:62)/30;
     
     % paths
     dataPath      = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
     figPath       = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/figures';
     unionClusPath = '/media/nasko/WD_BLACK4/unionClusters/';     

     % obtain all sessions in the ephys dataset
     sessions = bot.listSessions('ephys');
    
     % obtain all sessions with the desired session_type and then select the first one
%      desiredSessions = sessions(sessions.session_type == 'brain_observatory_1.1', :); % use Brain Observatory for now
%      desiredSessions = sessions(sessions.session_type == 'functional_connectivity', :); % use Brain Observatory for now
     desiredSessions = sessions;

     animalsList = desiredSessions.id;
    
     for loopAnimal = 1:length(animalsList)
        
        animalID = animalsList(loopAnimal);
         
        display(['animal: ' num2str(animalID)])
        
        if ~exist([figPath '/' num2str(animalID)],'dir')
            mkdir([figPath '/' num2str(animalID)])
        end
        
        load([dataPath ['neuronsMouse' num2str(animalsList(loopAnimal)) '.mat']])
            
        probeString = string(neurons.probeName);
        probeList   = unique(probeString);
        
        pairGroupStatsTable = table;

        %% optotag stim times
        
        if strcmp(desiredSessions.full_genotype{find(desiredSessions.id == animalID)},"Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt")
            geneStr = 'Sst';
        elseif strcmp(desiredSessions.full_genotype{find(desiredSessions.id == animalID)},"Vip-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt")
            geneStr = 'Vip';
        elseif strcmp(desiredSessions.full_genotype{find(desiredSessions.id == animalID)},"Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt")
            geneStr = 'Pvalb';
        else
            geneStr = '';
        end
        
        if ~isempty(geneStr)
            load([dataPath 'optoStimData_' geneStr '_' num2str(animalID) '.mat'])
            idx = find((optoStimData.duration > 0.009) & (optoStimData.duration < 0.02));

            optoStimTrials.start_time    = optoStimData.start_time(idx);
            optoStimTrials.stop_time     = optoStimData.stop_time(idx);
            optoStimTrials.condition     = optoStimData.condition(idx,:);
            optoStimTrials.level         = optoStimData.level(idx);
            optoStimTrials.duration      = optoStimData.duration(idx);
            optoStimTrials.stimulus_name = optoStimData.stimulus_name(idx,:);
            
        end
        
        %%
        
        for loopProbe = 1:length(probeList)
            
            display(['probe: ' num2str(probeList(loopProbe))])
            
            probeFilter  = strcmp(probeList(loopProbe),probeString);    
            neuronsProbe = neurons(probeFilter,:);

            brainRegionString = string(neuronsProbe.brainRegion);
            brainRegionList   = unique(brainRegionString);
        
            % load([unionClusPath 'validCluster_mouse'   num2str(animalID) '_' num2str(probeList(loopProbe)) '.mat']);
            % load([unionClusPath 'invalidCluster_mouse' num2str(animalID) '_' num2str(probeList(loopProbe)) '.mat']);
            % 
            % unionSpikeTimes = sort([validCluster.spikeTimes; invalidCluster.spikeTimes]);
            
            for loopBrainRegions = 1:length(brainRegionList)
                
                display(['region: ' brainRegionList{loopBrainRegions}])
                
                brainRegionFilter  = strcmp(brainRegionList(loopBrainRegions),brainRegionString); 
                neuronsBrainRegion = neuronsProbe(brainRegionFilter,:); 
                
                if nullFlag
                    dataJscale1 = load([dataPath 'BOTmouse' num2str(animalID) 'NULL' num2str(probeList(loopProbe)) ...
                                        'region' brainRegionList{loopBrainRegions} '_jscale1_alpha5_pairs.mat']);
%                     dataJscale5 = load([dataPath 'BOTmouse' num2str(animalID) 'NULL' num2str(probeList(loopProbe)) ...
%                                         'region' brainRegionList{loopBrainRegions} '_jscale5_alpha5_pairs.mat']);
                else
                    dataJscale1 = load([dataPath 'BOTmouse' num2str(animalID) num2str(probeList(loopProbe)) ...
                                        'region' brainRegionList{loopBrainRegions} '_jscale1_alpha5_pairs.mat']);
%                     dataJscale5 = load([dataPath 'BOTmouse' num2str(animalID) num2str(probeList(loopProbe)) ...
%                                         'region' brainRegionList{loopBrainRegions} '_jscale5_alpha5_pairs.mat']);
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

%                 [~,idxUnique,~] = unique(ExcPairs,'rows');
%                 ExcPairs(idxUnique,:) = [];
% 
%                 [~,idxUnique,~] = unique(InhPairs,'rows');
%                 InhPairs(idxUnique,:) = [];
% 
%                 [~,idxUnique,~] = unique(GapPairs,'rows');
%                 GapPairs(idxUnique,:) = [];

                allPairs = [ExcPairs; InhPairs; GapPairs];

                [~,idxUnique,~] = unique(allPairs,'rows');
                allPairs = allPairs(idxUnique,:);

                %%

                % Monitor specific plot settings.
                screensize = get(0,'screensize');
                % initiate figure
                hcomb = figure(102);

                res_type = 'QHD';
                pos = [70 230 1440 2560]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
                arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

                for loopPairs = 1:size(allPairs,1)

                    display([num2str(100*(loopPairs/size(allPairs,1))) '% done'])

                    idx1      = find(neuronsBrainRegion.unitID == allPairs(loopPairs,1));
                    idx2      = find(neuronsBrainRegion.unitID == allPairs(loopPairs,2));
                    
                    if nullFlag
                        res1    = dataJscale1.pairs.spikeTimes((dataJscale1.pairs.spikeInd) == allPairs(loopPairs,1));
                        res2    = dataJscale1.pairs.spikeTimes((dataJscale1.pairs.spikeInd) == allPairs(loopPairs,2));
                    else
                        res1    = neuronsBrainRegion.spikeTimes{idx1};
                        res2    = neuronsBrainRegion.spikeTimes{idx2};
                    end
                    
                    % refUnion  = setdiff(unionSpikeTimes,res1);
                    % tarUnion  = setdiff(unionSpikeTimes,res2);
                    
                    % refUnion  = refUnion(sort(randperm(length(refUnion),1e6)));
                    % tarUnion  = tarUnion(sort(randperm(length(tarUnion),1e6)));

                    cellID1   = allPairs(loopPairs,1);
                    cellID2   = allPairs(loopPairs,2);
                    
                    cell1type = neuronsBrainRegion.putativeCellType{idx1};
                    cell2type = neuronsBrainRegion.putativeCellType{idx2};
                    
                    % 2X cutoff optotagging
                    % if strcmp(neuronsBrainRegion.putativeOptotagType2X{idx1},'No opto label')
                    %     cell1opto = '';
                    % else 
                    %     cell1opto = [' (' neuronsBrainRegion.putativeOptotagType2X{idx1} ')'];
                    %     % res1optotagStim = extractOptotagStimSpikeTimes(res1,optoStimTrials);
                    % end
                    % if strcmp(neuronsBrainRegion.putativeOptotagType2X{idx2},'No opto label')
                    %     cell2opto = '';
                    % else
                    %     cell2opto = [' (' neuronsBrainRegion.putativeOptotagType2X{idx2} ')'];
                    %     % res2optotagStim = extractOptotagStimSpikeTimes(res2,optoStimTrials);
                    % end
                    
                    % 5X cutoff optotagging
                    if strcmp(neuronsBrainRegion.putativeOptotagType5X{idx1},'No opto label')
                        cell1opto = '';
                    else 
                        cell1opto = [' (' neuronsBrainRegion.putativeOptotagType2X{idx1} ')'];
                        res1optotagStim = extractOptotagStimSpikeTimes(res1,optoStimTrials);
                    end
                    if strcmp(neuronsBrainRegion.putativeOptotagType5X{idx2},'No opto label')
                        cell2opto = '';
                    else
                        cell2opto = [' (' neuronsBrainRegion.putativeOptotagType2X{idx2} ')'];
                        res2optotagStim = extractOptotagStimSpikeTimes(res2,optoStimTrials);
                    end

                    ch1       = double(neuronsBrainRegion.peakChan(idx1));
                    ch2       = double(neuronsBrainRegion.peakChan(idx2));

                    probe1    = char(neuronsBrainRegion.probeName{idx1});
                    probe2    = char(neuronsBrainRegion.probeName{idx2});
                    
                    y1 = double(neuronsBrainRegion.probeVertPos(idx1));
                    y2 = double(neuronsBrainRegion.probeVertPos(idx2));
                    
                    x1 = double(neuronsBrainRegion.probeHortPos(idx1));
                    x2 = double(neuronsBrainRegion.probeHortPos(idx2));
                    
                    d = sqrt((x2-x1)^2 + (y2-y1)^2);
                    
                    if d <= 55
                        continue 
                    end
                    
                    region1 = char(neuronsBrainRegion.brainRegion{idx1});
                    region2 = char(neuronsBrainRegion.brainRegion{idx2});

                    refAvgWaveMaxChan = neuronsBrainRegion.avgWaveforms{idx1}(ch1,:);
                    refAvgWaveOpsChan = neuronsBrainRegion.avgWaveforms{idx1}(ch2,:);

                    tarAvgWaveMaxChan = neuronsBrainRegion.avgWaveforms{idx2}(ch2,:);
                    tarAvgWaveOpsChan = neuronsBrainRegion.avgWaveforms{idx2}(ch1,:);

%                     if plotFlag
%                         if plotUnionCCG
%                             tiledlayout(5,1)
%                             nexttile
%                         else
%                             tiledlayout(3,1)
%                             nexttile
%                         end
%                     end
                    
                    if ~isempty(cell1opto) || ~isempty(cell2opto)
                        yyaxis left
                    end
                        
                    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
                              CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                        'plot_output', get(fig_use, 'Number'), ...
                                        'njitter', njitter, 'alpha', alpha,...
                                        'for_grant', for_grant,'plot_ACbands',true);

                    if ~(GSPExc(round(length(tR)/2)) == 1) 

                       close all

                       hcomb = figure(102);
                       arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

                       continue 
                    end

                    %% group stats data table
                    
                    pairGroupStatsTempTable = table;

                    pairGroupStatsTempTable.animalID              = animalID;
                    pairGroupStatsTempTable.probeID               = {probeList{loopProbe}};
                    pairGroupStatsTempTable.brainRegion           = {brainRegionList{loopBrainRegions}};
                    
                    pairGroupStatsTempTable.refNeuronID           = cellID1;
                    pairGroupStatsTempTable.tarNeuronID           = cellID2;

                    pairGroupStatsTempTable.refSpikeTimes         = {res1};
                    pairGroupStatsTempTable.tarSpikeTimes         = {res2};
                    
                    pairGroupStatsTempTable.refChannel            = ch1;
                    pairGroupStatsTempTable.tarChannel            = ch2;
                    
                    pairGroupStatsTempTable.refCellExplorerType   = {cell1type};
                    pairGroupStatsTempTable.tarCellExplorerType   = {cell2type};

                    pairGroupStatsTempTable.refOptoType2X         = {neuronsBrainRegion.putativeOptotagType2X{idx1}};
                    pairGroupStatsTempTable.tarOptoType2X         = {neuronsBrainRegion.putativeOptotagType2X{idx2}};

                    pairGroupStatsTempTable.refOptoType5X         = {neuronsBrainRegion.putativeOptotagType5X{idx1}};
                    pairGroupStatsTempTable.tarOptoType5X         = {neuronsBrainRegion.putativeOptotagType5X{idx2}};

                    pairGroupStatsTempTable.refWaveforms          = {neuronsBrainRegion.avgWaveforms{idx1}};
                    pairGroupStatsTempTable.tarWaveforms          = {neuronsBrainRegion.avgWaveforms{idx2}};
                    pairGroupStatsTempTable.tWave                 = {tWave};

                    pairGroupStatsTempTable.refSNR                = neuronsBrainRegion.snr(idx1);
                    pairGroupStatsTempTable.tarSNR                = neuronsBrainRegion.snr(idx2); 

                    pairGroupStatsTempTable.refFiringRate         = neuronsBrainRegion.firingRate(idx1);
                    pairGroupStatsTempTable.tarFiringRate         = neuronsBrainRegion.firingRate(idx2); 

                    pairGroupStatsTempTable.refTroughToPeakLength = neuronsBrainRegion.troughToPeakLength(idx1);
                    pairGroupStatsTempTable.tarTroughToPeakLength = neuronsBrainRegion.troughToPeakLength(idx2);

                    % pairGroupStatsTempTable.refUnion            = {refUnion};
                    % pairGroupStatsTempTable.tarUnion            = {tarUnion};
                    
                    pairGroupStatsTempTable.pairDistance          = d;

                    pairGroupStatsTempTable.pairRawCCG            = {ccgR(:,1,2)};
                    pairGroupStatsTempTable.refRawACG             = {ccgR(:,1,1)};
                    pairGroupStatsTempTable.tarRawACG             = {ccgR(:,2,2)};
                    pairGroupStatsTempTable.jitterMean            = {ccgjm};
                    pairGroupStatsTempTable.CCGbinLagTimes        = {tR};
                    pairGroupStatsTempTable.GSPExc                = {GSPExc};
                    pairGroupStatsTempTable.GSPInh                = {GSPInh};
                    pairGroupStatsTempTable.pvalE                 = {pvalE};
                    pairGroupStatsTempTable.pvalI                 = {pvalI};
                    pairGroupStatsTempTable.LSPExc                = {LSPExc};
                    pairGroupStatsTempTable.LSPInh                = {LSPInh};
                    
                    pairGroupStatsTable = [pairGroupStatsTable; pairGroupStatsTempTable];

                    %% extra plotting options

                    if plotFlag

                        pairStr = ['Mouse (GJ screen): ' num2str(animalID) ' - ' num2str(cellID1) cell1type cell1opto ' (ch: ' num2str(ch1) ', probe: ' probe1 ', region: ' region1 ')' ' v ' ...
                                                                     num2str(cellID2) cell2type cell2opto ' (ch: ' num2str(ch2) ', probe: ' probe2 ', region: ' region2 ')'   ...
                                                                     ' peak chan to peak chan distance: ' num2str(d) ' microns'];

                        title(pairStr,'FontSize',8)

                        ylims = get(gca,'ylim');

                        if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                        if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
                        
                        if plotOptoCCG
                            if ~isempty(cell1opto2X) || ~isempty(cell2opto2X)
                                yyaxis right

                                if ~isempty(cell1opto2X) && isempty(cell2opto2X)
                                    ccg = CCG([res1optotagStim;res2],[ones(size(res1optotagStim));2*ones(size(res2))], ...
                                                'binSize', binSize, 'duration', duration, 'Fs', 1/fs, 'norm', 'counts');
                                    ccg   = ccg(:,1,2);
                                    plot(tR*1000,ccg)
                                    legend('ar opto stim spike times v ref')
                                    ylabel('Count')
                                elseif ~isempty(cell1opto2X) && ~isempty(cell2opto2X)
                                    ccg = CCG([res1optotagStim;res2],[ones(size(res1optotagStim));2*ones(size(res2))], ...
                                                'binSize', binSize, 'duration', duration, 'Fs', 1/fs, 'norm', 'counts');
                                    ccg   = ccg(:,1,2);
                                    plot(tR*1000,ccg)  

                                    ccg = CCG([res1;res2optotagStim],[ones(size(res1));2*ones(size(res2optotagStim))], ...
                                                'binSize', binSize, 'duration', duration, 'Fs', 1/fs, 'norm', 'counts');
                                    ccg   = ccg(:,1,2);
                                    hold on
                                    plot(tR*1000,ccg)  
                                    hold off
                                    
                                    legend('tar opto stim spike times v ref','tar v ref opto stim spike times')
                                    ylabel('Count')
                                elseif isempty(cell1opto2X) && ~isempty(cell2opto2X)
                                    ccg = CCG([res1;res2optotagStim],[ones(size(res1));2*ones(size(res2optotagStim))], ...
                                                'binSize', binSize, 'duration', duration, 'Fs', 1/fs, 'norm', 'counts');
                                    ccg   = ccg(:,1,2);
                                    plot(tR*1000,ccg)
                                    legend('tar v ref opto stim spike times')
                                    ylabel('Count')
                                end
                            end
                        end
                        
                         % plot waveforms
%                         nexttile
%                         plot(tWave,refAvgWaveMaxChan/1000,'LineWidth',2)
%                         hold on
%                         plot(tWave,refAvgWaveOpsChan/1000,'LineWidth',2)
%                         hold off
%                         set(gca, 'YDir','reverse')
%                         xlabel('time [ms]')
%                         ylabel('revere order voltage [mV]')
%                         title(['Reference e-spike: ' num2str(cellID1) ' (ch: ' num2str(ch1) ', probe: ' probe1 ', region: ' region1 ')'])
%                         legend('ref. activity at ref. max chan','ref. activity at tar. max chan')
%                         xlim([-3.5 3.5])
%             %             xlim([-1 1])
% 
%                         nexttile
%                         plot(tWave,tarAvgWaveMaxChan/1000,'LineWidth',2)
%                         hold on
%                         plot(tWave,tarAvgWaveOpsChan/1000,'LineWidth',2)
%                         hold off
%                         set(gca, 'YDir','reverse')
%                         set(gca, 'XDir','reverse')
%                         xlabel('reverse order time [ms]')
%                         ylabel('revere order voltage [mV]')
%                         title(['Target e-spike: ' num2str(cellID2) ' (ch: ' num2str(ch2) ', probe: ' probe2 ', region: ' region2 ')'])
%                         legend('tar. activity at tar. max chan','tar. activity at ref. max chan')
%                         xlim([-3.5 3.5])
%             %             xlim([-1 1])
                        
                        % plot union CCGs
                        if plotUnionCCG

                            nexttile
                            [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                                  CCG_jitter(res1,tarUnion,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                            'plot_output', get(fig_use, 'Number'), ...
                                            'njitter', njitter, 'alpha', alpha,...
                                            'for_grant', for_grant,'plot_ACbands',true);
                            
                            title('reference unit vs shank union of all valid/invalid units (minus reference unit)')
                            ylims = get(gca,'ylim');
    
                            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
                            
                            nexttile
                            [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                                  CCG_jitter(refUnion,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                            'plot_output', get(fig_use, 'Number'), ...
                                            'njitter', njitter, 'alpha', alpha,...
                                            'for_grant', for_grant,'plot_ACbands',true);
                            
                            title('shank union of all valid/invalid units (minus target unit) vs target unit')
                            ylims = get(gca,'ylim');

                            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end   

                        end
                        
                        % save plots
                        save_file = fullfile([figPath '/' num2str(animalsList(loopAnimal))], [pairStr '.jpg']);
                        print(fig_use, save_file,'-djpeg',resolution_use);

                        close all

                        hcomb = figure(102);
                        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
                    end

                end
            end
        end
        
        if nullFlag
            save([dataPath 'groupStatsMousePutativeGJ' num2str(animalID) 'NULL.mat'],'pairGroupStatsTable')
        else
            save([dataPath 'groupStatsMousePutativeGJ' num2str(animalID) '.mat'],    'pairGroupStatsTable')
        end
        
     end
    
end

function optotagStimSpikeTimes = extractOptotagStimSpikeTimes(spikeTimes,optoStimTrials)

    startBuffer = 0.001;
    stopBuffer  = 0.009;

    optotagStimSpikeTimes = [];
    for loopTrials = 1:length(optoStimTrials.start_time)
        in_range = (spikeTimes > (optoStimTrials.start_time(loopTrials) + startBuffer)) & ...
                   (spikeTimes < (optoStimTrials.start_time(loopTrials) + stopBuffer));
               
        optotagStimSpikeTimes = [spikeTimes(in_range); optotagStimSpikeTimes];
    end

end