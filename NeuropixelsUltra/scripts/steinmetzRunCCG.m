function steinmetzRunCCG
    
%     parpool(6)

    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    tWave          = (-40:41)/30;

    figPath  = '/home/nasko/CUNY_Work_NPUltraWaveforms/figures';
    dataPath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';

    brainRegion = steinmetzBrainArea;
    brainRegionList = unique(brainRegion);
    
    pairGroupStatTable = table;
    
    
    for loopRegions = 1:length(brainRegionList)
        
        display([num2str((loopRegions/length(brainRegionList))*100) '% done'])
        display(['region: ' brainRegionList{loopRegions}])
        
        if strcmp(brainRegionList{loopRegions},'MOs2/3')
            brainRegionSaveStr = 'MOs2and3';
        elseif strcmp(brainRegionList{loopRegions},'ORBl2/3')
            brainRegionSaveStr = 'ORBl2and3';
        elseif strcmp(brainRegionList{loopRegions},'VISam2/3')
            brainRegionSaveStr = 'VISam2and3';
        elseif strcmp(brainRegionList{loopRegions},'VISpm2/3')
            brainRegionSaveStr = 'VISpm2and3';
        else
            brainRegionSaveStr = brainRegionList{loopRegions};
        end
        
        load([dataPath 'neuronsSteinmetzDataArea' brainRegionSaveStr '.mat']);
        
        dataJscale1 = load([dataPath 'SteinmetzMouseRegion' brainRegionSaveStr '_jscale1_alpha5_pairs.mat']);
        dataJscale5 = load([dataPath 'SteinmetzMouseRegion' brainRegionSaveStr '_jscale5_alpha5_pairs.mat']);

        if ~isempty(dataJscale1.pairs.ExcPairs) 
            ExcPairs = [dataJscale1.pairs.ExcPairs(:,1:2); dataJscale5.pairs.ExcPairs(:,1:2)];
        else
            ExcPairs = [];
        end 
        
        if ~isempty(dataJscale1.pairs.InhPairs)
            InhPairs = [dataJscale1.pairs.InhPairs(:,1:2); dataJscale5.pairs.InhPairs(:,1:2)];
        else
            InhPairs = [];
        end
            
        if ~isempty(dataJscale1.pairs.GapPairs)
            GapPairs = [dataJscale1.pairs.GapPairs(:,1:2); dataJscale5.pairs.GapPairs(:,1:2)];
        else
            GapPairs = [];
        end

        [~,idxUnique,~] = unique(ExcPairs,'rows');
        ExcPairs(idxUnique,:) = [];

        [~,idxUnique,~] = unique(InhPairs,'rows');
        InhPairs(idxUnique,:) = [];
    
        [~,idxUnique,~] = unique(GapPairs,'rows');
        GapPairs(idxUnique,:) = [];
    
        allPairs = [ExcPairs; InhPairs; GapPairs];

        [~,idxUnique,~] = unique(allPairs,'rows');
        allPairs = allPairs(idxUnique,:);

        %%

       
        % Monitor specific plot settings.
        screensize = get(0,'screensize');
        % initiate figure
        hcomb = figure(102);

        res_type = 'QHD';
        pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

        for i = 1:size(allPairs,1)

            display([num2str(100*(i/size(allPairs,1))) '% done'])

            idx1    = find(neurons.unitID == allPairs(i,1));
            idx2    = find(neurons.unitID == allPairs(i,2));

            res1    = neurons.spikeTimes{1,idx1};
            res2    = neurons.spikeTimes{1,idx2};
            
%             if length(res1) == length(res2)
%                continue 
%             end

            cellID1 = allPairs(i,1);
            cellID2 = allPairs(i,2);
            
            cell1type = neurons.putativeCellType{idx1};
            cell2type = neurons.putativeCellType{idx2};

            ch1     = double(neurons.peakChan(idx1));
            ch2     = double(neurons.peakChan(idx2));

            probe1  = num2str(neurons.probeID{idx1});
            probe2  = num2str(neurons.probeID{idx2});
            
            y1 = double(neurons.probeY(ch1));
            y2 = double(neurons.probeY(ch2));

            x1 = double(neurons.probeX(ch1));
            x2 = double(neurons.probeX(ch2));
            
            d = sqrt((x2-x1)^2 + (y2-y1)^2);
    
            region1 = char(neurons.brainRegion{idx1});
            region2 = char(neurons.brainRegion{idx2});

            refAvgWaveMaxChan = neurons.avgWaveforms{idx1}(:,ch1);
            refAvgWaveOpsChan = neurons.avgWaveforms{idx1}(:,ch2);

            tarAvgWaveMaxChan = neurons.avgWaveforms{idx2}(:,ch2);
            tarAvgWaveOpsChan = neurons.avgWaveforms{idx2}(:,ch1);

            if plotFlag
                tiledlayout(5,2)
                nexttile(1)
            end

            [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
                      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);

            pairCCGscreen.pair(i,:) = [cellID1,cellID2];
            pairCCGscreen.GSPExc{i} = GSPExc;
            pairCCGscreen.GSPInh{i} = GSPInh;

            if (sum(GSPExc) == 0) && (sum(GSPInh) == 0)

               close all

               hcomb = figure(102);
               arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

               continue 
            end

            if plotFlag

                pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', probe: ' probe1 ', region: ' region1 ')' ' v ' ...
                           num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', probe: ' probe2 ', region: ' region2 ')' ...
                          ' peak chan to peak chan distance: ' num2str(d) ' microns'];

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
                
                %%
                
                pairGroupStatsTempTable = table;
                
                pairGroupStatsTempTable.probeID       = num2str(neurons.probeID{idx1});
                pairGroupStatsTempTable.brainRegion   = {char(neurons.brainRegion{idx1})};
                
                pairGroupStatsTempTable.refNeuronID   = cellID1;
                pairGroupStatsTempTable.tarNeuronID   = cellID2;
                
                pairGroupStatsTempTable.refSpikeTimes = {res1};
                pairGroupStatsTempTable.tarSpikeTimes = {res2};
                
                pairGroupStatsTempTable.refChannel    = ch1;
                pairGroupStatsTempTable.tarChannel    = ch2;
                
                pairGroupStatsTempTable.refCellExplorerType = {cell1type};
                pairGroupStatsTempTable.tarCellExplorerType = {cell2type};
                
                pairGroupStatsTempTable.refWaveforms = {neurons.avgWaveforms{idx1}};
                pairGroupStatsTempTable.tarWaveforms = {neurons.avgWaveforms{idx2}};
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
                
                pairGroupStatTable = [pairGroupStatTable; pairGroupStatsTempTable];
                
                %%
                
                nexttile(3) 
                histogram(1000*diff(res1),'BinWidth',1); hold on; histogram(1000*diff(res2),'BinWidth',1); hold off
                xlim([0 150])
                xlabel('[ms]')
                
                nexttile(5)
                plot(1000*tR,ccgR(:,1,1),'LineWidth',2); hold on; plot(1000*tR,ccgR(:,2,2),'LineWidth',2); hold off
                
                nexttile(7)
                plot(tWave,refAvgWaveMaxChan/1000,'LineWidth',2)
                hold on
                plot(tWave,refAvgWaveOpsChan/1000,'LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                xlabel('time [ms]')
                ylabel('revere order voltage [mV]')
                title(['Reference e-spike: ' num2str(cellID1) ' (ch: ' num2str(ch1) ', probe: ' probe1 ', region: ' region1 ')'])
                legend('ref. activity at ref. max chan','ref. activity at tar. max chan')
                xlim([-3.5 3.5])
%                 xlim([-10 10])
    %             xlim([-1 1])

                nexttile(9)
                plot(tWave,tarAvgWaveMaxChan/1000,'LineWidth',2)
                hold on
                plot(tWave,tarAvgWaveOpsChan/1000,'LineWidth',2)
                hold off
                set(gca, 'YDir','reverse')
                set(gca, 'XDir','reverse')
                xlabel('reverse order time [ms]')
                ylabel('revere order voltage [mV]')
                title(['Target e-spike: ' num2str(cellID2) ' (ch: ' num2str(ch2) ', probe: ' probe2 ', region: ' region2 ')'])
                legend('tar. activity at tar. max chan','tar. activity at ref. max chan')
                xlim([-3.5 3.5])
%                 xlim([-10 10])
    %             xlim([-1 1])
                
                %% activity plot
                
                nexttile([5 1])
                
                waveform = neurons.avgWaveforms{idx1} + neurons.avgWaveforms{idx2};

                troughIdx = 43;

                chanLayout = flipud(reshape(1:384,[8,48])');

                x = 0:5:5*7;
                y = 0:5:5*47;

                for loopChans = 1:(48*8)

                    [a,b] = find(chanLayout == loopChans);
                    coordinates(loopChans,2) = y(a);
                    coordinates(loopChans,1) = x(b);

                    heat_values(loopChans) = waveform(troughIdx,loopChans);
                end

                % Create a 3D surface plot
                tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
                trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
                xlabel('x [\mum]');
                ylabel('y [\mum]');
                title('Activity Plot');

                axis equal

                % xlim([0 35])
                % ylim([0 235])

                % Add colorbar
                c = colorbar;
                c.Label.String = '[\muV]';
                colormap parula
                cmp = colormap;
                cmp = flipud(cmp);
                colormap(cmp);

                % Adjust view for better visualization
                view(0, 90);

                hold on
                for loopChans = 1:(48*8)
                    plot3(0.05*(1:size(waveform,1))  + coordinates(loopChans,1) - 2.5, ...
                          0.02*waveform(:,loopChans) + coordinates(loopChans,2), ...
                          max(waveform(troughIdx,:))*ones(size(waveform,1),1),'k','LineWidth',1)
                end
                hold off
                
                %%
    
                save_file = fullfile(figPath, pairStr);
                print(fig_use, save_file,'-djpeg',resolution_use);

                close all

                hcomb = figure(102);
                arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
            end

        end
    end
    
	save([dataPath 'groupStatsSteinmetzDataset.mat'],'pairGroupStatTable')
    
end