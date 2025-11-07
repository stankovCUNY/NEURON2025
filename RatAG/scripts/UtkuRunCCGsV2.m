function UtkuRunCCGsV2(dayNo,nullFlag)
    
    SampleRate = 30000; % SampleRate for Eran

    figpath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    if strcmp(dayNo,'1') || strcmp(dayNo,'1')
        UtkuData = loadUtkuV2(dayNo);
        datapath = [datapath 'AG_2019-12-23_NSD' '/'];
%         load([datapath 'AG_2019-12_23_NSD_adjustedSpikeTimes.mat'])
    elseif strcmp(dayNo,'2')
        UtkuData = loadUtkuV2(dayNo);
        datapath = [datapath 'AG_2019-12-27_NSD' '/'];
%         load([datapath 'AG_2019-12_27_NSD_adjustedSpikeTimes.mat'])
    elseif strcmp(dayNo,'B1')  
        datapath      = [datapath 'B1_2022-06-24_NSD_CA1_24hrs/2022-06-24_NSD_CA1_24hrs/2022-06-24_NSD_CA1_24hrscrs-merged_cleaned.GUI' '/'];
        spikeTimes    = double(readNPY([datapath 'spike_times.npy']))/SampleRate;
        spikeInd      = double(readNPY([datapath 'spike_clusters.npy']));
        clusterinfo   = ratB1clusterInfo;
        channelShanks = double(readNPY([datapath 'channel_shanks.npy'])) + 1;
    end   
    
    if nullFlag
        dataJscale1 = load(['UtkuDay' dayNo 'null_jscale1_alpha5_pairs.mat']);
        dataJscale5 = load(['UtkuDay' dayNo 'null_jscale5_alpha5_pairs.mat']);
    else
        dataJscale1 = load(['AG' dayNo '_jscale1_alpha5_pairs.mat']);
        dataJscale5 = load(['AG' dayNo '_jscale5_alpha5_pairs.mat']);
    end
        
    ExcPairs = [dataJscale1.pairs.ExcPairs(:,1:2); dataJscale5.pairs.ExcPairs(:,1:2)];
    if ~isempty(dataJscale1.pairs.InhPairs) % null set 
        InhPairs = [dataJscale1.pairs.InhPairs(:,1:2); dataJscale5.pairs.InhPairs(:,1:2)];
    else 
        InhPairs = [];
    end
    GapPairs = [dataJscale1.pairs.GapPairs(:,1:2); dataJscale5.pairs.GapPairs(:,1:2)];

    [~,idxUnique,~] = unique(ExcPairs,'rows');
    ExcPairs(idxUnique,:) = [];

    [~,idxUnique,~] = unique(InhPairs,'rows');
    InhPairs(idxUnique,:) = [];

    [~,idxUnique,~] = unique(GapPairs,'rows');
    GapPairs(idxUnique,:) = [];

    if ~isempty(dataJscale1.pairs.InhPairs) % null set 
        allPairs = [ExcPairs; GapPairs];
    else
        allPairs = [ExcPairs; InhPairs; GapPairs];
    end
    [~,idxUnique,~] = unique(allPairs,'rows');
    allPairs = allPairs(idxUnique,:);

    %%

    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = false;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    for i = 1:size(allPairs,1)
        
        tiledlayout(3,1)
        
        if strcmp(dayNo,'1') || strcmp(dayNo,'2')
            
            if nullFlag
                res1    = dataJscale1.pairs.spikeTimes((dataJscale1.pairs.spikeInd) == allPairs(i,1));
                res2    = dataJscale1.pairs.spikeTimes((dataJscale1.pairs.spikeInd) == allPairs(i,2));
            else
                res1    = UtkuData.cell_metrics.spikes.times{1,allPairs(i,1)};
                res2    = UtkuData.cell_metrics.spikes.times{1,allPairs(i,2)};
            end

            cellID1 = UtkuData.cell_metrics.UID(allPairs(i,1));
            cellID2 = UtkuData.cell_metrics.UID(allPairs(i,2));

            ch1     = UtkuData.cell_metrics.maxWaveformCh1(allPairs(i,1));
            ch2     = UtkuData.cell_metrics.maxWaveformCh1(allPairs(i,2));

            clu1    = UtkuData.cell_metrics.cluID(allPairs(i,1));
            clu2    = UtkuData.cell_metrics.cluID(allPairs(i,2));

            shank1  = UtkuData.cell_metrics.shankID(allPairs(i,1));
            shank2  = UtkuData.cell_metrics.shankID(allPairs(i,2));

            cell1type = UtkuData.cell_metrics.putativeCellType{1,allPairs(i,1)};
            cell2type = UtkuData.cell_metrics.putativeCellType{1,allPairs(i,2)};
            
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
                        
            noiseClusterRefShank = load([datapath 'Utkus_RatAG_day' dayNo '_shank' num2str(shank1) '_spikeTimesNoise.mat']);
            noiseClusterTarShank = load([datapath 'Utkus_RatAG_day' dayNo '_shank' num2str(shank2) '_spikeTimesNoise.mat']);

            disp([num2str(cellID1) cell1type ' v ' num2str(cellID2) cell2type ' - ' num2str(i/size(allPairs,1)*100) '% done'])
            
            %%
            
            nexttile
            
            [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);


            pairCCGscreen.pair(i,:)  = [cellID1,cellID2];
            pairCCGscreen.GSPExc{i}  = GSPExc;
            pairCCGscreen.GSPInh{i}  = GSPInh;
            pairCCGscreen.spikeTimes = dataJscale1.pairs.spikeTimes;
            pairCCGscreen.spikeInd   = dataJscale1.pairs.spikeInd;

            if plotFlag

                pairStr = ['AG' dayNo ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ')' ' v ' ...
                                            num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ')'];

                title(pairStr)
                ylims = get(gca,'ylim');
                if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
                
                %%
                
                nexttile
                
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,noiseClusterTarShank.spikeTimesNoise,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);


                pairStr = ['AG' dayNo ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ')' ' v ' 'noise/MUA cluster from sh: ' num2str(shank2)];

                title(pairStr)
                ylims = get(gca,'ylim');
                if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
                
                %%
                
                nexttile
                
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(noiseClusterRefShank.spikeTimesNoise,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);


                pairStr = ['AG' dayNo ' - ' 'noise/MUA cluster from sh: ' num2str(shank1) ' v ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ')'];

                title(pairStr)
                ylims = get(gca,'ylim');
                if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
                
                %%
                
                if (sum(GSPExc) ~= 0) && (sum(GSPInh) ~= 0)
                    save_file = fullfile(figpath, pairStr);
                    print(fig_use, save_file,'-djpeg',resolution_use);
                end
                
                close all

                hcomb = figure(102);
                arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
            end
            
        elseif strcmp(dayNo,'B1')    
            
             res1    = spikeTimes(spikeInd == allPairs(i,1));
             res2    = spikeTimes(spikeInd == allPairs(i,2));
             
             cellID1 = allPairs(i,1);
             cellID2 = allPairs(i,2);

             ch1     = clusterinfo.ch(find(clusterinfo.cluster_id == allPairs(i,1))) + 1;
             ch2     = clusterinfo.ch(find(clusterinfo.cluster_id == allPairs(i,2))) + 1;

             clu1    = allPairs(i,1);
             clu2    = allPairs(i,2);

             shank1  = channelShanks(ch1);
             shank2  = channelShanks(ch2);
             
             [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);


             pairCCGscreen.pair(i,:) = [cellID1,cellID2];
             pairCCGscreen.GSPExc{i} = GSPExc;
             pairCCGscreen.GSPInh{i} = GSPInh;
             
             if plotFlag

                pairStr = [dayNo ' (aged)'  ' - ' num2str(cellID1) ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ')' ' v ' ...
                                                  num2str(cellID2) ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ')'];

                title(pairStr)
                ylims = get(gca,'ylim');
                if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
                if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end

                save_file = fullfile(figpath, pairStr);
                print(fig_use, save_file,'-djpeg',resolution_use);

                close all

                hcomb = figure(102);
                arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
            end
             
        end
    end
    
    if strcmp(dayNo,'1') || strcmp(dayNo,'2')
        if nullFlag
            save([datapath '/AG' dayNo 'null_jitterScreen_pairs.mat'],'pairCCGscreen')
        else
            save([datapath '/AG' dayNo '_jitterScreen_pairs.mat'],'pairCCGscreen')
        end
    elseif strcmp(dayNo,'B1')
        save([datapath '/' dayNo '_jitterScreen_pairs.mat'],'pairCCGscreen')
    end 
    
end