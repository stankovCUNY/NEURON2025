function BapunRunCCG(RatID,nullFlag)
    
    data_dir = ['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/Rat' RatID '/'];

    if strcmp(RatID,'N')       
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
    elseif strcmp(RatID,'S')
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
    elseif strcmp(RatID,'U')
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');        
    end
    
    if nullFlag
        dataJscale1 = load([data_dir 'Bapuns_Rat' RatID 'null_jscale1_alpha5_pairs.mat']);
        dataJscale5 = load([data_dir 'Bapuns_Rat' RatID 'null_jscale5_alpha5_pairs.mat']);
    else
        dataJscale1 = load([data_dir 'Bapuns_Rat' RatID '_jscale1_alpha5_pairs.mat']);
        dataJscale5 = load([data_dir 'Bapuns_Rat' RatID '_jscale5_alpha5_pairs.mat']);
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

    figpath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    for loopPairs = 1:size(allPairs,1)
        
        disp([num2str(100*(loopPairs/size(allPairs,1))) '% done'])
        
        if nullFlag
            res1    = dataJscale1.pairs.spikeTimes((dataJscale1.pairs.spikeInd) == allPairs(loopPairs,1));
            res2    = dataJscale1.pairs.spikeTimes((dataJscale1.pairs.spikeInd) == allPairs(loopPairs,2));
        else
            res1    = neurons.spiketrains{1,allPairs(loopPairs,1)};
            res2    = neurons.spiketrains{1,allPairs(loopPairs,2)};
        end

        cellID1 = allPairs(loopPairs,1);
        cellID2 = allPairs(loopPairs,2);

        ch1     = double(neurons.peak_channels(allPairs(loopPairs,1))) + 1;
        ch2     = double(neurons.peak_channels(allPairs(loopPairs,2))) + 1;

        shank1  = double(neurons.shank_ids(allPairs(loopPairs,1))) + 1;
        shank2  = double(neurons.shank_ids(allPairs(loopPairs,2))) + 1;
        
        cell1type = neurons.neuron_type(allPairs(loopPairs,1),:);
        cell2type = neurons.neuron_type(allPairs(loopPairs,2),:);
        
        % reject MUAs at this stage
        if strcmp(cell1type,'mua  ') || strcmp(cell2type,'mua  ')
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
        
        if plotFlag
            if strcmp(RatID,'N')
                muaClusterRefShank    = load([data_dir 'Bapuns_Rat' RatID '_shank_' num2str(shank1) '_spikeTimesMUAs.mat']);
                noiseClusterRefShank  = load([data_dir 'Bapuns_Rat' RatID '_shank_' num2str(shank1) '_spikeTimesNoise.mat']);

                muaClusterTarShank    = load([data_dir 'Bapuns_Rat' RatID '_shank_' num2str(shank1) '_spikeTimesMUAs.mat']);
                noiseClusterTarShank  = load([data_dir 'Bapuns_Rat' RatID '_shank_' num2str(shank1) '_spikeTimesNoise.mat']);
            else
                muaCluster            = load([data_dir 'Bapuns_Rat' RatID '_spikeTimesMUAs.mat']);
                noiseCluster          = load([data_dir 'Bapuns_Rat' RatID '_spikeTimesNoise.mat']);
            end
            
            tiledlayout(5,1)
            nexttile
        end
        
        %% pair CCG
        
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant, 'plot_ACbands',true);


        pairCCGscreen.pair(loopPairs,:) = [cellID1,cellID2];
        pairCCGscreen.GSPExc{loopPairs} = GSPExc;
        pairCCGscreen.GSPInh{loopPairs} = GSPInh;
        pairCCGscreen.spikeTimes        = dataJscale1.pairs.spikeTimes;
        pairCCGscreen.spikeInd          = dataJscale1.pairs.spikeInd;

        if plotFlag

            pairStr = ['Rat' RatID ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' ...
                                         num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];

            title(pairStr)
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
            
            %% ref vs MUA cluster
            
            nexttile
            if strcmp(RatID,'N')
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,muaClusterTarShank.spikeTimesMUAs,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
            
                pairStr = ['Rat' RatID ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' 'MUA cluster from sh: ' num2str(shank2)];
                
            else
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,muaCluster.spikeTimesMUAs,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
            
                pairStr = ['Rat' RatID ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' 'collated MUAs cluster'];
            end
            
            title(pairStr)
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
            
            %% ref vs noise cluster
            
            nexttile
            if strcmp(RatID,'N')
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,noiseClusterTarShank.spikeTimesNoise ,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
            
                pairStr = ['Rat' RatID ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' 'noise cluster from sh: ' num2str(shank2)];
                
            else
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(res1,noiseCluster.spikeTimesNoise,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
            
                pairStr = ['Rat' RatID ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ')' ' v ' 'collated noise cluster'];
            end
            
            title(pairStr)
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
            
            %% MUA cluster vs tar
            
            nexttile
            if strcmp(RatID,'N')
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(muaClusterRefShank.spikeTimesMUAs,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
            
                pairStr = ['Rat' RatID ' - ' 'MUA cluster from sh: ' num2str(shank1) ' v ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
                
            else
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(muaCluster.spikeTimesMUAs,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
            
                pairStr = ['Rat' RatID ' - ' 'collated MUA cluster' ' v ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
            end
            
            title(pairStr)
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
            
            %% noise cluster vs tar
            
            nexttile
            if strcmp(RatID,'N')
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(noiseClusterRefShank.spikeTimesNoise,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
                            
                pairStr = ['Rat' RatID ' - ' 'noise cluster from sh: ' num2str(shank1) ' v ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
            else
                [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                      CCG_jitter(noiseCluster.spikeTimesNoise,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                'plot_output', get(fig_use, 'Number'), ...
                                'njitter', njitter, 'alpha', alpha,...
                                'for_grant', for_grant);
                            
                pairStr = ['Rat' RatID ' - ' 'collated noise cluster' ' v ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ')'];
            end
            
            title(pairStr)
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv'); end
            
            %%
            
%             save_file = fullfile(figpath, pairStr);
%             print(fig_use, save_file,'-djpeg',resolution_use);

            close all

            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
        end

    end
    
    if nullFlag
        save([data_dir '/Bapuns_Rat' RatID 'null_jitterScreen_pairs.mat'],'pairCCGscreen')
    else
        save([data_dir '/Bapuns_Rat' RatID '_jitterScreen_pairs.mat'],'pairCCGscreen')
    end
    
end