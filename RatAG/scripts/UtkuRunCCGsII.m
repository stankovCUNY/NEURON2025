function UtkuRunCCGsII(dayNo)

    UtkuData      = loadUtkuV2(dayNo); % NOTE: spike time data is in minutes!!!
    filteredPairs = ii_conv_list(dayNo);

    %%

    jscale         = 5;
    alpha_name     = 5;
    duration       = 1;
    fs             = 30000;
    binSize        = 1/1000;
%     binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    
    figpath     = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    UtkuPath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    day = str2double(dayNo);
    if day == 1
        datapath = [UtkuPath 'AG_2019-12-23_NSD' '/'];
    elseif day == 2
        datapath = [UtkuPath 'AG_2019-12-27_NSD' '/'];
    end 

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));


    for i = 1:size(filteredPairs,1)

        res1    = UtkuData.cell_metrics.spikes.times{1,filteredPairs(i,1)};
        res2    = UtkuData.cell_metrics.spikes.times{1,filteredPairs(i,2)};

        cellID1 = UtkuData.cell_metrics.UID(filteredPairs(i,1));
        cellID2 = UtkuData.cell_metrics.UID(filteredPairs(i,2));

        ch1     = UtkuData.cell_metrics.maxWaveformCh1(filteredPairs(i,1));
        ch2     = UtkuData.cell_metrics.maxWaveformCh1(filteredPairs(i,2));

        clu1    = UtkuData.cell_metrics.cluID(filteredPairs(i,1));
        clu2    = UtkuData.cell_metrics.cluID(filteredPairs(i,2));
        
        shank1  = UtkuData.cell_metrics.shankID(filteredPairs(i,1));
        shank2  = UtkuData.cell_metrics.shankID(filteredPairs(i,2));
        
        cell1type = UtkuData.cell_metrics.putativeCellType{1,filteredPairs(i,1)};
        cell2type = UtkuData.cell_metrics.putativeCellType{1,filteredPairs(i,2)};
        
        region1 = UtkuData.cell_metrics.brainRegion{1,filteredPairs(i,1)}; 
        region2 = UtkuData.cell_metrics.brainRegion{1,filteredPairs(i,2)};
        
        if strcmp(cell1type,'Pyramidal Cell')
            cell1type = 'p';
        elseif strcmp(cell1type,'Narrow Interneuron') || strcmp(cell1type,'Wide Interneuron')
            cell1type = 'i';
        end
        
        if strcmp(cell2type,'Pyramidal Cell')
            cell2type = 'p';
        elseif strcmp(cell2type,'Narrow Interneuron') || strcmp(cell2type,'Wide Interneuron')
            cell2type = 'i';
        end
                        
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);


        pairCCGscreen.pair(i,:) = [cellID1,cellID2];
        pairCCGscreen.GSPExc{i} = GSPExc;
        pairCCGscreen.GSPInh{i} = GSPInh;

        if plotFlag

            pairStr = ['AG II-pairs Day' dayNo ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ', region: ' region1 ')' ' v ' ...
                                        num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ', region: ' region2 ')'];

            title(pairStr)

            ylims = get(gca,'ylim');

            if any(GSPExc)
                hold on;
                plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'k^');
            end
            if any(GSPInh)
                hold on;
                plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'kv');
            end
                        
            save_file = fullfile(figpath, ['conv screening only - ' pairStr]);
            print(fig_use, save_file,'-djpeg',resolution_use);

            close all

            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
        end

    end
    
end