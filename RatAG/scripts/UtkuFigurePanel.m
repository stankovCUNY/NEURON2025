function utkuFigurePanel(pairAna,dayNo,plot_pointwiseBands,plotType,flagFlipOrder)

    UtkuData      = loadUtkuV2(dayNo); % NOTE: spike time data is in minutes!!!
    pairType      = 'exquisite';
    filteredPairs = UtkuScreenJitterV2(dayNo,pairType);

    %%

    jscale         = 5;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    fpass          = 300;
%     binSize        = 1/5000;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    Nwaveforms     = 100;
    filterFlag     = true;
    
    channelLayout = {[1:32],[33:64],[65:96],[97,128],[129:160],[161:192]}; 
    
    preLength  = 36;
    postLength = 36;
    
    % organizing and plotting variabls
    waveTimeShort   = UtkuData.cell_metrics.waveforms.time{1,1};     % 48 time points
    waveTimeLong    = UtkuData.cell_metrics.waveforms.time_all{1,1}; % 72 time points
    electrodeGroups = UtkuData.cell_metrics.general.electrodeGroups; 

    figpath     = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    UtkuPath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    day = str2double(dayNo);
    if day == 1
        datapath = [UtkuPath 'AG_2019-12-23_NSD' '/'];
        RawDataPath = ['/media/nasko/WD_BLACK/UtkuRawMultiArray/AG/' 'AG_2019-12_23_NSD' '/'];
    elseif day == 2
        datapath = [UtkuPath 'AG_2019-12-27_NSD' '/'];
        RawDataPath = ['/media/nasko/WD_BLACK/UtkuRawMultiArray/AG/' 'AG_2019-12_27_NSD' '/'];
    end 

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
                                   
    res1    = UtkuData.cell_metrics.spikes.times{1,pairAna(1)};
    res2    = UtkuData.cell_metrics.spikes.times{1,pairAna(2)};

    cellID1 = UtkuData.cell_metrics.UID(pairAna(1));
    cellID2 = UtkuData.cell_metrics.UID(pairAna(2));

    ch1     = UtkuData.cell_metrics.maxWaveformCh1(pairAna(1));
    ch2     = UtkuData.cell_metrics.maxWaveformCh1(pairAna(2));

    clu1    = UtkuData.cell_metrics.cluID(pairAna(1));
    clu2    = UtkuData.cell_metrics.cluID(pairAna(2));

    shank1  = UtkuData.cell_metrics.shankID(pairAna(1));
    shank2  = UtkuData.cell_metrics.shankID(pairAna(2));

    cell1type = UtkuData.cell_metrics.putativeCellType{1,pairAna(1)};
    cell2type = UtkuData.cell_metrics.putativeCellType{1,pairAna(2)};

    region1 = UtkuData.cell_metrics.brainRegion{1,pairAna(1)}; 
    region2 = UtkuData.cell_metrics.brainRegion{1,pairAna(2)};

    cell1waveform = UtkuData.cell_metrics.waveforms.filt_all{1,pairAna(1)};
    cell2waveform = UtkuData.cell_metrics.waveforms.filt_all{1,pairAna(2)};

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
    
    if strcmp(plotType,'ccgOnly')
        
        if flagFlipOrder
           
           temp1 = res1;
           temp2 = res2;
           
           res1 = temp2;
           res2 = temp1;
           
           clear temp1 temp2
           
        end 
        
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', true, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,'plot_pointwiseBands',plot_pointwiseBands); 
        
        ylims = get(gca,'ylim');
                    
        if any(GSPExc)
            hold on;
            plot(tR(GSPExc == 1)*1000, 0.95*ylims(2),'b^');
        end

        if any(GSPInh)
            hold on;
            plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'rv');
        end

        ylim(ylims)
        xlim([-1 1])
%         title(pairStr)
    
    elseif strcmp(plotType,'ccgRefWave')
        
        if flagFlipOrder
           
           temp1 = res1;
           temp2 = res2;
           
           res1 = temp2;
           res2 = temp1;
           
           clear temp1 temp2
           
        end 
        
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', true, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,'plot_pointwiseBands',plot_pointwiseBands); 
        
        ylims = get(gca,'ylim');
                    
        if any(GSPExc)
            hold on;
            plot(tR(GSPExc == 1)*1000, 0.95*ylims(2),'b^');
        end

        if any(GSPInh)
            hold on;
            plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'rv');
        end

        ylim(ylims)
        xlim([-1 1])
%         title(pairStr)
        
        % removing the synchronous spikes                
        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(211,1));

        res1sync = [];
        res2sync = [];

        for k = 1:211

            if (GSPExc(k) == 1) || (GSPInh(k) == 1)
                res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
                res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
            end

        end

        % remove sync spikes for waveforms
        res1nosync = setdiff(res1,res1sync);
        res2nosync = setdiff(res2,res2sync);

        if Nwaveforms < length(res1nosync)
            res1nosync = res1nosync(randperm(length(res1nosync),Nwaveforms));
        end
        if Nwaveforms < length(res2nosync)
            res2nosync = res2nosync(randperm(length(res2nosync),Nwaveforms));
        end

        spikeTimeIndxCell1 = round((res1nosync)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2nosync)*fs) + 1;                
        
    %             pairStr = ['AG' dayNo ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ', region: ' region1 ')' ' v ' ...
    %                                         num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ', region: ' region2 ')'];

        pairStr = [num2str(cellID1) cell1type ' v ' ...
                   num2str(cellID2) cell2type];

        title(pairStr)

        nexttile(3)
%             plot(waveTimeLong,cell1waveform(waveform1Idx,:)/1000,'LineWidth',2,'Color','#0072BD')
        chanDataAtCell1 = load([RawDataPath '/shank' num2str(shank1) '/rawDataCh' num2str(ch1) '.mat' ]);
        [chTimeSeriesCell1, chWaveformsCell1] = waveformAvg(double(chanDataAtCell1.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);

        chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(ch2) '.mat' ]);
        [chTimeSeriesCell2, chWaveformsCell2] = waveformAvg(double(chanDataAtCell2.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);

        plot(waveTimeLong,chTimeSeriesCell1,'LineWidth',3,'Color','#0072BD')
        hold on
        plot(waveTimeLong,chTimeSeriesCell2,'--','LineWidth',3,'Color','#0072BD')
        hold off

        xlim([-1 1])
        legend(['Ref e-spike at ref max chan (no sync spikes, ' num2str(Nwaveforms) ' random spikes)'], ...
               ['Ref e-spike at tar max chan (no sync spikes, ' num2str(Nwaveforms) ' random spikes)'])
        xlabel('Time [ms]')
        ylabel('voltage [mV]')
        title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ', region: ' region1 ') average'])
        
    else
    
        tiledlayout(9,1, 'Padding', 'none', 'TileSpacing', 'compact'); 

    %         nexttile([1 ,2])
        nexttile

        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);

        % removing the synchronous spikes                
        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(211,1));

        res1sync = [];
        res2sync = [];

        for k = 1:211

            if (GSPExc(k) == 1) || (GSPInh(k) == 1)
                res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
                res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
            end

        end

        % remove sync spikes for waveforms
        res1nosync = setdiff(res1,res1sync);
        res2nosync = setdiff(res2,res2sync);

        if Nwaveforms < length(res1nosync)
            res1nosync = res1nosync(randperm(length(res1nosync),Nwaveforms));
        end
        if Nwaveforms < length(res2nosync)
            res2nosync = res2nosync(randperm(length(res2nosync),Nwaveforms));
        end

        spikeTimeIndxCell1 = round((res1nosync)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2nosync)*fs) + 1;                

        if plotFlag

    %             pairStr = ['AG' dayNo ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ', region: ' region1 ')' ' v ' ...
    %                                         num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ', region: ' region2 ')'];

            pairStr = [num2str(cellID1) cell1type ' v ' ...
                       num2str(cellID2) cell2type];

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
            xlim([-(duration/2)*1000 (duration/2)*1000])

            nexttile
    %             plot(waveTimeLong,cell1waveform(waveform1Idx,:)/1000,'LineWidth',2,'Color','#0072BD')
            chanDataAtCell1 = load([RawDataPath '/shank' num2str(shank1) '/rawDataCh' num2str(ch1) '.mat' ]);
            [chTimeSeriesCell1, chWaveformsCell1] = waveformAvg(double(chanDataAtCell1.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);

            chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(ch2) '.mat' ]);
            [chTimeSeriesCell2, chWaveformsCell2] = waveformAvg(double(chanDataAtCell2.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);

            plot(waveTimeLong,chTimeSeriesCell1,'LineWidth',3,'Color','#0072BD')
            hold on
            plot(waveTimeLong,chTimeSeriesCell2,'--','LineWidth',3,'Color','#0072BD')
            hold off

            set(gca, 'YDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            legend(['Ref e-spike at ref max chan (no sync spikes, ' num2str(Nwaveforms) ' random spikes)'], ...
                   ['Ref e-spike at tar max chan (no sync spikes, ' num2str(Nwaveforms) ' random spikes)'])
            xlabel('Time [ms]')
            ylabel('revere order voltage [mV]')
            title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ', region: ' region1 ') average'])

            nexttile

            patchline(waveTimeLong,chWaveformsCell1(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
            hold on
            for k = 2:size(chWaveformsCell1,2)
                patchline(waveTimeLong,chWaveformsCell1(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
            end
            hold off

            set(gca, 'YDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            xlabel('Time [ms]')
            ylabel('revere order voltage [mV]')
            title(['Reference e-spike: ' num2str(cellID1) cell1type ' at ref max channel (' num2str(Nwaveforms) ' random spikes)'])

            nexttile
            if shank1 == shank2
    %                 hold on
    %                 plot(waveTimeLong,cell1waveform(waveform2Idx,:)/1000,'--','LineWidth',2,'Color','#0072BD')
    %                 hold off

                  for k = 1:32
    %                     chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(ch2) '.mat' ]);
                      chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(channelLayout{shank2}(k)) '.mat' ]);
                      [chTimeSeries(:,k), chWaveforms{k}] = waveformAvg(double(chanDataAtCell2.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
                  end
    %                 plot(waveTimeLong,chTimeSeries,'--','LineWidth',3,'Color','#0072BD')

                hold on
                for k = 1:size(chWaveformsCell1,2)
                    patchline(waveTimeLong,chWaveforms{find(channelLayout{shank2} == ch2)}(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
                end
                hold off
            elseif (shank1 ~= shank2)

                  for k = 1:32
    %                   chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(ch2) '.mat' ]);
                      chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(channelLayout{shank2}(k)) '.mat' ]);
                      [chTimeSeries(:,k), chWaveforms{k}] = waveformAvg(double(chanDataAtCell2.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
                  end
    %                 plot(waveTimeLong,chTimeSeries,'--','LineWidth',3,'Color','#0072BD')

                hold on
                for k = 1:size(chWaveformsCell1,2)
                    patchline(waveTimeLong,chWaveforms{find(channelLayout{shank2} == ch2)}(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
                end
                hold off
            end
            set(gca, 'YDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            xlabel('Time [ms]')
            ylabel('revere order voltage [mV]')
            title(['Reference e-spike: ' num2str(cellID1) cell1type ' at tar max channel (' num2str(Nwaveforms) ' random spikes)'])

            nexttile
            plot(waveTimeLong,chTimeSeries)
            set(gca, 'YDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            xlabel('Time [ms]')
            ylabel('revere order voltage [mV]')
            title(['Reference e-spike: ' num2str(cellID1) cell1type ' on all ' num2str(cellID2) cell2type '''s shank ' num2str(shank2) ' channels' ])
            legend(cellstr(num2str(channelLayout{shank2}', 'ch%-d')),'Location','northeast','NumColumns',4)

            nexttile
    %             plot(waveTimeLong,cell2waveform(waveform2Idx,:)/1000,'LineWidth',2,'Color','#D95319')

            chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/rawDataCh' num2str(ch2) '.mat' ]);
            [chTimeSeriesCell2, chWaveformsCell2] = waveformAvg(double(chanDataAtCell2.data),spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);

            chanDataAtCell1 = load([RawDataPath '/shank' num2str(shank1) '/rawDataCh' num2str(ch1) '.mat' ]);
            [chTimeSeriesCell1, chWaveformsCell1] = waveformAvg(double(chanDataAtCell1.data),spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);

            plot(waveTimeLong,chTimeSeriesCell2,'LineWidth',3,'Color','#D95319')
            hold on
            plot(waveTimeLong,chTimeSeriesCell1,'--','LineWidth',3,'Color','#D95319')
            hold off

            set(gca, 'YDir','reverse','XDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            legend(['Tar e-spike at tar max chan (no sync spikes, ' num2str(Nwaveforms) ' random spikes)'],...
                   ['Tar e-spike at ref max chan (no sync spikes, ' num2str(Nwaveforms) ' random spikes)'])
            xlabel('Time revere order [ms]')
            ylabel('revere order voltage [mV]')
            title(['Target e-spike: ' num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ', region: ' region2 ') average'])

            nexttile
            patchline(waveTimeLong,chWaveformsCell2(:,1),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
            hold on
            for k = 2:size(chWaveformsCell2,2)
                patchline(waveTimeLong,chWaveformsCell2(:,k),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
            end
            hold off

            set(gca, 'YDir','reverse','XDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            xlabel('Time revere order [ms]')
            ylabel('revere order voltage [mV]')
            title(['Target e-spike: ' num2str(cellID2) cell2type ' at tar max channel (' num2str(Nwaveforms) ' random spikes)'])

            nexttile
            if shank1 == shank2
    %                 hold on
    %                 plot(waveTimeLong,cell2waveform(waveform1Idx,:)/1000,'--','LineWidth',2,'Color','#D95319')
    %                 hold off

                for k = 1:32
                      chanDataAtCell1 = load([RawDataPath '/shank' num2str(shank1) '/rawDataCh' num2str(channelLayout{shank1}(k)) '.mat' ]);
                      [chTimeSeries(:,k), chWaveforms{k}] = waveformAvg(double(chanDataAtCell1.data),spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);
                end

                hold on
                for k = 1:size(chWaveformsCell2,2)
                    patchline(waveTimeLong,chWaveforms{find(channelLayout{shank1} == ch1)}(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
                end
                hold off
    %                 patchline(waveTimeLong,chWaveformsCell1(:,1),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
    %                 hold on
    %                 for k = 2:Nwaveforms
    %                     patchline(waveTimeLong,chWaveformsCell1(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
    %                 end
    %                 hold off    
            elseif (shank1 ~= shank2)

                for k = 1:32
                      chanDataAtCell1 = load([RawDataPath '/shank' num2str(shank1) '/rawDataCh' num2str(channelLayout{shank1}(k)) '.mat' ]);
                      [chTimeSeries(:,k), chWaveforms{k}] = waveformAvg(double(chanDataAtCell1.data),spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);
                end

                hold on
                for k = 1:size(chWaveformsCell2,2)
                    patchline(waveTimeLong,chWaveforms{find(channelLayout{shank1} == ch1)}(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
                end
                hold off
    %                 patchline(waveTimeLong,chWaveformsCell1(:,1),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
    %                 hold on
    %                 for k = 2:Nwaveforms
    %                     patchline(waveTimeLong,chWaveformsCell1(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
    %                 end
    %                 hold off                
            end
            set(gca, 'YDir','reverse','XDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            xlabel('Time revere order [ms]')
            ylabel('revere order voltage [mV]')
            title(['Target e-spike: ' num2str(cellID2) cell2type ' at ref max channel (' num2str(Nwaveforms) ' random spikes)'])

            nexttile
            plot(waveTimeLong,chTimeSeries)
            set(gca, 'YDir','reverse','XDir','reverse')
            xlim([-(duration/2)*1000 (duration/2)*1000])
            xlabel('Time revere order [ms]')
            ylabel('revere order voltage [mV]')
            title(['Target e-spike: ' num2str(cellID2) cell2type ' on all ' num2str(cellID1) cell1type '''s shank ' num2str(shank1) ' channels' ])
            legend(cellstr(num2str(channelLayout{shank1}', 'ch%-d')),'Location','northeast','NumColumns',4)


            save_file = fullfile(figpath, ['waveform ana - ' pairStr]);
            print(fig_use, save_file,'-djpeg',resolution_use);

            close all

            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

            toc
        end
    end
end