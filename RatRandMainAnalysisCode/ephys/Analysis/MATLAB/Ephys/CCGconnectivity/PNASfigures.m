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

figureCode = '5';

%% Figure 0 (graphical abstract)

if figureCode == '0'
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    set(gcf, 'Renderer', 'painters', 'Position', [1720 2562 2*144 256]);

    res_type = 'QHD';
    
    tiledlayout(2,1,'Padding','none','TileSpacing','tight')

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable

    % idx 23
    ref = pairGroupStatsTable.tarSpikeTimes{23,1};
    tar = pairGroupStatsTable.refSpikeTimes{23,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{23,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{23,1};
    d = pairGroupStatsTable.pairDistance(23);
    region = pairGroupStatsTable.brainRegion{23,1};

    [ccgR, tR] = CCG([ref;tar],[ones(size(ref));2*ones(size(tar))], ...
                'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                'norm', 'counts');
    
    nexttile(2)
    line(tR*1e3,ccgR(:,1,2)/length(ref),'color','k', 'LineWidth',1)

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    set(gca,'FontSize',5)
    set(gcf, 'Renderer', 'Painters');

    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    nexttile(1) 
    plot(pairGroupStatsTable.tWave{23,1}, ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),:) - ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),1),'LineWidth',1,'color','#0072BD')

    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

end

%% Figure 1 exquisite data
if figureCode == '1'
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 (3/5)*1920*0.3 1080*0.4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    tiledlayout(3,3,'Padding','none','TileSpacing','compact')
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable

    % idx 13
    ref = pairGroupStatsTable.refSpikeTimes{13,1};
    tar = pairGroupStatsTable.tarSpikeTimes{13,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{13,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{13,1};
    d = pairGroupStatsTable.pairDistance(13);
    region = pairGroupStatsTable.brainRegion{13,1};
    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % idx 18
    ref = pairGroupStatsTable.refSpikeTimes{18,1};
    tar = pairGroupStatsTable.tarSpikeTimes{18,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{18,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{18,1};
    d = pairGroupStatsTable.pairDistance(18);
    region = pairGroupStatsTable.brainRegion{18,1};
    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false,...
                                    'norm_flag', true);


    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % idx 2
    ref = pairGroupStatsTable.refSpikeTimes{2,1};
    tar = pairGroupStatsTable.tarSpikeTimes{2,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{2,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{2,1};
    d = pairGroupStatsTable.pairDistance(2);
    region = pairGroupStatsTable.brainRegion{2,1};
    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/groupStatsRatAGday2.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable

    % idx 128
    ref = pairGroupStatsTable.refSpikeTimes{128,1};
    tar = pairGroupStatsTable.tarSpikeTimes{128,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{128,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{128,1};
    d = pairGroupStatsTable.pairDistance(128);
    region = pairGroupStatsTable.refBrainRegion(128,:);
    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    title({'Rat AG';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse794812542.mat')
    % 
    % % idx 28
    % ref = pairGroupStatsTable.refSpikeTimes{28,1};
    % tar = pairGroupStatsTable.tarSpikeTimes{28,1}; 
    % nexttile
    % [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
    %           CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
    %                                 'plot_output', get(fig_use, 'Number'), ...
    %                                 'njitter', njitter, 'alpha', alpha,...
    %                                 'for_grant', for_grant, 'plot_pointwiseBands', false, ...
    %                                 'lineWidth', 4, 'axisLabels', false);
    %                             
    % 
    % xlim([-1,1])
    % ylims = get(gca,'ylim');
    % if any(GSPExc)
    %     hold on;
    %     plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
    % end
    % if any(GSPInh)
    %     hold on;
    %     plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
    % end
    % set(gca,'FontSize',12)
    % 
    % % idx 125
    % load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse839068429.mat')
    % 
    % ref = pairGroupStatsTable.refSpikeTimes{125,1};
    % tar = pairGroupStatsTable.tarSpikeTimes{125,1}; 
    % nexttile
    % [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
    %           CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
    %                                 'plot_output', get(fig_use, 'Number'), ...
    %                                 'njitter', njitter, 'alpha', alpha,...
    %                                 'for_grant', for_grant, 'plot_pointwiseBands', false, ...
    %                                 'lineWidth', 4, 'axisLabels', false);
    %                             
    % xlim([-1,1])
    % ylims = get(gca,'ylim');
    % if any(GSPExc)
    %     hold on;
    %     plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
    % end
    % if any(GSPInh)
    %     hold on;
    %     plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
    % end
    % set(gca,'FontSize',12)

%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse762120172.mat')
    load('/media/nasko/WD_BLACK31/BOTtemp/groupStatsMouse762120172.mat')

    % idx 41
    ref = pairGroupStatsTable.refSpikeTimes{41,1};
    tar = pairGroupStatsTable.tarSpikeTimes{41,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{41,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{41,1};
    d = pairGroupStatsTable.pairDistance(41);
    region = pairGroupStatsTable.brainRegion{41,1};
    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
    % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
%     title(['i' '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    title({'Mouse 762120172';['undef. ' region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse762602078.mat')
    load('/media/nasko/WD_BLACK31/BOTtemp/groupStatsMouse762602078.mat')
    
    % idx 41
    ref = pairGroupStatsTable.refSpikeTimes{41,1};
    tar = pairGroupStatsTable.tarSpikeTimes{41,1};
    refCelltype = pairGroupStatsTable.refCellExplorerType{41,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{41,1};
    d = pairGroupStatsTable.pairDistance(41);
    region = pairGroupStatsTable.brainRegion{41,1};
    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
    % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
%     title([refCelltype '-' 'i' ' pair ' num2str(d) '\mum ' region])
    title({'Mouse 762602078 ';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse768515987.mat')
    load('/media/nasko/WD_BLACK31/BOTtemp/groupStatsMouse768515987')
    
    % idx 58
    ref = pairGroupStatsTable.refSpikeTimes{58,1};
    tar = pairGroupStatsTable.tarSpikeTimes{58,1};
    refCelltype = pairGroupStatsTable.refCellExplorerType{58,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{58,1};
    d = pairGroupStatsTable.pairDistance(58);
    region = pairGroupStatsTable.brainRegion{58,1};

    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
    % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
%     title(['i' '-' 'i' ' pair ' num2str(d) '\mum ' region])
    title({'Mouse 768515987 ';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % idx 235
    ref = pairGroupStatsTable.refSpikeTimes{235,1};
    tar = pairGroupStatsTable.tarSpikeTimes{235,1};
    refCelltype = pairGroupStatsTable.refCellExplorerType{235,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{235,1};
    d = pairGroupStatsTable.pairDistance(235);
    region = pairGroupStatsTable.brainRegion{235,1};

    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
    % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
%     title(['p' '-' 'p' ' pair ' num2str(d) '\mum ' region])
    title({'Mouse 768515987';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse774875821.mat')
    load('/media/nasko/WD_BLACK31/BOTtemp/groupStatsMouse774875821.mat')
    
    % idx 102
    ref = pairGroupStatsTable.refSpikeTimes{102,1};
    tar = pairGroupStatsTable.tarSpikeTimes{102,1};
    refCelltype = pairGroupStatsTable.refCellExplorerType{102,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{102,1};
    d = pairGroupStatsTable.pairDistance(102);
    region = pairGroupStatsTable.brainRegion{102,1};

    nexttile
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);

    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    box off
    % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
%     title(['p' '-' 'p' ' pair ' num2str(d) '\mum ' region])
    title({'Mouse 774875821';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    
%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/groupStatsMouse793224716.mat')
    load('/media/nasko/WD_BLACK31/BOTtemp/groupStatsMouse793224716.mat')
    
%     % idx 30
%     ref = pairGroupStatsTable.refSpikeTimes{30,1};
%     tar = pairGroupStatsTable.tarSpikeTimes{30,1};
%     refCelltype = pairGroupStatsTable.refCellExplorerType{30,1};
%     tarCelltype = pairGroupStatsTable.tarCellExplorerType{30,1};
%     d = pairGroupStatsTable.pairDistance(30);
%     region = pairGroupStatsTable.brainRegion{30,1};
% 
%     nexttile
%    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%               CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
%                                     'plot_output', get(fig_use, 'Number'), ...
%                                     'njitter', njitter, 'alpha', alpha,...
%                                     'for_grant', for_grant, 'plot_pointwiseBands', false, ...
%                                     'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
%                                     'norm_flag', true);
% 
%     xlim([-1,1])
%     ylims = get(gca,'ylim');
%     ylabel('Spike Probability')
%     xlabel('[ms]')
%     box off
%     % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
% %     title(['p' '-' 'p' ' pair ' num2str(d) '\mum ' region])
% 
% %     if any(GSPExc)
% %         hold on;
% %         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
% %     end
% %     if any(GSPInh)
% %         hold on;
% %         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
% %     end
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
    
%     tiledlayout(1,2)
    
%     nexttile
    RoySummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsV2.mat','fracExq');
    AGsummaryStats  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsV2.mat','fracExq','brainRegions','screenPairType');
    NSUsummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsV2.mat','fracExq');
%     
%     bar(categorical(AGsummaryStats.brainRegions), ...
%         100*[[RoySummaryStats.fracExq{1,1} + AGsummaryStats.fracExq{1,1}(1) + NSUsummaryStats.fracExq{1,1}]/6 ...
%              [RoySummaryStats.fracExq{2,1} + AGsummaryStats.fracExq{2,1}(1) + NSUsummaryStats.fracExq{2,1}]/6 ...
%              [RoySummaryStats.fracExq{4,1} + AGsummaryStats.fracExq{4,1}(1) + NSUsummaryStats.fracExq{4,1}]/6; ...
%              [0 +                          AGsummaryStats.fracExq{1,1}(2) + 0]/2 ...
%              [0 +                          AGsummaryStats.fracExq{2,1}(2) + 0]/2 ...
%              [0 +                          AGsummaryStats.fracExq{4,1}(2) + 0]/2])
%     legend(AGsummaryStats.screenPairType{1}, AGsummaryStats.screenPairType{2}, ...
%         [AGsummaryStats.screenPairType{3} ' or ' AGsummaryStats.screenPairType{4}])
%     
%     nexttile
    BOTsummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStats.mat','fracExq','brainRegions','screenPairType');
    
     brainRegions  = ["CA1"   
                      "CA3" 
                      "APN"   
                      "CA1"   
                      "CA3"   
                      "DG"    
                      "Eth"   
                      "HPF"   
                      "IntG"  
                      "LGd"   
                      "LGv"   
                      "LP"    
                      "MB"    
                      "MGd"   
                      "MGm"   
                      "MGv"   
                      "MRN"   
                      "NOT"   
                      "PO"    
                      "POL"   
                      "ProS"  
                      "SCig"  
                      "SGN"   
                      "SUB"   
                      "TH"    
                      "VIS"   
                      "VISal" 
                      "VISam" 
                      "VISl"  
                      "VISli" 
                      "VISmma"
                      "VISmmp"
                      "VISp"  
                      "VISpm" 
                      "VISrl" 
                      "VL"    
                      "VPM"   
                      "undef. grey"  ];
    
    close 102
           
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920*0.04 (2/3)*1080/5]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    
%     tiledlayout(1,6,'Padding','none','TileSpacing','compact')
%     nexttile([1,4])
    
    t = tiledlayout(2,1,'Padding','none','TileSpacing','compact');
    bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
    bgAx.Layout.TileSpan = [2 1];
    
    ax1 = axes(t);
    b = bar(100*[[RoySummaryStats.fracExq{1,1} + AGsummaryStats.fracExq{1,1}(1) + NSUsummaryStats.fracExq{1,1}]/6 ...
                 [RoySummaryStats.fracExq{2,1} + AGsummaryStats.fracExq{2,1}(1) + NSUsummaryStats.fracExq{2,1}]/6 ...
                 [RoySummaryStats.fracExq{4,1} + AGsummaryStats.fracExq{4,1}(1) + NSUsummaryStats.fracExq{4,1}]/6; ...
                 [0 +                            AGsummaryStats.fracExq{1,1}(2) + 0]/2 ...
                 [0 +                            AGsummaryStats.fracExq{2,1}(2) + 0]/2 ...
                 [0 +                            AGsummaryStats.fracExq{4,1}(2) + 0]/2]);
    
    b(1).FaceColor = [0.8500 0.3250 0.0980];
    b(2).FaceColor = [0      0.4470 0.7410];
    
%     ylabel('Percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTick',[])
    ax1.Box = 'off';
    ax1.XAxis.Visible = 'off';
    ylim([7.5 8.5])
    yline(ax1,7.5);
    
%     legend(AGsummaryStats.screenPairType{1}, AGsummaryStats.screenPairType{2}, ...
%           [AGsummaryStats.screenPairType{3} ' or ' AGsummaryStats.screenPairType{4}])
%     set(ax1,'YTick',[]);
    
    ax2 = axes(t);
    ax2.Layout.Tile = 2;
    b = bar(categorical(brainRegions(1:2)),100*[[RoySummaryStats.fracExq{1,1} + AGsummaryStats.fracExq{1,1}(1) + NSUsummaryStats.fracExq{1,1}]/6 ...
                                         [RoySummaryStats.fracExq{2,1} + AGsummaryStats.fracExq{2,1}(1) + NSUsummaryStats.fracExq{2,1}]/6 ...
                                         [RoySummaryStats.fracExq{4,1} + AGsummaryStats.fracExq{4,1}(1) + NSUsummaryStats.fracExq{4,1}]/6; ...
                                         [0 +                            AGsummaryStats.fracExq{1,1}(2) + 0]/2 ...
                                         [0 +                            AGsummaryStats.fracExq{2,1}(2) + 0]/2 ...
                                         [0 +                            AGsummaryStats.fracExq{4,1}(2) + 0]/2]);
    b(1).FaceColor = [0.8500 0.3250 0.0980];
    b(2).FaceColor = [0      0.4470 0.7410];

%     ylabel('Percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    ylim([0 1])
    ax2.Box = 'off';
    ax1.XAxis.Visible = 'off';
    yline(ax2,1);
    
    close 102
    
    hcomb = figure(102);

%     tiledlayout(1,7,'Padding','none','TileSpacing','compact')
%     nexttile([1,4])
    
    res_type = 'QHD';
    pos = [70 230 1920*0.3 (2/3)*1080/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    maskNoPairs = find(~(((BOTsummaryStats.fracExq{1,1} < 0.05) & (BOTsummaryStats.fracExq{2,1} < 0.05)) & (BOTsummaryStats.fracExq{3,1} < 0.05)));
    [~,tempIdx] = sort(BOTsummaryStats.fracExq{1,1}(maskNoPairs),'descend');
    maskNoPairs = maskNoPairs(tempIdx);
    
    brainRegionsMice = brainRegions(3:end);
    brainRegionsMice = categorical(brainRegionsMice(maskNoPairs));
    brainRegionsMice = reordercats(brainRegionsMice,cellstr(brainRegionsMice)');
    
    b = bar(brainRegionsMice,100*[BOTsummaryStats.fracExq{1,1}(maskNoPairs)/58 BOTsummaryStats.fracExq{2,1}(maskNoPairs)/58 BOTsummaryStats.fracExq{4,1}(maskNoPairs)/58]);
    legend(AGsummaryStats.screenPairType{1}, AGsummaryStats.screenPairType{2}, ...
    [AGsummaryStats.screenPairType{3} ' or ' AGsummaryStats.screenPairType{4}])
    b(1).FaceColor = [0.8500 0.3250 0.0980];
    b(2).FaceColor = [0      0.4470 0.7410];
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
%     t = tiledlayout(2,1,'TileSpacing','compact');
%     bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
%     bgAx.Layout.TileSpan = [2 1];
%     
%     ax1 = axes(t);
%     b = bar(100*[BOTsummaryStats.fracExq{1,1}/58 BOTsummaryStats.fracExq{2,1}/58 BOTsummaryStats.fracExq{4,1}/58]);
%     b(1).FaceColor = [0.8500 0.3250 0.0980];
%     b(2).FaceColor = [0      0.4470 0.7410];
%     ylabel('percent')
%     set(gca,'FontSize',5)
%     set(gca,'XTick',[])
%     ax1.Box = 'off';
%     ax1.XAxis.Visible = 'off';
%     ylim([0.5 2])
%     yline(ax1,0.5);
%     legend(AGsummaryStats.screenPairType{1}, AGsummaryStats.screenPairType{2}, ...
%           [AGsummaryStats.screenPairType{3} ' or ' AGsummaryStats.screenPairType{4}])
%     set(ax1,'YTick',[]);
%     
%     ax2 = axes(t);
%     ax2.Layout.Tile = 2;
%     b = bar(categorical(brainRegions(3:end)),100*[BOTsummaryStats.fracExq{1,1}/58 BOTsummaryStats.fracExq{2,1}/58 BOTsummaryStats.fracExq{4,1}/58]);
%     b(1).FaceColor = [0.8500 0.3250 0.0980];
%     b(2).FaceColor = [0      0.4470 0.7410];
%     ylabel('percent')
%     set(gca,'FontSize',5)
%     ylim([0 0.5])
%     ax2.Box = 'off';
%     ax1.XAxis.Visible = 'off';
%     yline(ax2,1);
    
    %%
    
    animalsList = [715093703
                   719161530
                   721123822
                   732592105
                   737581020
                   739448407
                   742951821
                   743475441
                   744228101
                   746083955
                   750332458
                   750749662
                   751348571
                   754312389
                   754829445
                   755434585
                   756029989
                   757216464
                   757970808
                   758798717
                   759883607
                   760345702
                   760693773
                   761418226
                   762120172
                   762602078
                   763673393
                   766640955
                   767871931
                   768515987
                   771160300
                   771990200
                   773418906
                   774875821
                   778240327
                   778998620
                   779839471
                   781842082
                   786091066
                   787025148
                   789848216
                   791319847
                   793224716
                   794812542
                   797828357
                   798911424
                   799864342
                   816200189
                   819186360
                   819701982
                   821695405
                   829720705
                   831882777
                   835479236
                   839068429
                   839557629
                   840012044
                   847657808];

    datasetID   = 'AllenInstitute';
%     datapath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
    datapath    = '/media/nasko/WD_BLACK31/BOTtemp/';
    
    pairGroupStatTableExq = table;

    for loopAnimals = 1:58

        display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])

        load([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);

        pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
        
        pairGroupStatTableTemp = table;
        pairGroupStatTableTemp.refTroughToPeak     = cell2num(pairGroupStatTable.refTroughToPeakLength);
        pairGroupStatTableTemp.tarTroughToPeak     = cell2num(pairGroupStatTable.tarTroughToPeakLength);
        pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
        pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
        pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
        pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
        pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
        
    end
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStats.mat')
    pairGroupStatTable     = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = table;
    pairGroupStatTableTemp.refTroughToPeak     = pairGroupStatTable.refTroughToPeak;
    pairGroupStatTableTemp.tarTroughToPeak     = pairGroupStatTable.tarTroughToPeak;
    pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
    pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
    pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
    pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
    pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat')
    pairGroupStatTable     = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = table;
    pairGroupStatTableTemp.refTroughToPeak     = pairGroupStatTable.refTroughToPeak;
    pairGroupStatTableTemp.tarTroughToPeak     = pairGroupStatTable.tarTroughToPeak;
    pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
    pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
    pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
    pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
    pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStats.mat')
    pairGroupStatTable     = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = table;
    pairGroupStatTableTemp.refTroughToPeak     = pairGroupStatTable.refTroughToPeak;
    pairGroupStatTableTemp.tarTroughToPeak     = pairGroupStatTable.tarTroughToPeak;
    pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
    pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
    pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
    pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
    pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStats.mat')
    pairGroupStatTable     = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = table;
    pairGroupStatTableTemp.refTroughToPeak     = pairGroupStatTable.refTroughToPeak;
    pairGroupStatTableTemp.tarTroughToPeak     = pairGroupStatTable.tarTroughToPeak;
    pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
    pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
    pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
    pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
    pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')
    pairGroupStatTable     = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = table;
    pairGroupStatTableTemp.refTroughToPeak     = pairGroupStatTable.refTroughToPeak;
    pairGroupStatTableTemp.tarTroughToPeak     = pairGroupStatTable.tarTroughToPeak;
    pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
    pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
    pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
    pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
    pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStats.mat')
    pairGroupStatTable     = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = table;
    pairGroupStatTableTemp.refTroughToPeak     = pairGroupStatTable.refTroughToPeak;
    pairGroupStatTableTemp.tarTroughToPeak     = pairGroupStatTable.tarTroughToPeak;
    pairGroupStatTableTemp.refACGtauRise       = pairGroupStatTable.refACGtauRise;
    pairGroupStatTableTemp.tarACGtauRise       = pairGroupStatTable.tarACGtauRise;
    pairGroupStatTableTemp.refCellExplorerType = pairGroupStatTable.refCellExplorerType;
    pairGroupStatTableTemp.tarCellExplorerType = pairGroupStatTable.tarCellExplorerType;
    pairGroupStatTableExq                      = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
%     nexttile([1,3])

    hcomb = figure(102);
    
    res_type = 'QHD';
    pos = [70 230 1920*0.125 0.6*1080/3]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    plot([pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'p'));       ...
               pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'p'));        ...
               pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-wide'));  ...
               pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-wide'))], ... 
              [pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'p'));       ...
               pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'p'));        ...
               pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-wide'));  ...
               pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-wide'))], ...
              '^','LineWidth',1)
    hold on
    plot([pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'i'));       ...
               pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i'));        ...
               pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-narrow'));  ...
               pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-narrow'))], ... 
              [pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'i'));       ...
               pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i'));        ...
               pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-narrow'));  ...
               pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-narrow'))], ...
              'o','LineWidth',1)
    for loopPairs = 1:size(pairGroupStatTableExq,1)
        patchline([pairGroupStatTableExq.refTroughToPeak(loopPairs),pairGroupStatTableExq.tarTroughToPeak(loopPairs)], ... 
                   [pairGroupStatTableExq.refACGtauRise(loopPairs),pairGroupStatTableExq.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',0.1)
    end
    hold off
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
    legend('p','i','Location','NorthWest')
    ylim([0.05 10])
%     xlim([0.1 3])
    xlim([0.1 0.8])
    ylabel('ACG \tau_r_i_s_e')
    xlabel('Trough-to-peak [ms]')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    set(gcf, 'Renderer', 'Painters');
    
    hcomb = figure(103);
    
    res_type = 'QHD';
    pos = [70 230 1920*0.125 0.6*1080/3]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    hist3([[pairGroupStatTableExq.refTroughToPeak;   ...
            pairGroupStatTableExq.tarTroughToPeak], ...
           [pairGroupStatTableExq.refACGtauRise;     ...
            pairGroupStatTableExq.tarACGtauRise]],   ...
           'CDataMode','auto','FaceColor','interp','Edges',{0:1/30:1.5 0:1/5:10},'EdgeColor','none')
    view(2)
    
%     set(gca, 'YScale', 'log')
    xlim([0.1 0.8])
    ylim([0.05 10])
    clim([0 50])
    
    ylabel('ACG \tau_r_i_s_e')
    xlabel('Trough-to-peak [ms]')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    set(gcf, 'Renderer', 'Painters');
    
    colormap turbo
    c = colorbar('north','color','w');

    %%
    
    royTable = countPairsHiro;
    AGtable  = countPairsUtku;
    NSUtable = countPairsBapun;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/noPairCompTable.mat');
    BOTtable = noPairCompTable;
    %%
    
    iiTotal       = sum(royTable.nPairCompII) + ...
                    sum(AGtable.nPairCompII)  + ...
                    sum(NSUtable.nPairCompII) + ...
                    sum(BOTtable.nPairCompII);
    ppTotal       = sum(royTable.nPairCompPP) + ...
                    sum(AGtable.nPairCompPP)  + ...
                    sum(NSUtable.nPairCompPP) + ...
                    sum(BOTtable.nPairCompPP);
    ip_piTotal    = sum(royTable.nPairCompIPorPI) + ...
                    sum(AGtable.nPairCompIPorPI)  + ...
                    sum(NSUtable.nPairCompIPorPI) + ...
                    sum(BOTtable.nPairCompIPorPI);
    
    iiExqTotal    = 0;
    ppExqTotal    = 0;
    ip_piExqTotal = 0;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStats.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiExqTotal    = sum(noExq{1,1})    + iiExqTotal;
    ppExqTotal    = sum(noExq{2,1})    + ppExqTotal;
    ip_piExqTotal = 2*sum(noExq{3,1})  + ip_piExqTotal;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStats.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiExqTotal    = sum(noExq{1,1})    + iiExqTotal;
    ppExqTotal    = sum(noExq{2,1})    + ppExqTotal;
    ip_piExqTotal = 2*sum(noExq{3,1})  + ip_piExqTotal;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStats.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiExqTotal    = sum(noExq{1,1})    + iiExqTotal;
    ppExqTotal    = sum(noExq{2,1})    + ppExqTotal;
    ip_piExqTotal = 2*sum(noExq{3,1})  + ip_piExqTotal;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStats.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiExqTotal    = sum(noExq{1,1})    + iiExqTotal;
    ppExqTotal    = sum(noExq{2,1})    + ppExqTotal;
    ip_piExqTotal = 2*sum(noExq{3,1})  + ip_piExqTotal;
    
    (iiExqTotal/iiTotal)*100
    (ppExqTotal/ppTotal)*100
    (ip_piExqTotal/ip_piTotal)*100
    
    % (27/iiTotal)*100
    % (10/ppTotal)*100
    % (8/ip_piTotal)*100
    
    % ratioIiExqGJ    = iiExqTotal/27;
    % ratioPpExqGJ    = ppExqTotal/10;
    % ratioIp_PiExqGJ = ip_piExqTotal/8;
    
    iiGJtotal    = 0;
    ppGJtotal    = 0;
    ip_piGJtotal = 0;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiGJtotal    = sum(noGJ{1,1})    + iiGJtotal;
    ppGJtotal    = sum(noGJ{2,1})    + ppGJtotal;
    ip_piGJtotal = 2*sum(noGJ{3,1})  + ip_piGJtotal;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiGJtotal    = sum(noGJ{1,1})    + iiGJtotal;
    ppGJtotal    = sum(noGJ{2,1})    + ppGJtotal;
    ip_piGJtotal = 2*sum(noGJ{3,1})  + ip_piGJtotal;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiGJtotal    = sum(noGJ{1,1})    + iiGJtotal;
    ppGJtotal    = sum(noGJ{2,1})    + ppGJtotal;
    ip_piGJtotal = 2*sum(noGJ{3,1})  + ip_piGJtotal;
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ.mat')
%     iiTotal       = sum(noPair{1,1})   + iiTotal;
%     ppTotal       = sum(noPair{2,1})   + ppTotal;
%     ip_piTotal    = 2*sum(noPair{3,1}) + ip_piTotal;
    iiGJtotal    = sum(noGJ{1,1})    + iiGJtotal;
    ppGJtotal    = sum(noGJ{2,1})    + ppGJtotal;
    ip_piGJtotal = 2*sum(noGJ{3,1})  + ip_piGJtotal;
    
    (iiGJtotal/iiTotal)*100
    (ppGJtotal/ppTotal)*100
    (ip_piGJtotal/ip_piTotal)*100

    iiGJtotal/iiExqTotal
    ppGJtotal/ppExqTotal
    ip_piGJtotal/ip_piExqTotal

end

%% Figure 2 millisecond data
if figureCode == '2'
    
%    runPetersenAna
    runPetersenGJAna
   
    RoySummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ.mat','fracGJ');
    AGsummaryStats  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ.mat','fracGJ','brainRegions','screenPairType');
    NSUsummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ.mat','fracGJ');

    BOTsummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ.mat','fracGJ','brainRegions','screenPairType');

    brainRegions  = [ "CA1"   
                      "CA3" 
                      "APN"   
                      "CA1"   
                      "CA3"   
                      "DG"    
                      "Eth"   
                      "HPF"   
                      "IntG"  
                      "LGd"   
                      "LGv"   
                      "LP"    
                      "MB"    
                      "MGd"   
                      "MGm"   
                      "MGv"   
                      "MRN"   
                      "NOT"   
                      "PO"    
                      "POL"   
                      "ProS"  
                      "SCig"  
                      "SGN"   
                      "SUB"   
                      "TH"    
                      "VIS"   
                      "VISal" 
                      "VISam" 
                      "VISl"  
                      "VISli" 
                      "VISmma"
                      "VISmmp"
                      "VISp"  
                      "VISpm" 
                      "VISrl" 
                      "VL"    
                      "VPM"   
                      "undef. grey"  ];

    hcomb = figure(104);

    res_type = 'QHD';
    pos = [70 230 1920*0.04 (0.55)*1080/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    b = bar(categorical(brainRegions(1:2)),100*[[RoySummaryStats.fracGJ{1,1} + AGsummaryStats.fracGJ{1,1}(1) + NSUsummaryStats.fracGJ{1,1}]/6 ...
                                                [RoySummaryStats.fracGJ{2,1} + AGsummaryStats.fracGJ{2,1}(1) + NSUsummaryStats.fracGJ{2,1}]/6 ...
                                                [RoySummaryStats.fracGJ{4,1} + AGsummaryStats.fracGJ{4,1}(1) + NSUsummaryStats.fracGJ{4,1}]/6; ...
                                                [0 +                           AGsummaryStats.fracGJ{1,1}(2) + 0]/2 ...
                                                [0 +                           AGsummaryStats.fracGJ{2,1}(2) + 0]/2 ...
                                                [0 +                           AGsummaryStats.fracGJ{4,1}(2) + 0]/2]);
                                            
    b(1).FaceColor = [0.8500 0.3250 0.0980];
    b(2).FaceColor = [0      0.4470 0.7410];

    ylabel('Percent')
    ylim([0 20])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

%     close 104

    hcomb = figure(105);

    %     tiledlayout(1,7,'Padding','none','TileSpacing','compact')
    %     nexttile([1,4])

    res_type = 'QHD';
    pos = [70 230 1920*0.3 (2/3)*1080/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    maskNoPairs = find(~(((BOTsummaryStats.fracGJ{1,1} < 0.05) & (BOTsummaryStats.fracGJ{2,1} < 0.05)) & (BOTsummaryStats.fracGJ{3,1} < 0.05)));
    [~,tempIdx] = sort(BOTsummaryStats.fracGJ{1,1}(maskNoPairs),'descend');
    maskNoPairs = maskNoPairs(tempIdx);

    brainRegionsMice = brainRegions(3:end);
    brainRegionsMice = categorical(brainRegionsMice(maskNoPairs));
    brainRegionsMice = reordercats(brainRegionsMice,cellstr(brainRegionsMice)');

    b = bar(brainRegionsMice,100*[BOTsummaryStats.fracGJ{1,1}(maskNoPairs)/58 BOTsummaryStats.fracGJ{2,1}(maskNoPairs)/58 BOTsummaryStats.fracGJ{4,1}(maskNoPairs)/58]);
    legend(AGsummaryStats.screenPairType{1}, AGsummaryStats.screenPairType{2}, ...
    [AGsummaryStats.screenPairType{3} ' or ' AGsummaryStats.screenPairType{4}])

    b(1).FaceColor = [0.8500 0.3250 0.0980];
    b(2).FaceColor = [0      0.4470 0.7410];
    ylim([0 20])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
   
%     animalsList = [715093703
%                    719161530
%                    721123822
%                    732592105
%                    737581020
%                    739448407
%                    742951821
%                    743475441
%                    744228101
%                    746083955
%                    750332458
%                    750749662
%                    751348571
%                    754312389
%                    754829445
%                    755434585
%                    756029989
%                    757216464
%                    757970808
%                    758798717
%                    759883607
%                    760345702
%                    760693773
%                    761418226
%                    762120172
%                    762602078
%                    763673393
%                    766640955
%                    767871931
%                    768515987
%                    771160300
%                    771990200
%                    773418906
%                    774875821
%                    778240327
%                    778998620
%                    779839471
%                    781842082
%                    786091066
%                    787025148
%                    789848216
%                    791319847
%                    793224716
%                    794812542
%                    797828357
%                    798911424
%                    799864342
%                    816200189
%                    819186360
%                    819701982
%                    821695405
%                    829720705
%                    831882777
%                    835479236
%                    839068429
%                    839557629
%                    840012044
%                    847657808];
% 
%     datasetID   = 'AllenInstitute';
%     datapath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
% 
%     % Allen
% 
%     pairGroupStatTableExq = table;
% 
%     for loopAnimals = 1:58
% 
%         display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])
% 
%         load([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);
% 
%         pairGroupStatTableExq = [pairGroupStatTableExq; pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:)];
% 
%     end
% 
%     %
% 
%     figure('Renderer', 'painters', 'Position', [1720 2562 2*560 420])
% 
%     tiledlayout(1,2)
% 
%     nexttile
% 
%     plot(pairGroupStatTableExq.CCGbinLagTimes{1,1}*1000,...
%          mean([pairGroupStatTableExq.pairRawCCG{:,1}],2),...
%           'k','LineWidth',4);
%     hold on 
%     for loopPairs = 1:size(pairGroupStatTableExq,1)
%         patchline(pairGroupStatTableExq.CCGbinLagTimes{1,1}*1000, ... 
%                  [pairGroupStatTableExq.pairRawCCG{loopPairs,1}], ...
%                  'linewidth',2,'edgealpha',0.1)
%     end
%     hold off
%     xlim([-3.5 3.5])
%     ylabel('counts')
%     xlabel('[ms]')
%     set(gca,'FontSize',12)
% 
%     nexttile
%     semilogy(cell2num([pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'p'));       ...
%                        pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'p'));        ...
%                        pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-wide'));  ...
%                        pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-wide'))]), ... 
%                       [pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'p'));       ...
%                        pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'p'));        ...
%                        pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-wide'));  ...
%                        pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-wide'))], ...
%                      '^','LineWidth',4)
% 
%     hold on
%     semilogy(cell2num([pairGroupStatTableExq.refTroughToPeak(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-narrow'));  ...
%                        pairGroupStatTableExq.tarTroughToPeak(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-narrow'))]), ... 
%                       [pairGroupStatTableExq.refACGtauRise(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-narrow'));  ...
%                        pairGroupStatTableExq.tarACGtauRise(strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-narrow'))], ...
%                       'o','LineWidth',4)
%     for loopPairs = 1:size(pairGroupStatTableExq,1)
%         patchline(cell2num([pairGroupStatTableExq.refTroughToPeak(loopPairs),pairGroupStatTableExq.tarTroughToPeak(loopPairs)]), ... 
%                            [pairGroupStatTableExq.refACGtauRise(loopPairs),pairGroupStatTableExq.tarACGtauRise(loopPairs)], ...
%                           'linewidth',2,'edgealpha',0.5)
%     end
%     hold off
%     set(gca, 'YScale', 'log')
%     legend('p','i','Location','SouthEast')
%     ylim([0.05 30])
%     xlim([0.1 2])
%     ylabel('ACG \tau_r_i_s_e')
%     xlabel('trough-to-peak [ms]')
%     set(gca,'FontSize',12)
%     
%     set(gca,'FontSize',12)
%     set(gcf, 'Renderer', 'Painters');
%     
% 
%     ii_count = sum(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-narrow'));
%     pp_count = sum((strcmp(pairGroupStatTableExq.refCellExplorerType,'p') | strcmp(pairGroupStatTableExq.refCellExplorerType,'i-wide') ) & ...
%                    (strcmp(pairGroupStatTableExq.tarCellExplorerType,'p') | strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-wide')));
%     ip_count = sum(strcmp(pairGroupStatTableExq.refCellExplorerType,'i-narrow')  & ...
%                   (strcmp(pairGroupStatTableExq.tarCellExplorerType,'p') | strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-wide')));
%     pi_count = sum((strcmp(pairGroupStatTableExq.refCellExplorerType,'p') | strcmp(pairGroupStatTableExq.refCellExplorerType,'i-wide') ) & ...
%                     strcmp(pairGroupStatTableExq.tarCellExplorerType,'i-narrow'));
%     % box on
%     % boxplot([ii_count,pp_count,ip_count+pi_count])
%     % box off

    %%
    
    hcomb = figure(106);

    res_type = 'QHD';
    pos = [70 230 1920*0.4 (2/3)*1080/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    RoySummaryStatsGJ = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ_CA1.mat');
    AGsummaryStatsGJ  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ_CA1.mat');
    NSUsummaryStatsGJ = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ_CA1.mat');
    
    RoySummaryStatsExq = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStats_CA1.mat');
    AGsummaryStatsExq  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStats_CA1.mat');
    NSUsummaryStatsExq = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStats_CA1.mat');
    
    BOTsummaryStatsGJ  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ_CA1.mat');
    BOTsummaryStatsExq = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStats_CA1.mat');
    
    tiledlayout(1,6)
    
    tableTemp = table;
    
    tableTemp.category = [string(repmat('ms',58,1)); string(repmat('exq',58,1))];
    
    nexttile
    tableTemp.number = [100*BOTsummaryStatsGJ.fracGJperInd{1,1}{1,1}'; ...
                        100*BOTsummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1}'];
    boxchart(tableTemp.number,'GroupByColor',categorical(tableTemp.category))
    title({'Allen Int.','CA1 ii'})
    ylabel('percent')
    legend('location','southoutside')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTickLabel',[]);
    
    nexttile
    tableTemp.number = [100*BOTsummaryStatsGJ.fracGJperInd{2,1}{1,1}'; ...
                        100*BOTsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1}'];
    boxchart(tableTemp.number,'GroupByColor',categorical(tableTemp.category))
    title({'Allen Int.','CA1 pp'})
    ylabel('percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTickLabel',[]);
    
    nexttile
    tableTemp.number = [100*BOTsummaryStatsGJ.fracGJperInd{3,1}{1,1}'; ...
                        100*BOTsummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1}'];
    boxchart(tableTemp.number,'GroupByColor',categorical(tableTemp.category))
    title({'Allen Int.','CA1 ip or pi'})
    ylabel('percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTickLabel',[]);
    
    tableTemp = table;
    
    tableTemp.category = [string(repmat('ms',5,1)); string(repmat('exq',5,1))];
    
    nexttile
    tableTemp.number = [100*[     RoySummaryStatsGJ.fracGJperInd{1,1}{1,1}{1,1} ...
                   mean(AGsummaryStatsGJ.fracGJperInd{1,1}{1,1}{1,1}) ...
                        NSUsummaryStatsGJ.fracGJperInd{1,1}{1,1}{1,1}]';
                        100*[     RoySummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1} ...
                   mean(AGsummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1}) ...
                        NSUsummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1}]'];
    boxchart(tableTemp.number,'GroupByColor',categorical(tableTemp.category))
    title({'rats','CA1 ii'})
    ylabel('percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTickLabel',[]);
    
    nexttile
    tableTemp.number = [100*[     RoySummaryStatsGJ.fracGJperInd{2,1}{1,1}{1,1} ...
                   mean(AGsummaryStatsGJ.fracGJperInd{2,1}{1,1}{1,1}) ...
                        NSUsummaryStatsGJ.fracGJperInd{2,1}{1,1}{1,1}]';
                        100*[     RoySummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1} ...
                   mean(AGsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1}) ...
                        NSUsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1}]'];
    boxchart(tableTemp.number,'GroupByColor',categorical(tableTemp.category))
    title({'rats','CA1 pp'})
    ylabel('percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTickLabel',[]);
    
    nexttile
    tableTemp.number = [100*[     RoySummaryStatsGJ.fracGJperInd{3,1}{1,1}{1,1} ...
                   mean(AGsummaryStatsGJ.fracGJperInd{3,1}{1,1}{1,1}) ...
                        NSUsummaryStatsGJ.fracGJperInd{3,1}{1,1}{1,1}]';
                        100*[     RoySummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1} ...
                   mean(AGsummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1}) ...
                        NSUsummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1}]'];
    boxchart(tableTemp.number,'GroupByColor',categorical(tableTemp.category))
    title({'rats','CA1 ip or pi'})
    ylabel('percent')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    set(gca,'XTickLabel',[]);

end

%% Figure 3i common input simulations
if strcmp(figureCode,'3i')

%     figure('Renderer', 'painters', 'Position', [1720 2562 560 420*(2.5)])
    figure('Renderer', 'painters', 'Position', [1720 2562 2*144 0.5*256])
    
    tiledlayout(1,2,'Padding','none','TileSpacing','compact')
    
    timeMin = 12; % in seconds 
    DC      = 3.55;
    
    % % GJ common input
    % HodgkinHuxleyGJohmicCommonInput(timeMin,DC)

    % GJ common input
    HodgkinHuxleyGJohmic(timeMin,DC)
    
    DC      = 3.5;

    % E/E     
    refPSPvalence = 1;
    tarPSPvalence = 1;
    delayRef = 0; % ms
    delayTar = 0;% ms
    HodgkinHuxleyCommonInput(timeMin,refPSPvalence,tarPSPvalence,delayRef,delayTar,DC)
    
%     % I/I 
%     refPSPvalence = -1;
%     tarPSPvalence = -1;
%     delayRef = 0; % ms
%     delayTar = 0;% ms
%     HodgkinHuxleyCommonInput(timeMin,refPSPvalence,tarPSPvalence,delayRef,delayTar)
%   
%     % E/I 
%     refPSPvalence = 1;
%     tarPSPvalence = -1;
%     delayRef = 0; % ms
%     delayTar = 0;% ms
%     HodgkinHuxleyCommonInput(timeMin,refPSPvalence,tarPSPvalence,delayRef,delayTar)
    
    % Eph common input
%     inputType = 'datasetEspike';
%     HodgkinHuxleyEphapseCommonInput(timeMin,inputType)
    
%     inputType = 'deltaFunc';
%     HodgkinHuxleyEphapseCommonInput(timeMin,inputType)
    
end

%% Figure 3ii pair-wise simulations

if strcmp(figureCode,'3ii')
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    set(gcf, 'Renderer', 'painters', 'Position', [1720 2562 2*144 256]);

    res_type = 'QHD';
    
    tiledlayout(2,2,'Padding','none','TileSpacing','tight')
    
    timeSec = 12; % in seconds
    DC      = 3.5;
    
%     HodgkinHuxleyMonosynapse(timeSec)
    HodgkinHuxleyRectGJohmic(timeSec,DC)
    HodgkinHuxleyEpSohmic(timeSec,DC)
%     inputType = 'datasetEspike';
%     HodgkinHuxleyEphapse(timeSec,inputType)
%     inputType = 'deltaFunc';
%     HodgkinHuxleyEphapse(timeSec,inputType)
    
end

%% Figure 4i alignment figure 

if strcmp(figureCode,'4i')
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    set(gcf, 'Renderer', 'painters', 'Position', [1720 2562 3*144 2*256*(2/3)]);

    res_type = 'QHD';
    
    tiledlayout(3,4,'Padding','none','TileSpacing','tight')

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable
    
    % idx 23
    ref = pairGroupStatsTable.tarSpikeTimes{23,1};
    tar = pairGroupStatsTable.refSpikeTimes{23,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{23,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{23,1};
    d = pairGroupStatsTable.pairDistance(23);
    region = pairGroupStatsTable.brainRegion{23,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));
    
    nexttile(1)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    % text(-0.85,0.85*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz'];['tar ' num2str(tarFR,'%.1f') 'hz']},'FontSize',5)
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile(5)
    plot(pairGroupStatsTable.tWave{23,1}, ...
         pairGroupStatsTable.refWaveforms{23,1}(pairGroupStatsTable.refChannel(23),:) - ...
         pairGroupStatsTable.refWaveforms{23,1}(pairGroupStatsTable.refChannel(23),1),'LineWidth',1,'color',"#D95319")
    hold on 
    plot(pairGroupStatsTable.tWave{23,1}, ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.tarChannel(23),:) - ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.tarChannel(23),1),'LineWidth',1,'color','#0072BD')
    plot(pairGroupStatsTable.tWave{23,1}, ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),:) - ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),1),'-.','LineWidth',1,'color','#0072BD')
    
    hold off
    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Tar e-spike at tar max ch',...
%            'Ref e-spike at ref max ch',...
%            'Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

    nexttile(9)
    plot(pairGroupStatsTable.tWave{23,1}, ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),:) - ...
         pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),1),'-.','LineWidth',1,'color','#0072BD')

    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

    % idx 16
    ref = pairGroupStatsTable.tarSpikeTimes{16,1};
    tar = pairGroupStatsTable.refSpikeTimes{16,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{16,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{16,1};
    d = pairGroupStatsTable.pairDistance(16);
    region = pairGroupStatsTable.brainRegion{16,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(3)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    % text(-0.85,0.85*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz'];['tar ' num2str(tarFR,'%.1f') 'hz']},'FontSize',5)
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile(7)
    plot(pairGroupStatsTable.tWave{16,1}, ...
         pairGroupStatsTable.refWaveforms{16,1}(pairGroupStatsTable.refChannel(16),:) - ...
         pairGroupStatsTable.refWaveforms{16,1}(pairGroupStatsTable.refChannel(16),1),'LineWidth',1,'color',"#D95319")
    hold on 
    plot(pairGroupStatsTable.tWave{16,1}, ...
         pairGroupStatsTable.tarWaveforms{16,1}(pairGroupStatsTable.tarChannel(16),:) - ...
         pairGroupStatsTable.tarWaveforms{16,1}(pairGroupStatsTable.tarChannel(16),1),'LineWidth',1,'color','#0072BD')
    plot(pairGroupStatsTable.tWave{16,1}, ...
         pairGroupStatsTable.tarWaveforms{16,1}(pairGroupStatsTable.refChannel(16),:) - ...
         pairGroupStatsTable.tarWaveforms{16,1}(pairGroupStatsTable.refChannel(16),1),'-.','LineWidth',1,'color','#0072BD')
 
    hold off
    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Tar e-spike at tar max ch',...
%            'Ref e-spike at ref max ch',...
%            'Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

    nexttile(11)
    plot(pairGroupStatsTable.tWave{16,1}, ...
         pairGroupStatsTable.tarWaveforms{16,1}(pairGroupStatsTable.refChannel(16),:) - ...
         pairGroupStatsTable.tarWaveforms{16,1}(pairGroupStatsTable.refChannel(16),1),'-.','LineWidth',1,'color','#0072BD')

    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/groupStatsRatAGday2.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable
    
    % idx 21
    ref = pairGroupStatsTable.tarSpikeTimes{21,1};
    tar = pairGroupStatsTable.refSpikeTimes{21,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{21,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{21,1};
    d = pairGroupStatsTable.pairDistance(21);
    region = pairGroupStatsTable.tarBrainRegion(21,:);
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(2)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat AG';[region ' ' num2str(d,'%.0f') '\mum'];newline})
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    % text(-0.85,0.85*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz'];['tar ' num2str(tarFR,'%.1f') 'hz']},'FontSize',5)
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile(6)
    plot(pairGroupStatsTable.tWave{21,1}, ...
         pairGroupStatsTable.refWaveforms{21,1}(pairGroupStatsTable.refChannel(21),:) - ...
         pairGroupStatsTable.refWaveforms{21,1}(pairGroupStatsTable.refChannel(21),1),'LineWidth',1,'color',"#D95319")
    hold on 
    plot(pairGroupStatsTable.tWave{21,1}, ...
         pairGroupStatsTable.tarWaveforms{21,1}(pairGroupStatsTable.tarChannel(21),:) - ...
         pairGroupStatsTable.tarWaveforms{21,1}(pairGroupStatsTable.tarChannel(21),1),'LineWidth',1,'color','#0072BD')
    plot(pairGroupStatsTable.tWave{21,1}, ...
         pairGroupStatsTable.tarWaveforms{21,1}(pairGroupStatsTable.refChannel(21),:) - ...
         pairGroupStatsTable.tarWaveforms{21,1}(pairGroupStatsTable.refChannel(21),1),'-.','LineWidth',1,'color','#0072BD')
     
     
    hold off
    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Tar e-spike at tar max ch',...
%            'Ref e-spike at ref max ch',...
%            'Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

    nexttile(10)
    plot(pairGroupStatsTable.tWave{21,1}, ...
         pairGroupStatsTable.tarWaveforms{21,1}(pairGroupStatsTable.refChannel(21),:) - ...
         pairGroupStatsTable.tarWaveforms{21,1}(pairGroupStatsTable.refChannel(21),1),'-.','LineWidth',1,'color','#0072BD')
    
    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off
    
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/groupStatsRatU')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable
    
    % idx 73
    ref = pairGroupStatsTable.refSpikeTimes{76,1};
    tar = pairGroupStatsTable.tarSpikeTimes{76,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{76,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{76,1};
    d = pairGroupStatsTable.pairDistance(76);
    region = pairGroupStatsTable.brainRegion(76,:);
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(4)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat U';[char(region) ' ' num2str(d,'%.0f') '\mum'];newline})
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
    % text(-0.85,0.35*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz'];['tar ' num2str(tarFR,'%.1f') 'hz']},'FontSize',5)
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile(8)
    plot(pairGroupStatsTable.tWave{14,1}, ...
         pairGroupStatsTable.tarWaveforms{76,1}(pairGroupStatsTable.tarChannel(76),:) - ...
         pairGroupStatsTable.tarWaveforms{76,1}(pairGroupStatsTable.tarChannel(76),1),'LineWidth',1,'color',"#D95319")
    hold on 
    plot(pairGroupStatsTable.tWave{21,1}, ...
         pairGroupStatsTable.refWaveforms{76,1}(pairGroupStatsTable.refChannel(76),:) - ...
         pairGroupStatsTable.refWaveforms{76,1}(pairGroupStatsTable.refChannel(76),1),'LineWidth',1,'color','#0072BD')
    plot(pairGroupStatsTable.tWave{21,1}, ...
         pairGroupStatsTable.refWaveforms{76,1}(pairGroupStatsTable.tarChannel(76),:) - ...
         pairGroupStatsTable.refWaveforms{76,1}(pairGroupStatsTable.tarChannel(76),1),'-.','LineWidth',1,'color','#0072BD')
     
     
    hold off
    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Tar e-spike at tar max ch',...
%            'Ref e-spike at ref max ch',...
%            'Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

    nexttile(12)
    plot(pairGroupStatsTable.tWave{14,1}, ...
         pairGroupStatsTable.refWaveforms{76,1}(pairGroupStatsTable.tarChannel(76),:) - ...
         pairGroupStatsTable.refWaveforms{76,1}(pairGroupStatsTable.tarChannel(76),1),'-.','LineWidth',1,'color','#0072BD')
    
    set(gca, 'YDir','reverse')
    xlim([-1,1])
%     legend('Ref e-spike at tar max ch','Location','NorthWest')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    xlabel('[ms]');
    ylabel('[mV]');
    box off

%     load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/groupStatsRatN')
%     pairGroupStatsTable = pairGroupStatTable;
%     clear pairGroupStatTable
% 
%     % idx 14
%     ref = pairGroupStatsTable.refSpikeTimes{14,1};
%     tar = pairGroupStatsTable.tarSpikeTimes{14,1}; 
%     refCelltype = pairGroupStatsTable.refCellExplorerType{14,1};
%     tarCelltype = pairGroupStatsTable.tarCellExplorerType{14,1};
%     d = pairGroupStatsTable.pairDistance(14);
%     region = pairGroupStatsTable.brainRegion(14,:);
%     refFR = length(ref)/(ref(end)-ref(1));
%     tarFR = length(tar)/(tar(end)-tar(1));
% 
%     nexttile(4)
%     [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%               CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
%                                     'plot_output', get(fig_use, 'Number'), ...
%                                     'njitter', njitter, 'alpha', alpha,...
%                                     'for_grant', for_grant, 'plot_pointwiseBands', false, ...
%                                     'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
%                                     'norm_flag', true);
%     xlim([-1,1])
%     ylims = get(gca,'ylim');
%     ylabel('Spike Probability')
%     xlabel('[ms]')
%     title({'Rat N';[char(region) ' ' num2str(d,'%.0f') '\mum'];newline})
% %     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
% %     if any(GSPExc)
% %         hold on;
% %         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
% %     end
% %     if any(GSPInh)
% %         hold on;
% %         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
% %     end
%     text(-0.85,0.85*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz'];['tar ' num2str(tarFR,'%.1f') 'hz']},'FontSize',5)
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
% 
%     nexttile(8)
%     plot(pairGroupStatsTable.tWave{14,1}, ...
%          pairGroupStatsTable.tarWaveforms{14,1}(pairGroupStatsTable.tarChannel(14),:) - ...
%          pairGroupStatsTable.tarWaveforms{14,1}(pairGroupStatsTable.tarChannel(14),1),'LineWidth',1,'color',"#D95319")
%     hold on 
%     plot(pairGroupStatsTable.tWave{21,1}, ...
%          pairGroupStatsTable.refWaveforms{14,1}(pairGroupStatsTable.refChannel(14),:) - ...
%          pairGroupStatsTable.refWaveforms{14,1}(pairGroupStatsTable.refChannel(14),1),'LineWidth',1,'color','#0072BD')
%     plot(pairGroupStatsTable.tWave{21,1}, ...
%          pairGroupStatsTable.refWaveforms{14,1}(pairGroupStatsTable.refChannel(14),:) - ...
%          pairGroupStatsTable.refWaveforms{14,1}(pairGroupStatsTable.refChannel(14),1),'-.','LineWidth',1,'color','#0072BD')
% 
% 
%     hold off
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Tar e-spike at tar max ch',...
% %            'Ref e-spike at ref max ch',...
% %            'Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off
% 
%     nexttile(12)
%     plot(pairGroupStatsTable.tWave{14,1}, ...
%          pairGroupStatsTable.refWaveforms{14,1}(pairGroupStatsTable.tarChannel(14),:) - ...
%          pairGroupStatsTable.refWaveforms{14,1}(pairGroupStatsTable.tarChannel(14),1),'-.','LineWidth',1,'color','#0072BD')
% 
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off

    %%
    
%     refChan = pairGroupStatsTable.refChannel(14); 
%     tarChan = pairGroupStatsTable.tarChannel(14); 
%     chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(refChan) 'highpass300hz.mat']);
%     chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(tarChan) 'highpass300hz.mat']);
% 
%     lagLimit = duration/2;
%     % removing the synchronous spikes                
%     [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(ref',tar',binSize,lagLimit);
% 
%     % peak
%     refPeak = spikeTimesInBin{106}(:,1);
%     tarPeak = spikeTimesInBin{106}(:,2);
% 
%     refPeakLat   = round(refPeak*fs) + 1;
%     tarPeakLat   = round(tarPeak*fs) + 1;
% 
%     fpass      = 300;
% 
%     preLength  = 36;
%     postLength = 36;
% 
%     noDemeanFlag = true;
% 
%     % non-synch
%     [refPeakLatWaveMean,  refPeakLatWaveforms]   = waveformAvg(chanDataRef.data,refPeakLat,    preLength,postLength,fpass,fs,false,noDemeanFlag);
%     [tarPeakLatWaveMean,  tarPeakLatWaveforms]   = waveformAvg(chanDataTar.data,tarPeakLat,    preLength,postLength,fpass,fs,false,noDemeanFlag);
% 
%     nexttile(19)
%     plot(pairGroupStatsTable.tWave{14,1},refPeakLatWaveMean)
%     hold on 
%     plot(pairGroupStatsTable.tWave{14,1},tarPeakLatWaveMean)
%     hold off
% 
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off
% 
%     % idx 19
%     ref = pairGroupStatsTable.refSpikeTimes{19,1};
%     tar = pairGroupStatsTable.tarSpikeTimes{19,1}; 
%     refCelltype = pairGroupStatsTable.refCellExplorerType{19,1};
%     tarCelltype = pairGroupStatsTable.tarCellExplorerType{19,1};
%     d = pairGroupStatsTable.pairDistance(19);
%     region = pairGroupStatsTable.brainRegion(19,:);
%     refFR = length(ref)/(ref(end)-ref(1));
%     tarFR = length(tar)/(tar(end)-tar(1));
% 
%     nexttile(5)
%     [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%               CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
%                                     'plot_output', get(fig_use, 'Number'), ...
%                                     'njitter', njitter, 'alpha', alpha,...
%                                     'for_grant', for_grant, 'plot_pointwiseBands', false, ...
%                                     'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
%                                     'norm_flag', true);
%     xlim([-1,1])
%     ylims = get(gca,'ylim');
%     ylabel('Spike Probability')
%     xlabel('[ms]')
%     title({'Rat N';[char(region) ' ' num2str(d,'%.0f') '\mum'];newline})
% %     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
% %     if any(GSPExc)
% %         hold on;
% %         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
% %     end
% %     if any(GSPInh)
% %         hold on;
% %         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
% %     end
%     text(-0.85,0.85*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz'];['tar ' num2str(tarFR,'%.1f') 'hz']},'FontSize',5)
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
% 
%     nexttile(10)
%     plot(pairGroupStatsTable.tWave{19,1}, ...
%          pairGroupStatsTable.tarWaveforms{19,1}(pairGroupStatsTable.tarChannel(19),:) - ...
%          pairGroupStatsTable.tarWaveforms{19,1}(pairGroupStatsTable.tarChannel(19),1),'LineWidth',1,'color',"#D95319")
%     hold on 
%     plot(pairGroupStatsTable.tWave{21,1}, ...
%          pairGroupStatsTable.refWaveforms{19,1}(pairGroupStatsTable.refChannel(19),:) - ...
%          pairGroupStatsTable.refWaveforms{19,1}(pairGroupStatsTable.refChannel(19),1),'LineWidth',1,'color','#0072BD')
%     plot(pairGroupStatsTable.tWave{21,1}, ...
%          pairGroupStatsTable.refWaveforms{19,1}(pairGroupStatsTable.refChannel(19),:) - ...
%          pairGroupStatsTable.refWaveforms{19,1}(pairGroupStatsTable.refChannel(19),1),'-.','LineWidth',1,'color','#0072BD')
% 
% 
%     hold off
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Tar e-spike at tar max ch',...
% %            'Ref e-spike at ref max ch',...
% %            'Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off
% 
%     nexttile(15)
%     plot(pairGroupStatsTable.tWave{19,1}, ...
%          pairGroupStatsTable.refWaveforms{19,1}(pairGroupStatsTable.tarChannel(19),:) - ...
%          pairGroupStatsTable.refWaveforms{19,1}(pairGroupStatsTable.tarChannel(19),1),'-.','LineWidth',1,'color','#0072BD')
% 
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off
% 
%     %%
% 
%     refChan = pairGroupStatsTable.refChannel(19); 
%     tarChan = pairGroupStatsTable.tarChannel(19); 
%     chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(refChan) 'highpass300hz.mat']);
%     chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(tarChan) 'highpass300hz.mat']);
% 
%     lagLimit = duration/2;
%     % removing the synchronous spikes                
%     [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(ref',tar',binSize,lagLimit);
% 
%     % peak
%     refPeak = spikeTimesInBin{106}(:,1);
%     tarPeak = spikeTimesInBin{106}(:,2);
% 
%     refPeakLat   = round(refPeak*fs) + 1;
%     tarPeakLat   = round(tarPeak*fs) + 1;
% 
%     fpass      = 300;
% 
%     preLength  = 36;
%     postLength = 36;
% 
%     noDemeanFlag = true;
% 
%     % non-synch
%     [refPeakLatWaveMean,  refPeakLatWaveforms]   = waveformAvg(chanDataRef.data,refPeakLat,    preLength,postLength,fpass,fs,false,noDemeanFlag);
%     [tarPeakLatWaveMean,  tarPeakLatWaveforms]   = waveformAvg(chanDataTar.data,tarPeakLat,    preLength,postLength,fpass,fs,false,noDemeanFlag);
% 
%     nexttile(20)
%     plot(pairGroupStatsTable.tWave{19,1},refPeakLatWaveMean)
%     hold on 
%     plot(pairGroupStatsTable.tWave{19,1},tarPeakLatWaveMean)
%     hold off
% 
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off

    set(gcf, 'Renderer', 'Painters');
    
%     
%     nexttile(11)
%     plot(pairGroupStatsTable.tWave{23,1}(2:end)*1000, ...
%          diff(pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),:) - ...
%               pairGroupStatsTable.tarWaveforms{23,1}(pairGroupStatsTable.refChannel(23),1)), ...
%          'LineWidth',4,'color','#0072BD')
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     xlabel('Time Lag [\musec]');
% %     ylabel('voltage [mV]');
    
%     timeSec = 12; % in seconds
%     
% %     HodgkinHuxleyMonosynapse(timeSec)
% %     HodgkinHuxleyGJohmic(timeSec)
%     HodgkinHuxleyEpSohmic(timeSec)
%     inputType = 'datasetEspike';
%     HodgkinHuxleyEphapse(timeSec,inputType)
% %     inputType = 'deltaFunc';
% %     HodgkinHuxleyEphapse(timeSec,inputType)
    
end

%% Figure 4ii alignment figure 

if strcmp(figureCode,'4ii')
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    set(gcf, 'Renderer', 'painters', 'Position', [1720 2562 3*144 256]);

    res_type = 'QHD';
    
    tiledlayout(2,4,'Padding','none','TileSpacing','tight')

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable
    
    % idx 3
    ref = pairGroupStatsTable.tarSpikeTimes{3,1};
    tar = pairGroupStatsTable.refSpikeTimes{3,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{3,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{3,1};
    d = pairGroupStatsTable.pairDistance(3);
    region = pairGroupStatsTable.brainRegion{3,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));
    
    nexttile(1)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 19
    ref = pairGroupStatsTable.refSpikeTimes{19,1};
    tar = pairGroupStatsTable.tarSpikeTimes{19,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{19,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{19,1};
    d = pairGroupStatsTable.pairDistance(19);
    region = pairGroupStatsTable.brainRegion{19,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(2)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 42
    ref = pairGroupStatsTable.refSpikeTimes{42,1};
    tar = pairGroupStatsTable.tarSpikeTimes{42,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{42,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{42,1};
    d = pairGroupStatsTable.pairDistance(42);
    region = pairGroupStatsTable.brainRegion{42,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(3)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 6
    ref = pairGroupStatsTable.refSpikeTimes{6,1};
    tar = pairGroupStatsTable.tarSpikeTimes{6,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{6,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{6,1};
    d = pairGroupStatsTable.pairDistance(6);
    region = pairGroupStatsTable.brainRegion{6,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(4)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    %%

    % idx = 5
    
    % load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
    
    channelSet = {1:8, ...
                  9:16, ... 
                  17:24, ...
                  25:32, ...
                  33:40, ...
                  41:48, ...
                  48:56, ...
                  57:64};

    channelSetPlotIdx = [1,2,4,5]; % idx 3,19,42,6
    
    waveforms = pairGroupStatsTable.refWaveforms{5};
    
    shank       = pairGroupStatsTable.refShank(5);
    chan        = pairGroupStatsTable.refChannel(5);
    cellID      = pairGroupStatsTable.refNeuronID(5);
    cellType    = pairGroupStatsTable.refCellExplorerType{5};
    
    % titleStr    = {['rat ' animalNames{6} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    
    tWave = pairGroupStatsTable.tWave{1};
    
    for loopSets = 1:size(channelSetPlotIdx,2)
    
        nexttile(loopSets+4)
        plot(tWave,waveforms(channelSet{channelSetPlotIdx(loopSets)},:)');
        % 
        % if loopSets == 1
        %     title(titleStr)
        % else 
        %     title({[' sh. ' num2str(loopSets)],[]})
        % end
        
        xlim([-1,1])
        % ylim([-7 2])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        set(gca, 'YDir','reverse')
        box off
    
    end

end

%% Figure 4iii alignment figure 

if strcmp(figureCode,'4iii')
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    set(gcf, 'Renderer', 'painters', 'Position', [1720 2562 0.75*144 2.5*256]);

    res_type = 'QHD';
    
    tiledlayout(6,1,'Padding','none','TileSpacing','tight')

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable
    
    % idx 2
    ref = pairGroupStatsTable.tarSpikeTimes{2,1};
    tar = pairGroupStatsTable.refSpikeTimes{2,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{2,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{2,1};
    d = pairGroupStatsTable.pairDistance(2);
    region = pairGroupStatsTable.brainRegion{2,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));
    
    nexttile(1)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    % title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 22
    ref = pairGroupStatsTable.refSpikeTimes{22,1};
    tar = pairGroupStatsTable.tarSpikeTimes{22,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{22,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{22,1};
    d = pairGroupStatsTable.pairDistance(22);
    region = pairGroupStatsTable.brainRegion{22,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(2)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    % title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 38
    ref = pairGroupStatsTable.refSpikeTimes{38,1};
    tar = pairGroupStatsTable.tarSpikeTimes{38,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{38,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{38,1};
    d = pairGroupStatsTable.pairDistance(38);
    region = pairGroupStatsTable.brainRegion{38,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(3)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    % title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 13
    ref = pairGroupStatsTable.refSpikeTimes{13,1};
    tar = pairGroupStatsTable.tarSpikeTimes{13,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{13,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{13,1};
    d = pairGroupStatsTable.pairDistance(13);
    region = pairGroupStatsTable.brainRegion{13,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(5)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    % title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    % idx 12
    ref = pairGroupStatsTable.refSpikeTimes{12,1};
    tar = pairGroupStatsTable.tarSpikeTimes{12,1}; 
    refCelltype = pairGroupStatsTable.refCellExplorerType{12,1};
    tarCelltype = pairGroupStatsTable.tarCellExplorerType{12,1};
    d = pairGroupStatsTable.pairDistance(12);
    region = pairGroupStatsTable.brainRegion{12,1};
    refFR = length(ref)/(ref(end)-ref(1));
    tarFR = length(tar)/(tar(end)-tar(1));

    nexttile(4)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    % title({'Rat R';[region ' ' num2str(d,'%.0f') '\mum'];newline})
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    %%

    % idx = 5
    
    % load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
    
    channelSet = {1:8, ...
                  9:16, ... 
                  17:24, ...
                  25:32, ...
                  33:40, ...
                  41:48, ...
                  48:56, ...
                  57:64};

    channelSetPlotIdx = [1]; 
    
    waveforms = pairGroupStatsTable.refWaveforms{12};
    
    shank       = pairGroupStatsTable.refShank(12);
    chan        = pairGroupStatsTable.refChannel(12);
    cellID      = pairGroupStatsTable.refNeuronID(12);
    cellType    = pairGroupStatsTable.refCellExplorerType{12};
    
    % titleStr    = {['rat ' animalNames{6} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    
    tWave = pairGroupStatsTable.tWave{1};
    
    for loopSets = 1:size(channelSetPlotIdx,2)
    
        nexttile
        plot(tWave,waveforms(channelSet{channelSetPlotIdx(loopSets)},:)');
        % 
        % if loopSets == 1
        %     title(titleStr)
        % else 
        %     title({[' sh. ' num2str(loopSets)],[]})
        % end
        
        xlim([-1,1])
        % ylim([-7 2])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        set(gca, 'YDir','reverse')
        box off
    
    end

end

%% Figure 5 Steinmetz 
if figureCode == '5'
    load('/home/nasko/CUNY_Work_NPUltraWaveforms/data/groupStatsSteinmetzDataset.mat')

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

%     res_type = 'QHD';
%     pos = [1720 2562 3*560 0.75*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
%     arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    res_type = 'QHD';
    pos = [1720 2562 4*560*0.20 6*420*0.3]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    neuronIdx = [22,52,62,76];

    % tiledlayout(5,3,'Padding','compact','TileSpacing','compact')
    tiledlayout(6,6,'Padding','compact','TileSpacing','compact')
%     tiledlayout(6,4)


    refCellType = {'i','i','p','i'};  
    tarCelltype = {'p','i','p','i'}; 

    for i = 1:length(neuronIdx)-1

        nexttile([1,2])

        ref = pairGroupStatTable.refSpikeTimes{neuronIdx(i),1};
        tar = pairGroupStatTable.tarSpikeTimes{neuronIdx(i),1};

        region = pairGroupStatTable.brainRegion{neuronIdx(i),1};

        d = pairGroupStatTable.pairDistance(neuronIdx(i));

%         [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%               CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
%                                     'plot_output', get(fig_use, 'Number'), ...
%                                     'njitter', njitter, 'alpha', alpha,...
%                                     'for_grant', for_grant, 'plot_pointwiseBands', false, ...
%                                     'lineWidth', 4, 'axisLabels', false, 'microsecFlag', false);

        [ccgR, tR] = CCG([ref';tar'],[ones(size(ref'));2*ones(size(tar'))], ...
                'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                'norm', 'counts');
        line(tR*1e3,ccgR(:,1,2)/length(ref),'color','k', 'LineWidth',1)
        
        p_sum = sum(ccgR(76:136,1,2)/length(ref));
        
        xlim([-1,1])
        ylims = get(gca,'ylim');
        ylabel('Spike probability')
        xlabel('[ms]')
        set(gca,'FontSize',5)
        set(gcf, 'Renderer', 'Painters');
        % title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum'])
%         title([refCellType{i} '-' tarCelltype{i} ' pair ' num2str(d) '\mum ' region])
        title({['Mouse'],[region ' ' num2str(d,'%.0f') '\mum'],newline,['-1 to 1ms spike'],['probability = ' num2str(p_sum,'%.2f')],newline})

    %     if any(GSPExc)
    %         hold on;
    %         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
    %     end
    %     if any(GSPInh)
    %         hold on;
    %         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
    %     end
        set(gca,'FontSize',5)
    end
   
    for i = 1:length(neuronIdx)-1 % fourth example edited out Dec 2024
        
        %% prep title 
        ref = pairGroupStatTable.refSpikeTimes{neuronIdx(i),1};
        tar = pairGroupStatTable.tarSpikeTimes{neuronIdx(i),1};

        region = pairGroupStatTable.brainRegion{neuronIdx(i),1};

        d = pairGroupStatTable.pairDistance(neuronIdx(i));
        
        [ccgR, tR] = CCG([ref';tar'],[ones(size(ref'));2*ones(size(tar'))], ...
                'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
                'norm', 'counts');
%         line(tR*1e3,ccgR(:,1,2)/length(ref),'color','k', 'LineWidth',1)
        
        p_sum = sum(ccgR(76:136,1,2)/length(ref));
        
%         title({['Mouse'],[region ' ' num2str(d,'%.0f') '\mum'],['-1 to 1ms spike'],['probability = ' num2str(p_sum,'%.2f')],newline})

        %%
        
        waveformRef = pairGroupStatTable.refWaveforms{neuronIdx(i),1};
        waveformTar = pairGroupStatTable.tarWaveforms{neuronIdx(i),1};

        idx = 43; % time stamp for extracell trough

        chanLayout = flipud(reshape(1:384,[8,48])');

        x = 0:5:5*7;
        y = 0:5:5*47;

        for loopChans = 1:(48*8)

            [a,b] = find(chanLayout == loopChans);
            coordinates(loopChans,2) = y(a);
            coordinates(loopChans,1) = x(b);

            heat_values_ref(loopChans) = waveformRef(idx,loopChans);
            heat_values_tar(loopChans) = waveformTar(idx,loopChans);

        end
        
        %% ref

        % Create a 3D surface plot
        nexttile([5,1])
        tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
        trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values_ref, 'FaceColor', 'interp', 'EdgeColor', 'none');
        xlabel('[\mum]');
        ylabel('[\mum]');
    %     title('Activity Plot');

        axis equal

        % xlim([0 35])
        % ylim([0 235])

        % Add colorbar
        c = colorbar;
        c.Label.String = '[A.U.]';
        c.Location = 'SouthOutside';
        c.FontSize = 6;
        colormap parula
        cmp = colormap;
        cmp = flipud(cmp);
        colormap(cmp);

        % Adjust view for better visualization
        view(0, 90);

        hold on
        for loopChans = 1:(48*8)
            plot3(0.05*(1:size(waveformRef,1))  + coordinates(loopChans,1) - 2.5, ...
                  0.02*waveformRef(:,loopChans) + coordinates(loopChans,2), ...
                  max(waveformRef(idx,:))*ones(size(waveformRef,1),1),'k','LineWidth',0.5)
        end
        hold off
        % title({['Mouse'],[region ' ' num2str(d,'%.0f') '\mum'],['-1 to 1ms spike'],['probability = ' num2str(p_sum,'%.2f')],newline})
        set(gca,'FontSize',5)
        set(gcf, 'Renderer', 'Painters');
        title({'ref','template'})

        %% tar

        % Create a 3D surface plot
        nexttile([5,1])
        tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
        trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values_tar, 'FaceColor', 'interp', 'EdgeColor', 'none');
        xlabel('[\mum]');
        ylabel('[\mum]');
    %     title('Activity Plot');

        axis equal

        % xlim([0 35])
        % ylim([0 235])

        % Add colorbar
        c = colorbar;
        c.Label.String = '[A.U.]';
        c.Location = 'SouthOutside';
        c.FontSize = 6;
        colormap parula
        cmp = colormap;
        cmp = flipud(cmp);
        colormap(cmp);

        % Adjust view for better visualization
        view(0, 90);

        hold on
        for loopChans = 1:(48*8)
            plot3(0.05*(1:size(waveformTar,1))  + coordinates(loopChans,1) - 2.5, ...
                  0.02*waveformTar(:,loopChans) + coordinates(loopChans,2), ...
                  max(waveformTar(idx,:))*ones(size(waveformTar,1),1),'k','LineWidth',0.5)
        end
        hold off
        % title({['Mouse'],[region ' ' num2str(d,'%.0f') '\mum'],['-1 to 1ms spike'],['probability = ' num2str(p_sum,'%.2f')],newline})
        set(gca,'FontSize',5)
        set(gcf, 'Renderer', 'Painters');
        title({'tar','template'})

    end
end

%% Figure 6 global snapshot

if figureCode == '6'
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    
    animalNames = {'N','S','U','AG1','AG2','Roy'};
    
    % % idx = 23 tar
    % load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat')  
    % channelSet = {1:16, ...
    %               17:32, ... 
    %               33:48, ...
    %               49:64, ...
    %               65:80, ...
    %               81:96, ...
    %               97:112, ...
    %               113:128};
    % 
    % waveforms = pairGroupStatTable.tarWaveforms{23};
    % 
    % shank       = pairGroupStatTable.tarShank(23);
    % chan        = pairGroupStatTable.tarChannel(23);
    % cellID      = pairGroupStatTable.tarNeuronID(23);
    % cellType    = pairGroupStatTable.tarCellExplorerType{23};
    % 
    % titleStr    = {['rat ' animalNames{1} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    % 
    % tWave = pairGroupStatTable.tWave{1};
    % 
    % hcomb = figure(101);
    % 
    % res_type = 'QHD';
    % pos = [1720 2562 560 0.4*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    % arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    % 
    % tiledlayout(1,size(channelSet,2),'Padding','none','TileSpacing','tight')
    % for loopSets = 1:size(channelSet,2)
    % 
    %     nexttile(loopSets)
    %     plot(tWave,waveforms(channelSet{loopSets},:)');
    % 
    %     if loopSets == 1
    %         title(titleStr)
    %     else 
    %         title({[' sh. ' num2str(loopSets)],[]})
    %     end
    % 
    %     xlim([-1,1])
    %     set(gca,'FontSize',5)
    %     set(gca,'FontName','Arial')
    %     xlabel('[ms]');
    %     if loopSets == 1
    %         ylabel('[mV]');
    %     end
    %     set(gca, 'YDir','reverse')
    %     box off
    % 
    % end
    
    % idx = 9
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAna.mat') 
    
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
    
    waveforms = pairGroupStatTable.tarWaveforms{9};
    
    shank       = pairGroupStatTable.tarShank(9);
    chan        = pairGroupStatTable.tarChannel(9);
    cellID      = pairGroupStatTable.tarNeuronID(9);
    cellType    = pairGroupStatTable.tarCellExplorerType{9};
    
    titleStr    = {['rat ' animalNames{3} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    
    tWave = pairGroupStatTable.tWave{1};
    
    hcomb = figure(102);
    
    res_type = 'QHD';
    pos = [1720 2562 560 0.6*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    % tiledlayout(1,size(channelSet,2),'Padding','none','TileSpacing','tight')
    tiledlayout(2,8,'Padding','none','TileSpacing','tight')
    for loopSets = 1:8 % size(channelSet,2)
    
        nexttile(loopSets)
        plot(tWave,waveforms(channelSet{loopSets},:)');
        
        if loopSets == 1
            title(titleStr)
        else 
            title({[' sh. ' num2str(loopSets)],[]})
        end
        
        xlim([-1,1])
        ylim([-0.4 0.05])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        % set(gca, 'YDir','reverse')
        box off
    
    end
    
    for loopSets = 1:8 % size(channelSet,2)
    
        nexttile(loopSets+8)
        plot(tWave,waveforms(channelSet{loopSets},:)');
         
        xlim([-1,1])
        ylim([-0.05 0.05])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        % set(gca, 'YDir','reverse')
        box off
    
    end

    % idx = 17
    
    waveforms = pairGroupStatTable.tarWaveforms{17};
    
    shank       = pairGroupStatTable.tarShank(17);
    chan        = pairGroupStatTable.tarChannel(17);
    cellID      = pairGroupStatTable.tarNeuronID(17);
    cellType    = pairGroupStatTable.tarCellExplorerType{17};
    
    titleStr    = {['rat ' animalNames{3} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    
    tWave = pairGroupStatTable.tWave{1};
    
    hcomb = figure(103);
    
    res_type = 'QHD';
    pos = [1720 2562 560 0.6*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    % tiledlayout(1,size(channelSet,2),'Padding','none','TileSpacing','tight')
    tiledlayout(2,8,'Padding','none','TileSpacing','tight')
    for loopSets = 1:8 % size(channelSet,2)
    
        nexttile(loopSets)
        plot(tWave,waveforms(channelSet{loopSets},:)');
        
        if loopSets == 1
            title(titleStr)
        else 
            title({[' sh. ' num2str(loopSets)],[]})
        end
        
        xlim([-1,1])
        ylim([-0.8 0.5])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        % set(gca, 'YDir','reverse')
        box off
    
    end

    for loopSets = 1:8 % size(channelSet,2)
    
        nexttile(loopSets+8)
        plot(tWave,waveforms(channelSet{loopSets},:)');
        
        xlim([-1,1])
        ylim([-0.06 0.05])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        % set(gca, 'YDir','reverse')
        box off
    
    end
    
    % % idx = 24
    % 
    % waveforms = pairGroupStatTable.tarWaveforms{24};
    % 
    % shank       = pairGroupStatTable.tarShank(24);
    % chan        = pairGroupStatTable.tarChannel(24);
    % cellID      = pairGroupStatTable.tarNeuronID(24);
    % cellType    = pairGroupStatTable.tarCellExplorerType{24};
    % 
    % titleStr    = {['rat ' animalNames{3} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    % 
    % tWave = pairGroupStatTable.tWave{1};
    % 
    % hcomb = figure(104);
    % 
    % res_type = 'QHD';
    % pos = [1720 2562 560 0.8*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    % arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    % 
    % % tiledlayout(1,size(channelSet,2),'Padding','none','TileSpacing','tight')
    % tiledlayout(2,8,'Padding','none','TileSpacing','compact')
    % for loopSets = 1:8 % size(channelSet,2)
    % 
    %     nexttile(loopSets)
    %     plot(tWave,waveforms(channelSet{loopSets},:)');
    % 
    %     if loopSets == 1
    %         title(titleStr)
    %     else 
    %         title({[' sh. ' num2str(loopSets)],[]})
    %     end
    % 
    %     xlim([-1,1])
    %     set(gca,'FontSize',5)
    %     set(gca,'FontName','Arial')
    %     xlabel('[ms]');
    %     if loopSets == 1
    %         ylabel('[mV]');
    %     end
    %     set(gca, 'YDir','reverse')
    %     box off
    % 
    % end
    
    % idx = 4
    
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
    
    channelSet = {1:8, ...
                  9:16, ... 
                  17:24, ...
                  25:32, ...
                  33:40, ...
                  41:48, ...
                  48:56, ...
                  57:64};
    
    waveforms = pairGroupStatTable.tarWaveforms{4};
    
    shank       = pairGroupStatTable.tarShank(4);
    chan        = pairGroupStatTable.tarChannel(4);
    cellID      = pairGroupStatTable.tarNeuronID(4);
    cellType    = pairGroupStatTable.tarCellExplorerType{4};
    
    titleStr    = {['rat ' animalNames{6} ': ' num2str(cellID) cellType],['(sh: ' num2str(shank) ')'],[' sh. ' num2str(1)],[]};
    
    tWave = pairGroupStatTable.tWave{1};
    
    hcomb = figure(105);
    
    res_type = 'QHD';
    pos = [1720 2562 560 0.6*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    tiledlayout(2,size(channelSet,2),'Padding','none','TileSpacing','tight')
    for loopSets = 1:size(channelSet,2)
    
        nexttile(loopSets)
        plot(tWave,waveforms(channelSet{loopSets},:)');
        
        if loopSets == 1
            title(titleStr)
        else 
            title({[' sh. ' num2str(loopSets)],[]})
        end
        
        xlim([-1,1])
        ylim([-7 2])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        % set(gca, 'YDir','reverse')
        box off
    
    end
    
    for loopSets = 1:size(channelSet,2)
    
        nexttile(loopSets+8)
        plot(tWave,waveforms(channelSet{loopSets},:)');
        
        xlim([-1,1])
        ylim([-0.5 0.5])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        if loopSets == 1
            ylabel('[mV]');
        end
        % set(gca, 'YDir','reverse')
        box off
    
    end

end

%% Figure 7 spike sorting clusters
if figureCode == '7'
    
    traceType = 'raw'; % 'pca' or 'raw'
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    res_type = 'QHD';
%     pos = [1720 2562 2*560*0.4 2*420*0.4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    pos = [1720 2562 2*560*0.4 420*0.35]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
%     tiledlayout(2,2,'Padding','none','TileSpacing','compact')
    tiledlayout(1,3,'Padding','none','TileSpacing','compact')

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    pairGroupStatsTable = pairGroupStatTable;
    clear pairGroupStatTable

    ref = pairGroupStatsTable.tarSpikeTimes{2,1};
    tar = pairGroupStatsTable.refSpikeTimes{2,1}; 
    refCelltype = pairGroupStatsTable.tarCellExplorerType{2,1};
    tarCelltype = pairGroupStatsTable.refCellExplorerType{2,1};
    d = pairGroupStatsTable.pairDistance(2);
    region = pairGroupStatsTable.brainRegion{2,1};
    
    tarShank = 1;
    tarCh    = 1;
    
    nexttile(1)
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
              CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                                    'plot_output', get(fig_use, 'Number'), ...
                                    'njitter', njitter, 'alpha', alpha,...
                                    'for_grant', for_grant, 'plot_pointwiseBands', false, ...
                                    'lineWidth', 1, 'axisLabels', false, 'microsecFlag', false, ...
                                    'norm_flag', true);
    xlim([-1,1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    hold on
    scatter(tR(110)*1e3,ccgR(110,1,2)/length(ref),100,'ro','LineWidth',1)
    scatter(tR(114)*1e3,ccgR(114,1,2)/length(ref),100,'bo','LineWidth',1)
    hold off
    set(gcf, 'Renderer', 'Painters');
    box off

%     nexttile(3)
%     % plot(pairGroupStatsTable.tWave{2,1}*1000, ...
%     %      pairGroupStatsTable.tarWaveforms{2,1}(pairGroupStatsTable.refChannel(2),:) - ...
%     %      pairGroupStatsTable.tarWaveforms{2,1}(pairGroupStatsTable.refChannel(2),1),':','LineWidth',4,'color','#0072BD')
%     plot(pairGroupStatsTable.tWave{2,1}, ...
%          pairGroupStatsTable.tarWaveforms{2,1}(1,:) - ...
%          pairGroupStatsTable.tarWaveforms{2,1}(1,1),'-.','LineWidth',1,'color','#0072BD')
% 
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
% %     legend('Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     xlabel('[ms]');
%     ylabel('[mV]');
%     box off

    lagLimit = duration/2;

    % removing the synchronous spikes                
    [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(ref',tar',binSize,lagLimit);
    
    % peak
    refPeak = spikeTimesInBin{110}(:,1);
    tarPeak = spikeTimesInBin{110}(:,2);
    
    % trough
    refTrough = spikeTimesInBin{114}(:,1);
    tarTrough = spikeTimesInBin{114}(:,2);
    
%     tar = setdiff(tar,[tarPeak; tarTrough]);
    
    % inset
    % hAxes2 = axes('Position',[.2 .7 .1 .2]);
    % plot(tR*1e3,ccgR(:,1,2),'k')
    % hold on
    % scatter(tR(110)*1e3,ccgR(110,1,2),20,'ro','LineWidth',1)
    % scatter(tR(114)*1e3,ccgR(114,1,2),20,'bo','LineWidth',1)
    % hold off
    % set(hAxes2,'xtick',[])
    % set(hAxes2,'ytick',[])
    % set(hAxes2,'xlim',[-1 1])
    
    %% PCA traces
    
    if strcmp(traceType,'pca')
        
        fpass      = 300;

        preLength  = 36;
        postLength = 36;

        nPCsRecovery = 3;

        chanDataDir  = ['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/'];

        chanData = load([ chanDataDir 'ch' num2str(tarCh) 'highpass300hz.mat'],'data');

        noDemeanFlag = true;

        spikeTimeIndxCell = round((tar-chanData.data.onsetTime/1e6)*fs) + 1;
        [waveAvg,waveforms] = waveformAvg(double(chanData.data.channel), ... 
                                                  spikeTimeIndxCell, ... 
                                                  preLength,postLength,fpass,fs,false,noDemeanFlag);
        [coeff, score] = pca(waveforms,'NumComponents',nPCsRecovery);

        [~,idxPeak,~  ] = intersect(tar,tarPeak);
        [~,idxTrough,~] = intersect(tar,tarTrough);

        filterSpon   = score * coeff(randperm(size(coeff,1),10000),:)';
        filterPeak   = score * coeff(idxPeak,:)';
        filterTrough = score * coeff(idxTrough,:)';

    %     tiledlayout(1,2)

    %     nexttile; plot(filterSpon(:,1:100000),'k'); ylim([-6 2])
    %     nexttile; plot(filterPeak,'k');             ylim([-6 2])
    %     nexttile; plot(filterTrough,'k');           ylim([-6 2])

    %     nexttile
    %     patchline(pairGroupStatsTable.tWave{1,1}*1000, ... 
    %                  filterSpon(:,1), ...
    %                  'linewidth',2,'edgealpha',0.1)
    %     hold on 
    %     for i = 2:length(filterSpon)
    %         patchline(pairGroupStatsTable.tWave{1,1}*1000, ... 
    %                  filterSpon(:,i), ...
    %                  'linewidth',2,'edgealpha',0.1)
    %     end
    %     hold off
    %     ylim([-6 2])

%         nexttile(2)
        nexttile
        patchline(pairGroupStatsTable.tWave{1,1}, ... 
                     filterSpon(:,1), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:length(filterSpon)
            patchline(pairGroupStatsTable.tWave{1,1}, ... 
                     filterSpon(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:length(filterPeak)
            patchline(pairGroupStatsTable.tWave{1,1}, ... 
                     filterPeak(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor','r')
        end
        hold off
        ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

%         nexttile(4)
        nexttile
        patchline(pairGroupStatsTable.tWave{1,1}, ... 
                     filterSpon(:,1), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:length(filterSpon)
            patchline(pairGroupStatsTable.tWave{1,1}, ... 
                     filterSpon(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:length(filterTrough)
            patchline(pairGroupStatsTable.tWave{1,1}, ... 
                     filterTrough(:,i), ...
                     'linewidth',1,'edgealpha',0.2,'edgecolor','b') 
        end
        hold off
        ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        set(gcf, 'Renderer', 'Painters');
        box off

    %     legend('spon','peak','trough')
    
    %% snippets 
    elseif strcmp(traceType,'raw')
    
        % weird timeshift in Roy's data
        timeShift      = 514;

        chanData       = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarCh) 'highpass300hz.mat']);

    %     chanData       = chanData.data.channel;

        filterFlag     = false; 
        noDemeanFlag   = false;
        fpass          = 300;
        preLength      = 5*30;
        postLength     = 5*30;

        bx = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-behavior.mat');

    %     tarLat       = round(tar*fs)'       - bx.behavior.RoyMaze1.datFrame(2,1) - timeShift;
    %     tarPeakLat   = round(tarPeak*fs)'   - bx.behavior.RoyMaze1.datFrame(2,1) - timeShift;
    %     tarTroughLat = round(tarTrough*fs)' - bx.behavior.RoyMaze1.datFrame(2,1) - timeShift;

        tarLat       = round((tar      -chanData.data.onsetTime/1e6)*fs) + 1;
        tarPeakLat   = round((tarPeak  -chanData.data.onsetTime/1e6)*fs) + 1;
        tarTroughLat = round((tarTrough-chanData.data.onsetTime/1e6)*fs) + 1;

        tarLatRandPerm = sort(tarLat(randperm(length(tarLat),10000)));

        [tarSponWaveMean,     tarSponWaveforms]      = waveformAvg(chanData.data.channel,tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
        [tarPeakLatWaveMean,  tarPeakLatWaveforms]   = waveformAvg(chanData.data.channel,tarPeakLat,    preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
        [tarTroughLatWaveMean,tarTroughLatWaveforms] = waveformAvg(chanData.data.channel,tarTroughLat,  preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

%         nexttile(2)
        nexttile

        patchline((-preLength:postLength-1)*(30/1000), ... 
                         tarSponWaveforms(:,1), ...
                         'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(tarSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,i), ...
                     'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
        end
        for i = 1:100 %size(tarPeakLatWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPeakLatWaveforms(:,i), ...
                     'linewidth',0.01,'edgealpha',0.1,'edgecolor','r')
        end
        hold off
        ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off

%         nexttile(4)
        nexttile

        patchline((-preLength:postLength-1)*(30/1000), ... 
                         tarSponWaveforms(:,1), ...
                         'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(tarSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,i), ...
                     'linewidth',0.01,'edgealpha',0.1,'edgecolor',[.7 .7 .7])
        end
        for i = 1:100 %size(tarTroughLatWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarTroughLatWaveforms(:,i), ...
                     'linewidth',0.01,'edgealpha',0.1,'edgecolor','b')
        end
        hold off
        ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        box off
        set(gcf, 'Renderer', 'Painters');
    
    end
    
    %% cluster plots
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(103);

    res_type = 'QHD';
    pos = [1720 2562 3*560*0.33 2*420*0.5]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact')
    
    % peak
    tarUnitSpikeCountAllSess = 517090;
    plotParam = 'ind';
    synchSpikesFlag = true;
    refUnitSpikeTimes = refPeak;
    binOfInterestBackup = [];

    curateClustersForPlots(tarShank, ...
                           tarCh, ...
                           tarUnitSpikeCountAllSess, ...
                           plotParam, ...
                           synchSpikesFlag, ...
                           refUnitSpikeTimes, ...
                           'r',...
                           binOfInterestBackup)
    
    % trough                   
    tarShank = 1;
    tarCh    = 1;
    tarUnitSpikeCountAllSess = 517090;
    plotParam = 'ind';
    synchSpikesFlag = true;
    refUnitSpikeTimes = refTrough;
    binOfInterestBackup = [];

    curateClustersForPlots(tarShank, ...
                           tarCh, ...
                           tarUnitSpikeCountAllSess, ...
                           plotParam, ...
                           synchSpikesFlag, ...
                           refUnitSpikeTimes,...
                           'b',...
                           binOfInterestBackup)
                       
    set(gcf, 'Renderer', 'Painters');
end

%% Figure 8 Allen Institute cluster quality
if figureCode == '8'
    
    datapath = '/media/nasko/WD_BLACK31/BOTtemp/';

    connType = 'exq';

    animalsList = [715093703
                   719161530
                   721123822
                   732592105
                   737581020
                   739448407
                   742951821
                   743475441
                   744228101
                   746083955
                   750332458
                   750749662
                   751348571
                   754312389
                   754829445
                   755434585
                   756029989
                   757216464
                   757970808
                   758798717
                   759883607
                   760345702
                   760693773
                   761418226
                   762120172
                   762602078
                   763673393
                   766640955
                   767871931
                   768515987
                   771160300
                   771990200
                   773418906
                   774875821
                   778240327
                   778998620
                   779839471
                   781842082
                   786091066
                   787025148
                   789848216
                   791319847
                   793224716
                   794812542
                   797828357
                   798911424
                   799864342
                   816200189
                   819186360
                   819701982
                   821695405
                   829720705
                   831882777
                   835479236
                   839068429
                   839557629
                   840012044
                   847657808];

    %% Allen SNR all units
    
    refIDexqList                = [];
    tarIDexqList                = [];

    refSNRexqList               = [];
    tarSNRexqList               = [];

    refLratioExqList            = [];
    tarLratioExqList            = [];

    refDprimeExqList            = [];
    tarDprimeExqList            = [];

    refIsolationDistanceExqList = [];
    tarIsolationDistanceExqList = [];

    refIDeqxList                = [];
    tarIDeqxList                = [];

    refTypeEqxList              = {};
    tarTypeEqxList              = {};

    for loopAnimals = 1:58

        display(['animal: '  num2str(animalsList(loopAnimals))])

        %%
        if strcmp(connType,'exq')
            load([datapath 'groupStatsMouse' num2str(animalsList(loopAnimals))])
        elseif strcmp(connType,'GJ')
            load([datapath 'groupStatsMousePutativeGJ' num2str(animalsList(loopAnimals))])
        end

        pairGroupStatTable = pairGroupStatsTable;
        clear pairGroupStatsTable

         % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);

        load(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/neuronsMouse' num2str(animalsList(loopAnimals)) 'qualityMetrics'])

        for loopPairs = 1:size(pairGroupStatTable,1)

            pairGroupStatTable.refLratio(loopPairs)            = neurons.Lratio(neurons.unitID == pairGroupStatTable.refNeuronID(loopPairs));
            pairGroupStatTable.tarLratio(loopPairs)            = neurons.Lratio(neurons.unitID == pairGroupStatTable.tarNeuronID(loopPairs));

            pairGroupStatTable.refDprime(loopPairs)            = neurons.Dprime(neurons.unitID == pairGroupStatTable.refNeuronID(loopPairs));
            pairGroupStatTable.tarDprime(loopPairs)            = neurons.Dprime(neurons.unitID == pairGroupStatTable.tarNeuronID(loopPairs));

            pairGroupStatTable.refIsolationDistance(loopPairs) = neurons.isolationDistance(neurons.unitID == pairGroupStatTable.refNeuronID(loopPairs));
            pairGroupStatTable.tarIsolationDistance(loopPairs) = neurons.isolationDistance(neurons.unitID == pairGroupStatTable.tarNeuronID(loopPairs));

        end
        
        refIDexqList                = [refIDexqList; pairGroupStatTable.refNeuronID];
        tarIDexqList                = [tarIDexqList; pairGroupStatTable.tarNeuronID];

        refSNRexqList               = [refSNRexqList; pairGroupStatTable.refSNR];
        tarSNRexqList               = [tarSNRexqList; pairGroupStatTable.tarSNR];

        refLratioExqList            = [refLratioExqList; pairGroupStatTable.refLratio];
        tarLratioExqList            = [tarLratioExqList; pairGroupStatTable.tarLratio];

        refDprimeExqList            = [refDprimeExqList; pairGroupStatTable.refDprime];
        tarDprimeExqList            = [tarDprimeExqList; pairGroupStatTable.tarDprime];

        refIsolationDistanceExqList = [refIsolationDistanceExqList; pairGroupStatTable.refIsolationDistance];
        tarIsolationDistanceExqList = [tarIsolationDistanceExqList; pairGroupStatTable.tarIsolationDistance];

        refIDeqxList                = [refIDeqxList; pairGroupStatTable.refNeuronID];
        tarIDeqxList                = [tarIDeqxList; pairGroupStatTable.tarNeuronID];

        refTypeEqxList              = [refTypeEqxList; pairGroupStatTable.refCellExplorerType];
        tarTypeEqxList              = [tarTypeEqxList; pairGroupStatTable.tarCellExplorerType];

    end
    
    IDexqList                = [refIDexqList;                tarIDexqList];
    SNRexqList               = [refSNRexqList;               tarSNRexqList];
    LratioExqList            = [refLratioExqList;            tarLratioExqList];
    DprimeExqList            = [refDprimeExqList;            tarDprimeExqList];
    isolationDistanceExqList = [refIsolationDistanceExqList; tarIsolationDistanceExqList];
    cellTypeEqxList          = [refTypeEqxList;              tarTypeEqxList];

    [~,idxTemp] = unique([refIDeqxList; tarIDeqxList]);
    
    IDexqList                = IDexqList(idxTemp);
    SNRexqList               = SNRexqList(idxTemp);
    LratioExqList            = LratioExqList(idxTemp);
    DprimeExqList            = DprimeExqList(idxTemp);
    isolationDistanceExqList = isolationDistanceExqList(idxTemp); 
    cellTypeEqxList          = cellTypeEqxList(idxTemp); 

    % (sum(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'))/length(cellTypeEqxList))*100
    % (sum(strcmp(cellTypeEqxList,'i-narrow'))/length(cellTypeEqxList))*100

    %% Allen SNR exq units

    datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';

    IDlist                = [];
    SNRlist               = [];
    LratioList            = [];
    DprimeList            = [];
    isolationDistanceList = [];
    cellTypeList          = {};

    for loopAnimals = 1:58

        display(['animal: '  num2str(animalsList(loopAnimals))])

        load([datapath 'neuronsMouse' num2str(animalsList(loopAnimals)) 'qualityMetrics']);

        IDlist                = [IDlist;                   neurons.unitID];
        SNRlist               = [SNRlist;                  neurons.snr];
        LratioList            = [LratioList;               neurons.Lratio];
        DprimeList            = [DprimeList;               neurons.Dprime];
        isolationDistanceList = [isolationDistanceList;    neurons.isolationDistance];
        cellTypeList          = [cellTypeList;             neurons.putativeCellType];

    end
    
    % filter out exquisite units from all units
    [~,idxFilt] = setdiff(double(IDlist),IDexqList);
    
    IDlist                = double(IDlist(idxFilt));
    SNRlist               = SNRlist(idxFilt);
    LratioList            = LratioList(idxFilt);
    DprimeList            = DprimeList(idxFilt);
    isolationDistanceList = isolationDistanceList(idxFilt); 
    cellTypeList          = cellTypeList(idxFilt);

    % (sum(strcmp(cellTypeList,'p') | strcmp(cellTypeList,'i-wide'))/length(cellTypeList))*100
    % (sum(strcmp(cellTypeList,'i-narrow'))/length(cellTypeList))*100

    tiledlayout(4,2)

    %% Allen SNR

    % pyramids
    allUnits = SNRlist(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
    exqUnits = SNRexqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
    nexttile
    histogram(allUnits,'Normalization','probability','BinWidth',0.1)
    hold on 
    histogram(exqUnits,'Normalization','probability','BinWidth',0.1)
    hold off
    d = (mean(allUnits)-mean(exqUnits))/std(allUnits);
    xlim([0 10])
    ylabel('probability')
    title(['Allen Int. pyramids SNR, d = ' num2str(d)])
    legend('non-exq','exq')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % interneuron
    allUnits = SNRlist(   strcmp(cellTypeList,   'i-narrow'));
    exqUnits = SNRexqList(strcmp(cellTypeEqxList,'i-narrow'));
    nexttile
    histogram(allUnits,'Normalization','probability','BinWidth',0.1)
    hold on 
    histogram(exqUnits,'Normalization','probability','BinWidth',0.1)
    hold off
    d = (mean(allUnits)-mean(exqUnits))/std(allUnits);
    xlim([0 10])
    ylabel('probability')
    title(['Allen Int. interneurons SNR, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % [P,H] = ranksum(SNRlist,SNRexqList) 

    %% Allen Lratio

    % pyramids
    allUnits = LratioList(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
    exqUnits = LratioExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
    nexttile
    histogram(allUnits,'Normalization','probability','BinWidth',0.001)
    hold on 
    histogram(exqUnits,'Normalization','probability','BinWidth',0.001)
    hold off
    d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
    xlim([0 0.05])
    ylabel('probability')
    title(['Allen Int. pyramids L ratio, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % interneuron
    allUnits = LratioList(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
    exqUnits = LratioExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
    nexttile
    histogram(LratioList(   strcmp(cellTypeList,   'i-narrow')),'Normalization','probability','BinWidth',0.001)
    hold on 
    histogram(LratioExqList(strcmp(cellTypeEqxList,'i-narrow')),'Normalization','probability','BinWidth',0.001)
    hold off
    d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
    xlim([0 0.05])
    ylabel('probability')
    title(['Allen Int. interneurons L ratio, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')


    %% Allen Dprime

    % pyramids
    allUnits = DprimeList(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
    exqUnits = DprimeExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
    nexttile
    histogram(allUnits,'Normalization','probability','BinWidth',0.5)
    hold on 
    histogram(exqUnits,'Normalization','probability','BinWidth',0.5)
    hold off
    d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
    xlim([0 20])
    ylabel('probability')
    title(['Allen Int. pyramids D prime, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % interneuron
    allUnits = DprimeList(   strcmp(cellTypeList,   'i-narrow'));
    exqUnits = DprimeExqList(strcmp(cellTypeEqxList,'i-narrow'));
    nexttile
    histogram(allUnits,'Normalization','probability','BinWidth',0.5)
    hold on 
    histogram(exqUnits,'Normalization','probability','BinWidth',0.5)
    hold off
    d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
    xlim([0 20])
    ylabel('probability')
    title(['Allen Int. interneurons D prime, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    %% Allen isolation distance

    % pyramids
    allUnits = isolationDistanceList(   strcmp(cellTypeList,'p') |    strcmp(cellTypeList,'i-wide'));
    exqUnits = isolationDistanceExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
    nexttile
    histogram(allUnits,'Normalization','probability','BinWidth',10)
    hold on 
    histogram(exqUnits,'Normalization','probability','BinWidth',10)
    hold off
    d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
    xlim([0 500])
    ylabel('probability')
    title(['Allen Int. pyramids iso. distance, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    % interneuron
    allUnits = isolationDistanceList(   strcmp(cellTypeList,   'i-narrow'));
    exqUnits = isolationDistanceExqList(strcmp(cellTypeEqxList,'i-narrow'));
    nexttile
    histogram(isolationDistanceList(   strcmp(cellTypeList,   'i-narrow')),'Normalization','probability','BinWidth',10)
    hold on 
    histogram(isolationDistanceExqList(strcmp(cellTypeEqxList,'i-narrow')),'Normalization','probability','BinWidth',10)
    hold off
    d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
    xlim([0 500])
    ylabel('probability')
    title(['Allen Int. interneurons iso. distance, d = ' num2str(d)])
    % legend('non-exq units','exq units')
    box off
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

end

%% Figure 9 quartile split analysis
if figureCode == '9'

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);
    
    res_type = 'QHD';
%     pos = [1720 2562 2*560*0.4 2*420*0.4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    pos = [1720 2562 1.15*560 (1.65)*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    testType = 'noMinus';

%     tiledlayout(2,2,'Padding','none','TileSpacing','compact')
    tiledlayout(4,4,'Padding','compact','TileSpacing','compact')

    UtkuRunSpWaveAnaMedianSplit('troughAmplitude','1',8, testType)  % idx 8

    UtkuRunSpWaveAnaMedianSplit('troughAmplitude','2',71,testType) % idx 71

    UtkuRunSpWaveAnaMedianSplit('troughAmplitude','1',20,testType) % idx 20

    UtkuRunSpWaveAnaMedianSplit('troughAmplitude','1',24,testType) % idx 24
    
    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(103);
    
    res_type = 'QHD';
%     pos = [1720 2562 2*560*0.4 2*420*0.4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    pos = [1720 2562 1.15*560 (1.65)*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
%     tiledlayout(2,2,'Padding','none','TileSpacing','compact')
    tiledlayout(4,4,'Padding','compact','TileSpacing','compact')

    BapunRunSpWaveAnaMedianSplit('S','troughAmplitude',86,testType)

    BapunRunSpWaveAnaMedianSplit('S','troughAmplitude',90,testType)

    BapunRunSpWaveAnaMedianSplit('U','troughAmplitude',18,testType)

    BapunRunSpWaveAnaMedianSplit('U','troughAmplitude',40,testType)

    % runSpWaveAnaMedianSplit('troughAmplitude','allStates',2)
     
    % runSpWaveAnaMedianSplit('troughAmplitude','allStates',23) % idx 23

    % runSpWaveAnaMedianSplit('troughAmplitude','allStates',28)

    % UtkuRunSpWaveAnaMedianSplit('halfWidth','1',20)
    % 
    % UtkuRunSpWaveAnaMedianSplit('halfWidth','1',24)
    % 
    % runSpWaveAnaMedianSplit('halfWidth','allStates',2)
    
end

%% Figure S1 Allen Optogenetics
if strcmp(figureCode,'S1')
    BOTprocessOptotaggingData
end

%% simulation noise
if strcmp(figureCode,'S2')
    noiseAnaCCG
end

%% supplement axo-Axonic 
if strcmp(figureCode,'S3')

    caseNo = 1;
    axoAxonicAna(caseNo)

    caseNo = 2;
    axoAxonicAna(caseNo)

    caseNo = 3;
    axoAxonicAna(caseNo)

end 

%% Figure 4 alignment figure 2 (obsolete)

% if figureName == '4'
% 
%     % Monitor specific plot settings.
%     screensize = get(0,'screensize');
%     % initiate figure
%     hcomb = figure(102);
% 
%     res_type = 'QHD';
%     pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
%     arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% 
%     tiledlayout(3,1)
% 
%     load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
%     pairGroupStatsTable = pairGroupStatTable;
%     clear pairGroupStatTable
% 
%     % idx 34
%     ref = pairGroupStatsTable.refSpikeTimes{34,1};
%     tar = pairGroupStatsTable.tarSpikeTimes{34,1}; 
%     refCelltype = pairGroupStatsTable.tarCellExplorerType{34,1};
%     tarCelltype = pairGroupStatsTable.refCellExplorerType{34,1};
%     d = pairGroupStatsTable.pairDistance(34);
%     region = pairGroupStatsTable.brainRegion{34,1};
% 
%     nexttile
%     [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%               CCG_jitter(ref,tar,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
%                                     'plot_output', get(fig_use, 'Number'), ...
%                                     'njitter', njitter, 'alpha', alpha,...
%                                     'for_grant', for_grant, 'plot_pointwiseBands', false, ...
%                                     'lineWidth', 4, 'axisLabels', false, 'microsecFlag', false);
% 
%     xlim([-1,1])
%     ylims = get(gca,'ylim');
%     ylabel('Counts')
%     xlabel('Time [\musec]')
%     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
%     if any(GSPExc)
%         hold on;
%         plot(tR(GSPExc == 1)*1e6, 1*ylims(2), 'r^', 'lineWidth', 4);
%     end
%     if any(GSPInh)
%         hold on;
%         plot(tR(GSPInh == 1)*1e6, 1*ylims(2),'bv', 'lineWidth', 4);
%     end
%     set(gca,'FontSize',12)
% 
%     nexttile
%     plot(pairGroupStatsTable.tWave{34,1}*1000, ...
%          pairGroupStatsTable.refWaveforms{34,1}(pairGroupStatsTable.refChannel(34),:) - ...
%          pairGroupStatsTable.refWaveforms{34,1}(pairGroupStatsTable.refChannel(34),1),'LineWidth',4,'color',"#D95319")
%     hold on 
%     plot(pairGroupStatsTable.tWave{34,1}*1000, ...
%          pairGroupStatsTable.tarWaveforms{34,1}(pairGroupStatsTable.tarChannel(34),:) - ...
%          pairGroupStatsTable.tarWaveforms{34,1}(pairGroupStatsTable.tarChannel(34),1),'LineWidth',4,'color','#0072BD')
%     plot(pairGroupStatsTable.tWave{34,1}*1000, ...
%          pairGroupStatsTable.tarWaveforms{34,1}(pairGroupStatsTable.refChannel(34),:) - ...
%          pairGroupStatsTable.tarWaveforms{34,1}(pairGroupStatsTable.refChannel(34),1),':','LineWidth',4,'color','#0072BD')
% 
%     hold off
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
%     legend('Tar e-spike at tar max ch',...
%            'Ref e-spike at ref max ch',...
%            'Ref e-spike at tar max ch','Location','NorthWest')
%     set(gca,'FontSize',12)
%     xlabel('Time Lag [\musec]');
%     ylabel('voltage [mV]');
% 
%     nexttile 
%     plot(pairGroupStatsTable.tWave{34,1}*1000, ...
%          1000*(pairGroupStatsTable.refWaveAtTarShank{34,1} - ...
%          pairGroupStatsTable.refWaveAtTarShank{34,1}(:,1))*1,'LineWidth',4) 
%     set(gca, 'YDir','reverse')
%     xlim([-1,1])
%     set(gca,'FontSize',12)
%     xlabel('Time Lag [\musec]');
%     ylabel('voltage [mV]');
% 
% end

