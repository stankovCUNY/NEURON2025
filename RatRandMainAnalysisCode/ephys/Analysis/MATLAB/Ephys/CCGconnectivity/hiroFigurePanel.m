function hiroFigurePanel(pairAna,plot_pointwiseBands,plotType,flagFlipOrder)

    % constants
    session_name   = 'RoyMaze1';
    conn_type      = {'GapPairs','ExcPairs','InhPairs'};
    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    fpass          = 300;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    filterFlag     = false;
    Nwaveforms     = 100;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

    preLength  = 36;
    postLength = 36;

    shankChanList = {1:8  , ...
                     9:16 , ...
                     17:24, ...
                     25:32, ...
                     33:40, ...
                     41:48, ...
                     49:56, ...
                     57:64};

%     % Monitor specific plot settings.
%     screensize = get(0,'screensize');
%     % initiate figure
%     hcomb = figure(102);
% 
%     res_type = 'QHD';
%     pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
%     arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    % pair data 
    datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";

    [jit_var_out,~] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

    % get all pairs
    cell_pairs = {jit_var_out.cell_pair}';

    % convert to matrix
    cell_pairs_mat = [];
    for i = 1:size(cell_pairs,1)
        cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
    end

    %% get session maze time markers
    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    [data_dir, name, ~] = fileparts(spike_data_fullpath);
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');

    %% get channels from waveform files. Not computaitonal;y savvy but it's fastest solution
    extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
    RoySession1 = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);

    tWave = (-26:27)*(1/30);

    %%

    waveDataSpCounts = [];
    for i = 1:size(RoySession1.times,2)
        waveDataSpCounts(i) = size(RoySession1.times{i},1);
    end

    current_pair = pairAna;

    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));

    n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;

    n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
    n2shank = jit_var_out(current_pair_indices(2)).cell2shank;

    numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                     length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                     length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

    numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                     length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                     length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

    waveIdxCell1 = find(numSpikesCell1 == waveDataSpCounts);
    waveIdxCell2 = find(numSpikesCell2 == waveDataSpCounts);

    chCell1 = RoySession1.maxWaveformCh(waveIdxCell1);
    chCell2 = RoySession1.maxWaveformCh(waveIdxCell2);

    cell1type = jit_var_out(current_pair_indices(2)).cell1type;
    cell2type = jit_var_out(current_pair_indices(2)).cell2type;

    cellID1 = current_pair(1);
    cellID2 = current_pair(2);

    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;

    [chanDataAtCell1, onsetTime] = HiroLoadRawNSC(session_name,chCell1);
    [chanDataAtCell2, onsetTime] = HiroLoadRawNSC(session_name,chCell2);

%     [chanDataAtRandom, onsetTime] = HiroLoadRawNSC(40);

    spikeTime = (1/30)*(-preLength+1:postLength);

    if strcmp(plotType,'ccgOnly')
       
       if flagFlipOrder
           
           temp1 = res1;
           temp2 = res2;
           
           res1 = temp2;
           res2 = temp1;
           
           clear temp1 temp2
           
           pairStr = [session_name ' - '  num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')' ' v ' ...
                                          num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'];
           
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

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;

        % randomly selected spikes
        spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;

        [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
        [spikeAvgMaxChanCell2, waveformsMaxCell2] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
        [spikeAvgOpsChanCell2, waveformsOpsCell2] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);


        clear chanDataAtCell1
        clear chanDataAtCell2

    %     spikeAvgOpsChanRandomCell1nosync = waveformAvg(chanDataAtRandom,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
    %     spikeAvgOpsChanRandomCell2nosync = waveformAvg(chanDataAtRandom,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        allChansAtCell1 = [shankChanList{n1shank}];
        allChansAtCell2 = [shankChanList{n2shank}];

        spikeAvgOpsShankCell1nosync = zeros(8,preLength+postLength);
        spikeAvgOpsShankCell2nosync = zeros(8,preLength+postLength);

        for k = 1:8

            [chanDataAtCell2temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell2(k));
            spikeAvgOpsShankCell1nosync(k,:) = waveformAvg(chanDataAtCell2temp,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);

        end
        clear chanDataAtCell2temp

        for k = 1:8

            [chanDataAtCell1temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell1(k));
            spikeAvgOpsShankCell2nosync(k,:) = waveformAvg(chanDataAtCell1temp,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        end
        clear chanDataAtCell1temp
        
        nexttile(3)
        plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',2,'Color','#0072BD')
        hold on
        plot(spikeTime,spikeAvgOpsChanCell1,      '--','LineWidth',2,'Color','#0072BD')
    %     plot(spikeTime,spikeAvgOpsChanCell1nosync,':', 'LineWidth',2,'Color','#0072BD')
    %     plot(spikeTime,spikeAvgOpsChanRandomCell1nosync,   'LineWidth',2,'Color','k')
        hold off
        set(gca, 'YDir','reverse')
        legend('Ref e-spike at ref max chan', ...
               'Ref e-spike at tar max chan', ...
               'Ref e-spike at tar max chan (no sync spikes)', ...
               'Ref e-spike at chan 40 shank 5')
        xlabel('Time [ms]')
        ylabel('revere order voltage [mV]')
        title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
        xlim([-3.5 3.5])
        
    else
        
        pairStr = [session_name ' - '  num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                                       num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];
        
        tiledlayout(9,1, 'Padding', 'none', 'TileSpacing', 'compact'); 

        nexttile

        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', true, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,'plot_pointwiseBands',plot_pointwiseBands);

        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(211,1));
        
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
        title(pairStr)

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

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;

        % randomly selected spikes
        spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;

        [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
        [spikeAvgMaxChanCell2, waveformsMaxCell2] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
        [spikeAvgOpsChanCell2, waveformsOpsCell2] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);


        clear chanDataAtCell1
        clear chanDataAtCell2

    %     spikeAvgOpsChanRandomCell1nosync = waveformAvg(chanDataAtRandom,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
    %     spikeAvgOpsChanRandomCell2nosync = waveformAvg(chanDataAtRandom,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        allChansAtCell1 = [shankChanList{n1shank}];
        allChansAtCell2 = [shankChanList{n2shank}];

        spikeAvgOpsShankCell1nosync = zeros(8,preLength+postLength);
        spikeAvgOpsShankCell2nosync = zeros(8,preLength+postLength);

        for k = 1:8

            [chanDataAtCell2temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell2(k));
            spikeAvgOpsShankCell1nosync(k,:) = waveformAvg(chanDataAtCell2temp,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);

        end
        clear chanDataAtCell2temp

        for k = 1:8

            [chanDataAtCell1temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell1(k));
            spikeAvgOpsShankCell2nosync(k,:) = waveformAvg(chanDataAtCell1temp,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);

        end
        clear chanDataAtCell1temp

        nexttile
        plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',2,'Color','#0072BD')
        hold on
        plot(spikeTime,spikeAvgOpsChanCell1,      '--','LineWidth',2,'Color','#0072BD')
    %     plot(spikeTime,spikeAvgOpsChanCell1nosync,':', 'LineWidth',2,'Color','#0072BD')
    %     plot(spikeTime,spikeAvgOpsChanRandomCell1nosync,   'LineWidth',2,'Color','k')
        hold off
        set(gca, 'YDir','reverse')
        legend('Ref e-spike at ref max chan', ...
               'Ref e-spike at tar max chan', ...
               'Ref e-spike at tar max chan (no sync spikes)', ...
               'Ref e-spike at chan 40 shank 5')
        xlabel('Time [ms]')
        ylabel('revere order voltage [mV]')
        title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
        xlim([-3.5 3.5])

        nexttile
        patchline(spikeTime,waveformsMaxCell1(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
        hold on
        for k = 2:Nwaveforms
            patchline(spikeTime,waveformsMaxCell1(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'EdgeAlpha',0.25);
        end
        hold off
        set(gca, 'YDir','reverse')
        xlabel('Time [ms]')
        ylabel('revere order voltage [mV]')
        xlim([-3.5 3.5])

        nexttile
        patchline(spikeTime,waveformsOpsCell1(:,1),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
        hold on
        for k = 2:Nwaveforms
            patchline(spikeTime,waveformsOpsCell1(:,k),'EdgeColor','#0072BD','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
        end
        hold off
        set(gca, 'YDir','reverse')
        xlabel('Time [ms]')
        ylabel('revere order voltage [mV]')
        xlim([-3.5 3.5])

        nexttile
    %     plot(spikeTime,spikeAvgOpsChanCell1nosync + flipud(spikeAvgOpsChanCell2nosync),'LineWidth',2,'Color','#0072BD')
        plot(spikeTime,spikeAvgOpsShankCell1nosync')
        set(gca, 'YDir','reverse')
        legend(cellstr(num2str(shankChanList{n2shank}', 'ch%-d')))
        xlabel('Time [ms]')
        ylabel('revere order voltage [mV]')
        xlim([-3.5 3.5])
        title(['Reference e-spike: ' num2str(cellID1) cell1type ' on all ' num2str(cellID2) cell1type '''s shank ' num2str(n2shank) ' channels'])

        nexttile
        plot(spikeTime,spikeAvgMaxChanCell2,'LineWidth',2,'Color','#D95319')
        hold on
        plot(spikeTime,spikeAvgOpsChanCell2,      '--','LineWidth',2,'Color','#D95319')
    %     plot(spikeTime,spikeAvgOpsChanCell2nosync,':', 'LineWidth',2,'Color','#D95319')
    %     plot(spikeTime,spikeAvgOpsChanRandomCell2nosync,   'LineWidth',2,'Color','k')
        hold off
        set(gca, 'YDir','reverse','XDir','reverse')
        legend('Tar e-spike at tar max chan', ...
               'Tar e-spike at ref max chan', ...
               'Tar e-spike at ref max chan (no sync spikes)', ...
               'Tar e-spike at chan 40 shank 5')
        xlabel('Time revere order [ms]')
        ylabel('revere order voltage [mV]')
        title(['Target e-spike: ' num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'])
        xlim([-3.5 3.5])

        nexttile
        patchline(spikeTime,waveformsMaxCell2(:,1),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
        hold on
        for k = 2:Nwaveforms
            patchline(spikeTime,waveformsMaxCell2(:,k),'EdgeColor','#D95319','LineWidth',0.5,'EdgeAlpha',0.25);
        end
        hold off
        set(gca, 'YDir','reverse','XDir','reverse')
        xlabel('Time revere order [ms]')
        ylabel('revere order voltage [mV]')
        xlim([-3.5 3.5])

        nexttile
        patchline(spikeTime,waveformsOpsCell2(:,1),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
        hold on
        for k = 2:Nwaveforms
            patchline(spikeTime,waveformsOpsCell2(:,k),'EdgeColor','#D95319','LineWidth',0.5,'LineStyle','--','EdgeAlpha',0.25);
        end
        hold off
        set(gca, 'YDir','reverse','XDir','reverse')
        xlabel('Time revere order [ms]')
        ylabel('revere order voltage [mV]')
        xlim([-3.5 3.5])

        nexttile
        plot(spikeTime,spikeAvgOpsShankCell2nosync')
        set(gca, 'YDir','reverse','XDir','reverse')
        legend(cellstr(num2str(shankChanList{n1shank}', 'ch%-d')))
        ylabel('revere order voltage [mV]')
        xlim([-3.5 3.5])
        title(['Target e-spike: '    num2str(cellID2) cell2type ' on all ' num2str(cellID1) cell1type '''s shank ' num2str(n1shank) ' channels'])
    end
        
%     save_file = fullfile(data_dir, ['waveformAnaV2_' session_name '_jscale' num2str(jscale) '_' ,...
%         num2str(cellID1) cell1type '_v_' num2str(cellID2) cell2type]);
% %     set(fig_use,'PaperOrientation','landscape');
% %     print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
%     print(fig_use, save_file,'-dpng',resolution_use);
% 
%     % reset plot
%     close 102
% 
%     hcomb = figure(102);
%     arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

end

                         
                         
                         