close all; clear all; clc; 

%% hardcoding pairs of interest for waveform analysis

% constants
session_name   = 'RoyMaze1';
% session_name   = 'KevinMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1;
alpha_name     = 5;
duration       = 0.01;
fpass          = 300;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
filterFlag     = false;
Nwaveforms     = 1000;
resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using

preLength  = 36;
postLength = 36;

if strcmp(session_name,'RoyMaze1')
    
    fs        = 30000;
    binSize   = 1/fs;
    spikeTime = (1/30)*(-preLength+1:postLength);
    pairsAna = [3   44;             20  45;             20  79;             34  79;
                79  91;             79  118;            91  118;            44  79;
                91  109;            73  91;             119 91;             45  4;
                45  19;             79  38;             3   45;             4   79;
                18  45;             32  79;             79  38;             3   38;
                91  101;            45  12;             82  91;             3   25;
                70  91;             20  34;             20  44;             40  79;
                79  58;             3   32;             4   25;             20  38;
                46  79;             79  87;             87  119;            98  119;
                119 122;            45  15;             17  45;             79  19;
                91  68;             79  95;             91  110;            45  79];

elseif strcmp(session_name,'KevinMaze1')
    
    fs        = 32000;
    binSize  = 1/fs;
    spikeTime = (1/32)*(-preLength+1:postLength);
    pairsAna = [70  22;             70  54;             70  50;             3   22];
    
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

res_type = 'QHD';
% pos = [70 230 2660/4 1860/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
% pos = [70 230 2660/4 1860*(3/4)]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
% pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
pos = [70 230 1860 2660]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';

arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

% pair data 
datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

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
load(fullfile(data_dir, 'wake-spikes.mat'));


%% get channels from waveform files. Not computaitonal;y savvy but it's fastest solution
extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';

if strcmp(session_name,'RoyMaze1')
    Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);
elseif strcmp(session_name,'KevinMaze1')
    Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Kevin-maze1'],'noPrompts',true);
end

%%

waveDataSpCounts = [];
for i = 1:size(Session.times,2)
    waveDataSpCounts(i) = size(Session.times{i},1);
end

for j = 1:size(pairsAna,1)
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = pairsAna(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
   
    n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
    n2shank = jit_var_out(current_pair_indices(2)).cell2shank;
    
    if strcmp(session_name,'RoyMaze1')
        
        numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

        numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell2_spike_times);
                     
    elseif strcmp(session_name,'KevinMaze1')
        
        numSpikesCell1 = length(spikes.KevinMaze1(current_pair(1)).time);
        numSpikesCell2 = length(spikes.KevinMaze1(current_pair(2)).time);
        
    end
       
    waveIdxCell1 = find(numSpikesCell1 == waveDataSpCounts);
    waveIdxCell2 = find(numSpikesCell2 == waveDataSpCounts);
    
%     chCell1 = Session.maxWaveformCh(waveIdxCell1) + 1;
%     chCell2 = Session.maxWaveformCh(waveIdxCell2) + 1;
    
    cellID1 = current_pair(1);
    cellID2 = current_pair(2);

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')
    
    refIdx = find(unitsGlobalSnapshots.neuronID == cellID1);
    tarIdx = find(unitsGlobalSnapshots.neuronID == cellID2);

    refWave = unitsGlobalSnapshots.waveforms{refIdx};
    tarWave = unitsGlobalSnapshots.waveforms{tarIdx};
    
    [~,chCell1] = min(min(refWave'));
    [~,chCell2] = min(min(tarWave'));
    
    cell1type = jit_var_out(current_pair_indices(2)).cell1type;
    cell2type = jit_var_out(current_pair_indices(2)).cell2type;
    
    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
%     res2 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
%     res1 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    
%     [chanDataAtCell1, onsetTime] = HiroLoadRawNSC(session_name,chCell1);
%     [chanDataAtCell2, onsetTime] = HiroLoadRawNSC(session_name,chCell2);
    
    [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,chCell1);
    [chanDataAtCell2, onsetTime] = HiroLoad300hz(session_name,chCell2);
    
    %%
    
%     if strcmp(session_name,'KevinMaze1') % adjust for difference in sampling freq
%         res1 = ( (res1 - res1(1))*(32e3/30e3) ) + res1(1);
%         res2 = ( (res2 - res2(2))*(32e3/30e3) ) + res2(2);
%     end
    
%     [chanDataAtRandom, onsetTime] = HiroLoadRawNSC(40);
        
    pairStr = [session_name ' - '  num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                                   num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];

    
    tiledlayout(9,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
%      tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact'); 

    nexttile
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant,  'plot_pointwiseBands',false);
    
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(length(GSPExc),1),duration);
    
    res1sync = [];
    res2sync = [];
    
    for k = 1:length(GSPExc)
       
        if (GSPExc(k) == 1) || (GSPInh(k) == 1)
            res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
            res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
        end
            
    end

%     res1sync = SyncSpBinAll{147}(:,1);
%     res2sync = SyncSpBinAll{147}(:,2);
    
    % remove sync spikes for waveforms
    res1nosync = setdiff(res1,res1sync);
    res2nosync = setdiff(res2,res2sync);
    
    if Nwaveforms < length(res1nosync)
        res1nosync = res1nosync(randperm(length(res1nosync),Nwaveforms));
    end
    if Nwaveforms < length(res2nosync)
        res2nosync = res2nosync(randperm(length(res2nosync),Nwaveforms));
    end
    
%     if strcmp(session_name,'RoyMaze1')
%         
%        % get waveforms
%        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
%        spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;
%        
%        % randomly selected spikes
%        spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
%        spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;
%                      
%     elseif strcmp(session_name,'KevinMaze1')
%         
%        % get waveforms
%        spikeTimeIndxCell1 = round(res1*fs) + 1; 
%        spikeTimeIndxCell2 = round(res2*fs) + 1;
%        
%        % randomly selected spikes
%        spikeTimeIndxCell1nosync = round((res1nosync)*fs) + 1; 
%        spikeTimeIndxCell2nosync = round((res2nosync)*fs) + 1;
%         
%     end
   
%     onsetTime = behavior.KevinMaze1.time(1,2);
    
    % get waveforms
    spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
    spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;

    % randomly selected spikes
    spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
    spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;
    
    spikeTimeIndxCell1sync = round((res1sync-onsetTime/1e6)*fs) + 1; 
    spikeTimeIndxCell2sync = round((res2sync-onsetTime/1e6)*fs) + 1;
    
    [~, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
    [~, waveformsMaxCell2] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);

    [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
    [spikeAvgMaxChanCell2, waveformsMaxCell2] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
    
    [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
    [spikeAvgOpsChanCell2, waveformsOpsCell2] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
    
%     
    
%     [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
%     [spikeAvgMaxChanCell2, waveformsMaxCell2] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);
%     
%     [spikeAvgOpsChanCell1, waveformsOpsCell1] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
%     [spikeAvgOpsChanCell2, waveformsOpsCell2] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);
%     
%     [spikeAvgOpsChanCell1nosync, waveformsOpsCell1nosync] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
%     [spikeAvgOpsChanCell2nosync, waveformsOpsCell2nosync] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
    
    clear chanDataAtCell1
    clear chanDataAtCell2
    
%     spikeAvgOpsChanRandomCell1nosync = waveformAvg(chanDataAtRandom,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
%     spikeAvgOpsChanRandomCell2nosync = waveformAvg(chanDataAtRandom,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
    
    allChansAtCell1 = [shankChanList{n1shank}];
    allChansAtCell2 = [shankChanList{n2shank}];
    
    spikeAvgOpsShankCell1nosync = zeros(8,preLength+postLength);
    spikeAvgOpsShankCell2nosync = zeros(8,preLength+postLength);
    
    for k = 1:8
        
%         [chanDataAtCell2temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell2(k));
        [chanDataAtCell2temp, onsetTime] = HiroLoad300hz(session_name,allChansAtCell2(k));
        spikeAvgOpsShankCell1nosync(k,:) = waveformAvg(chanDataAtCell2temp,spikeTimeIndxCell1nosync,preLength,postLength,fpass,fs,filterFlag);
        clear chanDataAtCell2temp
        
    end
    
    for k = 1:8
        
%         [chanDataAtCell1temp, onsetTime] = HiroLoadRawNSC(session_name,allChansAtCell1(k));
        [chanDataAtCell1temp, onsetTime] = HiroLoad300hz(session_name,allChansAtCell1(k));
        spikeAvgOpsShankCell2nosync(k,:) = waveformAvg(chanDataAtCell1temp,spikeTimeIndxCell2nosync,preLength,postLength,fpass,fs,filterFlag);
        clear chanDataAtCell1temp
        
    end
    
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

%     set(gca,'FontSize',16)
    
    nexttile
    plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',2,'Color','#0072BD')
%     plot(spikeTime,spikeAvgOpsChanCell1,'LineWidth',2,'Color','#0072BD')
    hold on
    plot(spikeTime,spikeAvgOpsChanCell1,'--','LineWidth',2,'Color','#0072BD')
%     plot(spikeTime,spikeAvgMaxChanCell1,'--','LineWidth',2,'Color','#0072BD')
%     plot(spikeTime,spikeAvgOpsChanCell1nosync,':', 'LineWidth',2,'Color','#0072BD')
%     plot(spikeTime,spikeAvgOpsChanRandomCell1nosync,   'LineWidth',2,'Color','k')
    hold off
    set(gca, 'YDir','reverse')
    legend('Ref e-spike at ref max chan', ...79i v 20i (idx 3 reversed)
           'Ref e-spike at tar max chan', ...
           'Ref e-spike at tar max chan (no sync spikes)', ...
           'Ref e-spike at chan 40 shank 5')
    xlabel('Time [ms]')
    ylabel('revere order voltage [mV]')
    title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
    xlim([-3.5 3.5])
%     xlim([-1 1])
    
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
%     xlim([-1 1])
    
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
%     xlim([-1 1])
    
    nexttile
%     plot(spikeTime,spikeAvgOpsChanCell1nosync + flipud(spikeAvgOpsChanCell2nosync),'LineWidth',2,'Color','#0072BD')
    plot(spikeTime,spikeAvgOpsShankCell1nosync','LineWidth',2)
    set(gca, 'YDir','reverse')
    legend(cellstr(num2str(shankChanList{n2shank}', 'ch%-d')))
    xlabel('Time [ms]')
    ylabel('revere order voltage [mV]')
    xlim([-3.5 3.5])
%     xlim([-1 1])
    title(['Reference e-spike: ' num2str(cellID1) cell1type ' on all ' num2str(cellID2) cell1type '''s shank ' num2str(n2shank) ' channels'])
    
    nexttile
    plot(spikeTime,spikeAvgMaxChanCell2,'LineWidth',2,'Color','#D95319')
    hold on
    plot(spikeTime,spikeAvgOpsChanCell2,'--','LineWidth',2,'Color','#D95319')
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
%     xlim([-1 1])
    
    
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
%     xlim([-1 1])
    
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
%     xlim([-1 1])
    
    nexttile
    plot(spikeTime,fliplr(spikeAvgOpsShankCell2nosync'),'LineWidth',2)
    set(gca, 'YDir','reverse','XDir','reverse')
    set(gca,'YDir','reverse')
    legend(cellstr(num2str(shankChanList{n1shank}', 'ch%-d')))
    ylabel('revere order voltage [mV]')
     xlabel('Time [ms]')
    xlim([-3.5 3.5])
%     xlim([-1 1])
    title(['Target e-spike: '    num2str(cellID2) cell2type ' on all ' num2str(cellID1) cell1type '''s shank ' num2str(n1shank) ' channels'])
    
    save_file = fullfile(data_dir, ['waveformAnaV2_' session_name '_jscale' num2str(jscale) '_' ,...
        num2str(cellID1) cell1type '_v_' num2str(cellID2) cell2type]);
%     set(fig_use,'PaperOrientation','landscape');
%     print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
    print(fig_use, save_file,'-dpng',resolution_use);
            
    % reset plot
    close 102

    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

close 102


                         
                         
                         