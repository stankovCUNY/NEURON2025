close all; clear all; clc; 

%% hardcoding pairs of interest for waveform analysis

pairsAna = [82 91];
       
% constants
session_name   = 'RoyMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1;
alpha_name     = 5;
duration       = 0.007;
fs             = 30000;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
numWaves       = 100;
resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

preLength  = 36;
postLength = 36;

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
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

%% get channels from waveform files. Not computaitonal;y savvy but it's fastest solution
extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';
RoySession1 = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);

tWave = (-26:27)*(1/30);

%%

waveDataSpCounts = [];
for i = 1:size(RoySession1.times,2)
    waveDataSpCounts(i) = size(RoySession1.times{i},1);
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

    [chanDataAtCell1, onsetTime] = HiroLoadRawNSC(chCell1);
    [chanDataAtCell2, onsetTime] = HiroLoadRawNSC(chCell2);
    
%     [chanDataAtRandom, onsetTime] = HiroLoadRawNSC(40);
    
    spikeTime = (1/30)*(-preLength+1:postLength);
    
    pairStr = [num2str(cellID1) cell1type ' v ' ...
               num2str(cellID2) cell2type ];

    
    tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
        
    nexttile
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res2,res1,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant);
    
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
    
    % get waveforms
    spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
    spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;
    
    spikeTimeIndxCell1nosync = round((res1nosync-onsetTime/1e6)*fs) + 1; 
    spikeTimeIndxCell2nosync = round((res2nosync-onsetTime/1e6)*fs) + 1;
    
    spikeAvgMaxChanCell1 = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength)/1000;
    spikeAvgMaxChanCell2 = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2,preLength,postLength)/1000;
    
    spikeAvgOpsChanCell1 = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1,preLength,postLength)/1000;
    spikeAvgOpsChanCell2 = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2,preLength,postLength)/1000;
    
    spikeAvgOpsChanCell1nosync = waveformAvg(chanDataAtCell2,spikeTimeIndxCell1nosync,preLength,postLength)/1000;
    spikeAvgOpsChanCell2nosync = waveformAvg(chanDataAtCell1,spikeTimeIndxCell2nosync,preLength,postLength)/1000;
    
    clear chanDataAtCell1
    clear chanDataAtCell2
           
    ylims = get(gca,'ylim');
    
    if any(GSPExc)
        hold on;
        plot(tR(GSPExc == 1)*1000, 0.95*ylims(2),'k^');
    end

    if any(GSPInh)
        hold on;
        plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'kv');
    end
    
    ylim(ylims)
    xlim([-1.5 1.5])
    title(pairStr)

%     set(gca,'FontSize',16)
    
    nexttile
    plot(spikeTime,spikeAvgMaxChanCell2,'LineWidth',2,'Color','#0072BD')
    hold on
%     plot(spikeTime,spikeAvgOpsChanCell2,      '--','LineWidth',2,'Color','#0072BD')
    plot(spikeTime,spikeAvgOpsChanCell2nosync,':', 'LineWidth',2,'Color','#0072BD')
%     plot(spikeTime,spikeAvgOpsChanRandomCell1nosync,   'LineWidth',2,'Color','k')
    hold off
    set(gca, 'YDir','reverse')
    legend('Ref e-spike at ref max chan', ...
           'Ref e-spike at tar max chan (no sync spikes)')
    xlabel('Time [ms]')
    ylabel('revere order voltage [mV]')
    title(['Reference e-spike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
    xlim([-1.5 1.5])
    
%     nexttile
%     plot(spikeTime,spikeAvgMaxChanCell2,'LineWidth',2,'Color','#D95319')
%     hold on
%     plot(spikeTime,spikeAvgOpsChanCell2,      '--','LineWidth',2,'Color','#D95319')
%     plot(spikeTime,spikeAvgOpsChanCell2nosync,':', 'LineWidth',2,'Color','#D95319')
% %     plot(spikeTime,spikeAvgOpsChanRandomCell2nosync,   'LineWidth',2,'Color','k')
%     hold off
%     set(gca, 'YDir','reverse','XDir','reverse')
%     legend('Tar e-spike at tar max chan', ...
%            'Tar e-spike at ref max chan', ...
%            'Tar e-spike at ref max chan (no sync spikes)', ...
%            'Tar e-spike at chan 40 shank 5')
%     xlabel('Time revere order [ms]')
%     ylabel('revere order voltage [mV]')
%     title(['Target e-spike: ' num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'])
%     xlim([-3.5 3.5])
    
    save_file = fullfile(data_dir, ['waveformAnaV2_' session_name '_jscale' num2str(jscale) '_' ,...
        num2str(cellID1) cell1type '_v_' num2str(cellID2) cell2type]);
    set(fig_use,'PaperOrientation','landscape');
%     print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
    print(fig_use, save_file,'-dpng',resolution_use);
            
    % reset plot
    close 102

    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

close 102


                         
                         
                         