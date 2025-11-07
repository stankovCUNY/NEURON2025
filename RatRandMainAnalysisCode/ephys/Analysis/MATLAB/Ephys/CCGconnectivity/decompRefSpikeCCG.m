pairsAna = [45 3];

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

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 2660/4 1860/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
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
   
    A  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    B  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    LOI_maze = jit_var_out(current_pair_indices(2)).GSPExc;
    
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(A,B,LOI_maze);
    synchNames = sprintf("pair%d%sv%d%s", ... 
                 jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                 jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type);
    syncData.(synchNames).MazeSyncSpikes    = SyncSp;
    syncData.(synchNames).MazeSyncSpikesAll = SyncSpBinAll;
    syncData.(synchNames).MazeSyncCCG       = SyncCCG;
    
    if isempty(SyncSp)
        continue
    end
    
    AvBcondA = SyncSp(:,1); % sync reference spikes
    AvBcondB = SyncSp(:,2);
    
    % plot regular CCG 
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(A,B,fs,binSize,duration,'jscale',jscale, 'plot_flag', true, ...
                'plot_output', get(fig_use, 'Number'), 'subfig', 1, ...
                'subplot_size', [1, 1], 'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
    
    % sync Ref spikes
    [~,~,~,~,ccgRsync,~,~,~,~,~] = ...
      CCG_jitter(AvBcondA,B,fs,binSize,duration,'jscale',jscale, 'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), 'subfig', 1, ...
                'subplot_size', [1, 1], 'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
            
    % non-sync Ref spikes
    [~,~,~,~,ccgRnonsync,~,~,~,~,~] = ...
      CCG_jitter(setdiff(A,AvBcondA),B,fs,binSize,duration,'jscale',jscale, 'plot_flag', false, ...
                'plot_output', get(fig_use, 'Number'), 'subfig', 3, ...
                'subplot_size', [3, 1], 'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);        
            
    tRms = tR*1000;        
    
    hold on
%     line(tRms,ccgR(:,2,1),'color','k','LineWidth',3)
    line(tRms,ccgRsync(:,1,2),'color','r','LineWidth',2)
    line(tRms,ccgRnonsync(:,1,2),'color','b','LineWidth',2)
    ylabel('Count')
    xlim([tRms(1),tRms(end)])
    
    title("45 v 3")
    legend("","","","","Sync. ref. spikes","Non-sync. ref. spikes")
    set(gca,'FontSize',12)
    
    hold off
%     title(sprintf('%s - ACG %d%s (sh %d) - %s', ...
%                 jit_var_out(current_pair_indices(2)).epoch, ...
%                 jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
%                 n1shank, ...
%                 conn_type{conn_type_idx(current_pair_indices(2))}))
    
    
    save_file = fullfile(datapath, ['Exq_pairs_ACG_CCG_refDecomp' session_name]);
    set(fig_use,'PaperOrientation','landscape');
    print(fig_use, save_file,'-dpdf',resolution_use);
            
%     printNK(['waveformAna_' session_name '_jscale' num2str(jscale)],...
%                     'SuFiSyn', 'hfig', fig_use, 'append', true);
            
    reset plot
    close 102

%     hcomb = figure(102);
%     arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

    

