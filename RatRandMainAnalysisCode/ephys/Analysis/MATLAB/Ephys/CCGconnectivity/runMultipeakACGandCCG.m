close all; clear all; clc; 

%% hardcoding pairs of interest for waveform analysis

pairsAna = [
            3   44;
            20  45;
            20  79;
            34  79;
            79  91;
            79  118;
            91  118;
            44  79;
            91  109;
            73  91;
            119 91;
            45  4;
            45  19;
            79  38
            ];
       
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
                 
    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,9*duration,'jscale',jscale, 'plot_flag', true, ...
                'plot_output', get(fig_use, 'Number'), 'subfig', 1, ...
                'subplot_size', [3, 1], 'njitter', njitter, 'alpha', alpha,...
                'for_grant', for_grant);
    
    tRms = tR*1000;        
    
    tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact');
    
    nexttile       
    line(tRms,ccgR(:,1,1),'color','k','LineWidth',1)
    ylabel('Count')
    xlim([tRms(1),tRms(end)])
    title(sprintf('%s - ACG %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                n1shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
    set(gca,'FontSize',16)
    
    nexttile
    line(tRms,ccgR(:,1,2),'color','k','LineWidth',1)
    xlim([tRms(1),tRms(end)])
    title(sprintf('%s - CCG %d%s (sh %d) v %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                n1shank, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                n2shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
    set(gca,'FontSize',16)
    
    nexttile
    line(tRms,ccgR(:,2,1),'color','k','LineWidth',1)
    ylabel('Count')
    xlabel('Time Lag [ms]')
    xlim([tRms(1),tRms(end)])
    title(sprintf('%s - CCG %d%s (sh %d) v %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                n2shank, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                n1shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
    set(gca,'FontSize',16)
    
    nexttile
    line(tRms,ccgR(:,2,2),'color','k','LineWidth',1)
    xlabel('Time Lag [ms]')
    xlim([tRms(1),tRms(end)])        
    title(sprintf('%s - ACG %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                n2shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
    set(gca,'FontSize',16)
    
    save_file = fullfile(datapath, ['Exq_pairs_ACG_CCG_' session_name]);
    set(fig_use,'PaperOrientation','landscape');
    print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
            
%     printNK(['waveformAna_' session_name '_jscale' num2str(jscale)],...
%                     'SuFiSyn', 'hfig', fig_use, 'append', true);
            
    % reset plot
    close 102

    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

close 102


                         
                         
                         