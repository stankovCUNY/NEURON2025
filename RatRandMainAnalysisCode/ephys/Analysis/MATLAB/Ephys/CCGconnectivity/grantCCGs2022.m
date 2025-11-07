close all; clear all; clc; 

%% hardcoding pairs of interest for waveform analysis

pairsAna = [91 118; 16 38];


% constants
session_name   = 'RoyMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1; % in samples
alpha_name     = 5;
duration       = 0.025;
fs             = 30000;
binSize        = (1/fs)*5;
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

%% get session maze time markers
spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
[data_dir, name, ~] = fileparts(spike_data_fullpath);
load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');


for j = 1:size(pairsAna,1)
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = pairsAna(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
   
    n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
    n2shank = jit_var_out(current_pair_indices(2)).cell2shank;
                 
    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    
%     tiledlayout(3,1, 'Padding', 'none', 'TileSpacing', 'compact');
%     nexttile
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',5, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant);
            
    ylims = get(gca,'ylim');
    
    if any(GSPExc)
        hold on;
        plot(tR(GSPExc == 1)*1000, 0.99*ylims(2),'r^');
    end

    if any(GSPInh)
        hold on;
        plot(tR(GSPInh == 1)*1000, 0.99*ylims(2),'bv');
    end
    
    ylim(ylims)
    
    title(sprintf('%s - %d%s (sh %d) v %d%s (sh %d) - %s', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                n1shank, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                n2shank, ...
                conn_type{conn_type_idx(current_pair_indices(2))}))
    set(gca,'FontSize',16)
    
    save_file = fullfile(data_dir, ['grantCCGs_Feb2022_' session_name '_jscale' '5']);
    set(fig_use,'PaperOrientation','landscape');
    print(fig_use, save_file,'-dpsc',resolution_use,'-append','-bestfit');
            
    % reset plot
    close 102

    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
end

close 102
       