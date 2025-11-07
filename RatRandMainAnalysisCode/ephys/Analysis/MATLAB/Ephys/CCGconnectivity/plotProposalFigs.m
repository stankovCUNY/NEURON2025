close all; clear all; clc; 

%% hardcoding pairs for proposal plots
pairsAna = [45 3
            20 45
            91 119
            79 26];
        
ylimList = {[0   3500]
            [350 2000]
            [200 1100]
            [0   80]};

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
datapath   = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";
figurepath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/figures/";

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

% get all pairs
cell_pairs = {jit_var_out.cell_pair}';

% convert to matrix
cell_pairs_mat = [];
for i = 1:size(cell_pairs,1)
    cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
end

tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 

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
        
    %% make CCGs
    nexttile
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant);
                
%     ylims = get(gca,'ylim');
    ylims = ylimList{j};
    ylim(ylims)
    if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
    if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
            
    title(sprintf('%d%s v %d%s', ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type))
    
    if (j == 1 || j == 2); xlabel(''); end
    if (j == 2 || j == 4); ylabel(''); end
            
    make_figure_pretty(102)
%     set(gca,'FontSize',16)
               
end

save_file = fullfile(figurepath, 'Fig1');
set(fig_use,'PaperOrientation','landscape');
print(fig_use, save_file,'-dpdf',resolution_use,'-bestfit');

close 102

%% (NOTE: the indexing work here is a mess and redundant but it works!)

%maze only
idxMaze = cellfun(@(a) strmatch(a,{jit_var_out.epoch}),{'Maze'},'uniform',false);

jit_var_outMaze = jit_var_out(idxMaze{1});

matIdx = [];
for i = 1:size(jit_var_outMaze,2)
    matIdx = [matIdx; jit_var_outMaze(i).cell_pair];
end

% unique pairs only
[~,~,IC] = unique(sort(matIdx,2),'row','first');

matIdxNew = [];
for i = 1:max(IC)
    tmp = find(IC == i);
    matIdxNew = [matIdxNew; matIdx(tmp(1),:)];
end
    
matIdx = matIdxNew;

excHist = zeros(211,1);
inhHist = zeros(211,1);

pairCCHstack        = [];
pairCCHstackIDs     = {};
pairCCHstackIDschar = {};

for j = 1:size(matIdx,1)
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = matIdx(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    excHist = excHist + jit_var_out(current_pair_indices(2)).GSPExc;
    inhHist = inhHist + jit_var_out(current_pair_indices(2)).GSPInh;
    
    stack = zeros(211,1);
    stack(find(jit_var_out(current_pair_indices(2)).GSPExc)) =  1;
    stack(find(jit_var_out(current_pair_indices(2)).GSPInh)) = -1;
    
    if sum(jit_var_out(current_pair_indices(2)).GSPExc + ...
           jit_var_out(current_pair_indices(2)).GSPInh) ~= 0
        pairCCHstack = [pairCCHstack stack];
        pairCCHstackIDs{j} = current_pair;
    end
end
pairCCHstackIDs = pairCCHstackIDs(~cellfun('isempty',pairCCHstackIDs));

for i = 1:size(pairCCHstackIDs,2)
    matTemp = pairCCHstackIDs{i};
    pairCCHstackIDschar{i} = sprintf("%d %d",matTemp(1),matTemp(2));
end

excHist = excHist + fliplr(excHist);
inhHist = inhHist + fliplr(inhHist);

tR = jit_var_out(1).tR*1000;

[~,sortIdx] = sort(sum(abs(pairCCHstack)),'descend');
pairCCHstack = pairCCHstack(:,sortIdx);
pairCCHstackIDschar = pairCCHstackIDschar(sortIdx);

subMilliInh   = sum(tR(106:136).*inhHist(106:136))/sum(inhHist(106:136));
superMilliInh = sum(tR(136:end).*inhHist(136:end))/sum(inhHist(136:end));

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'compact'); 

nexttile

imagesc(pairCCHstack')
colormap bluewhitered
ylabel('Pair','FontSize',16)
xlabel('lag [ms]','FontSize',16)
set(gca, 'XTick', 1:30:211, 'XTickLabel', tR(1:30:211)        ,'FontSize',16) 
% set(gca, 'YTick', 1:46,     'YTickLabel', pairCCHstackIDschar ,'FontSize',10)

nexttile

bar(tR(106:end),excHist(106:end),'FaceColor',[0.6350 0.0780 0.1840])
hold on
bar(tR(106:end),inhHist(106:end),'FaceColor',[0 0.4470 0.7410])
hxl1 = xline(subMilliInh,  '--',sprintf('Sub Inh. Avg. %0.2f ms',subMilliInh));
hxl1.FontSize  = 14;
hxl1.LineWidth = 3;
hxl2 = xline(superMilliInh,'--',sprintf('Super Inh. Avg. %0.2f ms',superMilliInh));
hxl2.FontSize = 14;
hxl2.LineWidth = 3;
hold off
ylabel('Aggregate bin count')
xlabel('absolute lag [ms]')
set(gca,'FontSize',16)

% make_figure_pretty(1)

save_file = fullfile(figurepath, 'Fig2');
set(102,'PaperOrientation','landscape');
print(102, save_file,'-dpdf',resolution_use,'-bestfit');

%%

decompRefSpikeCCG