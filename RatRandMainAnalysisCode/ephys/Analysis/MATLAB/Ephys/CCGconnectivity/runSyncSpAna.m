clc; clear all;

%% run informed neuron analysis
%  This will have to be universal code among datasets

% constants
session_name = "RoyMaze1";
conn_type    = {'GapPairs','ExcPairs','InhPairs'};
jscale_name  = 1;
alpha_name   = 5;
duration     = 0.007;

% pair data 
datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale_name,alpha_name,datapath);

% constants
nepochs_plot = 3;
nrows        = 4;
njitter      = 100;
alpha        = 0.05;
for_grant    = false;
fs           = 30000;
jscale       = jscale_name;
fig_use      = 102;

% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

% get all pairs
cell_pairs = {jit_var_out.cell_pair}';

% convert to matrix
cell_pairs_mat = [];
for i = 1:size(cell_pairs,1)
    cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
end
cell_pairs_mat_unique = unique(cell_pairs_mat,'rows');

% individual neurons
[pre_neurons, maze_neurons, post_neurons, cell_type] = saveAllNeuronData;


%% CCGs

syncData = struct;

% loop pairs
for j = 1:size(cell_pairs_mat_unique,1)
    
    % reset plot for new pair 
    close 102    
    hcomb = figure(102);
    
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    idx = j; % pair to run, right now we are working with pre-maze time period
    
    current_pair = cell_pairs_mat_unique(idx,:);
    
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    % create pair informed spikes for each session (pre, maze, post)
    % constants
    lagReal = jit_var_out(1).tR;
    lagCCGrun = -105:105;
%     lagCCGrun = (-20:20)*15;
    binSize = 1/fs; % one sample for fine jitter
%     binSize = 15; % one sample for fine jitter
    binSizeSyncSpk = 1;

    % prep pair spike times
    n1pre   = jit_var_out(current_pair_indices(1)).cell1_spike_times;
    n2pre   = jit_var_out(current_pair_indices(1)).cell2_spike_times;
%     LOI_pre  = jit_var_out(current_pair_indices(1)).GSPExc | ...
%                jit_var_out(current_pair_indices(1)).GSPInh;
    LOI_pre  = jit_var_out(current_pair_indices(1)).GSPExc;

    n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
%     LOI_maze = jit_var_out(current_pair_indices(2)).GSPExc | ...
%                jit_var_out(current_pair_indices(2)).GSPInh;
    LOI_maze = jit_var_out(current_pair_indices(2)).GSPExc;

    n1post  = jit_var_out(current_pair_indices(3)).cell1_spike_times;
    n2post  = jit_var_out(current_pair_indices(3)).cell2_spike_times;
%     LOI_post = jit_var_out(current_pair_indices(3)).GSPExc | ...
%                jit_var_out(current_pair_indices(3)).GSPInh;
    LOI_post = jit_var_out(current_pair_indices(3)).GSPExc;

%     LOImidbins = zeros(211,1);
%     LOImidbins([104 105 106 107 108]) = 1;
%     LOI_pre  = LOImidbins & LOI_pre;
%     LOI_maze = LOImidbins & LOI_maze;
%     LOI_post = LOImidbins & LOI_post;
    
    tic
    
    % extract informed pair
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(n1pre, n2pre, LOI_pre);
    synchNames = sprintf("pair%d%sv%d%s", ... 
                 jit_var_out(current_pair_indices(1)).cell_pair(1),jit_var_out(current_pair_indices(1)).cell1type, ...
                 jit_var_out(current_pair_indices(1)).cell_pair(2),jit_var_out(current_pair_indices(1)).cell2type);
    syncData.(synchNames).PreSyncSpikes    = SyncSp;
    syncData.(synchNames).PreSyncSpikesAll = SyncSpBinAll;
    syncData.(synchNames).PreSyncCCG       = SyncCCG;
    SyncSp = min(SyncSp,[],2);
    res1_pre  = SyncSp;  
    
%     plot(jit_var_out(current_pair_indices(1)).tR,jit_var_out(current_pair_indices(1)).ccgR(:,1,2))
%     hold on
%     plot(jit_var_out(current_pair_indices(1)).tR,SyncCCG,'r')
%     hold off
    
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(n1maze,n2maze,LOI_maze);
    synchNames = sprintf("pair%d%sv%d%s", ... 
                 jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                 jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type);
    syncData.(synchNames).MazeSyncSpikes    = SyncSp;
    syncData.(synchNames).MazeSyncSpikesAll = SyncSpBinAll;
    syncData.(synchNames).MazeSyncCCG       = SyncCCG;
    SyncSp = min(SyncSp,[],2);
    res1_maze = SyncSp;
    
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(n1post,n2post,LOI_post);
    synchNames = sprintf("pair%d%sv%d%s", ... 
                 jit_var_out(current_pair_indices(3)).cell_pair(1),jit_var_out(current_pair_indices(3)).cell1type, ...
                 jit_var_out(current_pair_indices(3)).cell_pair(2),jit_var_out(current_pair_indices(3)).cell2type);
    syncData.(synchNames).PostSyncSpikes    = SyncSp;
    syncData.(synchNames).PostSyncSpikesAll = SyncSpBinAll;
    syncData.(synchNames).PostSyncCCG       = SyncCCG;
    SyncSp = min(SyncSp,[],2);
    res1_post = SyncSp;
    
    toc
    
    if (isempty(res1_pre) && isempty(res1_maze)) && isempty(res1_post)
        continue
    end
    
    
    synchNameFilename = sprintf("%s_%s_informed_spikes_CCGs_of_pair_%d%sv%d%s_data_jscale%d_alpha%d", ... 
                                session_name,conn_type{conn_type_idx(current_pair_indices(1))}, ...
                                jit_var_out(current_pair_indices(1)).cell_pair(1),jit_var_out(current_pair_indices(1)).cell1type, ...
                                jit_var_out(current_pair_indices(1)).cell_pair(2),jit_var_out(current_pair_indices(1)).cell2type, ...
                                jscale_name,alpha_name);
    
    if ~isfile(synchNameFilename)
        save([datapath + synchNameFilename],'syncData')   
    end
    
    cell1ShankID = jit_var_out(current_pair_indices(1)).cell1shank(1);
    cell2ShankID = jit_var_out(current_pair_indices(1)).cell2shank(1);
    cellSweepShankID = pre_neurons.shankID';
    
    cellSweep = find((cellSweepShankID ~= cell1ShankID) & (cellSweepShankID ~= cell2ShankID));
    
    sub_plot_row_idx = repmat(1:nrows,1,ceil(length(cellSweep)/nrows));
    sub_plot_row_idx = sub_plot_row_idx(1:length(cellSweep));
    
    % loop individual neurons
    for ii = 1:length(cellSweep)
        
        % lazy work around
        i = cellSweep(ii);
        
        % loop maze sessions
        for k = 1:3
            
            if k == 1
                res1 = res1_pre;
                res2 = pre_neurons.spikeTimes{i}/1000; % convert from seconds to milliseconds
            elseif k == 2
                res1 = res1_maze;
                res2 = maze_neurons.spikeTimes{i}/1000; 
            elseif k == 3
                res1 = res1_post;
                res2 = post_neurons.spikeTimes{i}/1000; 
            end
            
            [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                            'plot_output', get(fig_use, 'Number'), 'subfig', k + (sub_plot_row_idx(ii)-1)*nepochs_plot, ...
                            'subplot_size', [nrows, 3], 'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
            
            ylims = get(gca,'ylim');
            
            if any(GSPExc)
                hold on;
                plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^');
            end
            
            if any(GSPInh)
                hold on;
                plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'bv');
            end
                        
            title(sprintf('%s (%d%s sh %d v %d%s sh %d) v %d%s sh %d - %s', ...
                jit_var_out(current_pair_indices(k)).epoch, ...
                jit_var_out(current_pair_indices(k)).cell_pair(1),jit_var_out(current_pair_indices(k)).cell1type, ...
                cell1ShankID, ...
                jit_var_out(current_pair_indices(k)).cell_pair(2),jit_var_out(current_pair_indices(k)).cell2type, ...
                cell2ShankID, ...
                pre_neurons.UID(i),cell_type(i), ...
                cellSweepShankID(i), ...
                conn_type{conn_type_idx(current_pair_indices(k))}))
        end
                
        if (k == 3) && (sub_plot_row_idx(ii) == 4)
            
%             make_figure_pretty(102);
            
            saveFilename = sprintf("%s_%s_informed_spikes_CCGs_of_pair_%d%sv%d%s_v_allNeurons_jscale%d_alpha%d", ... 
                                    session_name,conn_type{conn_type_idx(current_pair_indices(k))}, ...
                                    jit_var_out(current_pair_indices(k)).cell_pair(1),jit_var_out(current_pair_indices(k)).cell1type, ...
                                    jit_var_out(current_pair_indices(k)).cell_pair(2),jit_var_out(current_pair_indices(k)).cell2type, ...
                                    jscale_name,alpha_name);
            saveFilename = char(saveFilename);
            
            printNK(saveFilename,'SuFiSyn','hfig',hcomb,'append',true);
            
            % reset plot
            close 102
            
            hcomb = figure(102);
            arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
        end
        
    end
end
% 
% %% ACGs
% 
% % sub_plot_row_idx = repmat(1:nrows,1,ceil(length(cell_pairs_mat_unique)/nrows));
% % sub_plot_row_idx = sub_plot_row_idx(1:length(cell_pairs_mat_unique));
% % 
% % %new plot
% % close 102    
% % hcomb = figure(102);
% % arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% % 
% % % loop pairs
% % for j = 1:size(cell_pairs_mat_unique,1)
% %        
% %     idx = j; % pair to run, right now we are working with pre-maze time period
% %     
% %     current_pair = cell_pairs_mat(idx,:);
% %     
% %     current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
% %                                 (current_pair(2) == cell_pairs_mat(:,2)));
% %     
% %     % create pair informed spikes for each session (pre, maze, post)
% %     % constants
% %     lagReal = jit_var_out(1).tR;
% %     lagCCGrun = -105:105;
% %     binSize = 1; % one sample for fine jitter
% % 
% %     % prep pair spike times
% %     n1_pre   = jit_var_out(current_pair_indices(1)).cell1_spike_times;
% %     n2_pre   = jit_var_out(current_pair_indices(1)).cell2_spike_times;
% %     LOI_pre  = jit_var_out(current_pair_indices(1)).GSPExc | ...
% %                jit_var_out(current_pair_indices(1)).GSPInh;
% % 
% %     n1_maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
% %     n2_maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
% %     LOI_maze = jit_var_out(current_pair_indices(2)).GSPExc | ...
% %                jit_var_out(current_pair_indices(2)).GSPInh;
% % 
% %     n1_post  = jit_var_out(current_pair_indices(3)).cell1_spike_times;
% %     n2_post  = jit_var_out(current_pair_indices(3)).cell2_spike_times;
% %     LOI_post = jit_var_out(current_pair_indices(3)).GSPExc | ...
% %                jit_var_out(current_pair_indices(3)).GSPInh;
% %     
% %     tic
% %     % extract informed pair
% %     [SyncSp,SyncCCG] = SyncSpikes(n1_pre, n2_pre, LOI_pre, binSize,lagCCGrun);
% %     res1_pre  = SyncSp;
% %     
% % %     plot(jit_var_out(idx).tR,jit_var_out(idx).ccgR(:,1,2))
% % %     hold on
% % %     plot(jit_var_out(idx).tR,SyncCCG,'r')
% % %     hold off
% %     
% %     [SyncSp,SyncCCG] = SyncSpikes(n1_maze,n2_maze,LOI_maze,binSize,lagCCGrun);
% %     res1_maze = SyncSp;
% %     [SyncSp,SyncCCG] = SyncSpikes(n1_post,n2_post,LOI_post,binSize,lagCCGrun);
% %     res1_post = SyncSp;
% %     toc
% %         
% %     % loop maze sessions
% %     for k = 1:3
% % 
% %         if k == 1
% %             res1 = res1_pre;
% %         elseif k == 2
% %             res1 = res1_maze;
% %         elseif k == 3
% %             res1 = res1_post;
% %         end
% % 
% %         duration = 0.007;
% % 
% %         [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
% %                 CCG_jitter(res1,res1,fs,binSize,duration,'jscale',jscale, ...
% %                           'plot_output', get(fig_use, 'Number'), 'subfig', k + (sub_plot_row_idx(idx)-1)*nepochs_plot, ...
% %                           'subplot_size', [nrows, 3], 'njitter', njitter, 'alpha', alpha,...
% %                           'for_grant', for_grant);
% % 
% %         title(sprintf('%s (%d%s v %d%s) v (%d%s v %d%s) - %s', ...
% %             jit_var_out(current_pair_indices(k)).epoch, ...
% %             jit_var_out(current_pair_indices(k)).cell_pair(1),jit_var_out(current_pair_indices(k)).cell1type, ...
% %             jit_var_out(current_pair_indices(k)).cell_pair(2),jit_var_out(current_pair_indices(k)).cell2type, ...
% %             jit_var_out(current_pair_indices(k)).cell_pair(1),jit_var_out(current_pair_indices(k)).cell1type, ...
% %             jit_var_out(current_pair_indices(k)).cell_pair(2),jit_var_out(current_pair_indices(k)).cell2type, ...
% %             conn_type{conn_type_idx(current_pair_indices(k))}))
% %         
% %         if (k == 3) && (sub_plot_row_idx(idx) == 4)
% % 
% %             make_figure_pretty(102);
% % 
% %             saveFilename = sprintf("%s_%s_informed_spikes_ACGs_jscale%d_alpha%d", ... 
% %                                     session_name,conn_type{conn_type_idx(current_pair_indices(k))}, ...
% %                                     jscale_name,alpha_name);
% %             saveFilename = char(saveFilename);
% % 
% %             printNK(saveFilename,'SuFiSyn','hfig',hcomb,'append',true);
% % 
% %             % reset plot
% %             close 102
% % 
% %             hcomb = figure(102);
% %             arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% %         end
% %     end
% % 
% %     
% % end
%     
% %%
% 
% % % [SyncSp,SyncCCG] = SyncSpikes(res1,res2,ones(211,1),1/30000,tR,ccgR(:,2,1));
% % 
% % ts1 = jit_var_out(idx).cell1_spike_times*1e3;
% % ts2 = jit_var_out(idx).cell2_spike_times*1e3;
% % 
% % [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(ts1, ts2, [-10 10]);
% % hist(tsOffsets, 211);
% % 
% % ts1bin = zeros(int64(max(ts2)),1);
% % ts2bin = zeros(int64(max(ts2)),1);
% % 
% % ts1bin(int64(ts1)) = 1;
% % ts2bin(int64(ts2)) = 1;
% % 
% % ts1bin = sparse(ts1bin);
% % ts2bin = sparse(ts2bin);
% % 
% % count = 0;
% % 
% % for i = 1:length(ts1bin)
% %     
% %     if ts1bin(i) == ts2bin(i)
% %         
% %         count = count + 1;
% %         
% %     end
% %     
% % end
% 
% %%
% 
% % binSize_Fs = 1;
% % halfBins   = 105;
% % 
% % fs = 20000;
% % 
% % times2  = round([n1; n2]*fs);
% % groups2 = [ones(length(n1),1); 2*ones(length(n2),1)];
% % 
% % [times2,sortidx] = sort(times2);
% % groups2 = groups2(sortidx);
% % 
% % counts = double(CCGHeart(times2,uint32(groups2),binSize_Fs,uint32(halfBins)));
% % 
% % plot(counts)








