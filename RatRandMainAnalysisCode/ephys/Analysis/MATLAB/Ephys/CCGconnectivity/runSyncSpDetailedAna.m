% all based around Maze session
% A is the reference (R)
% B is the target (T)
% AvB are synchronous spikes of A and B
% C is the possible input


clc; clear all;

%% run informed neuron analysis
%  This will have to be universal code among datasets

% constants
session_name = "RoyMaze1";
conn_type    = {'GapPairs','ExcPairs','InhPairs'};
jscale_name  = 1;
alpha_name   = 5;
duration     = 0.007;
nepochs_plot = 3;
nrows        = 4;
njitter      = 100;
alpha        = 0.05;
for_grant    = false;
fs           = 30000;
jscale       = jscale_name;
fig_use      = 102;

% pair data 
datapath = "/home/nasko/CUNY_Work_Han_Kanram_Nat/data/wake_new/";

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale_name,alpha_name,datapath);

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
    
    current_pair         = cell_pairs_mat_unique(idx,:);
    current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                (current_pair(2) == cell_pairs_mat(:,2)));
    
    % create pair informed spikes for each session (pre, maze, post)
    % constants
    lagReal = jit_var_out(1).tR;
    lagCCGrun = -105:105;
    binSize = 1/fs; % one sample for fine jitter
    binSizeSyncSpk = 1;
    
    LOI_select = zeros(211,1);
    LOI_select(106) = 1;
    
    A  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    B  = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    LOI_maze = jit_var_out(current_pair_indices(2)).GSPExc & LOI_select;
%     LOI_maze = jit_var_out(current_pair_indices(2)).GSPInh;
    
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(A,B,LOI_maze);
    synchNames = sprintf("pair%d%sv%d%s", ... 
                 jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                 jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type);
    syncData.(synchNames).MazeSyncSpikes    = SyncSp;
    syncData.(synchNames).MazeSyncSpikesAll = SyncSpBinAll;
    syncData.(synchNames).MazeSyncCCG       = SyncCCG;
    SyncSp = min(SyncSp,[],2);
    AvB = SyncSp;
    
    if isempty(AvB)
        continue
    end
    
    AshankID = jit_var_out(current_pair_indices(1)).cell1shank(1);
    BshankID = jit_var_out(current_pair_indices(1)).cell2shank(1);
    cellSweepShankID = pre_neurons.shankID';
    
    cellSweep = find((cellSweepShankID ~= AshankID) & (cellSweepShankID ~= BshankID));
    
    subPlotIdx = [1 5 9 13];
    subPlotIdx = repmat(subPlotIdx,1,20);
    subPlotIdx = subPlotIdx(1:length(cellSweep));
    
    % loop individual neurons
    for ii = 1:length(cellSweep)
        
        % lazy work around
        i = cellSweep(ii);
            
            % plot AvB
            [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
                  CCG_jitter(A,B,fs,binSize,duration,'jscale',jscale, ...
                            'plot_output', get(fig_use, 'Number'), 'subfig', subPlotIdx(ii), ...
                            'subplot_size', [nrows, 4], 'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);   
                        
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
            
            title(sprintf('%s - %s - %d%s sh %d v %d%s sh %d', ...
                jit_var_out(current_pair_indices(2)).epoch, ...
                conn_type{conn_type_idx(current_pair_indices(2))}, ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                AshankID, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                BshankID))
            
            % plot Cv(AvB)
            C = maze_neurons.spikeTimes{i}/1000; 

            [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
                  CCG_jitter(C,AvB,fs,binSize,duration,'jscale',jscale, ...
                            'plot_output', get(fig_use, 'Number'), 'subfig', subPlotIdx(ii) + 1, ...
                            'subplot_size', [nrows, 4], 'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
            
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
            
            title(sprintf('%d%s sh %d v (%d%s sh %d v %d%s sh %d)', ...
                pre_neurons.UID(i),cell_type(i), ...
                cellSweepShankID(i), ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                AshankID, ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                BshankID))
            
            % plot CvA

            [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
                  CCG_jitter(C,A,fs,binSize,duration,'jscale',jscale, ...
                            'plot_output', get(fig_use, 'Number'), 'subfig', subPlotIdx(ii) + 2, ...
                            'subplot_size', [nrows, 4], 'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
            
            ylims = get(gca,'ylim');
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
            
            title(sprintf('%d%s sh %d v %d%s sh %d', ...
                pre_neurons.UID(i),cell_type(i), ...
                cellSweepShankID(i), ...
                jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                AshankID))
            
            % plot CvB

            [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
                  CCG_jitter(C,B,fs,binSize,duration,'jscale',jscale, ...
                            'plot_output', get(fig_use, 'Number'), 'subfig', subPlotIdx(ii) + 3, ...
                            'subplot_size', [nrows, 4], 'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
            
            ylims = get(gca,'ylim');
            
            if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
            if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
            
            title(sprintf('%d%s sh %d v %d%s sh %d', ...
                pre_neurons.UID(i),cell_type(i), ...
                cellSweepShankID(i), ...
                jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                BshankID))
             
        if subPlotIdx(ii) == 13
                        
            saveFilename = sprintf("%s_%s_INH_sync_CCG_spikes_detailed_analysis_of_pair_%d%sv%d%s_v_allNeurons_jscale%d_alpha%d", ... 
                                    session_name,conn_type{conn_type_idx(current_pair_indices(2))}, ...
                                    jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                                    jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
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

