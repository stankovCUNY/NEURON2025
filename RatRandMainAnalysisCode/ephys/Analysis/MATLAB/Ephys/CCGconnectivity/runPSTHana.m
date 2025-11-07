clc; clear all

data_dir = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new";
session  = "RoyMaze1";

load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
load(fullfile(data_dir, 'wake-ripple.mat'));

mazeMarkers   = behavior.(session).time;
rippleMarkers = ripple.(session).time;

rippleMarkers = rippleMarkers((rippleMarkers(:,1) > mazeMarkers(2,1)) & (rippleMarkers(:,2) < mazeMarkers(2,2)),:);
rippleMarkers = rippleMarkers/(1e6); % raw data is saved in microseconds???

trialLength = 0.5;

%%

% all based around Maze session
% A is the reference (R)
% B is the target (T)
% AvB are synchronous spikes of A and B
% C is the possible input

%% run informed neuron analysis
%  This will have to be universal code among datasets

% constants
session_name = "RoyMaze1";
conn_type    = {'GapPairs','ExcPairs','InhPairs'};
jscale_name  = 1;
alpha_name   = 5;
duration     = 0.007;
nepochs_plot = 3;
nrows        = 3;
njitter      = 100;
alpha        = 0.05;
for_grant    = false;
fs           = 30000;
jscale       = jscale_name;
fig_use      = 102;

% pair data 
datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";

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

subPlotRow = 0;

%% CCGs

syncData = struct;

% loop pairs
for j = 5:size(cell_pairs_mat_unique,1)
    
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
    
    AvBcondA = SyncSp(:,1);
    AvBcondB = SyncSp(:,2);
    
    subplot(2,2,1)
    [N_A,Ntrial_A] = PSTH(A,rippleMarkers,trialLength,binSize,'k');
    hold on
    [N_AvBcondA,Ntrial_AvBcondA] = PSTH(AvBcondA,rippleMarkers,trialLength,binSize,'g');
    hold off
    
    subplot(2,2,2)
    [N_B,Ntrial_B] = PSTH(B,rippleMarkers,trialLength,binSize,'k');
    hold on
    [N_AvBcondB,Ntrial_AvBcondB] = PSTH(AvBcondB,rippleMarkers,trialLength,binSize,'g');
    hold off
        
%     tic
%     JPSTH_AvB = JPSTH(Ntrial_A,Ntrial_B);
%     toc
    
    subplot(2,2,3)
    imagesc(flipud(JPSTH_AvB))
    
    k = -105:105;
    for l = 1:length(k)
        
        diagCCG(l) = sum(diag(JPSTH_AvB,k(l)));
        
    end
    
    subplot(2,2,4)
    plot(k*1/(3e4),diagCCG)
    
    
%     AshankID = jit_var_out(current_pair_indices(1)).cell1shank(1);
%     BshankID = jit_var_out(current_pair_indices(1)).cell2shank(1);
%     cellSweepShankID = maze_neurons.shankID';
%     
%     cellSweep = find((cellSweepShankID ~= AshankID) & (cellSweepShankID ~= BshankID));
%    
%     % initiate neuron group
%     C_group = [];
%     C_cellCount = 0;
%     
%     % loop individual neurons to create group C
%     for ii = 1:length(cellSweep)
%         
%         CvAidx_GSPExc = [];
%         CvAidx_GSPInh = [];
%         CvBidx_GSPExc = [];
%         CvBidx_GSPInh = [];
%         
%         % lazy work around
%         i = cellSweep(ii);
%         C = maze_neurons.spikeTimes{i}/1000;
%                 
%         % A
%         CvAidx = find((maze_neurons.UID(i) == cell_pairs_mat(:,1)) & ...
%                       (current_pair(1)     == cell_pairs_mat(:,2)));
%         if ~isempty(CvAidx)
%             CvAidx_GSPExc = jit_var_out(CvAidx(2)).GSPExc;
%             CvAidx_GSPInh = jit_var_out(CvAidx(2)).GSPInh;
%         end
%         
%         % B
%         CvBidx = find((maze_neurons.UID(i) == cell_pairs_mat(:,1)) & ...
%                       (current_pair(2)     == cell_pairs_mat(:,2)));
%         if ~isempty(CvBidx)
%             CvBidx_GSPExc = jit_var_out(CvBidx(2)).GSPExc;
%             CvBidx_GSPInh = jit_var_out(CvBidx(2)).GSPInh;
%         end
%                 
%         if ((sum(CvAidx_GSPExc) ~= 0) || (sum(CvAidx_GSPInh) ~= 0)) || ((sum(CvBidx_GSPExc) ~= 0) || (sum(CvBidx_GSPInh) ~= 0))
%             % group significant third neurons 
%             C_group = [C_group; C];
%             C_cellCount = C_cellCount + 1;
%         end
%     end
%     
%     if isempty(C_group)
%         continue
%     end
%     
%     C_group = sort(C_group);
%         
%     % plot AvB
%     [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
%           CCG_jitter(A,B,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_output', get(fig_use, 'Number'), 'subfig', subPlotRow*5 + 1, 'plot_flag', true,...
%                     'subplot_size', [nrows, 5], 'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant);   
%                         
%     ylims = get(gca,'ylim');
%     if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
%     if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
%             
%     title(sprintf('%s - %s - %d%s sh %d v %d%s sh %d', ...
%         jit_var_out(current_pair_indices(2)).epoch, ...
%         conn_type{conn_type_idx(current_pair_indices(2))}, ...
%         jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
%         AshankID, ...
%         jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
%         BshankID))
%             
%     % plot Cv(AvBcondA)
%     [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
%           CCG_jitter(C_group,AvBcondA,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_output', get(fig_use, 'Number'), 'subfig', subPlotRow*5 + 2, 'plot_flag', true, ...
%                     'subplot_size', [nrows, 5], 'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant);
%             
%     ylims = get(gca,'ylim');
%     if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
%     if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
%             
%     title(sprintf('C (cell count = %d) v %d%s sh %d sync spikes', ...
%         C_cellCount, ...
%         jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
%         AshankID))
%             
%     % plot Cv(AvBcondB)
%     [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
%           CCG_jitter(C_group,AvBcondB,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_output', get(fig_use, 'Number'), 'subfig', subPlotRow*5 + 3, 'plot_flag', true, ...
%                     'subplot_size', [nrows, 5], 'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant);
%             
%     ylims = get(gca,'ylim');
%     if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
%     if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
%             
%     title(sprintf('C (cell count = %d) v %d%s sh %d sync spikes', ...
%         C_cellCount, ...
%         jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
%         BshankID))
%             
%     % plot CvA
%     [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
%           CCG_jitter(C_group,A,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_output', get(fig_use, 'Number'), 'subfig', subPlotRow*5 + 4, 'plot_flag', true, ...
%                     'subplot_size', [nrows, 5], 'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant);
%             
%     ylims = get(gca,'ylim');
%     if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
%     if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
%             
%     title(sprintf('C (cell count = %d) v %d%s sh %d', ...
%         C_cellCount, ...
%         jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
%         AshankID))
%             
%     % plot CvB
%     [GSPExc,GSPInh,~,~,~,tR,~,~,~,~] = ...
%           CCG_jitter(C_group,B,fs,binSize,duration,'jscale',jscale, ...
%                     'plot_output', get(fig_use, 'Number'), 'subfig', subPlotRow*5 + 5, 'plot_flag', true, ...
%                     'subplot_size', [nrows, 5], 'njitter', njitter, 'alpha', alpha,...
%                     'for_grant', for_grant);
%             
%     ylims = get(gca,'ylim');
%     if any(GSPExc); hold on; plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'r^'); end
%     if any(GSPInh); hold on; plot(tR(GSPInh == 1)*1000, 0.95*ylims(2), 'bv'); end
%             
%     title(sprintf('C (cell count = %d) v %d%s sh %d', ...
%         C_cellCount, ...
%         jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
%         BshankID))
%        
%     if (subPlotRow == 2) || (j == 123)
% 
%         saveFilename = sprintf("%s_sync_CCG_spikes_detailed_analysis_sigNetworkC_jscale%d_alpha%d", ... 
%                                 session_name,jscale_name,alpha_name);
%         saveFilename = char(saveFilename);
% 
%         printNK(saveFilename,'SuFiSyn','hfig',hcomb,'append',true);
% 
%         % reset plot
%         close 102
% 
%         hcomb = figure(102);
%         arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
%         
%         subPlotRow = 0;
%     else
%         subPlotRow = subPlotRow + 1;
%     end
    
    

end
