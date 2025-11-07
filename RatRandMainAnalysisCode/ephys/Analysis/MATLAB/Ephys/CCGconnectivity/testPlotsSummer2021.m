% tic
% % extract informed pair
% [InfSp,InfoCCG] = InformedSpikes(n1_pre, n2_pre, LOI_pre, binSizeInfSpk,lagCCGrun);
% res1_pre  = InfSp;
% 
% %     plot(jit_var_out(idx).tR,jit_var_out(idx).ccgR(:,1,2))
% %     hold on
% %     plot(jit_var_out(idx).tR,InfoCCG,'r')
% %     hold off
% 
% [InfSp,InfoCCG] = InformedSpikes(n1_maze,n2_maze,LOI_maze,binSizeInfSpk,lagCCGrun);
% res1_maze = InfSp;
% [InfSp,InfoCCG] = InformedSpikes(n1_post,n2_post,LOI_post,binSizeInfSpk,lagCCGrun);
% res1_post = InfSp;
% toc

% constants
nepochs_plot = 3;
nrows        = 1;
njitter      = 100;
alpha        = 0.05;
for_grant    = false;
fs           = 30000;
jscale       = 1/fs;
fig_use      = 102;

binSize = 1/fs; % one sample for fine jitter
duration     = 0.007;

load("RoyMaze1_GapPairs_informed_spikes_CCGs_of_pair_20iv45i_data_jscale1_alpha5")
pair20iv45i = synchData.pair20iv45i;
load("RoyMaze1_GapPairs_informed_spikes_CCGs_of_pair_44iv79i_data_jscale1_alpha5")
pair44iv79i = synchData.pair44iv79i;

clear synchData

screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

% res_type = 'QHD';
% pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
% arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

[~,~,~,~,ccgR,tR,~,~,~,~] = ...
CCG_jitter(pair20iv45i.PreSyncSpikes,pair44iv79i.PreSyncSpikes,fs,binSize,duration,'jscale',jscale, ...
            'plot_output', get(fig_use, 'Number'), 'subfig', 1, ...
            'subplot_size', [nrows, 3], 'njitter', njitter, 'alpha', alpha,...
            'for_grant', for_grant);
title("Pre (20i v 45i) v (44i v 79i)")
        
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
CCG_jitter(pair20iv45i.MazeSyncSpikes,pair44iv79i.MazeSyncSpikes,fs,binSize,duration,'jscale',jscale, ...
            'plot_output', get(fig_use, 'Number'), 'subfig', 2, ...
            'subplot_size', [nrows, 3], 'njitter', njitter, 'alpha', alpha,...
            'for_grant', for_grant);
title("Maze (20i v 45i) v (44i v 79i)")
        
[~,~,~,~,ccgR,tR,~,~,~,~] = ...
CCG_jitter(pair20iv45i.PostSyncSpikes,pair44iv79i.PostSyncSpikes,fs,binSize,duration,'jscale',jscale, ...
            'plot_output', get(fig_use, 'Number'), 'subfig', 3, ...
            'subplot_size', [nrows, 3], 'njitter', njitter, 'alpha', alpha,...
            'for_grant', for_grant);
title("Post (20i v 45i) v (44i v 79i)")
        
 make_figure_pretty(102);