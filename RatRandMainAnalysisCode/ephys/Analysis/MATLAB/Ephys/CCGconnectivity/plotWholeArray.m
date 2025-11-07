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
                91  68;             79  95;             91  110];

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

for j = 37:size(pairsAna,1)
    
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
    
    chCell1 = Session.maxWaveformCh(waveIdxCell1) + 1;
    chCell2 = Session.maxWaveformCh(waveIdxCell2) + 1;
    
    cell1type = jit_var_out(current_pair_indices(2)).cell1type;
    cell2type = jit_var_out(current_pair_indices(2)).cell2type;
    
    cellID1 = current_pair(1);
    cellID2 = current_pair(2);
    
    % plot CCG 
    res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
    res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;
    
    loopChans = 1:64;
    
    %%
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', false, ...
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
    
    % remove sync spikes for waveforms
    res1nosync = setdiff(res1,res1sync);
    res2nosync = setdiff(res2,res2sync);
    
    %%
    
    for k = 1:64
        k
        tic 
        [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,k);
        [chanDataAtCell2, onsetTime] = HiroLoad300hz(session_name,k);
        
%         spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;
        
%         spikeTimeIndxCell1 = round((res1sync-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2 = round((res2sync-onsetTime/1e6)*fs) + 1;

        spikeTimeIndxCell1 = round((res1nosync-onsetTime/1e6)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2nosync-onsetTime/1e6)*fs) + 1;
        
        [spikeAvgMaxChanCell1(k,:), ~] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
        [spikeAvgMaxChanCell2(k,:), ~] = waveformAvg(chanDataAtCell2,spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);
        toc
    end

end

close 102


                         
                         
                         