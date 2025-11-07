clear all

%% 

metric = 'halfWidth';
state  = 'theta';

session_name   = 'RoyMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1;
alpha_name     = 5;
duration       = 0.007;
fs             = 30000;
fpass          = 300;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
filterFlag     = false;

preLength  = 36;
postLength = 36;

tWave     = (-26:27)*(1/30);
spikeTime = (1/30)*(-preLength+1:postLength);

% pair data 
datapath  = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
phasepath = ['/media/nasko/WD_BLACK3/HiroPhaseBinSpikeTimes/maze/' metric '/'];

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);
jit_var_out = addFlipped_Jit_var_out(jit_var_out);

% get all pairs
cell_pairs = {jit_var_out.cell_pair}';

% convert to matrix
cell_pairs_mat = [];
for i = 1:size(cell_pairs,1)
    cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
end

%% get channels from waveform files. Not computationally savvy but it's fastest solution
extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';

Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);
waveDataSpCounts = [];
for i = 1:size(Session.times,2)
    waveDataSpCounts(i) = size(Session.times{i},1);
end

%%

current_pair = [91 82];
    
current_pair_indices_jit_var = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));

resRef = jit_var_out(current_pair_indices_jit_var(2)).cell1_spike_times;
resTar = jit_var_out(current_pair_indices_jit_var(2)).cell2_spike_times;                        

numSpikesRef = length(jit_var_out(current_pair_indices_jit_var(1)).cell1_spike_times) + ...
               length(jit_var_out(current_pair_indices_jit_var(2)).cell1_spike_times) + ...
               length(jit_var_out(current_pair_indices_jit_var(3)).cell1_spike_times);

numSpikesTar = length(jit_var_out(current_pair_indices_jit_var(1)).cell2_spike_times) + ...
               length(jit_var_out(current_pair_indices_jit_var(2)).cell2_spike_times) + ...
               length(jit_var_out(current_pair_indices_jit_var(3)).cell2_spike_times);

waveIdxRef = find(numSpikesRef == waveDataSpCounts);
waveIdxTar = find(numSpikesTar == waveDataSpCounts);

chRef = Session.maxWaveformCh(waveIdxRef);
chTar = Session.maxWaveformCh(waveIdxTar);

[chanDataAtRef, onsetTime] = HiroLoad300hz(session_name,chRef);
[chanDataAtTar, onsetTime] = HiroLoad300hz(session_name,chTar);

%%

load('RoyMaze1GroupStats.mat')
current_pair_indices_group_stats = find((current_pair(1) == groupStat.mainRef) & ...
                                        (current_pair(2) == groupStat.mainTar));

[chanDataAtMaxCorrCh, onsetTime] = HiroLoad300hz(session_name,groupStat.mainMaxChan(current_pair_indices_group_stats));                                    
                                    
peakBinTime = groupStat.mainFeatCCGlatAll(current_pair_indices_group_stats);

figure(102)
[GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
        CCG_jitter(resRef,resTar,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', false, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant);
CCG = ccgR(:,1,2);

CCGbinsOfInterest = zeros(length(tR),1);
CCGbinsOfInterest((tR*1000) == peakBinTime) = 1;
[SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(resRef, resTar, CCGbinsOfInterest,duration);

spikeTimeIndxRef = round((SyncSp(:,1)-onsetTime/1e6)*fs) + 1;
[spikeAvgMaxChanRef, waveformsMaxRef]    = waveformAvg(chanDataAtRef,      spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);
[spikeAvgMaxCorrCh,  waveformsMaxCorrCh] = waveformAvg(chanDataAtMaxCorrCh,spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);

spikeTimeIndxTar = round((SyncSp(:,2)-onsetTime/1e6)*fs) + 1;
[spikeAvgMaxChanTar, waveformsMaxTar] = waveformAvg(chanDataAtTar,spikeTimeIndxTar,preLength,postLength,fpass,fs,filterFlag);

refFlucsProximal = waveformsMaxRef - spikeAvgMaxChanRef;
refFlucsDistal   = waveformsMaxCorrCh - spikeAvgMaxCorrCh;
tarFlucs = waveformsMaxTar - spikeAvgMaxChanTar;

matCorr = corr(refFlucsProximal',refFlucsDistal');
Rsq = matCorr.^2;

v = max(Rsq(:));           % value
[ii,jj] = find(Rsq == v);  % row and column number

scatter(refFlucs(ii,:),tarFlucs(jj,:))

spikeTime(ii)

tiledlayout(2,1)

nexttile
k = -7;
plot(spikeTime((-k/2)+1:end+(k/2)),diag(Rsq,k))

nexttile
plot(spikeTime,spikeAvgMaxCorrCh)
hold on
plot(spikeTime,spikeAvgMaxChanRef)
plot(spikeTime,spikeAvgMaxChanTar)
hold off

%%

spikeTimeIndxRef = round((SyncSp(:,1)-onsetTime/1e6)*fs) + 1;

[spikeAvgRefAtRefChSync, waveformsMaxRefAtRefChSync] = waveformAvg(chanDataAtRef, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);
[spikeAvgRefAtTarChSync, waveformsMaxRefAtTarChSync] = waveformAvg(chanDataAtTar, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);


resRefNoSync = setdiff(resRef,SyncSp(:,1));
spikeTimeIndxRef = round((resRefNoSync(:,1)-onsetTime/1e6)*fs) + 1;

[spikeAvgRefAtRefChNoSync, waveformsMaxRefAtRefChNoSync] = waveformAvg(chanDataAtRef, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);
[spikeAvgRefAtTarChNoSync, waveformsMaxRefAtTarChNoSync] = waveformAvg(chanDataAtTar, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);

%%

spikeTimeIndxRef = round((SyncSp(:,1)-onsetTime/1e6)*fs) + 1;
[spikeAvgRefAtMaxCorrChSync, waveformsMaxRefAtMaxCorrChSync] = waveformAvg(chanDataAtMaxCorrCh, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);

resRefNoSync = setdiff(resRef,SyncSp(:,1));
spikeTimeIndxRef = round((resRefNoSync(:,1)-onsetTime/1e6)*fs) + 1;
[spikeAvgRefAtMaxCorrChNoSync, waveformsMaxRefAtMaxCorrChNoSync] = waveformAvg(chanDataAtMaxCorrCh, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);

%%

offset = 0;
% offset = -0.000533333;

spikeTimeIndxRef = round(((SyncSp(:,2) + offset)-onsetTime/1e6)*fs) + 1;
[spikeAvgTarAtTarChSync, waveformsMaxTarAtTarChSync] = waveformAvg(chanDataAtTar, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);

resTarNoSync = setdiff(resTar,SyncSp(:,2)) + offset;
spikeTimeIndxRef = round((resTarNoSync(:,1)-onsetTime/1e6)*fs) + 1;
[spikeAvgTarAtTarChNoSync, waveformsMaxTarAtTarChNoSync] = waveformAvg(chanDataAtTar, spikeTimeIndxRef,preLength,postLength,fpass,fs,filterFlag);


%%

tiledlayout(6,1)

ylimits = [-1 0.5];

nexttile
plot(spikeTime,spikeAvgTarAtTarChNoSync)
% ylim(ylimits)
set(gca,'YDir','reverse')

nexttile
plot(spikeTime,spikeAvgTarAtTarChSync)
% ylim(ylimits)
set(gca,'YDir','reverse')

nexttile
plot(spikeTime,spikeAvgTarAtTarChSync - spikeAvgTarAtTarChNoSync)
% ylim(ylimits)
set(gca,'YDir','reverse')

nexttile
plot(spikeTime,spikeAvgRefAtTarChNoSync)
% ylim(ylimits])
set(gca,'YDir','reverse')

nexttile
plot(spikeTime,spikeAvgRefAtMaxCorrChNoSync)
% ylim(ylimits)
set(gca,'YDir','reverse')

nexttile
plot(tR*1000,CCG)
xlim([-1.5 1.5])





