%% phase/CCG/waveform plots

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

% pair data 
datapath  = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
phasepath = ['/media/nasko/WD_BLACK3/HiroPhaseBinSpikeTimes/maze/' metric '/'];

[jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

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

current_pair = [20 45];
    
current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                            (current_pair(2) == cell_pairs_mat(:,2)));

resRef = jit_var_out(current_pair_indices(2)).cell1_spike_times;
resTar = jit_var_out(current_pair_indices(2)).cell2_spike_times;                        

numSpikesRef = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
               length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
               length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

numSpikesTar = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
               length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
               length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

waveIdxRef = find(numSpikesRef == waveDataSpCounts);
waveIdxTar = find(numSpikesTar == waveDataSpCounts);

chRef = Session.maxWaveformCh(waveIdxRef);
chTar = Session.maxWaveformCh(waveIdxTar);

[chanDataAtRef, onsetTime] = HiroLoad300hz(session_name,chRef);
[chanDataAtTar, onsetTime] = HiroLoad300hz(session_name,chTar);

%%

refAboveMed = load([phasepath 'RoyMaze1' num2str(current_pair(1)) 'v' num2str(current_pair(2)) state 'PhaseBinsSpikeTimesForRefAboveMed.mat']);
refBelowMed = load([phasepath 'RoyMaze1' num2str(current_pair(1)) 'v' num2str(current_pair(2)) state 'PhaseBinsSpikeTimesForRefBelowMed.mat']);
tarAboveMed = load([phasepath 'RoyMaze1' num2str(current_pair(1)) 'v' num2str(current_pair(2)) state 'PhaseBinsSpikeTimesForTarAboveMed.mat']);
tarBelowMed = load([phasepath 'RoyMaze1' num2str(current_pair(1)) 'v' num2str(current_pair(2)) state 'PhaseBinsSpikeTimesForTarBelowMed.mat']);

figure(102)

CCGrefAboveMed = zeros(8,211);
CCGrefBelowMed = zeros(8,211);
CCGtarAboveMed = zeros(8,211);
CCGtarBelowMed = zeros(8,211);

spikeAvgMaxChanRefBelowMedAll = zeros(8,72);
spikeAvgMaxChanRefAboveMedAll = zeros(8,72);
spikeAvgMaxChanTarBelowMedAll = zeros(8,72);
spikeAvgMaxChanTarAboveMedAll = zeros(8,72);

for i = 1:8
    
    i
    tic
    
    %%
%     nexttile
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter([refAboveMed.phaseBinSpikeTimes{1,i}{1,1}],resTar,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefAboveMed(i,:) = ccgR(:,1,2);
    
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(refBelowMed.phaseBinSpikeTimes{1,i}{1,1},resTar,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefBelowMed(i,:) = ccgR(:,1,2);
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter([refAboveMed.phaseBinSpikeTimes{1,i}{1,1}; refBelowMed.phaseBinSpikeTimes{1,i}{1,1}], ...
                        resTar,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefAll(i,:) = ccgR(:,1,2);       
    
    spikeTimeIndxRefBelowMed = round((refBelowMed.phaseBinSpikeTimes{1,i}{1,1}-onsetTime/1e6)*fs) + 1; 
    [spikeAvgMaxChanRefBelowMed, waveformsMaxRefBelowMed] = waveformAvg(chanDataAtRef,spikeTimeIndxRefBelowMed,preLength,postLength,fpass,fs,filterFlag);
    spikeAvgMaxChanRefBelowMedAll(i,:) = spikeAvgMaxChanRefBelowMed;
    waveformsMaxRefBelowMedAll{i} =  waveformsMaxRefBelowMed;
    
    spikeTimeIndxRefAboveMed = round((refAboveMed.phaseBinSpikeTimes{1,i}{1,1}-onsetTime/1e6)*fs) + 1; 
    [spikeAvgMaxChanRefAboveMed, waveformsMaxRefAboveMed] = waveformAvg(chanDataAtRef,spikeTimeIndxRefAboveMed,preLength,postLength,fpass,fs,filterFlag);
    spikeAvgMaxChanRefAboveMedAll(i,:) = spikeAvgMaxChanRefAboveMed;
    waveformsMaxRefAboveMedAll{i} =  waveformsMaxRefAboveMed;
    
%     plot(tR,CCGrefAboveMed)
%     hold on
%     plot(tR,CCGrefBelowMed)
%     hold off
%     legend('CCGrefAboveMed','CCGrefBelowMed')
    
    %%
%     nexttile
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(resRef,tarAboveMed.phaseBinSpikeTimes{1,i}{1,1},fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGtarAboveMed(i,:) = ccgR(:,1,2);
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(resRef,tarBelowMed.phaseBinSpikeTimes{1,i}{1,1},fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGtarBelowMed(i,:) = ccgR(:,1,2);
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(resRef,[tarAboveMed.phaseBinSpikeTimes{1,i}{1,1}; tarBelowMed.phaseBinSpikeTimes{1,i}{1,1}], ...
                        fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGtarAll(i,:) = ccgR(:,1,2);
    
    spikeTimeIndxTarBelowMed = round((tarBelowMed.phaseBinSpikeTimes{1,i}{1,1}-onsetTime/1e6)*fs) + 1; 
    [spikeAvgMaxChanTarBelowMed, waveformsMaxTarBelowMed] = waveformAvg(chanDataAtTar,spikeTimeIndxTarBelowMed,preLength,postLength,fpass,fs,filterFlag);
    spikeAvgMaxChanTarBelowMedAll(i,:) = spikeAvgMaxChanTarBelowMed;
    waveformsMaxTarBelowMedAll{i} =  waveformsMaxTarBelowMed;
    
    spikeTimeIndxTarAboveMed = round((tarAboveMed.phaseBinSpikeTimes{1,i}{1,1}-onsetTime/1e6)*fs) + 1; 
    [spikeAvgMaxChanTarAboveMed, waveformsMaxTarAboveMed] = waveformAvg(chanDataAtTar,spikeTimeIndxTarAboveMed,preLength,postLength,fpass,fs,filterFlag);
    spikeAvgMaxChanTarAboveMedAll(i,:) = spikeAvgMaxChanTarAboveMed;
    waveformsMaxTarAboveMedAll{i} =  waveformsMaxTarAboveMed;
    
%     plot(tR,CCGtarAboveMed)
%     hold on
%     plot(tR,CCGtarBelowMed)
%     hold off
%     legend('CCGrefAboveMed','CCGrefBelowMed')

    toc
    
end

binLabels = {'-\pi to -3\pi/4','-3\pi/4 to -\pi/2','-\pi/2 to -\pi/4','-\pi/4 to 0', ...
             '0 to \pi/4',     '\pi/4 to \pi/2' ,  '\pi/2 to 3\pi/4', '3\pi/4 to \pi'};

tiledlayout(8,4)
for i = 1:8
    
    nexttile
    plot(tR*1000,CCGrefAll(i,:)/2,'LineWidth',2)
    hold on 
    plot(tR*1000,CCGrefAboveMed(i,:),'LineWidth',2)
    plot(tR*1000,CCGrefBelowMed(i,:),'LineWidth',2)
    hold off
    ylim([0 200])
    xlim([-3.5 3.5])
    title([num2str(current_pair(1)) ' v ' num2str(current_pair(2)) ': ' binLabels{i}])
    legend([num2str(current_pair(1)) ' ref all'], ...
           [num2str(current_pair(1)) ' ref above med'], ...
           [num2str(current_pair(1)) ' ref below med'])
       
    nexttile
    plot(spikeAvgMaxChanRefAboveMedAll(i,:))
%     plot(std(waveformsMaxRefAboveMedAll{i},0,2))
    hold on
    plot(spikeAvgMaxChanRefBelowMedAll(i,:))
%     plot(std(waveformsMaxRefBelowMedAll{i},0,2))
    hold off
    
    nexttile
    plot(tR*1000,CCGtarAll(i,:)/2,'LineWidth',2)
    hold on 
    plot(tR*1000,CCGtarAboveMed(i,:),'LineWidth',2)
    plot(tR*1000,CCGtarBelowMed(i,:),'LineWidth',2)
    hold off
    ylim([0 100])
    xlim([-3.5 3.5])
    title([num2str(current_pair(1)) ' v ' num2str(current_pair(2)) ': ' binLabels{i}])
    legend([num2str(current_pair(2)) ' tar all'], ...
           [num2str(current_pair(2)) ' tar above med'], ...
           [num2str(current_pair(2)) ' tar below med'])
       
    nexttile
    plot(spikeAvgMaxChanTarAboveMedAll(i,:))
%     plot(std(waveformsMaxTarAboveMedAll{i},0,2))
    hold on
    plot(spikeAvgMaxChanTarBelowMedAll(i,:))
%     plot(std(waveformsMaxTarBelowMedAll{i},0,2))
    hold off
    
end

%%

% tiledlayout(2,1)
% 
% load('res1_20i_aboveMed.mat')
% load('res1_20i_belowMed.mat')
% 
% nexttile
% histogram(diff(resAboveMed),100000,'Normalization','probability')
% hold on
% histogram(diff(resBelowMed),100000,'Normalization','probability')
% xlim([0 0.01])
% hold off
% 
% load('res1_45i_aboveMed.mat')
% load('res1_45i_belowMed.mat')
% 
% nexttile
% histogram(diff(resAboveMed),100000,'Normalization','probability')
% hold on
% histogram(diff(resBelowMed),100000,'Normalization','probability')
% xlim([0 0.01])
% hold off

%%

% session_name   = 'RoyMaze1';
% conn_type      = {'GapPairs','ExcPairs','InhPairs'};
% jscale         = 1;
% alpha_name     = 5;
% duration       = 0.02;
% fs             = 30000;
% fpass          = 300;
% binSize        = 1/fs;
% fig_use        = 102;
% njitter        = 500;
% alpha          = 0.05;
% for_grant      = false;
% filterFlag     = true;
% Nwaveforms     = 10000;
% 
% figure(102)
% tiledlayout(2,1)
% 
% load('res1_20i_aboveMed.mat')
% load('res1_20i_belowMed.mat')
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%             CCG_jitter(resAboveMed,resAboveMed,fs,binSize,duration,'jscale',jscale, ...
%                         'plot_flag', false, ...
%                         'plot_output', get(fig_use, 'Number'), ...
%                         'njitter', njitter, 'alpha', alpha,...
%                         'for_grant', for_grant);
%                     
% ACGabove = ccgR(:,1,2);
% ACGabove(301) = 0;
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%             CCG_jitter(resBelowMed,resBelowMed,fs,binSize,duration,'jscale',jscale, ...
%                         'plot_flag', false, ...
%                         'plot_output', get(fig_use, 'Number'), ...
%                         'njitter', njitter, 'alpha', alpha,...
%                         'for_grant', for_grant);
%                     
% ACGbelow = ccgR(:,1,2);
% ACGbelow(301) = 0;
% 
% nexttile
% plot(tR*1000,ACGabove/sum(ACGabove),'LineWidth',1.5)
% hold on
% plot(tR*1000,ACGbelow/sum(ACGbelow),'LineWidth',1.5)
% hold off
% legend('above','below')
% 
% load('res1_45i_aboveMed.mat')
% load('res1_45i_belowMed.mat')
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%             CCG_jitter(resAboveMed,resAboveMed,fs,binSize,duration,'jscale',jscale, ...
%                         'plot_flag', false, ...
%                         'plot_output', get(fig_use, 'Number'), ...
%                         'njitter', njitter, 'alpha', alpha,...
%                         'for_grant', for_grant);
%                     
% ACGabove = ccgR(:,1,2);
% ACGabove(301) = 0;
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%             CCG_jitter(resBelowMed,resBelowMed,fs,binSize,duration,'jscale',jscale, ...
%                         'plot_flag', false, ...
%                         'plot_output', get(fig_use, 'Number'), ...
%                         'njitter', njitter, 'alpha', alpha,...
%                         'for_grant', for_grant);
%                     
% ACGbelow = ccgR(:,1,2);
% ACGbelow(301) = 0;
% 
% nexttile
% plot(tR*1000,ACGabove/sum(ACGabove),'LineWidth',1.5)
% hold on
% plot(tR*1000,ACGbelow/sum(ACGbelow),'LineWidth',1.5)
% hold off
% legend('above','below')

%%

% Some data
% T.Country = string(('A':'Z')');
% T.Population = randi([1000 10000], size(T.Country));
% T = struct2table(T);
% n = height(T); 
% dtheta = 360/n;         % sector for each country
% gap = 0.1;              % gaps betwen adjacent contries
% npoints = 21;           % npoints for each sector
% figure; hold on;
% for i=1:n
%     theta = linspace((i-1)*dtheta, (i-gap)*dtheta, npoints);
%     patch([0 T.Population(i)*cosd(theta)], [0 T.Population(i)*sind(theta)], 'b');
%     rtext = T.Population(i)+300;
%     text(rtext*cosd((i-.5)*dtheta), rtext*sind((i-.5)*dtheta), T.Country(i));
% end
% axis equal
% axis off


