%% pair ISIs

subplot(1,3,1)
histogram(diff(res1),'BinWidth',1/100,'BinLimits',[0 .5],'Normalization','probability')
set(gca, 'XScale', 'log')
subplot(1,3,2)
histogram(diff(res2),'BinWidth',1/100,'BinLimits',[0 .5],'Normalization','probability')
set(gca, 'XScale', 'log')
subplot(1,3,3)
histogram(diff(res1(refSyncTriggers)),'BinWidth',1/100,'BinLimits',[0 0.5],'Normalization','probability')
set(gca, 'XScale', 'log')

%%

meanRefSyncWaveform   = mean(waveDataCell1Bin1,1)/1000;
meanTarSyncWaveform   = mean(waveDataCell2Bin1,1)/1000;

meanRefNoSyncWaveform = mean(waveDataCell1NoSync,1)/1000;
meanTarNoSyncWaveform = mean(waveDataCell2NoSync,1)/1000;

stanDevRefSyncWaveform   = std(waveDataCell1Bin1,0,1)/1000;
stanDevTarSyncWaveform   = std(waveDataCell2Bin1,0,1)/1000;

stanDevRefNoSyncWaveform = std(waveDataCell1NoSync,0,1)/1000;
stanDevTarNoSyncWaveform = std(waveDataCell2NoSync,0,1)/1000;

%%

Async   = meanRefSyncWaveform   + stanDevRefSyncWaveform;
Bsync   = meanRefSyncWaveform   - stanDevRefSyncWaveform;

AnoSync = meanRefNoSyncWaveform + stanDevRefNoSyncWaveform;
BnoSync = meanRefNoSyncWaveform - stanDevRefNoSyncWaveform;

plot(tWave,meanRefSyncWaveform,'--', 'color','r','LineWidth',2)
hold on
plot(tWave,meanRefNoSyncWaveform,'-', 'color','b','LineWidth',2)
patch([tWave'; flipud(tWave')],[Async';   flipud(Bsync')],   'r','EdgeAlpha',0, 'FaceAlpha', 0.1); 
patch([tWave'; flipud(tWave')],[AnoSync'; flipud(BnoSync')], 'b','EdgeAlpha',0, 'FaceAlpha', 0.1); 
hold off

%%

Async   = meanTarSyncWaveform   + stanDevTarSyncWaveform;
Bsync   = meanTarSyncWaveform   - stanDevTarSyncWaveform;

AnoSync = meanTarNoSyncWaveform + stanDevTarNoSyncWaveform;
BnoSync = meanTarNoSyncWaveform - stanDevTarNoSyncWaveform;

plot(tWave,meanTarSyncWaveform,'--', 'color','r','LineWidth',2)
hold on
plot(tWave,meanTarNoSyncWaveform,'-', 'color','b','LineWidth',2)
patch([tWave'; flipud(tWave')],[Async';   flipud(Bsync')],   'r','EdgeAlpha',0, 'FaceAlpha', 0.1); 
patch([tWave'; flipud(tWave')],[AnoSync'; flipud(BnoSync')], 'b','EdgeAlpha',0, 'FaceAlpha', 0.1); 
hold off

%% left over scatter code

%     plot(tWave,meanRefWaveform,'-', 'color', [0      0.4470 0.7410],'LineWidth',2)
%     plot(tWave,meanTarWaveform,'--', 'color',[0      0.4470 0.7410],'LineWidth',2)
    
%     current_pair
%     [~,idx] = max(meanRefWaveform(1:27));   idx
%     [~,idx] = max(meanRefWaveform(28:end)); idx + 27
%     [~,idx] = max(meanTarWaveform(1:27));   idx
%     [~,idx] = max(meanTarWaveform(28:end)); idx + 27

nexttile
trials = [(res1(refSyncTriggers) - PSTHpreTrialLen) (res1(refSyncTriggers) + PSTHtrialLen)];
PSTH(res2,trials,PSTHtrialLen,PSTHpreTrialLen,PSTHbinSize,'k')

nexttile
trials = [(res2(tarSyncTriggers) - PSTHpreTrialLen) (res2(tarSyncTriggers) + PSTHtrialLen)];
PSTH(res1,trials,PSTHtrialLen,PSTHpreTrialLen,PSTHbinSize,'k')

%%

STref = res1; 
STtar = res2;
mappedNeuron = 'tar';
fileVideoName = sprintf('%s - %d%s (sh %d) v %d%s (sh %d) - %', ...
                        jit_var_out(current_pair_indices(2)).epoch, ...
                        jit_var_out(current_pair_indices(2)).cell_pair(1),jit_var_out(current_pair_indices(2)).cell1type, ...
                        n1shank, ...
                        jit_var_out(current_pair_indices(2)).cell_pair(2),jit_var_out(current_pair_indices(2)).cell2type, ...
                        n2shank,mappedNeuron);
syncCCGanaFlag = false;
syncSTpeak = SyncSpIdxBin{CCGpeaks(j)};
syncSTref  = STref(syncSTpeak(:,1));
syncSTtar  = STtar(syncSTpeak(:,2));


arrayVoltageMap(STref,STtar,tWave,fileVideoName,syncCCGanaFlag,syncSTref,syncSTtar,mappedNeuron,CCGparams,CCGpeaks(j))

