close all
clear 

tiledlayout(5,3,'Padding','none','TileSpacing','tight')

PvalbList   = [721123822
               746083955
               760345702
               773418906
               797828357
               829720705
               839557629
               840012044];

VipList     = [751348571
               755434585
               762120172
               791319847
               798911424
               816200189
               819701982
               835479236];

SstList     = [715093703
               719161530
               756029989
               758798717
               760693773
               762602078
               786091066
               787025148
               789848216
               794812542
               831882777
               839068429];
           
datasetID   = 'AllenInstitute';
% datapath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
datapath    = '/media/nasko/WD_BLACK31/BOTtemp/';

SNRthresh = 2;

monoInhBinNum = 0;

%% Pvalb

pairGroupStatTablePvalb = table;

for loopAnimals = 1:length(PvalbList)

    display(['running dataset: ' datasetID ', animal: '  num2str(PvalbList(loopAnimals))])

    load([datapath datasetID num2str(PvalbList(loopAnimals)) 'processedGroupStats.mat']);
    
    pairGroupStatTablePvalb = [pairGroupStatTablePvalb; 
                               pairGroupStatTable((strcmp(pairGroupStatTable.refOptoType2X(1:size(pairGroupStatTable,1)/2),'Pvalb') & ...
                                                         (pairGroupStatTable.refSNR(1:size(pairGroupStatTable,1)/2) > SNRthresh)) | ...
                                                  (strcmp(pairGroupStatTable.tarOptoType2X(1:size(pairGroupStatTable,1)/2),'Pvalb') & ...
                                                         (pairGroupStatTable.tarSNR(1:size(pairGroupStatTable,1)/2) > SNRthresh)) ...
                                                         ,:)];
    
end

refMonosyp = [];
tarMonosyp = [];

for loopPairs = 1:size(pairGroupStatTablePvalb,1)
    refMonosyp(loopPairs) = sum(pairGroupStatTablePvalb.GSPInh{loopPairs}(136:167)) >= monoInhBinNum;
    tarMonosyp(loopPairs) = sum(pairGroupStatTablePvalb.GSPInh{loopPairs}(47:76))   >= monoInhBinNum;
end

refMonosyp = refMonosyp';
tarMonosyp = tarMonosyp';

nexttile(1)
scatter(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength(strcmp(pairGroupStatTablePvalb.refCellExplorerType,'p')); ...
                  pairGroupStatTablePvalb.tarTroughToPeakLength(strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'p'))]), ... 
                 [pairGroupStatTablePvalb.refACGtauRise(strcmp(pairGroupStatTablePvalb.refCellExplorerType,'p')); ...
                  pairGroupStatTablePvalb.tarACGtauRise(strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'p'))], ...
        'o','LineWidth',1)

hold on
text(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength((strcmp(pairGroupStatTablePvalb.refCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarTroughToPeakLength((strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)]), ... 
              [pairGroupStatTablePvalb.refACGtauRise(        (strcmp(pairGroupStatTablePvalb.refCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarACGtauRise(        (strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)], ...
              [pairGroupStatTablePvalb.refOptoType2X(        (strcmp(pairGroupStatTablePvalb.refCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarOptoType2X(        (strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)], ...
              'FontSize',5)
 
scatter(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength(strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-narrow')); pairGroupStatTablePvalb.tarTroughToPeakLength(strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-narrow'))]), ... 
                 [pairGroupStatTablePvalb.refACGtauRise(strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-narrow')); pairGroupStatTablePvalb.tarACGtauRise(strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-narrow'))], ...
         'o','LineWidth',1)
text(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength((strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarTroughToPeakLength((strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)]), ... 
              [pairGroupStatTablePvalb.refACGtauRise(        (strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarACGtauRise(        (strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)], ...
              [pairGroupStatTablePvalb.refOptoType2X(        (strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarOptoType2X(        (strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)], ...
              'FontSize',5)
          
scatter(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength(strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-wide')); pairGroupStatTablePvalb.tarTroughToPeakLength(strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-wide'))]), ... 
                 [pairGroupStatTablePvalb.refACGtauRise(strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-wide')); pairGroupStatTablePvalb.tarACGtauRise(strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-wide'))], ...
         'o','LineWidth',1)
text(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength((strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarTroughToPeakLength((strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)]), ... 
              [pairGroupStatTablePvalb.refACGtauRise(        (strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarACGtauRise(        (strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)], ...
              [pairGroupStatTablePvalb.refOptoType2X(        (strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb')) & refMonosyp); ...
               pairGroupStatTablePvalb.tarOptoType2X(        (strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')) & tarMonosyp)], ...
              'FontSize',5)
  
for loopPairs = 1:size(pairGroupStatTablePvalb,1)
    patchline(cell2num([pairGroupStatTablePvalb.refTroughToPeakLength(loopPairs),pairGroupStatTablePvalb.tarTroughToPeakLength(loopPairs)]), ... 
              [pairGroupStatTablePvalb.refACGtauRise(loopPairs),pairGroupStatTablePvalb.tarACGtauRise(loopPairs)], ...
              'linewidth',0.01,'edgealpha',0.1)
end
hold off
set(gca, 'YScale', 'log')
legend('p','i-narrow','i-wide')
xlim([0 2])
ylabel('ACG \tau_r_i_s_e')
xlabel('trough-to-peak [ms]')
title(['Pvalb SNR threshold: ' num2str(SNRthresh)])
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

idxPvalbExq = find((strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb') & (refMonosyp & cell2num(pairGroupStatTablePvalb.flagExq)) | ...
                   (strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb') & (tarMonosyp & cell2num(pairGroupStatTablePvalb.flagExq)))));

idxRef = find(strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb'));
idxTar = find(strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb'));

nexttile(4)
violin([cell2num(pairGroupStatTablePvalb.refTroughToPeakLength(idxRef)); ...
        cell2num(pairGroupStatTablePvalb.tarTroughToPeakLength(idxTar))],'xlabel',{'Pvalb'});
ylim([0 1])
ylabel('Trough-to-peak [ms]')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

nexttile(7)
ISIlist = [];
for loopPairs = 1:length(idxRef)
    ISIlist = [ISIlist log10(mean(diff(pairGroupStatTablePvalb.refSpikeTimes{idxRef(loopPairs)}))*1000)];
end
for loopPairs = 1:length(idxTar)
    ISIlist = [ISIlist log10(mean(diff(pairGroupStatTablePvalb.tarSpikeTimes{idxTar(loopPairs)}))*1000)];
end

violin(ISIlist','xlabel',{'Pvalb'})
ylim([0 5])
ylabel('Avg. ISI log([ms])')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

% pairGroupStatTablePvalb(idx,:);

nexttile(10)
hold on
for loopPairs = 1:length(idxRef)
    patchline(pairGroupStatTablePvalb.tWave{1,1}, ... 
             [pairGroupStatTablePvalb.refWaveforms{idxRef(loopPairs),1}(pairGroupStatTablePvalb.refChannel(idxRef(loopPairs)),:)]/1e6, ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.25 .25 .25])
end
for loopPairs = 1:length(idxTar)
    patchline(pairGroupStatTablePvalb.tWave{1,1}, ... 
             [pairGroupStatTablePvalb.tarWaveforms{idxTar(loopPairs),1}(pairGroupStatTablePvalb.tarChannel(idxTar(loopPairs)),:)]/1e6, ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.25 .25 .25])
end
hold off
xlim([pairGroupStatTablePvalb.tWave{1,1}(1) pairGroupStatTablePvalb.tWave{1,1}(end)])
xlabel('[ms]')
ylabel('[mV]')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

nexttile(13)
plot(pairGroupStatTablePvalb.CCGbinLagTimes{1,1}*1000, ... 
             [pairGroupStatTablePvalb.pairRawCCG{idxPvalbExq(1),1}]/sum([pairGroupStatTablePvalb.pairRawCCG{idxPvalbExq(1),1}]),'k', ...
             'linewidth',1)
hold on
xline(-1,'--r');
xline(1, '--r');
for loopPairs = 2:length(idxPvalbExq)
    stkplt
    plot(pairGroupStatTablePvalb.CCGbinLagTimes{1,1}*1000, ... 
             [pairGroupStatTablePvalb.pairRawCCG{idxPvalbExq(loopPairs),1}]/sum([pairGroupStatTablePvalb.pairRawCCG{idxPvalbExq(loopPairs),1}]),'k', ...
             'linewidth',1)
%     drawnow
%     pause
end
hold off
if mod(length(idxPvalbExq),2)
    stkplt
end
xlim([-3.5 3.5])

%% Vip

pairGroupStatTableVip = table;

for loopAnimals = 1:length(VipList)

    display(['running dataset: ' datasetID ', animal: '  num2str(VipList(loopAnimals))])

    load([datapath datasetID num2str(VipList(loopAnimals)) 'processedGroupStats.mat']);
    
    pairGroupStatTableVip = [pairGroupStatTableVip; 
                             pairGroupStatTable((strcmp(pairGroupStatTable.refOptoType2X(1:size(pairGroupStatTable,1)/2),'Vip') & ...
                                                       (pairGroupStatTable.refSNR(1:size(pairGroupStatTable,1)/2) > SNRthresh)) | ...
                                                (strcmp(pairGroupStatTable.tarOptoType2X(1:size(pairGroupStatTable,1)/2),'Vip') & ...
                                                       (pairGroupStatTable.tarSNR(1:size(pairGroupStatTable,1)/2) > SNRthresh)) ...
                                                       ,:)];
    
end

refMonosyp = [];
tarMonosyp = [];

for loopPairs = 1:size(pairGroupStatTableVip,1)
    refMonosyp(loopPairs) = sum(pairGroupStatTableVip.GSPInh{loopPairs}(136:167)) >= monoInhBinNum;
    tarMonosyp(loopPairs) = sum(pairGroupStatTableVip.GSPInh{loopPairs}(47:76))   >= monoInhBinNum;
end

refMonosyp = refMonosyp';
tarMonosyp = tarMonosyp';

nexttile(2)
scatter(cell2num([pairGroupStatTableVip.refTroughToPeakLength(strcmp(pairGroupStatTableVip.refCellExplorerType,'p')); ...
                  pairGroupStatTableVip.tarTroughToPeakLength(strcmp(pairGroupStatTableVip.tarCellExplorerType,'p'))]), ... 
                 [pairGroupStatTableVip.refACGtauRise(strcmp(pairGroupStatTableVip.refCellExplorerType,'p')); ...
                  pairGroupStatTableVip.tarACGtauRise(strcmp(pairGroupStatTableVip.tarCellExplorerType,'p'))], ...
        'o','LineWidth',1)

hold on
text(cell2num([pairGroupStatTableVip.refTroughToPeakLength((strcmp(pairGroupStatTableVip.refCellExplorerType,'p') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarTroughToPeakLength((strcmp(pairGroupStatTableVip.tarCellExplorerType,'p') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)]), ... 
              [pairGroupStatTableVip.refACGtauRise(        (strcmp(pairGroupStatTableVip.refCellExplorerType,'p') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarACGtauRise(        (strcmp(pairGroupStatTableVip.tarCellExplorerType,'p') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)], ...
              [pairGroupStatTableVip.refOptoType2X(        (strcmp(pairGroupStatTableVip.refCellExplorerType,'p') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarOptoType2X(        (strcmp(pairGroupStatTableVip.tarCellExplorerType,'p') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)], ...
              'FontSize',5)
 
scatter(cell2num([pairGroupStatTableVip.refTroughToPeakLength(strcmp(pairGroupStatTableVip.refCellExplorerType,'i-narrow')); pairGroupStatTableVip.tarTroughToPeakLength(strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-narrow'))]), ... 
                 [pairGroupStatTableVip.refACGtauRise(strcmp(pairGroupStatTableVip.refCellExplorerType,'i-narrow')); pairGroupStatTableVip.tarACGtauRise(strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-narrow'))], ...
         'o','LineWidth',1)
text(cell2num([pairGroupStatTableVip.refTroughToPeakLength((strcmp(pairGroupStatTableVip.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarTroughToPeakLength((strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)]), ... 
              [pairGroupStatTableVip.refACGtauRise(        (strcmp(pairGroupStatTableVip.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarACGtauRise(        (strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)], ...
              [pairGroupStatTableVip.refOptoType2X(        (strcmp(pairGroupStatTableVip.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarOptoType2X(        (strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)], ...
              'FontSize',5)
  
scatter(cell2num([pairGroupStatTableVip.refTroughToPeakLength(strcmp(pairGroupStatTableVip.refCellExplorerType,'i-wide')); pairGroupStatTableVip.tarTroughToPeakLength(strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-wide'))]), ... 
                 [pairGroupStatTableVip.refACGtauRise(strcmp(pairGroupStatTableVip.refCellExplorerType,'i-wide')); pairGroupStatTableVip.tarACGtauRise(strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-wide'))], ...
         'o','LineWidth',1)
text(cell2num([pairGroupStatTableVip.refTroughToPeakLength((strcmp(pairGroupStatTableVip.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarTroughToPeakLength((strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)]), ... 
              [pairGroupStatTableVip.refACGtauRise(        (strcmp(pairGroupStatTableVip.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarACGtauRise(        (strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)], ...
              [pairGroupStatTableVip.refOptoType2X(        (strcmp(pairGroupStatTableVip.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip')) & refMonosyp); ...
               pairGroupStatTableVip.tarOptoType2X(        (strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')) & tarMonosyp)], ...
              'FontSize',5)
          
for loopPairs = 1:size(pairGroupStatTableVip,1)
    patchline(cell2num([pairGroupStatTableVip.refTroughToPeakLength(loopPairs),pairGroupStatTableVip.tarTroughToPeakLength(loopPairs)]), ... 
                       [pairGroupStatTableVip.refACGtauRise(loopPairs),pairGroupStatTableVip.tarACGtauRise(loopPairs)], ...
                       'linewidth',0.01,'edgealpha',0.1)
end
hold off
set(gca, 'YScale', 'log')
legend('p','i-narrow','i-wide')
xlim([0 2])
ylabel('ACG \tau_r_i_s_e')
xlabel('trough-to-peak [ms]')
title(['Vip SNR threshold: ' num2str(SNRthresh)])
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

idxVipExq = find((strcmp(pairGroupStatTableVip.refOptoType2X,'Vip') & (refMonosyp & cell2num(pairGroupStatTableVip.flagExq)) | ...
                 (strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip') & (tarMonosyp & cell2num(pairGroupStatTableVip.flagExq)))));

idxRef = find(strcmp(pairGroupStatTableVip.refOptoType2X,'Vip'));
idxTar = find(strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip'));

nexttile(5)
violin([cell2num(pairGroupStatTableVip.refTroughToPeakLength(idxRef)); ...
        cell2num(pairGroupStatTableVip.tarTroughToPeakLength(idxTar))],'xlabel',{'Vip'});
ylim([0 1])
ylabel('Trough-to-peak [ms]')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

ISIlist = [];
for loopPairs = 1:length(idxRef)
    ISIlist = [ISIlist log10(mean(diff(pairGroupStatTableVip.refSpikeTimes{idxRef(loopPairs)}))*1000)];
end
for loopPairs = 1:length(idxTar)
    ISIlist = [ISIlist log10(mean(diff(pairGroupStatTableVip.tarSpikeTimes{idxTar(loopPairs)}))*1000)];
end

nexttile(8)
violin(ISIlist','xlabel',{'Vip'})
ylim([0 5])
ylabel('Avg. ISI log([ms])')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

% pairGroupStatTableVip(idx,:);

nexttile(11)
hold on
for loopPairs = 1:length(idxRef)
    patchline(pairGroupStatTableVip.tWave{1,1}, ... 
             [pairGroupStatTableVip.refWaveforms{idxRef(loopPairs),1}(pairGroupStatTableVip.refChannel(idxRef(loopPairs)),:)]/1e6, ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.25 .25 .25])
end
for loopPairs = 1:length(idxTar)
    patchline(pairGroupStatTableVip.tWave{1,1}, ... 
             [pairGroupStatTableVip.tarWaveforms{idxTar(loopPairs),1}(pairGroupStatTableVip.tarChannel(idxTar(loopPairs)),:)]/1e6, ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.25 .25 .25])
end
hold off
xlim([pairGroupStatTableVip.tWave{1,1}(1) pairGroupStatTableVip.tWave{1,1}(end)])
xlabel('[ms]')
ylabel('[mV]')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

nexttile(14)
plot(pairGroupStatTableVip.CCGbinLagTimes{1,1}*1000, ... 
             [pairGroupStatTableVip.pairRawCCG{idxVipExq(1),1}]/sum([pairGroupStatTableVip.pairRawCCG{idxVipExq(1),1}]),'k', ...
             'linewidth',1)
hold on
xline(-1,'--r');
xline(1, '--r');
for loopPairs = 2:length(idxVipExq)
    stkplt
    plot(pairGroupStatTableVip.CCGbinLagTimes{1,1}*1000, ... 
             [pairGroupStatTableVip.pairRawCCG{idxVipExq(loopPairs),1}]/sum([pairGroupStatTableVip.pairRawCCG{idxVipExq(loopPairs),1}]),'k', ...
             'linewidth',1)
%     drawnow
%     pause
end
hold off
if mod(length(idxVipExq),2)
    stkplt
end
xlim([-3.5 3.5])

%% Sst

pairGroupStatTableSst = table;

for loopAnimals = 1:length(SstList)

    display(['running dataset: ' datasetID ', animal: '  num2str(SstList(loopAnimals))])

    load([datapath datasetID num2str(SstList(loopAnimals)) 'processedGroupStats.mat']);
    
    pairGroupStatTableSst = [pairGroupStatTableSst; 
                             pairGroupStatTable((strcmp(pairGroupStatTable.refOptoType2X(1:size(pairGroupStatTable,1)/2),'Sst') & ...
                                                       (pairGroupStatTable.refSNR(1:size(pairGroupStatTable,1)/2) > SNRthresh)) | ...
                                                (strcmp(pairGroupStatTable.tarOptoType2X(1:size(pairGroupStatTable,1)/2),'Sst') & ...
                                                       (pairGroupStatTable.tarSNR(1:size(pairGroupStatTable,1)/2) > SNRthresh)) ...
                                                       ,:)];
end

refMonosyp = [];
tarMonosyp = [];

for loopPairs = 1:size(pairGroupStatTableSst,1)
    refMonosyp(loopPairs) = sum(pairGroupStatTableSst.GSPInh{loopPairs}(136:167)) >= monoInhBinNum;
    tarMonosyp(loopPairs) = sum(pairGroupStatTableSst.GSPInh{loopPairs}(47:76))   >= monoInhBinNum;
end

refMonosyp = refMonosyp';
tarMonosyp = tarMonosyp';

nexttile(3)
scatter(cell2num([pairGroupStatTableSst.refTroughToPeakLength(strcmp(pairGroupStatTableSst.refCellExplorerType,'p')); ...
                  pairGroupStatTableSst.tarTroughToPeakLength(strcmp(pairGroupStatTableSst.tarCellExplorerType,'p'))]), ... 
                 [pairGroupStatTableSst.refACGtauRise(strcmp(pairGroupStatTableSst.refCellExplorerType,'p')); ...
                  pairGroupStatTableSst.tarACGtauRise(strcmp(pairGroupStatTableSst.tarCellExplorerType,'p'))], ...
        'o','LineWidth',1)

hold on
text(cell2num([pairGroupStatTableSst.refTroughToPeakLength((strcmp(pairGroupStatTableSst.refCellExplorerType,'p') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarTroughToPeakLength((strcmp(pairGroupStatTableSst.tarCellExplorerType,'p') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)]), ... 
              [pairGroupStatTableSst.refACGtauRise(        (strcmp(pairGroupStatTableSst.refCellExplorerType,'p') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarACGtauRise(        (strcmp(pairGroupStatTableSst.tarCellExplorerType,'p') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)], ...
              [pairGroupStatTableSst.refOptoType2X(        (strcmp(pairGroupStatTableSst.refCellExplorerType,'p') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarOptoType2X(        (strcmp(pairGroupStatTableSst.tarCellExplorerType,'p') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)], ...
              'FontSize',5)
          
scatter(cell2num([pairGroupStatTableSst.refTroughToPeakLength(strcmp(pairGroupStatTableSst.refCellExplorerType,'i-narrow')); pairGroupStatTableSst.tarTroughToPeakLength(strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-narrow'))]), ... 
                 [pairGroupStatTableSst.refACGtauRise(strcmp(pairGroupStatTableSst.refCellExplorerType,'i-narrow')); pairGroupStatTableSst.tarACGtauRise(strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-narrow'))], ...
         'o','LineWidth',1)
text(cell2num([pairGroupStatTableSst.refTroughToPeakLength((strcmp(pairGroupStatTableSst.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarTroughToPeakLength((strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)]), ... 
              [pairGroupStatTableSst.refACGtauRise(        (strcmp(pairGroupStatTableSst.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarACGtauRise(        (strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)], ...
              [pairGroupStatTableSst.refOptoType2X(        (strcmp(pairGroupStatTableSst.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarOptoType2X(        (strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)], ...
              'FontSize',5)
          
scatter(cell2num([pairGroupStatTableSst.refTroughToPeakLength(strcmp(pairGroupStatTableSst.refCellExplorerType,'i-wide')); pairGroupStatTableSst.tarTroughToPeakLength(strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-wide'))]), ... 
                 [pairGroupStatTableSst.refACGtauRise(strcmp(pairGroupStatTableSst.refCellExplorerType,'i-wide')); pairGroupStatTableSst.tarACGtauRise(strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-wide'))], ...
         'o','LineWidth',1)
text(cell2num([pairGroupStatTableSst.refTroughToPeakLength((strcmp(pairGroupStatTableSst.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarTroughToPeakLength((strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)]), ... 
              [pairGroupStatTableSst.refACGtauRise(        (strcmp(pairGroupStatTableSst.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarACGtauRise(        (strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)], ...
              [pairGroupStatTableSst.refOptoType2X(        (strcmp(pairGroupStatTableSst.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst')) & refMonosyp); ...
               pairGroupStatTableSst.tarOptoType2X(        (strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')) & tarMonosyp)], ...
              'FontSize',5)
          
for loopPairs = 1:size(pairGroupStatTableSst,1)
    patchline(cell2num([pairGroupStatTableSst.refTroughToPeakLength(loopPairs),pairGroupStatTableSst.tarTroughToPeakLength(loopPairs)]), ... 
                       [pairGroupStatTableSst.refACGtauRise(loopPairs),pairGroupStatTableSst.tarACGtauRise(loopPairs)], ...
                       'linewidth',0.01,'edgealpha',0.1)
end
hold off
set(gca, 'YScale', 'log')
legend('p','i-narrow','i-wide')
xlim([0 2])
ylabel('ACG \tau_r_i_s_e')
xlabel('trough-to-peak [ms]')
title(['Sst SNR threshold: ' num2str(SNRthresh)])
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

idxSstExq = find((strcmp(pairGroupStatTableSst.refOptoType2X,'Sst') & (refMonosyp & cell2num(pairGroupStatTableSst.flagExq)) | ...
                 (strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst') & (tarMonosyp & cell2num(pairGroupStatTableSst.flagExq)))));

idxRef = find(strcmp(pairGroupStatTableSst.refOptoType2X,'Sst'));
idxTar = find(strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst'));

nexttile(6)
violin([cell2num(pairGroupStatTableSst.refTroughToPeakLength(idxRef)); ...
        cell2num(pairGroupStatTableSst.tarTroughToPeakLength(idxTar))],'xlabel',{'Sst'});
ylim([0 1])
ylabel('Trough-to-peak [ms]')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

ISIlist = [];
for loopPairs = 1:length(idxRef)
    ISIlist = [ISIlist log10(mean(diff(pairGroupStatTableSst.refSpikeTimes{idxRef(loopPairs)}))*1000)];
end
for loopPairs = 1:length(idxTar)
    ISIlist = [ISIlist log10(mean(diff(pairGroupStatTableSst.tarSpikeTimes{idxTar(loopPairs)}))*1000)];
end

nexttile(9)
violin(ISIlist','xlabel',{'Sst'})
ylim([0 5])
ylabel('Avg. ISI log([ms])')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

% pairGroupStatTableSst(idx,:);

nexttile(12)
hold on
for loopPairs = 1:length(idxRef)
    patchline(pairGroupStatTableSst.tWave{1,1}, ... 
             [pairGroupStatTableSst.refWaveforms{idxRef(loopPairs),1}(pairGroupStatTableSst.refChannel(idxRef(loopPairs)),:)]/1e6, ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.25 .25 .25])
end
for loopPairs = 1:length(idxTar)
    patchline(pairGroupStatTableSst.tWave{1,1}, ... 
             [pairGroupStatTableSst.tarWaveforms{idxTar(loopPairs),1}(pairGroupStatTableSst.tarChannel(idxTar(loopPairs)),:)]/1e6, ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.25 .25 .25])
end
hold off
hold off
xlim([pairGroupStatTableSst.tWave{1,1}(1) pairGroupStatTableSst.tWave{1,1}(end)])
xlabel('[ms]')
ylabel('[mV]')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')

nexttile(15)
plot(pairGroupStatTableSst.CCGbinLagTimes{1,1}*1000, ... 
             [pairGroupStatTableSst.pairRawCCG{idxSstExq(1),1}]/sum([pairGroupStatTableSst.pairRawCCG{idxSstExq(1),1}]),'k', ...
             'linewidth',1)
hold on
xline(-1,'--r');
xline(1, '--r');
for loopPairs = 2:length(idxSstExq)
    stkplt
    plot(pairGroupStatTableSst.CCGbinLagTimes{1,1}*1000, ... 
             [pairGroupStatTableSst.pairRawCCG{idxSstExq(loopPairs),1}]/sum([pairGroupStatTableSst.pairRawCCG{idxSstExq(loopPairs),1}]),'k', ...
             'linewidth',1)
%     drawnow
%     pause
end
hold off
if mod(length(idxSstExq),2)
    stkplt
end
xlim([-3.5 3.5])

set(gcf, 'Renderer', 'Painters')

%%
% Pvalb
sum([strcmp(pairGroupStatTablePvalb.refCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb'); ...
     strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'p') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')])
sum([strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb'); ...
     strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')])
sum([strcmp(pairGroupStatTablePvalb.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.refOptoType2X,'Pvalb'); ...
     strcmp(pairGroupStatTablePvalb.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTablePvalb.tarOptoType2X,'Pvalb')])

% Vip
sum([strcmp(pairGroupStatTableVip.refCellExplorerType,'p') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip'); ...
     strcmp(pairGroupStatTableVip.tarCellExplorerType,'p') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')])
sum([strcmp(pairGroupStatTableVip.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip'); ...
     strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')])
sum([strcmp(pairGroupStatTableVip.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.refOptoType2X,'Vip'); ...
     strcmp(pairGroupStatTableVip.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableVip.tarOptoType2X,'Vip')])
 
% Sst 
sum([strcmp(pairGroupStatTableSst.refCellExplorerType,'p') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst'); ...
     strcmp(pairGroupStatTableSst.tarCellExplorerType,'p') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')])
sum([strcmp(pairGroupStatTableSst.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst'); ...
     strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-narrow') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')])
sum([strcmp(pairGroupStatTableSst.refCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.refOptoType2X,'Sst'); ...
     strcmp(pairGroupStatTableSst.tarCellExplorerType,'i-wide') & strcmp(pairGroupStatTableSst.tarOptoType2X,'Sst')])
