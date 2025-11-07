load('medSplitSets_halfWidth_allStates.mat')

savepath = '/media/nasko/WD_BLACK3/HiroPhaseBinSpikeTimes/maze/halfWidth/';

fs           = 30000;
session_name = 'RoyMaze1';
freqRange    = {'agmon','theta','gamma'};

for i = 1:3
    for j = 1:size(medSplit.pairs,1)
        
        tic 
        
        chRef = medSplit.refCh(j);
        chTar = medSplit.tarCh(j);
        
        idRef = medSplit.pairs(j,1);
        idTar = medSplit.pairs(j,2);

        %% ref
        resRefAbove = medSplit.refAboveMed{j};
        resRefBelow = medSplit.refBelowMed{j};
        
        [phaseData,phaseBinSpikeTimes]                  = phaseAnalysis(resRefAbove,chRef,freqRange{i},session_name,fs);
        phaseBinCountsRefAboveMed.(freqRange{i}){j}     = phaseData;
%         phaseSpikeTimesRefAboveMed.(freqRange{i}){j}    = phaseBinSpikeTimes;
        save([savepath session_name num2str(idRef) 'v' num2str(idTar) freqRange{i} 'PhaseBinsSpikeTimesForRefAboveMed.mat'],'phaseBinSpikeTimes');
        
        [phaseData,phaseBinSpikeTimes]                  = phaseAnalysis(resRefBelow,chRef,freqRange{i},session_name,fs);
        phaseBinCountsRefBelowMed.(freqRange{i}){j}     = phaseData;
%         phaseSpikeTimesRefBelowMed.(freqRange{i}){j}    = phaseBinSpikeTimes;
        save([savepath session_name num2str(idRef) 'v' num2str(idTar) freqRange{i} 'PhaseBinsSpikeTimesForRefBelowMed.mat'],'phaseBinSpikeTimes');
        
        %% tar
        resTarAbove = medSplit.tarAboveMed{j};
        resTarBelow = medSplit.tarBelowMed{j};
        
        [phaseData,phaseBinSpikeTimes]                  = phaseAnalysis(resTarAbove,chRef,freqRange{i},session_name,fs);
        phaseBinCountsTarAboveMed.(freqRange{i}){j}     = phaseData;
%         phaseSpikeTimesTarAboveMed.(freqRange{i}){j}    = phaseBinSpikeTimes;
        save([savepath session_name num2str(idRef) 'v' num2str(idTar) freqRange{i} 'PhaseBinsSpikeTimesForTarAboveMed.mat'],'phaseBinSpikeTimes');
        
        [phaseData,phaseBinSpikeTimes]                  = phaseAnalysis(resTarBelow,chRef,freqRange{i},session_name,fs);
        phaseBinCountsTarBelowMed.(freqRange{i}){j}     = phaseData;
%         phaseSpikeTimesTarBelowMed.(freqRange{i}){j}    = phaseBinSpikeTimes;
        save([savepath session_name num2str(idRef) 'v' num2str(idTar) freqRange{i} 'PhaseBinsSpikeTimesForTarBelowMed.mat'],'phaseBinSpikeTimes');
        
        toc
        
        disp([num2str(((size(medSplit.pairs,1)*(i-1)+j)/(size(medSplit.pairs,1)*3))*100) ' % done!'])
        
    end
        
end

% j = 44;
% 
% tiledlayout(3,2)
% 
% xTickBinLabels = {'-\pi to -3\pi/4','-3\pi/4 to -\pi/2','-\pi/2 to -\pi/4','-\pi/4 to 0', ...
%                   '0 to \pi/4',     '\pi/4 to \pi/2' ,  '\pi/2 to 3\pi/4', '3\pi/4 to \pi'};
%               
% nexttile
% plot(phaseDataRefAboveMed.(freqRange{1}){j},'LineWidth',2)
% hold on
% plot(phaseDataRefBelowMed.(freqRange{1}){j},'LineWidth',2)
% hold off
% xticklabels(xTickBinLabels)
% ylabel('mean bin spike rate/mean spike rate [a.u.]')
% xlabel('phase bin [rad]')
% title(['halfWidth med split phase preference 350-450Hz ref neuron ' num2str(medSplit.pairs(j,1))])
% legend('above median','below median')
% 
% nexttile
% plot(phaseDataTarAboveMed.(freqRange{1}){j},'LineWidth',2)
% hold on
% plot(phaseDataTarBelowMed.(freqRange{1}){j},'LineWidth',2)
% hold off
% xticklabels(xTickBinLabels)
% ylabel('mean bin spike rate/mean spike rate [a.u.]')
% xlabel('phase bin [rad]')
% title(['halfWidth med split phase preference 350-450Hz tar neuron ' num2str(medSplit.pairs(j,2))])
% legend('above median','below median')
% 
% 
% nexttile
% plot(phaseDataRefAboveMed.(freqRange{2}){j},'LineWidth',2)
% hold on
% plot(phaseDataRefBelowMed.(freqRange{2}){j},'LineWidth',2)
% hold off
% xticklabels(xTickBinLabels)
% ylabel('mean bin spike rate/mean spike rate [a.u.]')
% xlabel('phase bin [rad]')
% title(['halfWidth med split phase preference theta ref neuron ' num2str(medSplit.pairs(j,1))])
% legend('above median','below median')
%         
% nexttile
% plot(phaseDataTarAboveMed.(freqRange{2}){j},'LineWidth',2)
% hold on
% plot(phaseDataTarBelowMed.(freqRange{2}){j},'LineWidth',2)
% hold off
% xticklabels(xTickBinLabels)
% ylabel('mean bin spike rate/mean spike rate [a.u.]')
% xlabel('phase bin [rad]')
% title(['halfWidth med split phase preference theta tar neuron ' num2str(medSplit.pairs(j,2))])
% legend('above median','below median')
% 
% nexttile
% plot(phaseDataRefAboveMed.(freqRange{3}){j},'LineWidth',2)
% hold on
% plot(phaseDataRefBelowMed.(freqRange{3}){j},'LineWidth',2)
% hold off
% xticklabels(xTickBinLabels)
% ylabel('mean bin spike rate/mean spike rate [a.u.]')
% xlabel('phase bin [rad]')
% title(['halfWidth med split phase preference gamma ref neuron ' num2str(medSplit.pairs(j,1))])
% legend('above median','below median')
%         
% nexttile
% plot(phaseDataTarAboveMed.(freqRange{3}){j},'LineWidth',2)
% hold on
% plot(phaseDataTarBelowMed.(freqRange{3}){j},'LineWidth',2)
% hold off
% xticklabels(xTickBinLabels)
% ylabel('mean bin spike rate/mean spike rate [a.u.]')
% xlabel('phase bin [rad]')
% title(['halfWidth med split phase preference gamma tar neuron ' num2str(medSplit.pairs(j,2))])
% legend('above median','below median')

