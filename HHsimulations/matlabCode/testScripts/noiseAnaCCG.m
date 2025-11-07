tiledlayout(4,2,'Padding','none','TileSpacing','tight')

nexttile(1)
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR1.mat','ccg','spikeTimesRef','tR')
plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'LineWidth',1)
hold on
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR5.mat','ccg','spikeTimesRef','tR')
plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'LineWidth',1)
hold off
% xlim([0 2])
legend('SNR1','SNR5','location','northwest')
title('H-H fast EpH sim CCG')
xlabel('[ms]')
ylabel('Probability')
%     title('R v T')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
box off

nexttile(3)
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR1.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',1)
hold on 
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR5.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',1)
hold off
% xlim([0 2])
title('H-H fast EpH evoked potential')
xlabel('[ms]')
ylabel('Vm [mV]')
% 	title('pulse')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
box off

nexttile(5)
%     yyaxis left
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR1.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1)
hold on
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR5.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1)
hold off
%     yyaxis right
%     plot((0:143)*dt*1000,waveform,'LineWidth',4) 
% xlim([0 2])
title('H-H fast EpH evoked potential derivative')
xlabel('[ms]')
ylabel('dVm/dt [mV/ms]')    
%     title('d(pulse)/dt')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
box off

nexttile(7)
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR1.mat','V')
histogram(V(3e4*60*2:end),'facealpha',.7,'edgecolor','none')
hold on
load('/media/nasko/WD_BLACK31/SimTemp/HHephSNR5.mat','V')
histogram(V(3e4*60*2:end),'facealpha',.7,'edgecolor','none')
hold off
title('H-H fast EpH membrane distribution')
xlabel('Vm [mV]')
ylabel('Counts')
% 	title('pulse')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
box off

%% 

nexttile(2)
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR1.mat','ccg','spikeTimesRef','tR')
plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'LineWidth',1)
hold on
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR5.mat','ccg','spikeTimesRef','tR')
plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'LineWidth',1)
hold off
legend('SNR1','SNR5','location','northwest')
title('EIF fast EpH sim CCG')
% xlim([0 2])
xlabel('[ms]')
ylabel('Probability')
%     title('R v T')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
box off

nexttile(4)
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR1.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',1)
hold on 
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR5.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',1)
hold off
% xlim([0 2])
title('EIF fast EpH evoked potential')
xlabel('[ms]')
ylabel('Vm [mV]')
% 	title('pulse')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
box off

nexttile(6)
%     yyaxis left
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR1.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1)
hold on
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR5.mat','preLength','postLength','dt','evokedSpikelet')
plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1)
hold off
%     yyaxis right
%     plot((0:143)*dt*1000,waveform,'LineWidth',4) 
% xlim([0 2])
title('EIF fast EpH evoked potential derivative')
xlabel('[ms]')
ylabel('dVm/dt [mV/ms]')    
%     title('d(pulse)/dt')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
box off

nexttile(8)
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR1.mat','V_trace')
histogram(V_trace(3e4*60*2:end),'facealpha',.7,'edgecolor','none')
hold on
load('/media/nasko/WD_BLACK31/SimTemp/EIFephSNR5.mat','V_trace')
histogram(V_trace(3e4*60*2:end),'facealpha',.7,'edgecolor','none')
hold off
title('EIF fast EpH membrane distribution')
xlabel('Vm [mV]')
ylabel('Counts')
% 	title('pulse')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
box off
