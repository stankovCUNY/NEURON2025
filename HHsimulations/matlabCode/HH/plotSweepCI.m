for loopSims = 1:3
    
    dcLevels = [5 7.5 10];
    
    tiledlayout(3,1)

    %%
    
    nexttile
    
    load(['/media/nasko/WD_BLACK31/SimTemp/HHciGJdc' num2str(dcLevels(1)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgCiT(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    hold on

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciGJdc' num2str(dcLevels(2)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgCiT(:,1,2)/length(spikeTimesRef),'k--','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciGJdc' num2str(dcLevels(3)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgCiT(:,1,2)/length(spikeTimesRef),'k:','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    hold off
    
    legend('DC 5','DC 7.5','DC 10')
    
    %%
    
    nexttile
    
    load(['/media/nasko/WD_BLACK31/SimTemp/HHciGJdc' num2str(dcLevels(1)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    hold on

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciGJdc' num2str(dcLevels(2)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k--','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciGJdc' num2str(dcLevels(3)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k:','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    hold off
    
    %%
    
    nexttile

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciMSdc' num2str(dcLevels(1)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    hold on

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciMSdc' num2str(dcLevels(2)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k--','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    load(['/media/nasko/WD_BLACK31/SimTemp/HHciMSdc' num2str(dcLevels(3)) '.mat'],'ccgCiT','ccgRT','spikeTimesCI','spikeTimesRef','spikeTimesTar','tR')
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k:','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

    hold off
   
end