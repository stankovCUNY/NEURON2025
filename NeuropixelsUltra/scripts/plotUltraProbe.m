dataPath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';
load([dataPath 'neuronsSteinmetzDataArea' brainRegionSaveStr '.mat']);

chanLayout = flipud(reshape(1:384,[8,48])');

tWave = (-40:41)/30;

% refID = 2541;
% tarID = 2509;

refID = 1361;
tarID = 1335;

refIdx = find(neurons.unitID == refID);
tarIdx = find(neurons.unitID == tarID);

refWave = neurons.avgWaveforms{refIdx};
tarWave = neurons.avgWaveforms{tarIdx};

ymin = min(min([refWave tarWave]));
ymax = max(max([refWave tarWave]));

figure
tiledlayout(48,8)

for i = 1:384
   
    nexttile(i)
    plot(tWave,refWave(:,i),'LineWidth',2)
    hold on
    plot(tWave+2.7,tarWave(:,i),'LineWidth',2)
    hold off 
    
%     xlim([-1 1])
    xlim([tWave(1) tWave(end)+2.7])
    ylim([ymin ymax])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
end
