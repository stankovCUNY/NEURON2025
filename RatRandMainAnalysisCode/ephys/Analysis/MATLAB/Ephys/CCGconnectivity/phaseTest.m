clear

load('medSplitSets_peakToTrough_allStates.mat')

session_name   = 'RoyMaze1';
conn_type      = {'GapPairs','ExcPairs','InhPairs'};
jscale         = 1;
alpha_name     = 5;
duration       = 0.01;
fs             = 30000;
fpass          = 300;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
filterFlag     = true;
Nwaveforms     = 10000;

winSize        = 0;
% winSize        = 30*10;
freqRange      = 'agmon';
medianSet      = 'sync';

pairs = {[45,20],...
         [91,79],...
         [79,91],...
         [91,118],...
         [79,44]};

refCh = [1 , 28, 18, 28, 18];
tarCh = [10,18, 28, 33, 12];

PairIdx        = 1;
phase          = phaseExtraction(tarCh(PairIdx),freqRange);

t = (0:1/fs:(1/fs)*length(phase))';

if strcmp(medianSet,'all')
    data = load(['/home/nasko/res_' num2str(pairs{PairIdx}(2)) 'i.mat']);
    
    try
        startTime = data.res1(1)*fs;
        resTar  = data.res1;
        resTarOrig  = data.res1;
    catch
        startTime = data.res2(1)*fs;
        resTar  = data.res2;
        resTarOrig  = data.res2;
    end
        
    resTar  = resTar*fs;
    resTar  = resTar-startTime;
    resTar  = round(resTar) + 1;
    
elseif strcmp(medianSet,'sync')
    
    data = load(['/home/nasko/res_' num2str(pairs{PairIdx}(2)) 'i.mat']);
    
    try
        startTime = data.res1(1)*fs;
    catch
        startTime = data.res2(1)*fs;
    end
    
    data = load(['/home/nasko/res_' num2str(pairs{PairIdx}(2)) 'i_sync.mat']);
    
    resTar  = data.res1sync;
    resTarOrig  = data.res1sync;
    resTar  = resTar*fs;
    resTar  = resTar-startTime;
    resTar  = round(resTar) + 1;
    
elseif  strcmp(medianSet,'below')
%     data = load(['res_' num2str(pairs{PairIdx}(2)) 'i_' medianSet 'Med.mat']);
%     
%     startTime = data.resBelowMed(1)*fs;
% 
%     resTar  = data.resBelowMed;
%     resTarOrig  = data.resBelowMed;
%     resTar  = resTar*fs;
%     resTar  = resTar-startTime;
%     resTar  = round(resTar) + 1;
    
    data = medSplit.tarBelowMed{2};
    
    startTime = data(1)*fs;
    
    resTar = data;
    resTarOrig = data;
    resTar  = resTar*fs;
    resTar  = resTar-startTime;
    resTar  = round(resTar) + 1;
    
elseif strcmp(medianSet,'above')
%     data = load(['res_' num2str(pairs{PairIdx}(2)) 'i_' medianSet 'Med.mat']);
%     
%     startTime = data.resAboveMed(1)*fs;
% 
%     resTar  = data.resAboveMed;
%     resTarOrig  = data.resAboveMed;
%     resTar  = resTar*fs;
%     resTar  = resTar-startTime;
%     resTar  = round(resTar) + 1;
    
    data = medSplit.tarAboveMed{2};
    
    startTime = data(1)*fs;
    
    resTar = data;
    resTarOrig = data;
    resTar  = resTar*fs;
    resTar  = resTar-startTime;
    resTar  = round(resTar) + 1;
    
end

data = load(['/home/nasko/res_' num2str(pairs{PairIdx}(1)) 'i.mat']);

try
    resRef  = data.res1;
catch
    resRef  = data.res2;
end

clear data

noSamples = length(t);

spikeBinTar = zeros(noSamples,1);
spikeBinTar(resTar) = 1;

if winSize == 0
    spikeCount = spikeBinTar;
else
    spikeCount = movsum(spikeBinTar,winSize); % how to normalize here?
end

phaseLabels = zeros(noSamples,1);

figure;

%%

phaseBins = -pi:pi/4:pi;
binCounts = zeros(length(phaseBins)-1,1);

for i = 1:length(phase)

    if ((phase(i) > phaseBins(1)) && (phase(i) < phaseBins(2)))
        phaseLabels(i) = 1;
    elseif ((phase(i) > phaseBins(2)) && (phase(i) < phaseBins(3)))
        phaseLabels(i) = 2;
    elseif ((phase(i) > phaseBins(3)) && (phase(i) < phaseBins(4)))
        phaseLabels(i) = 3;
    elseif ((phase(i) > phaseBins(4)) && (phase(i) < phaseBins(5)))
        phaseLabels(i) = 4;
    elseif ((phase(i) > phaseBins(5)) && (phase(i) < phaseBins(6)))
        phaseLabels(i) = 5;
    elseif ((phase(i) > phaseBins(6)) && (phase(i) < phaseBins(7)))
        phaseLabels(i) = 6;
    elseif ((phase(i) > phaseBins(7)) && (phase(i) < phaseBins(8)))
        phaseLabels(i) = 7;
    elseif ((phase(i) > phaseBins(8)) && (phase(i) < phaseBins(9)))
        phaseLabels(i) = 8;
    end

end


%%

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 1) = 1;
% FR20iGroupPhase1 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase1 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 2) = 1;
% FR20iGroupPhase2 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase2 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 3) = 1;
% FR20iGroupPhase3 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase3 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 4) = 1;
% FR20iGroupPhase4 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase4 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 5) = 1;
% FR20iGroupPhase5 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase5 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 6) = 1;
% FR20iGroupPhase6 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase6 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 7) = 1;
% FR20iGroupPhase7 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase7 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseLabelsTemp = zeros(noSamples,1);
phaseLabelsTemp(phaseLabels == 8) = 1;
% FR20iGroupPhase8 = ((mean(spikeCount(find(phaseLabelsTemp))) - mean(spikeCount))/mean(spikeCount))*100;
FR20iGroupPhase8 = ((mean(spikeCount(find(phaseLabelsTemp))))/mean(spikeCount))*100;

phaseData = [FR20iGroupPhase1 FR20iGroupPhase2 FR20iGroupPhase4 FR20iGroupPhase4 ...
             FR20iGroupPhase5 FR20iGroupPhase6 FR20iGroupPhase7 FR20iGroupPhase8];

plot(phaseData,'o','LineWidth',2)


  
%%

n = length(phaseData); 
dtheta = 360/n;         % sector for each country
gap = 0;              % gaps betwen adjacent contries
npoints = 21;           % npoints for each sector
figure; 
hold on;
for i=1:n
    theta = linspace((i-1)*dtheta, (i-gap)*dtheta, npoints);
    patch([0 (phaseData(i)+1)*cosd(theta)], [0 (phaseData(i)+1)*sind(theta)], 'b');
    rtext = phaseData(i)+1;
    text(rtext*cosd((i-.5)*dtheta), rtext*sind((i-.5)*dtheta), num2str(mean(phaseBins(i:i+1))));
end
axis equal
axis off
  
% for i = 1:8
%     
%     patch([0 sind(rad2deg(mean(phaseBins(i:1+i))))],[0 cosd(rad2deg(mean(phaseBins(i:1+i))))],'b')
%     
% end
  
%%

% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 1) = 1;
% res20iGroupPhase1 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 2) = 1;
% res20iGroupPhase2 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 3) = 1;
% res20iGroupPhase3 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 4) = 1;
% res20iGroupPhase4 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 5) = 1;
% res20iGroupPhase5 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 6) = 1;
% res20iGroupPhase6 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 7) = 1;
% res20iGroupPhase7 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% phaseLabelsTemp = zeros(noSamples,1);
% phaseLabelsTemp(phaseLabels == 8) = 1;
% res20iGroupPhase8 = find(spikeBinTar & phaseLabelsTemp) + startTime;
% 
% res20iGroupStr{1} = res20iGroupPhase1/fs;
% res20iGroupStr{2} = res20iGroupPhase2/fs;
% res20iGroupStr{3} = res20iGroupPhase3/fs;
% res20iGroupStr{4} = res20iGroupPhase4/fs;
% res20iGroupStr{5} = res20iGroupPhase5/fs;
% res20iGroupStr{6} = res20iGroupPhase6/fs;
% res20iGroupStr{7} = res20iGroupPhase7/fs;
% res20iGroupStr{8} = res20iGroupPhase8/fs;
% 
% [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%             CCG_jitter(resTarOrig,res45i,fs,binSize,duration,'jscale',jscale, ...
%                         'plot_flag', false, ...
%                         'plot_output', get(fig_use, 'Number'), ...
%                         'njitter', njitter, 'alpha', alpha,...
%                         'for_grant', for_grant);
%                     
% CCGall = ccgR(:,1,2)/8;
% 
% for i = 1:8
%     [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
%             CCG_jitter(res20iGroupStr{i},res45i,fs,binSize,duration,'jscale',jscale, ...
%                         'plot_flag', false, ...
%                         'plot_output', get(fig_use, 'Number'), ...
%                         'njitter', njitter, 'alpha', alpha,...
%                         'for_grant', for_grant);
%     
%     CCGphase(i,:) = ccgR(:,1,2);
%     
% end
% 
% tiledlayout(4,2)
% nexttile
% plot(tR*1000,[CCGall CCGphase(1,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(1)) ' to ' num2str(phaseBins(2))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(2,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(2)) ' to ' num2str(phaseBins(3))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(3,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(3)) ' to ' num2str(phaseBins(4))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(4,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(4)) ' to ' num2str(phaseBins(5))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(5,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(5)) ' to ' num2str(phaseBins(6))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(6,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(6)) ' to ' num2str(phaseBins(7))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(7,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(7)) ' to ' num2str(phaseBins(8))])
% nexttile
% plot(tR*1000,[CCGall CCGphase(8,:)'],'LineWidth',1.5)
% title(['phase ' num2str(phaseBins(8)) ' to ' num2str(phaseBins(9))])



