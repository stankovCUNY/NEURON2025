load('RoyMaze1GroupStats.mat')

distancePlot = [groupStat.preDistance        groupStat.mainDistance       groupStat.postDistance];
CCGplot      = [groupStat.preFeatCCGall      groupStat.mainFeatCCGall     groupStat.postFeatCCGall]; 
CCGplotNorm  = [groupStat.preFeatCCGnormAll  groupStat.mainFeatCCGnormAll groupStat.postFeatCCGnormAll];
CCGplotAgmon = [groupStat.preFeatCCGagmonIdx groupStat.mainPkAgmonIdx     groupStat.postFeatCCGagmonIdx];
WAVEplot     = [groupStat.preFeatWAVEall     groupStat.mainFeatWAVEall    groupStat.postFeatWAVEall];
CCGplotLat   = [groupStat.preFeatCCGlatAll   groupStat.mainFeatCCGlatAll  groupStat.postFeatCCGlatAll];
WAVEplotLat  = [groupStat.preFeatWAVElatAll  groupStat.mainFeatWAVElatAll groupStat.postFeatWAVElatAll];

tiledlayout(2,3)

%
nexttile
scatter(WAVEplotLat,CCGplotLat,'LineWidth',2)
xlabel('Waveform feature latency [microsec.]')
ylabel('CCG feature latency [microsec.]')

%
nexttile
% scatter(groupStat.mainDistance,groupStat.mainFeatCCGnormAll,'LineWidth',2)
scatter(groupStat.mainDistance,groupStat.mainPkAgmonIdx,'LineWidth',2)
xlabel('Pair distance - max ref channel to max correlation channel [micron]')
ylabel('Normalized CCG peak amp. [a.u.]')

%
nexttile
% scatter([groupStat.preDistance groupStat.postDistance],[groupStat.preFeatCCGnormAll groupStat.postFeatCCGnormAll],'LineWidth',2)
scatter([groupStat.preDistance groupStat.postDistance],[groupStat.preFeatCCGagmonIdx groupStat.postFeatCCGagmonIdx],'LineWidth',2)
xlabel('Pair distance - max ref channel to max correlation channel [micron]')
ylabel('Normalized CCG trough amp. [a.u.]')

%
nexttile
scatter([groupStat.preDistance groupStat.postDistance],[groupStat.preFeatWAVEall groupStat.postFeatWAVEall],'LineWidth',2)
xlabel('Pair distance - max ref channel to max correlation channel [micron]')
ylabel('Highpass filt. waveform peak amp. [mV]')

%
nexttile
scatter(groupStat.mainDistance,groupStat.mainFeatWAVEall,'LineWidth',2)
xlabel('Pair distance - max ref channel to max correlation channel [micron]')
ylabel('Highpass filt. waveform trough amp. [mV]')

%
nexttile
scatter(WAVEplot,CCGplotNorm,'LineWidth',2)
xlabel('Waveform feature amp. [mV]')
ylabel('Normalized CCG feature amp. [a.u.]')

%%

