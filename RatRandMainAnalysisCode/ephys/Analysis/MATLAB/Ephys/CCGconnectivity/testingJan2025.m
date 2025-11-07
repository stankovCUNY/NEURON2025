clear

%% Hiro's isolation distance for Roy

% spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
% dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
% 
% [data_dir, name, ~] = fileparts(spike_data_fullpath);
% 
% load(spike_data_fullpath, 'spikes')
% load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
% load(fullfile(data_dir, 'wake-basics.mat'),   'basics');
% 
% load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')
% 
% cellIDs         = 1:size(spikes.RoyMaze1,2);
% neuronTemp      = cell2num({spikes.RoyMaze1.id}');
% cellTypeList = strings(size(spikes.RoyMaze1,2),1);
% cellTypeList([spikes.RoyMaze1.quality] == 8) = "i";
% cellTypeList(neuronTemp(:,2) == 1)           = "noise";
% cellTypeList(strcmp(cellTypeList,""))        = "p"; % double check!!
% 
% isolationDistanceList = [spikes.RoyMaze1.isoDist];
% 
% % get exq lists
% load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
% 
% exqDlist = unique([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID]);
% 
% [~,idxTemp] = intersect(cellIDs,exqDlist);
% 
% isolationDistanceExqList = isolationDistanceList(idxTemp);
% cellTypeExqList          = cellTypeList(idxTemp);
% 
% % pyramids
% allUnits = isolationDistanceList(strcmp(cellTypeList,'p'));
% exqUnits = isolationDistanceExqList(strcmp(cellTypeExqList,'p'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',10)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',10)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 500])
% ylabel('probability')
% title(['Rat Roy pyramids iso. distance (from Hiro), d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% % 
% % % interneuron
% allUnits = isolationDistanceList(strcmp(cellTypeList,'i'));
% exqUnits = isolationDistanceExqList(strcmp(cellTypeExqList,'i'));
% nexttile
% histogram(isolationDistanceList(   strcmp(cellTypeList,   'i')),'Normalization','probability','BinWidth',10)
% hold on 
% histogram(isolationDistanceExqList(strcmp(cellTypeExqList,'i')),'Normalization','probability','BinWidth',10)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 500])
% ylabel('probability')
% title(['Rat Roy interneurons iso. distance (from Hiro), d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')

%% Allen Int cluster quality ana

% datapath = '/media/nasko/WD_BLACK31/BOTtemp/';
% 
% connType = 'exq';
% 
% animalsList = [715093703
%                719161530
%                721123822
%                732592105
%                737581020
%                739448407
%                742951821
%                743475441
%                744228101
%                746083955
%                750332458
%                750749662
%                751348571
%                754312389
%                754829445
%                755434585
%                756029989
%                757216464
%                757970808
%                758798717
%                759883607
%                760345702
%                760693773
%                761418226
%                762120172
%                762602078
%                763673393
%                766640955
%                767871931
%                768515987
%                771160300
%                771990200
%                773418906
%                774875821
%                778240327
%                778998620
%                779839471
%                781842082
%                786091066
%                787025148
%                789848216
%                791319847
%                793224716
%                794812542
%                797828357
%                798911424
%                799864342
%                816200189
%                819186360
%                819701982
%                821695405
%                829720705
%                831882777
%                835479236
%                839068429
%                839557629
%                840012044
%                847657808];
% 
% %% Allen SNR all units
%            
% refSNRexqList               = [];
% tarSNRexqList               = [];
% 
% refLratioExqList            = [];
% tarLratioExqList            = [];
% 
% refDprimeExqList            = [];
% tarDprimeExqList            = [];
% 
% refIsolationDistanceExqList = [];
% tarIsolationDistanceExqList = [];
% 
% refIDeqxList                = [];
% tarIDeqxList                = [];
% 
% refTypeEqxList = {};
% tarTypeEqxList = {};
% 
% for loopAnimals = 1:58
%     
%     display(['animal: '  num2str(animalsList(loopAnimals))])
% 
%     %%
%     if strcmp(connType,'exq')
%         load([datapath 'groupStatsMouse' num2str(animalsList(loopAnimals))])
%     elseif strcmp(connType,'GJ')
%         load([datapath 'groupStatsMousePutativeGJ' num2str(animalsList(loopAnimals))])
%     end
%     
%     pairGroupStatTable = pairGroupStatsTable;
%     clear pairGroupStatsTable
% 
%      % remove duplicates
%     [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
%     pairGroupStatTable = pairGroupStatTable(idx,:);
%     
%     load(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/neuronsMouse' num2str(animalsList(loopAnimals)) 'qualityMetrics'])
%     
%     for loopPairs = 1:size(pairGroupStatTable,1)
% 
%         pairGroupStatTable.refLratio(loopPairs)            = neurons.Lratio(neurons.unitID == pairGroupStatTable.refNeuronID(loopPairs));
%         pairGroupStatTable.tarLratio(loopPairs)            = neurons.Lratio(neurons.unitID == pairGroupStatTable.tarNeuronID(loopPairs));
% 
%         pairGroupStatTable.refDprime(loopPairs)            = neurons.Dprime(neurons.unitID == pairGroupStatTable.refNeuronID(loopPairs));
%         pairGroupStatTable.tarDprime(loopPairs)            = neurons.Dprime(neurons.unitID == pairGroupStatTable.tarNeuronID(loopPairs));
% 
%         pairGroupStatTable.refIsolationDistance(loopPairs) = neurons.isolationDistance(neurons.unitID == pairGroupStatTable.refNeuronID(loopPairs));
%         pairGroupStatTable.tarIsolationDistance(loopPairs) = neurons.isolationDistance(neurons.unitID == pairGroupStatTable.tarNeuronID(loopPairs));
%         
%     end
%     
%     refSNRexqList               = [refSNRexqList; pairGroupStatTable.refSNR];
%     tarSNRexqList               = [tarSNRexqList; pairGroupStatTable.tarSNR];
%     
%     refLratioExqList            = [refLratioExqList; pairGroupStatTable.refLratio];
%     tarLratioExqList            = [tarLratioExqList; pairGroupStatTable.tarLratio];
% 
%     refDprimeExqList            = [refDprimeExqList; pairGroupStatTable.refDprime];
%     tarDprimeExqList            = [tarDprimeExqList; pairGroupStatTable.tarDprime];
% 
%     refIsolationDistanceExqList = [refIsolationDistanceExqList; pairGroupStatTable.refIsolationDistance];
%     tarIsolationDistanceExqList = [tarIsolationDistanceExqList; pairGroupStatTable.tarIsolationDistance];
%     
%     refIDeqxList                = [refIDeqxList; pairGroupStatTable.refNeuronID];
%     tarIDeqxList                = [tarIDeqxList; pairGroupStatTable.tarNeuronID];
%     
%     refTypeEqxList              = [refTypeEqxList; pairGroupStatTable.refCellExplorerType];
%     tarTypeEqxList              = [tarTypeEqxList; pairGroupStatTable.tarCellExplorerType];
% 
% end
% 
% SNRexqList               = [refSNRexqList;  tarSNRexqList];
% LratioExqList            = [refLratioExqList; tarLratioExqList];
% DprimeExqList            = [refDprimeExqList; tarDprimeExqList];
% isolationDistanceExqList = [refIsolationDistanceExqList; tarIsolationDistanceExqList];
% cellTypeEqxList          = [refTypeEqxList; tarTypeEqxList];
% 
% [~,idxTemp] = unique([refIDeqxList; tarIDeqxList]);
% 
% SNRexqList               = SNRexqList(idxTemp);
% LratioExqList            = LratioExqList(idxTemp);
% DprimeExqList            = DprimeExqList(idxTemp);
% isolationDistanceExqList = isolationDistanceExqList(idxTemp); 
% cellTypeEqxList          = cellTypeEqxList(idxTemp); 
% 
% (sum(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'))/length(cellTypeEqxList))*100
% (sum(strcmp(cellTypeEqxList,'i-narrow'))/length(cellTypeEqxList))*100
% 
% %% Allen SNR exq units
% 
% datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
% 
% SNRlist               = [];
% LratioList            = [];
% DprimeList            = [];
% isolationDistanceList = [];
% cellTypeList          = {};
% 
% for loopAnimals = 1:58
%     
%     display(['animal: '  num2str(animalsList(loopAnimals))])
%         
%     load([datapath 'neuronsMouse' num2str(animalsList(loopAnimals)) 'qualityMetrics']);
%     
%     SNRlist               = [SNRlist;                  neurons.snr];
%     LratioList            = [LratioList;               neurons.Lratio];
%     DprimeList            = [DprimeList;               neurons.Dprime];
%     isolationDistanceList = [isolationDistanceList;    neurons.isolationDistance];
%     cellTypeList          = [cellTypeList;             neurons.putativeCellType];
%     
% end
% 
% (sum(strcmp(cellTypeList,'p') | strcmp(cellTypeList,'i-wide'))/length(cellTypeList))*100
% (sum(strcmp(cellTypeList,'i-narrow'))/length(cellTypeList))*100
% 
% tiledlayout(4,2)
% 
% %% Allen SNR
% 
% % pyramids
% allUnits = SNRlist(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
% exqUnits = SNRexqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',0.1)
% hold off
% d = (mean(allUnits)-mean(exqUnits))/std(allUnits);
% xlim([0 10])
% ylabel('probability')
% title(['Allen Int. pyramids SNR, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% % interneuron
% allUnits = SNRlist(   strcmp(cellTypeList,   'i-narrow'));
% exqUnits = SNRexqList(strcmp(cellTypeEqxList,'i-narrow'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',0.1)
% hold off
% d = (mean(allUnits)-mean(exqUnits))/std(allUnits);
% xlim([0 10])
% ylabel('probability')
% title(['Allen Int. interneurons SNR, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% % [P,H] = ranksum(SNRlist,SNRexqList) 
% 
% %% Allen Lratio
% 
% % pyramids
% allUnits = LratioList(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
% exqUnits = LratioExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',0.001)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',0.001)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 0.05])
% ylabel('probability')
% title(['Allen Int. pyramids L ratio, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% % interneuron
% allUnits = LratioList(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
% exqUnits = LratioExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
% nexttile
% histogram(LratioList(   strcmp(cellTypeList,   'i-narrow')),'Normalization','probability','BinWidth',0.001)
% hold on 
% histogram(LratioExqList(strcmp(cellTypeEqxList,'i-narrow')),'Normalization','probability','BinWidth',0.001)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 0.05])
% ylabel('probability')
% title(['Allen Int. interneurons L ratio, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% 
% %% Allen Dprime
% 
% % pyramids
% allUnits = DprimeList(   strcmp(cellTypeList,'p')    | strcmp(cellTypeList,'i-wide'));
% exqUnits = DprimeExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',0.5)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',0.5)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 20])
% ylabel('probability')
% title(['Allen Int. pyramids D prime, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% % interneuron
% allUnits = DprimeList(   strcmp(cellTypeList,   'i-narrow'));
% exqUnits = DprimeExqList(strcmp(cellTypeEqxList,'i-narrow'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',0.5)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',0.5)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 20])
% ylabel('probability')
% title(['Allen Int. interneurons D prime, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% %% Allen isolation distance
% 
% % pyramids
% allUnits = isolationDistanceList(   strcmp(cellTypeList,'p') |    strcmp(cellTypeList,'i-wide'));
% exqUnits = isolationDistanceExqList(strcmp(cellTypeEqxList,'p') | strcmp(cellTypeEqxList,'i-wide'));
% nexttile
% histogram(allUnits,'Normalization','probability','BinWidth',10)
% hold on 
% histogram(exqUnits,'Normalization','probability','BinWidth',10)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 500])
% ylabel('probability')
% title(['Allen Int. pyramids iso. distance, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% 
% % interneuron
% allUnits = isolationDistanceList(   strcmp(cellTypeList,   'i-narrow'));
% exqUnits = isolationDistanceExqList(strcmp(cellTypeEqxList,'i-narrow'));
% nexttile
% histogram(isolationDistanceList(   strcmp(cellTypeList,   'i-narrow')),'Normalization','probability','BinWidth',10)
% hold on 
% histogram(isolationDistanceExqList(strcmp(cellTypeEqxList,'i-narrow')),'Normalization','probability','BinWidth',10)
% hold off
% d = (mean(allUnits,'omitnan')-mean(exqUnits,'omitnan'))/std(allUnits,'omitnan');
% xlim([0 500])
% ylabel('probability')
% title(['Allen Int. interneurons iso. distance, d = ' num2str(d)])
% legend('all units','exq units')
% box off
% set(gca,'FontSize',5)
% set(gca,'FontName','Arial')
% % 
% % 
% % %% phy purity for rats AG,N,S,U units
% % if loopAnimal == 1
% %     load('RatNprocessedGroupStats.mat')
% %     clusterInfo = readtable(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(1) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t');
% %     for loopShanks = 2:8
% %         clusterInfo  = [clusterInfo; readtable(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t')];
% %     end
% % elseif loopAnimal == 2
% %     load('RatSprocessedGroupStats.mat')
% %     clusterInfo      = readtable('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/cluster_info.tsv', "FileType","text",'Delimiter', '\t');
% % elseif loopAnimal == 3
% %     load('RatUprocessedGroupStats.mat')
% %     clusterInfo      = readtable('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/cluster_info.tsv', "FileType","text",'Delimiter', '\t');
% % elseif loopAnimal == 4
% %     load('RatAGday1processedGroupStats.mat')
% %     clusterInfo  = readtable(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(1) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t');
% %     for loopShanks = 2:6
% %         clusterInfo  = [clusterInfo; readtable(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t')];
% %     end
% % elseif loopAnimal == 5
% %     load('RatAGday2processedGroupStats.mat')
% %     clusterInfo  = readtable(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(1) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t');
% %     for loopShanks = 2:5
% %         clusterInfo  = [clusterInfo; readtable(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t')];
% %     end
% % end
% % 
% % clusterInfo = clusterInfo(strcmp(clusterInfo.group,'good'),:);
% % 
% % phyPurity = clusterInfo.purity;
% % phyPurity(isnan(phyPurity)) = [];
% % 
% % histogram(phyPurity,'Normalization','probability','BinWidth',0.05)
% % 
% % %%
% % 
% % load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy-maze1/Roy-maze1.clusterQuality.mat')

%% prep data for collecting PC vals data

% animalNames = {'N','S','U','AG1','AG2','Roy'};
% 
% fs         = 30000;
% 
% loopAnimal = 1;
% 
% if loopAnimal == 1
% 
%     neuronsN    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
% 
%     probegroupN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
% 
%     channelID = double(cell2num(probegroupN.channel_id));
% 
%     channelIDfilt = channelID(cell2num(probegroupN.connected));
% 
%     shankIDfilt = probegroupN.shank_id(cell2num(probegroupN.connected));
% 
%     cellTypeList     =         neuronsN.neuron_type(  ~(sum(neuronsN.neuron_type(:,:)=='mua  ',2) == 5),:);
%     neuronChannel    =  double(neuronsN.peak_channels(~(sum(neuronsN.neuron_type(:,:)=='mua  ',2) == 5))) + 1;
%     neuronShank      =  double(neuronsN.shank_ids(    ~(sum(neuronsN.neuron_type(:,:)=='mua  ',2) == 5))) + 1;
%     neuronSpikeTimes = {neuronsN.spiketrains{1,       ~(sum(neuronsN.neuron_type(:,:)=='mua  ',2) == 5)}};
% 
%     cellIDs          =  double(neuronsN.neuron_ids) + 1;
%     cellTypeList     =         neuronsN.neuron_type;
%     neuronChannel    =  double(neuronsN.peak_channels) + 1;
%     neuronShank      =  double(neuronsN.shank_ids) + 1;
%     neuronSpikeTimes =         neuronsN.spiketrains;
%     for loopNeurons = 1:size(neuronSpikeTimes,2)
%         neuronNoSpikes(loopNeurons)    = length(neuronSpikeTimes{loopNeurons});
%     end
%     noiseClusterSpikeLat = {};
% 
%     channel_map_xml = {[2   1   0   12  3   11  4   10  5   9   6   8   7], ... 
%                        [29  18  30  17  31  16  28  19  27  20  26  21  25  22  24  23], ...
%                        [45  34  46  33  47  32  44  35  43  36  42  37  41  38  40  39], ...
%                        [61  50  62  49  48  60  51  59  52  58  53  57  54  56  55], ...
%                        [77  66  78  65  79  64  76  67  75  68  74  69  73  70  72  71], ...
%                        [93  82  94  81  95  80  92  83  91  84  90  85  89  86  88  87], ...
%                        [109 98  110 97  111 96  108 99  107 100 106 101 105 102 104 103] ,...
%                        [125 114 126 113 127 112 124 115 123 116 122 117 121 118 120 119]};
% 
%     channel_map_adjusted = channel_map_xml;
% 
%     for loopShanks = 1:8        
%         allSpikeLat{loopShanks}      = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/spike_times.npy']));
%         spikeTemplates{loopShanks}   = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/spike_templates.npy']));
%         templateChannels{loopShanks} = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/pc_feature_ind.npy']));
%         channel_map{loopShanks}      = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/channel_map.npy']));
% 
%         cluster_info     = readtable(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t');
%         noiseClusterIdx  = find(strcmp(cluster_info.group,'noise')); 
% 
%         templateList     = unique(spikeTemplates{loopShanks});
%         spikeNoTemplates = histc(spikeTemplates{loopShanks},templateList);
% 
%         noiseClusterType = string(cluster_info.group(noiseClusterIdx));
%         noiseClusterChan = cluster_info.ch(noiseClusterIdx) + 1;
% %         noiseClusterChan = channel_map_xml{loopShanks};
%         noiseClusterShank = ones(length(noiseClusterChan),1)*loopShanks;
% 
%         noiseClusterID   = [];
% 
%         lostNoiseUnitsIdx = [];
%         for loopNoiseCluster = 1:length(noiseClusterIdx)
%             tempIdx = find(cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster)) == spikeNoTemplates);
% 
%             if isempty(tempIdx) %|| (cluster_info.n_spikes(loopNoiseCluster) < min(spikeNoTemplates))
% 
%                 lostNoiseUnitsIdx = [lostNoiseUnitsIdx loopNoiseCluster];
% 
%                 continue
%             end
%             noiseClusterID       = [noiseClusterID str2double(['999' num2str(templateList(tempIdx))])]; % code noise units with 999 prefix
%             noiseClusterSpikeLat = [noiseClusterSpikeLat allSpikeLat{loopShanks}(    templateList(tempIdx) == spikeTemplates{loopShanks})];
%         end
% 
%         noiseClusterIdx(lostNoiseUnitsIdx)   = [];
%         noiseClusterType(lostNoiseUnitsIdx)  = [];
%         noiseClusterChan(lostNoiseUnitsIdx)  = [];
%         noiseClusterShank(lostNoiseUnitsIdx) = [];
% 
%         if isempty(intersect(neuronNoSpikes,cluster_info.n_spikes(noiseClusterIdx)))
%             display("Sanity check works!!")
%         end
% 
%         zerosIdx      = templateChannels{loopShanks} == 0;
%         zerosIdx(:,1) = 0;
%         templateChannels{loopShanks}(zerosIdx) = NaN;
% 
%         cellIDs          =  [cellIDs               noiseClusterID];
%         cellTypeList     =  [string(cellTypeList);  noiseClusterType];
%         neuronChannel    =  [neuronChannel         noiseClusterChan'];
%         neuronShank      =  [neuronShank           noiseClusterShank'];
% 
% %         % adjust channel_map
% %         channel_map_adjusted{1} = channel_map{1};
% %         if loopShanks > 1
% %             adjust(loopShanks-1) = length(channel_map{loopShanks-1}); 
% %             channel_map_adjusted{loopShanks} = channel_map{loopShanks} + sum(adjust(loopShanks-1));
% %         end
% % 
% %         channel_map_adjusted{loopShanks} = sort(channelIDfilt((loopShanks - 1) == shankIDfilt));
% 
%     end
% 
% elseif loopAnimal == 2
% 
%     neuronsS    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
% 
%     probegroupS = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');
% 
%     cellTypeList     =         neuronsS.neuron_type(  ~(sum(neuronsS.neuron_type(:,:)=='mua  ',2) == 5),:);
%     neuronChannel    =  double(neuronsS.peak_channels(~(sum(neuronsS.neuron_type(:,:)=='mua  ',2) == 5))) + 1;
%     neuronShank      =  double(neuronsS.shank_ids(    ~(sum(neuronsS.neuron_type(:,:)=='mua  ',2) == 5))) + 1;
%     neuronSpikeTimes = {neuronsS.spiketrains{1,       ~(sum(neuronsS.neuron_type(:,:)=='mua  ',2) == 5)}};
% 
%     cellIDs          =  double(neuronsS.neuron_ids) + 1;
%     cellTypeList     =         neuronsS.neuron_type;
%     neuronChannel    =  double(neuronsS.peak_channels) + 1;
%     neuronShank      =  double(neuronsS.shank_ids) + 1;
%     neuronSpikeTimes =         neuronsS.spiketrains;
%     for loopNeurons = 1:size(neuronSpikeTimes,2)
%         neuronNoSpikes(loopNeurons)    = length(neuronSpikeTimes{loopNeurons});
%     end
% 
%     allSpikeLat      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/spike_times.npy'));
%     spikeTemplates   = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/spike_templates.npy'));
%     templateChannels = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/pc_feature_ind.npy'));
%     channel_map      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/channel_map.npy'));
% 
%     cluster_info     = readtable("/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/cluster_info.tsv", "FileType","text",'Delimiter', '\t');
%     noiseClusterIdx  = find(strcmp(cluster_info.group,'noise')); 
% 
%     templateList     = unique(spikeTemplates);
%     spikeNoTemplates = histc(spikeTemplates,templateList);
% 
%     noiseClusterType = string(cluster_info.group(noiseClusterIdx));
%     noiseClusterChan = cluster_info.ch(noiseClusterIdx) + 1;
% 
%     noiseClusterID   = [];
% 
%     lostNoiseUnitsIdx = [];
%     for loopNoiseCluster = 1:length(noiseClusterIdx)
%         tempIdx = find(cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster)) == spikeNoTemplates);
% 
%         if isempty(tempIdx) %|| (cluster_info.n_spikes(loopNoiseCluster) < min(spikeNoTemplates))
% 
%             lostNoiseUnitsIdx = [lostNoiseUnitsIdx loopNoiseCluster];
% 
%             continue
%         end
%         noiseClusterID                         = [noiseClusterID str2double(['999' num2str(templateList(tempIdx))])]; % code noise units with 999 prefix
%         noiseClusterSpikeLat{loopNoiseCluster} = allSpikeLat(    templateList(tempIdx) == spikeTemplates);
%     end
% 
%     noiseClusterIdx(lostNoiseUnitsIdx)  = [];
%     noiseClusterType(lostNoiseUnitsIdx) = [];
%     noiseClusterChan(lostNoiseUnitsIdx) = [];
% 
%     if isempty(intersect(neuronNoSpikes,cluster_info.n_spikes(noiseClusterIdx)))
%         display("Sanity check works!!")
%     end
% %     noiseSpikeTimes  = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/Bapuns_RatS_spikeTimesNoise.mat');
% 
%     cellIDs          =  [cellIDs               noiseClusterID];
%     cellTypeList     =  [string(cellTypeList); noiseClusterType];
%     neuronChannel    =  [neuronChannel';       noiseClusterChan];
%     neuronShank      =  [neuronShank';        (probegroupS.shank_id(noiseClusterChan) + 1)'];
% 
%     zerosIdx      = templateChannels == 0;
%     zerosIdx(:,1) = 0;
%     templateChannels(zerosIdx) = NaN;
% 
% elseif loopAnimal == 3
% 
%     neuronsU    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
% 
%     probegroupU = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');
% 
%     cellTypeList     =         neuronsU.neuron_type(  ~(sum(neuronsU.neuron_type(:,:)=='mua  ',2) == 5),:);
%     neuronChannel    =  double(neuronsU.peak_channels(~(sum(neuronsU.neuron_type(:,:)=='mua  ',2) == 5))) + 1;
%     neuronShank      =  double(neuronsU.shank_ids(    ~(sum(neuronsU.neuron_type(:,:)=='mua  ',2) == 5))) + 1;
%     neuronSpikeTimes = {neuronsU.spiketrains{1,       ~(sum(neuronsU.neuron_type(:,:)=='mua  ',2) == 5)}};
% 
%     cellIDs          =  double(neuronsU.neuron_ids) + 1;
%     cellTypeList     =         neuronsU.neuron_type;
%     neuronChannel    =  double(neuronsU.peak_channels) + 1;
%     neuronShank      =  double(neuronsU.shank_ids) + 1;
%     neuronSpikeTimes =         neuronsU.spiketrains;
%     for loopNeurons = 1:size(neuronSpikeTimes,2)
%         neuronNoSpikes(loopNeurons)    = length(neuronSpikeTimes{loopNeurons});
%     end
%     noiseClusterSpikeLat = {};
% 
%     allSpikeLat      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/spike_times.npy'));
%     spikeTemplates   = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/spike_templates.npy'));
%     templateChannels = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/pc_feature_ind.npy'));
%     channel_map      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/channel_map.npy'));
% 
%     cluster_info     = readtable("/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/cluster_info.tsv", "FileType","text",'Delimiter', '\t');
%     noiseClusterIdx  = find(strcmp(cluster_info.group,'noise'));
% 
%     templateList     = unique(spikeTemplates);
%     spikeNoTemplates  = histc(spikeTemplates,templateList);
% 
%     noiseClusterType = string(cluster_info.group(noiseClusterIdx));
%     noiseClusterChan = cluster_info.ch(noiseClusterIdx) + 1;
% 
%     noiseClusterID   = [];
% 
%     lostNoiseUnitsIdx = [];
%     for loopNoiseCluster = 1:length(noiseClusterIdx)
%         tempIdx = find(cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster)) == spikeNoTemplates);
% 
%         if isempty(tempIdx) || ~isempty(intersect(neuronNoSpikes,cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster))))
% 
%             lostNoiseUnitsIdx = [lostNoiseUnitsIdx loopNoiseCluster];
% 
%             continue
%         end
% 
%         tempIdx = tempIdx(1);
% 
%         noiseClusterID                         = [noiseClusterID str2double(['999' num2str(templateList(tempIdx))])]; % code noise units with 999 prefix
%         noiseClusterSpikeLat{loopNoiseCluster} = allSpikeLat(    templateList(tempIdx) == spikeTemplates);
%     end
% 
%     noiseClusterIdx(lostNoiseUnitsIdx)  = [];
%     noiseClusterType(lostNoiseUnitsIdx) = [];
%     noiseClusterChan(lostNoiseUnitsIdx) = [];
% 
%     if isempty(intersect(neuronNoSpikes,cluster_info.n_spikes(noiseClusterIdx)))
%         display("Sanity check works!!")   
%     end
% 
% %     noiseSpikeTimes  = load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/Bapuns_RatU_spikeTimesNoise.mat');
% 
%     cellIDs          =  [cellIDs               noiseClusterID];
%     cellTypeList     =  [string(cellTypeList); noiseClusterType];
%     neuronChannel    =  [neuronChannel';       noiseClusterChan];
%     neuronShank      =  [neuronShank';        (probegroupU.shank_id(noiseClusterChan) + 1)'];
% 
%     zerosIdx      = templateChannels == 0;
%     zerosIdx(:,1) = 0;
%     templateChannels(zerosIdx) = NaN;
% 
% elseif loopAnimal == 4
% 
%     UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
%     load([UtkuPath 'AG_2019-12-23_NSD' '/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);
% 
%     cellIDs          = cell_metrics.UID;
%     cellTypeList     = cell_metrics.putativeCellType;
%     neuronChannel    = cell_metrics.maxWaveformCh + 1;
%     neuronShank      = cell_metrics.shankID;
%     neuronSpikeTimes = cell_metrics.spikes.times;
%     for loopNeurons = 1:size(neuronSpikeTimes,2)
%         neuronNoSpikes(loopNeurons)    = length(neuronSpikeTimes{loopNeurons});
%     end
%     noiseClusterSpikeLat = {};
% 
% %     channel_map_noise = [18 12 6 0 30 24 19 13 7 1 31 25 20 14 8 2 32 26 21 15 9 3 33 27 22 16 10 4 34 29 23 17 11 5 35 29];
%     channel_map_noise = [18 19 20 21 22 23 12 13 14 15 16 17 6 7 8 9 10 11 0 1 2 3 4 5 30 31 32 33 34 35 24 25 26 27 28 29];
% %     channel_map_noise = 0:31;
% 
%     for loopShanks = 1:6
%         allSpikeLat{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/spike_times.npy']));
%         spikeTemplates{loopShanks}   = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/spike_templates.npy']));
%         templateChannels{loopShanks} = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/pc_feature_ind.npy']));
%         channel_map{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/channel_map.npy']));
% 
%         % adjust channel_map
%         channel_map_adjusted{1} = channel_map{1};
%         channel_map_noise_adjusted{1} = channel_map_noise;
%         if loopShanks > 1
%             adjust(loopShanks-1) = length(channel_map{loopShanks-1}); 
%             channel_map_adjusted{loopShanks} = channel_map{loopShanks} + sum(adjust(1:loopShanks-1));
%             channel_map_noise_adjusted{loopShanks} = channel_map_noise + sum(adjust(1:loopShanks-1));
%         end
% 
%         cluster_info     = readtable(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t');
%         noiseClusterIdx  = find(strcmp(cluster_info.group,'noise'));
% 
%         templateList     = unique(spikeTemplates{loopShanks});
%         spikeNoTemplates = histc(spikeTemplates{loopShanks},templateList);
% 
%         noiseClusterType = string(cluster_info.group(noiseClusterIdx));
%         noiseClusterChan = channel_map_noise_adjusted{loopShanks}(cluster_info.ch(noiseClusterIdx) + 1) + 1;
%         noiseClusterChan = cluster_info.ch(noiseClusterIdx) + 1;
%         noiseClusterShank = ones(length(noiseClusterChan),1)*loopShanks;
% 
%         noiseClusterID   = [];
%         lostNoiseUnitsIdx = [];
%         for loopNoiseCluster = 1:length(noiseClusterIdx)
%             tempIdx = find(cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster)) == spikeNoTemplates);
% 
%             % some noise units not recoverable
%             if isempty(tempIdx) || ~isempty(intersect(neuronNoSpikes,cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster))))
% 
%                 lostNoiseUnitsIdx = [lostNoiseUnitsIdx loopNoiseCluster];
% 
%                 continue
%             end
% 
%             tempIdx = tempIdx(1);
% 
%             noiseClusterID                         = [noiseClusterID str2double(['999' num2str(templateList(tempIdx))])]; % code noise units with 999 prefix
%             noiseClusterSpikeLat = [noiseClusterSpikeLat allSpikeLat{loopShanks}(    templateList(tempIdx) == spikeTemplates{loopShanks})];
%         end
% 
%         noiseClusterIdx(lostNoiseUnitsIdx)   = [];
%         noiseClusterType(lostNoiseUnitsIdx)  = [];
%         noiseClusterChan(lostNoiseUnitsIdx)  = [];
%         noiseClusterShank(lostNoiseUnitsIdx) = [];
% 
%         zerosIdx      = templateChannels{loopShanks} == 0;
%         zerosIdx(:,1) = 0;
%         templateChannels{loopShanks}(zerosIdx) = NaN;
% 
%         cellIDs          =  [cellIDs               noiseClusterID];
%         cellTypeList     =  [string(cellTypeList)  noiseClusterType'];
%         neuronChannel    =  [neuronChannel         noiseClusterChan'];
%         neuronShank      =  [neuronShank           noiseClusterShank'];
% 
%     end
% 
% elseif loopAnimal == 5
% 
%     UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
%     load([UtkuPath 'AG_2019-12-27_NSD' '/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);
% 
%     cellIDs          = cell_metrics.UID;
%     cellTypeList     = cell_metrics.putativeCellType;
%     neuronChannel    = cell_metrics.maxWaveformCh + 1;
%     neuronShank      = cell_metrics.shankID;
%     neuronSpikeTimes = cell_metrics.spikes.times;
%     for loopNeurons = 1:size(neuronSpikeTimes,2)
%         neuronNoSpikes(loopNeurons)    = length(neuronSpikeTimes{loopNeurons});
%     end
%     noiseClusterSpikeLat = {};
% 
%     channel_map_noise = [18 19 20 21 22 23 12 13 14 15 16 17 6 7 8 9 10 11 0 1 2 3 4 5 30 31 32 33 34 35 24 25 26 27 28 29];
% 
%     for loopShanks = 1:5
%         allSpikeLat{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/spike_times.npy']));
%         spikeTemplates{loopShanks}   = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/spike_templates.npy']));
%         templateChannels{loopShanks} = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/pc_feature_ind.npy']));
%         channel_map{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/channel_map.npy']));
% 
%         %adjust channel_map
%         channel_map_adjusted{1} = channel_map{1};
%         channel_map_noise_adjusted{1} = channel_map_noise;
%         if loopShanks > 1
%             adjust(loopShanks-1) = length(channel_map{loopShanks-1}); 
%             channel_map_adjusted{loopShanks} = channel_map{loopShanks} + sum(adjust(1:loopShanks-1));
%             channel_map_noise_adjusted{loopShanks} = channel_map_noise + sum(adjust(1:loopShanks-1));
%         end
% 
%         cluster_info     = readtable(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/cluster_info.tsv'], "FileType","text",'Delimiter', '\t');
%         noiseClusterIdx  = find(strcmp(cluster_info.group,'noise'));
% 
%         templateList     = unique(spikeTemplates{loopShanks});
%         spikeNoTemplates = histc(spikeTemplates{loopShanks},templateList);
% 
%         noiseClusterType = string(cluster_info.group(noiseClusterIdx));
%         noiseClusterChan = channel_map_noise_adjusted{loopShanks}(cluster_info.ch(noiseClusterIdx) + 1) + 1;
%         noiseClusterShank = ones(length(noiseClusterChan),1)*loopShanks;
% 
%         noiseClusterID   = [];
%         lostNoiseUnitsIdx = [];
%         for loopNoiseCluster = 1:length(noiseClusterIdx)
%             tempIdx = find(cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster)) == spikeNoTemplates);
% 
%             % some noise units not recoverable
%             if isempty(tempIdx) || ~isempty(intersect(neuronNoSpikes,cluster_info.n_spikes(noiseClusterIdx(loopNoiseCluster))))
% 
%                 lostNoiseUnitsIdx = [lostNoiseUnitsIdx loopNoiseCluster];
% 
%                 continue
%             end
% 
%             tempIdx = tempIdx(1);
% 
%             noiseClusterID       = [noiseClusterID str2double(['999' num2str(templateList(tempIdx))])]; % code noise units with 999 prefix
%             noiseClusterSpikeLat = [noiseClusterSpikeLat allSpikeLat{loopShanks}(    templateList(tempIdx) == spikeTemplates{loopShanks})];
%         end
% 
%         noiseClusterIdx(lostNoiseUnitsIdx)   = [];
%         noiseClusterType(lostNoiseUnitsIdx)  = [];
%         noiseClusterChan(lostNoiseUnitsIdx)  = [];
%         noiseClusterShank(lostNoiseUnitsIdx) = [];
% 
%         zerosIdx      = templateChannels{loopShanks} == 0;
%         zerosIdx(:,1) = 0;
%         templateChannels{loopShanks}(zerosIdx) = NaN;
% 
%         cellIDs          =  [cellIDs               noiseClusterID];
%         cellTypeList     =  [string(cellTypeList)  noiseClusterType'];
%         neuronChannel    =  [neuronChannel         noiseClusterChan];
%         neuronShank      =  [neuronShank           noiseClusterShank'];
%     end
% 
% elseif loopAnimal == 6
% 
%     % HiroCurateNoiseClusters 
% 
%     spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
%     dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
% 
%     [data_dir, name, ~] = fileparts(spike_data_fullpath);
% 
%     load(spike_data_fullpath, 'spikes')
%     load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
%     load(fullfile(data_dir, 'wake-basics.mat'),   'basics');
% 
%     load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')
% 
%     cellIDs         = 1:size(spikes.RoyMaze1,2);
%     neuronTemp      = cell2num({spikes.RoyMaze1.id}');
%     cellTypeList = strings(size(spikes.RoyMaze1,2),1);
%     cellTypeList([spikes.RoyMaze1.quality] == 8) = "i";
%     cellTypeList(neuronTemp(:,2) == 1)           = "noise";
%     cellTypeList(strcmp(cellTypeList,""))        = "p"; % double check!!
% 
%     neuronShank      = neuronTemp(:,1);
%     neuronSpikeTimes = {spikes.RoyMaze1.time}';
% 
%     for loopNeurons = 1:size(neuronSpikeTimes,1)
% 
%         % noise 
%         if strcmp(cellTypeList(loopNeurons),'noise')
%             cellIDs(loopNeurons) = str2double(['999' num2str(loopNeurons)]);
%         end
% 
%         temp = neuronSpikeTimes{loopNeurons};
%         temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
%         neuronSpikeTimes{loopNeurons} = temp/1e6; % convert to seconds 
% 
%         neuronNoSpikes(loopNeurons)    = length(neuronSpikeTimes{loopNeurons});
% 
%         if ~strcmp(cellTypeList{loopNeurons},"noise")
%             wave = unitsGlobalSnapshots.waveforms{cellIDs(loopNeurons) == unitsGlobalSnapshots.neuronID};
%             [~,chTemp] = min(min(wave'));
% 
%             neuronChannel(loopNeurons) = chTemp; 
%         else
%             neuronChannel(loopNeurons) = NaN; 
%         end
% 
%     end
% 
% end
% 
% shankList = sort(unique(neuronShank));
% 
% for loopShanks = 1:length(shankList)
% 
%     % define PC sets for channels - for Roy PCAs 
%             PCsets = [1:3; ...
%                       4:6;...
%                       7:9;...
%                       10:12;...
%                       13:15;...
%                       16:18;...
%                       19:21;...
%                       22:24];
% 
%     clusterValNeuron         = table;
%     % clusterValNeuron.PCAvals = {};
% 
%     disp(['Shank loop:  ' num2str((loopShanks/length(shankList))*100) '% done'])
% 
%     skippedLoops = [];
% 
%     for loopNeurons = 1:length(cellIDs)
% 
%         disp(['Neuron loop:  ' num2str((loopNeurons/length(cellIDs))*100) '% done'])
% 
%         cellID        = cellIDs(loopNeurons);
%         channel       = neuronChannel(loopNeurons);
%         shank         = neuronShank(loopNeurons);
%         cellType      = cellTypeList{loopNeurons};
% 
%         if shank ~= shankList(loopShanks)
% 
%             skippedLoops = [skippedLoops loopNeurons];
% 
%             continue
%         end
% 
%         if loopNeurons <= length(neuronSpikeTimes)
%             spikeTimes    = neuronSpikeTimes{loopNeurons};
%             spikeTimesLat = round(spikeTimes*fs)';
%         else % else is for animal loops 1-5 which includes N, U, S, and AG
%             spikeTimesLat = noiseClusterSpikeLat{loopNeurons - length(neuronSpikeTimes)};
%         end
% 
%         if sum(loopAnimal == [1,4,5]) == 1
% 
%             [~,idx,     ~] = intersect(allSpikeLat{shank},spikeTimesLat);
% 
%             if ~strcmp(cellType,'noise')
%                 templateList   = unique(spikeTemplates{shank}(idx));
% 
%                 spikeTemplatesTemp = spikeTemplates{shank}(idx);
% 
%                 chanIdx = mode(spikeTemplatesTemp);
% 
%                 % max chan sanity check
%                 clusIdx = find( (find((channel - 1) == channel_map_adjusted{shank}) - 1) == (templateChannels{shank}(chanIdx + 1 ,:)) );
% 
%                 if isempty(clusIdx)
%                    clusIdx = NaN;
%                 end
%             else 
%                 spikeTemplatesTemp = {};
%             end
% 
%             clusterLoopMax = size(templateChannels{loopShanks},2);
% 
%         elseif sum(loopAnimal == [2,3]) == 1
% 
%             [~,idx,     ~] = intersect(allSpikeLat,spikeTimesLat);
% 
%             if isempty(idx)
%                 continue
%             end
% 
%             if ~strcmp(cellType,'noise')
%                 templateList   = unique(spikeTemplates(idx));
% 
%                 spikeTemplatesTemp = spikeTemplates(idx);
% 
%                 chanIdx = mode(spikeTemplatesTemp);
% 
%                 % max chan sanity check
%                 clusIdx = find( (find((channel - 1) == channel_map) - 1) == (templateChannels(chanIdx + 1 ,:)) );
% 
%                 if isempty(clusIdx)
%                    clusIdx = NaN;
%                 end
%             else 
%                 spikeTemplatesTemp = {};
%             end
% 
%             clusterLoopMax = size(templateChannels,2);
% 
%         elseif loopAnimal == 6
% 
%             % weird timeshift in Roy's data
%             adjustFactor = 514;
% 
%             fs = 30000;
% 
%             % load session markers
%             bx = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-behavior.mat');
% 
%             % adjust unit 
%             load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy-maze1/Roy-maze1.clusterQuality.mat')
%             for i = 1:size({clusterQuality.clus},2)
%                 clusterCount(i) = size(clusterQuality(i).clus,1);
%             end
% 
%             clusterCounts = cumsum(clusterCount);
%             clusterCounts = [0 clusterCounts];
% 
%         %     unitIdx = unit-clusterCounts(shank)-2*(shank);
% 
%             % adjust channel index
%             % chIdx = tarCh-(tarShank-1)*8;
%             % 
%             % load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy-maze1/Roy-maze1.clusterQuality.mat')
% 
%             % load specific shank
%             fetTable = LoadRoyFet(['Roy-maze1.fet.' num2str(loopShanks)]);
%             cluTable = table(load(['Roy-maze1.clu.' num2str(loopShanks)])); 
%             cluTable = cluTable{:,:};
%             cluTable = cluTable(2:end);
%             spikeTimes = table(load(['Roy-maze1.res.' num2str(loopShanks)]));
%             spikeTimes = spikeTimes{:,:};
% 
%             for loopPCsets = 1:size(PCsets,1)
% 
%                 [PC1,PC2,PC3] = ...
%                 filterPCsforChannel(fetTable, ...
%                                     cluTable, ...
%                                     spikeTimes, ...
%                                     bx, ...
%                                     PCsets, ...
%                                     loopPCsets, ...
%                                     length(spikes.RoyMaze1(loopNeurons).time), ...
%                                     adjustFactor); 
% 
%                 PC1all{loopPCsets} = PC1;
%                 PC2all{loopPCsets} = PC2;
%                 PC3all{loopPCsets} = PC3;
%             end
% 
%             tempCounts = [];
%             for loopPtemp = 1:size(PC1,2)
%                  tempCounts1(loopPtemp) = size(PC1{loopPtemp},1) - 3;
%                  tempCounts2(loopPtemp) = size(PC1{loopPtemp},1) - 2;
%                  tempCounts3(loopPtemp) = size(PC1{loopPtemp},1) - 1;
%                  tempCounts4(loopPtemp) = size(PC1{loopPtemp},1) + 0;
%                  tempCounts5(loopPtemp) = size(PC1{loopPtemp},1) + 1;
%                  tempCounts6(loopPtemp) = size(PC1{loopPtemp},1) + 2;
%                  tempCounts7(loopPtemp) = size(PC1{loopPtemp},1) + 3;
%             end
% 
%             tempIdx1 = neuronNoSpikes(loopNeurons) == tempCounts1;
%             tempIdx2 = neuronNoSpikes(loopNeurons) == tempCounts2;
%             tempIdx3 = neuronNoSpikes(loopNeurons) == tempCounts3;
%             tempIdx4 = neuronNoSpikes(loopNeurons) == tempCounts4;
%             tempIdx5 = neuronNoSpikes(loopNeurons) == tempCounts5;
%             tempIdx6 = neuronNoSpikes(loopNeurons) == tempCounts6;
%             tempIdx7 = neuronNoSpikes(loopNeurons) == tempCounts7;
% 
%             tempIdx = (((tempIdx1 | tempIdx2) | (tempIdx3 | tempIdx4))  | (tempIdx5 | tempIdx6)) | tempIdx7;
% 
%             clusterValAllChans = [];
%             for loopPCsets = 1:size(PCsets,1)
%                  clusterValAllChans{loopPCsets} = [PC1all{loopPCsets}{tempIdx} ...
%                                                    PC2all{loopPCsets}{tempIdx} ...
%                                                    PC3all{loopPCsets}{tempIdx}];
%             end
% 
%             spikeTemplatesTemp = {};
%         end
% 
%         if loopAnimal ~= 6
% 
%             clusterValAllChans = [];
% 
%             for loopClusterChannels = 1:clusterLoopMax %length(~isnan(templateChannels(chanIdx + 1 ,:)))
% 
%                 if loopAnimal == 1
%                     clusterVal = readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(shank) '/pc_features_PCchan' num2str(loopClusterChannels - 1) '.npy']);
%                 elseif loopAnimal == 2
%                     clusterVal = readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/pc_features_PCchan' num2str(loopClusterChannels - 1) '.npy']);
%                 elseif loopAnimal == 3
%                     clusterVal = readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/pc_features_PCchan' num2str(loopClusterChannels - 1) '.npy']);
%                 elseif loopAnimal == 4
%                     clusterVal = readNPY(['/media/nasko/WD_BLACK/UtkuRawMultiArrayUnpackedNPY/AG_2019-12-23_NSD/Units/s' num2str(shank) '/pc_features_PCchan' num2str(loopClusterChannels - 1) '.npy']);
%                 elseif loopAnimal == 5
%                     clusterVal = readNPY(['/media/nasko/WD_BLACK/UtkuRawMultiArrayUnpackedNPY/AG_2019-12-27_NSD/Units/s' num2str(shank) '/pc_features_PCchan' num2str(loopClusterChannels - 1) '.npy']);
%                 end
% 
%                 clusterValAllChans{loopClusterChannels} = clusterVal(idx,:);
% 
%             end
%         end
% 
%         clusterValNeuron.spikeTemplates{loopNeurons}   = spikeTemplatesTemp;
%         clusterValNeuron.clusIdx(loopNeurons)          = clusIdx;
%         clusterValNeuron.cellID(loopNeurons)           = cellID;
%         % clusterValNeuron.PCAvals{loopNeurons}          = clusterValAllChans;
%         clusterValNeuron.channel(loopNeurons)          = channel;
%         clusterValNeuron.shank(loopNeurons)            = shank;
% %         clusterValNeuron.noiseSpikeTimes(loopNeurons)  = noiseSpikeTimes;
%         if sum(loopAnimal == [1,4,5]) == 1
%             clusterValNeuron.templateChannels{loopNeurons} = templateChannels{shank}(chanIdx + 1 ,:);
%         elseif sum(loopAnimal == [2,3]) == 1
%             clusterValNeuron.templateChannels{loopNeurons} = templateChannels(chanIdx + 1 ,:);
%         elseif loopAnimal == 6
%             clusterValNeuron.templateChannels{loopNeurons} = 1:8;
%         end
% 
%         save(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_rat' animalNames{loopAnimal} '_shank' num2str(loopShanks) '_loopNeuron' num2str(loopNeurons)],'clusterValAllChans')
% 
%     end
% 
%     clusterValNeuron(skippedLoops(skippedLoops <= size(clusterValNeuron,1)),:) = [];
% 
%     save(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_rat' animalNames{loopAnimal} '_shank' num2str(loopShanks)],'clusterValNeuron')
% 
% end

%% prep ClusterAna

% for loopRat = 1:6
% 
%     if loopRat == 1
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStats.mat') % 8 shanks
%         pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
%     elseif loopRat == 2
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat') % 14 shanks
%         pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
%     elseif loopRat == 3
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStats.mat') % 12 shanks
%         pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
%     elseif loopRat == 4
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat') % 6 shanks
%         pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
%     elseif loopRat == 5
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStats.mat') % 5 shanks
%         pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
%     elseif loopRat == 6
%         load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStats.mat')
%         pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
%     end
% 
%     fs             = 30000;
%     binSize        = 1/fs;
% 
%     duration       = 0.007;
%     lagLimit       = duration/2;
% 
%     for loopPairs = 1:size(pairGroupStatTable,1)
% 
%         R      = pairGroupStatTable.refSpikeTimes{loopPairs,1};
%         T      = pairGroupStatTable.tarSpikeTimes{loopPairs,1};
% 
%         refLat = round(R*fs)';
%         tarLat = round(T*fs)';
% 
%         tarChannel = pairGroupStatTable.tarChannel(loopPairs,1);
% 
%         GSPExc = pairGroupStatTable.GSPExc{loopPairs,1};
% 
%         [CCG,spikeTimesInBin] = CCGtrackedSpikeTimes(R,T,binSize,lagLimit);
% 
%         Rsync = [];
%         Tsync = [];
%         for k = 1:211
%             if (GSPExc(k) == 1) 
%                 if isempty(spikeTimesInBin{k})
%                     continue
%                 else
%                     Rsync = [Rsync; spikeTimesInBin{k}(:,1)];
%                     Tsync = [Tsync; spikeTimesInBin{k}(:,2)];
%                 end
%             end
%         end
% 
%         pairGroupStatTable.refSynch{loopPairs} = Rsync;
%         pairGroupStatTable.tarSynch{loopPairs} = Tsync;
% 
%         pairGroupStatTable.refSynchLat{loopPairs} = round(Rsync/(1/3e4))';
%         pairGroupStatTable.tarSynchLat{loopPairs} = round(Tsync/(1/3e4))';
% 
%         if loopRat == 1
%             save('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat','pairGroupStatTable') % 8 shanks
%         elseif loopRat == 2
%             save('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStatsClusterAna.mat','pairGroupStatTable') % 14 shanks
%         elseif loopRat == 3
%             save('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAnamat','pairGroupStatTable') % 12 shanks
%         elseif loopRat == 4
%             save('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStatsClusterAna.mat','pairGroupStatTable') % 6 shanks
%         elseif loopRat == 5
%             save('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStatsClusterAna.mat','pairGroupStatTable') % 5 shanks
%         elseif loopRat == 6
%             save(('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat'))
%         end
% 
%     end
% end


%% colate data for Lratio and iso distance calc

% % vars for cellExplorer scripts for Lratio and iso distance
% cellIDall                = [];
% isolation_distance       = [];
% L_ratio                  = [];
% ratLoopNo                = [];
% cellTypeAll              = {};
% cellTypeExq              = {};
% 
% isolation_distance_exq   = [];
% L_ratio_exq              = [];
% 
% isolation_distance_synch = [];
% L_ratio_synch            = [];
% ratLoopNoSynch           = [];
% 
% Nrandspk           = 1000;
% minNochans         = [8,8,8,8,8,8]; 
% 
% fs                 = 30000;
% 
% noShanks           = [8,14,12,6,5,8];
% 
% for loopRat = 6:6
% 
%     if loopRat == 1
%         neuronsN    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
%     elseif loopRat == 2
%         neuronsS    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
%     elseif loopRat == 3
%         neuronsU    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
%     elseif loopRat == 4
%         UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
%         load([UtkuPath 'AG_2019-12-23_NSD' '/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);
%     elseif loopRat == 5
%         UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
%         load([UtkuPath 'AG_2019-12-27_NSD' '/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);
%     elseif loopRat == 6
%         spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
%         load(spike_data_fullpath, 'spikes')
% 
%         cellIDs         = 1:size(spikes.RoyMaze1,2);
%         neuronTemp      = cell2num({spikes.RoyMaze1.id}');
% 
%         cellTypeList = strings(size(spikes.RoyMaze1,2),1);
%         cellTypeList([spikes.RoyMaze1.quality] == 8) = "i";
%         cellTypeList(neuronTemp(:,2) == 1)           = "noise";
%         cellTypeList(strcmp(cellTypeList,""))        = "p"; % double check!!
%     end
% 
%     for loopShanks = 1:noShanks(loopRat)
% 
%         tic 
% 
%         if loopRat == 1
%             load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'   num2str(loopShanks) '.mat']) % 8 shanks
%             load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat')  
%         elseif loopRat == 2
%             load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'   num2str(loopShanks) '.mat']) % 14 shanks
%             load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStatsClusterAna.mat') 
%         elseif loopRat == 3
%             load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'   num2str(loopShanks) '.mat']) % 12 shanks
%             load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAna.mat') 
%         elseif loopRat == 4
%             load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(loopShanks) '.mat']) % 6 shanks
%             load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStatsClusterAna.mat')
%         elseif loopRat == 5
%             load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(loopShanks) '.mat']) % 5 shanks
%             load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStatsClusterAna.mat')
%         elseif loopRat == 6
%             load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(loopShanks) '.mat']) % 8 shanks
%             load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
%         end
% 
%         %% synch spike times from exquisite CCG peaks
%         % NOTE: channel numbers here from rats N and S are from not probefile 
%         pairGroupStatTableRefSynch = pairGroupStatTable(pairGroupStatTable.refShank == loopShanks,:);
%         pairGroupStatTableTarSynch = pairGroupStatTable(pairGroupStatTable.tarShank == loopShanks,:);
% 
%         if ~isempty(pairGroupStatTableRefSynch)
%             for loopNeurons = 1:size(pairGroupStatTableRefSynch,1) 
%                 [~,idxRef,~] = intersect(pairGroupStatTableRefSynch.refSpikeTimes{loopNeurons}, ...
%                                          pairGroupStatTableRefSynch.refSynch{loopNeurons});
%                 idxRefSynch{loopNeurons} = idxRef;
%             end
%         else
%             idxRefSynch = [];
%         end
%         if ~isempty(pairGroupStatTableTarSynch)
%             for loopNeurons = 1:size(pairGroupStatTableTarSynch,1)
%                 [~,idxTar,~] = intersect(pairGroupStatTableTarSynch.tarSpikeTimes{loopNeurons}, ...
%                                          pairGroupStatTableTarSynch.tarSynch{loopNeurons});
%                 idxTarSynch{loopNeurons} = idxTar;
%             end
%         else
%             idxTarSynch = [];
%         end
% 
%         synchUnitIDs         = [pairGroupStatTableRefSynch.refNeuronID; pairGroupStatTableTarSynch.tarNeuronID];
%         synchUnitIdx4PCAvals = [idxRefSynch';                           idxTarSynch'                          ];
% 
%         for loopNeurons1 = 1:size(clusterValNeuron,1)
% 
%             % skip units lacking 4 channels for PCA ana
%             if sum(~isnan(clusterValNeuron.templateChannels{loopNeurons1})) < minNochans(loopRat)
%                 display(['Unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                 continue
%             end
% 
%             % skip noise as main 
%             if contains(string(clusterValNeuron.cellID(loopNeurons1)),'999');
%                 display(['Noise unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                 continue
%             end
% 
%             %% tag cell type
%             if loopRat == 1
%                 cellType = neuronsN.neuron_type((neuronsN.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1),1:5);
%             elseif loopRat == 2
%                 cellType = neuronsS.neuron_type((neuronsS.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1),1:5);
%             elseif loopRat == 3
%                 cellType = neuronsU.neuron_type((neuronsU.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1),1:5);
%             elseif (loopRat == 4) || (loopRat == 5)
%                 cellType = cell_metrics.putativeCellType{cell_metrics.cellID == clusterValNeuron.cellID(loopNeurons1)};
%             elseif loopRat == 6
%                 cellType = cellTypeList(cellIDs == clusterValNeuron.cellID(loopNeurons1));
%             end
% 
%             if loopRat == 1
%                 if strcmp(neuronsN.neuron_type((neuronsN.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1),1:5),'mua  ')
%                     display(['MUA unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                     continue
%                 end
%             elseif loopRat == 2
%                 if strcmp(neuronsS.neuron_type((neuronsS.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1),1:5),'mua  ')
%                     display(['MUA unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                     continue
%                 end
%             elseif loopRat == 3
%                 if strcmp(neuronsU.neuron_type((neuronsU.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1),1:5),'mua  ')
%                     display(['MUA unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                     continue
%                 end
%             elseif(loopRat == 4) || (loopRat == 5) % skip mua's in Utku rat
%                 if strcmp(cell_metrics.label{cell_metrics.cellID == clusterValNeuron.cellID(loopNeurons1)},'mua')
%                     display(['MUA unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                     continue
%                 end
%             end
% 
%             % remove NaN from templateChannels
%             templateChannels = clusterValNeuron.templateChannels{loopNeurons1};
%             templateChannels(isnan(templateChannels)) = [];
% 
%             % random chosen channels for ana
%             randIdx   = sort(randperm(length(templateChannels),minNochans(loopRat)));   
%             randChans = templateChannels(randIdx);
% 
%             idxAnaPartnerUnits = [];
%             for loopNeurons2 = 1:size(clusterValNeuron,1)
%                idxAnaPartnerUnits(loopNeurons2)  = length(intersect(randChans,clusterValNeuron.templateChannels{loopNeurons2})) >= minNochans(loopRat);
%             end
% 
%             clusterValNeuronAna = clusterValNeuron(find(idxAnaPartnerUnits),:);
% 
%             %% noise only as metric partners
% 
%             noiseIdxBinary = contains(string(clusterValNeuronAna.cellID),'999');
%             noiseIdx       = find(noiseIdxBinary);
% 
%             % any noise units left
%             if isempty(noiseIdx)
%                 display(['Unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                 continue
%             end
% 
%             unitOfInterestIdx = find(clusterValNeuronAna.cellID == clusterValNeuron.cellID(loopNeurons1));
% 
%             clusterValNeuronAna = clusterValNeuronAna([unitOfInterestIdx; noiseIdx],:);
% 
%             % limit ana partners to 10 minus 1
%             if loopRat ~= 6
%                 if size(clusterValNeuronAna,1) >= 10
%                     clusterValNeuronAna = clusterValNeuronAna(1:10,:);
%                 end
%             end
% 
%             % only one unit for analysis? If so skip. 
%             if size(clusterValNeuronAna,1) <= 1
%                 display(['Unit skipped: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
%                 continue
%             end
% 
%             features      = [];
%             cluster_index = [];
%             shank_list    = [];
% 
%             features_synch      = {};
%             cluster_index_synch = {};
%             shank_list_synch    = {};
% 
%             cellID = [];
% 
%             % find if neuron is part of exquisite pairs
%             if sum(clusterValNeuron.cellID(loopNeurons1) == synchUnitIDs) > 0
%                 synchUnitIDsIdx = find(clusterValNeuron.cellID(loopNeurons1) == synchUnitIDs);
%                 synchUnitIdx4PCAvalsTemp = synchUnitIdx4PCAvals(synchUnitIDsIdx);
%             else 
%                 synchUnitIDsIdx = [];
%             end
% 
%             %% creating analysis groups  
% 
%             commonChannels = randChans;
% 
%             % colating the PCA vals
%             for loopNeurons2 = 1:length(clusterValNeuronAna.cellID)
% 
%                 [~,~,idxCommon] = intersect(commonChannels,clusterValNeuronAna.templateChannels{loopNeurons2});
% 
%                 PCAvals = [];
%     %             for loopClusterChannels = 1:size(clusterValNeuron.templateChannels{loopNeurons},2)
%     %                 PCAvals    = [PCAvals clusterValNeuron.PCAvals{loopNeurons}{1,loopClusterChannels}];
%     %             end
%                 for loopClusterChannels = 1:length(idxCommon)
%                     % PCAvals    = [PCAvals clusterValNeuronAna.PCAvals{loopNeurons2}{1,idxCommon(loopClusterChannels)}];
% 
%                     if loopRat == 1
%                         loopIdx = find((neuronsN.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1));
%                         load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(loopShanks) '_loopNeuron' num2str(loopIdx) '.mat']);
%                     elseif loopRat == 2
%                         loopIdx = find((neuronsS.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1));
%                         load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(loopShanks) '_loopNeuron' num2str(loopIdx) '.mat']);
%                     elseif loopRat == 3
%                         loopIdx = find((neuronsU.neuron_ids + 1) == clusterValNeuron.cellID(loopNeurons1));
%                         load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(loopShanks) '_loopNeuron' num2str(loopIdx) '.mat']);
%                     elseif loopRat == 4
%                         loopIdx = find(cell_metrics.cellID == clusterValNeuron.cellID(loopNeurons1));
%                         load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(loopShanks) '_loopNeuron' num2str(loopIdx) '.mat']);
%                     elseif loopRat == 5
%                         loopIdx = find(cell_metrics.cellID == clusterValNeuron.cellID(loopNeurons1));
%                         load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(loopShanks) '_loopNeuron' num2str(loopIdx) '.mat']);
%                     elseif loopRat == 6
%                         loopIdx = find(cellIDs == clusterValNeuron.cellID(loopNeurons1));
%                         load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(loopShanks) '_loopNeuron' num2str(loopIdx) '.mat'])
%                     end
% 
%                     PCAvals    = [PCAvals clusterValAllChans{1,idxCommon(loopClusterChannels)}];
%                 end
% 
%                 PCAsynch = {};
%                 % synch spikes (make extra PCA variables)
%                 if (clusterValNeuronAna.cellID(loopNeurons2)             == clusterValNeuron.cellID(loopNeurons1)) && ...
%                    (size(clusterValAllChans{1},1) == size(PCAvals,1)                      )
% 
%                     SynchClusterIdx = cluster_index;
% 
%                     if length(synchUnitIDsIdx) == 1
%                         PCAsynch{1} = PCAvals(synchUnitIdx4PCAvalsTemp{1}(synchUnitIdx4PCAvalsTemp{1} < size(PCAvals,1)),:);
%                     elseif length(synchUnitIDsIdx) > 1
%                         for loopExq = 1:size(synchUnitIdx4PCAvalsTemp,1)
%                             PCAsynch{loopExq} = PCAvals(synchUnitIdx4PCAvalsTemp{loopExq}(synchUnitIdx4PCAvalsTemp{loopExq} < size(PCAvals,1)),:);
%                         end
%                     end
%                 end
% 
%                 if size(PCAvals,1) > Nrandspk
%                     idxRand = randperm(size(PCAvals,1),Nrandspk);
%                     PCAvals = PCAvals(idxRand,:);
%                 end
% 
%                 if ~isempty(PCAsynch)
%                     for loopExq = 1:length(PCAsynch)
% 
%                         if size(PCAsynch{loopExq},1) > Nrandspk
%                             idxRand = randperm(size(PCAsynch{loopExq},1),Nrandspk);
%                             PCAsynch{loopExq} = PCAsynch{loopExq}(idxRand,:);
%                         end
% 
%                         features_synch{loopExq}      = PCAsynch{loopExq};
%                         cluster_index_synch{loopExq} = loopNeurons2*ones(size(PCAsynch{loopExq},1),1);
%                         shank_list_synch{loopExq}    = clusterValNeuronAna.shank(loopNeurons2)*ones(size(PCAsynch{loopExq},1),1);
%                     end
%                 end
% 
%                 features      = [features;      PCAvals                                                        ];
%                 cluster_index = [cluster_index; loopNeurons2*ones(size(PCAvals,1),1)                           ];
%                 shank_list    = [shank_list;    clusterValNeuronAna.shank(loopNeurons2)*ones(size(PCAvals,1),1)];
% 
%                 if ~isempty(features_synch) & isempty(PCAsynch)
%                     for loopExq = 1:length(features_synch)
% 
%                         features_synch{loopExq}      = [features_synch{loopExq};      PCAvals                                                        ];
%                         cluster_index_synch{loopExq} = [cluster_index_synch{loopExq}; loopNeurons2*ones(size(PCAvals,1),1)                           ];
%                         shank_list_synch{loopExq}    = [shank_list_synch{loopExq};    clusterValNeuronAna.shank(loopNeurons2)*ones(size(PCAvals,1),1)];
% 
%                     end
%                 end
% 
%                 cellID    = [cellID;    clusterValNeuronAna.cellID(loopNeurons2)];
%             end
% 
%             tempIdx = find(clusterValNeuron.cellID(loopNeurons1) == cellID);
% 
%             cellIDall = [cellIDall; cellID(tempIdx)];
% 
%             tempIsolation_distance = IsolationDistance_calc(features,cluster_index);
%             tempL_ratio            = L_ratio_calc(          features,cluster_index);
% 
%             isolation_distance     = [isolation_distance; tempIsolation_distance(tempIdx)];
%             L_ratio                = [L_ratio;            tempL_ratio(tempIdx)];
% 
%             ratLoopNo              = [ratLoopNo loopRat*ones(length(tempIdx),1)];
%             cellTypeAll            = [cellTypeAll cellType]; 
% 
%             for loopExq = 1:length(features_synch)
% 
%                 % any spike time
%                 isolation_distance_exq   = [isolation_distance_exq; tempIsolation_distance(tempIdx)];
%                 L_ratio_exq              = [L_ratio_exq;            tempL_ratio(tempIdx)           ];
% 
%                 % synch spike times
%                 tempIsolation_distance   = IsolationDistance_calc(features_synch{loopExq},cluster_index_synch{loopExq});
%                 tempL_ratio              = L_ratio_calc(          features_synch{loopExq},cluster_index_synch{loopExq});
% 
%                 isolation_distance_synch = [isolation_distance_synch; tempIsolation_distance(tempIdx)];
%                 L_ratio_synch            = [L_ratio_synch;            tempL_ratio(tempIdx)];
% 
%                 ratLoopNoSynch           = [ratLoopNoSynch loopRat*ones(length(tempIdx),1)];
%                 cellTypeExq              = [cellTypeExq cellType]; 
% 
%             end
% 
%             display(['Unit: ' num2str((loopNeurons1/size(clusterValNeuron,1))*100) '% done'])
% 
%         end
% 
%         toc
% 
%         display(['Rat loop ' num2str(loopRat) ': '  num2str((loopShanks/noShanks(loopRat))*100) '% done'])
% 
%     end 
% 
% end


%%

% tiledlayout(1,2)
% nexttile
% histogram(isolation_distance,'Normalization','probability','BinWidth',1)
% xlim([0 100])
% nexttile
% histogram(L_ratio,           'Normalization','probability','BinWidth',0.1)
% 
% % noiseIdxBinary = contains(string(cellIDall),'999');
% % 
% tiledlayout(2,2)
% 
% % idxAll = strcmp(cellTypeAll,'inter');
% % idxExq = strcmp(cellTypeExq,'inter');
% 
% idxAll = strcmp(cellTypeAll,'Narrow Interneuron');
% idxExq = strcmp(cellTypeExq,'Narrow Interneuron');
% 
% nexttile
% histogram(isolation_distance(idxAll),    'Normalization','probability','BinWidth',1)
% hold on 
% histogram(isolation_distance_exq(idxExq),'Normalization','probability','BinWidth',1)
% hold off
% title('interneuron isolation distance')
% legend('all','exquisite')
% 
% nexttile
% histogram(L_ratio(idxAll),    'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(L_ratio_exq(idxExq),'Normalization','probability','BinWidth',0.1)
% hold off
% title('interneuron L ratio')
% legend('all','exquisite')
% 
% nexttile
% histogram(isolation_distance_exq(idxExq),  'Normalization','probability','BinWidth',1,'FaceColor',"#D95319")
% hold on 
% histogram(isolation_distance_synch,'Normalization','probability','BinWidth',1,'FaceColor',"#7E2F8E")
% hold off
% title('interneuron isolation distance')
% legend('exquisite','exquisite synch')
% 
% nexttile
% histogram(L_ratio_exq(idxExq),  'Normalization','probability','BinWidth',0.1,'FaceColor',"#D95319")
% hold on 
% histogram(L_ratio_synch,'Normalization','probability','BinWidth',0.1,'FaceColor',"#7E2F8E")
% hold off
% title('interneuron L ratio')
% legend('exquisite','exquisite synch')

%%

% tiledlayout(2,2)
% 
% % idxAll = strcmp(cellTypeAll,'pyr  ');
% % idxExq = strcmp(cellTypeExq,'pyr  ');
% 
% idxAll = strcmp(cellTypeAll,'Pyramidal Cell') | strcmp(cellTypeAll,'Wide Interneuron');
% idxExq = strcmp(cellTypeExq,'Pyramidal Cell') | strcmp(cellTypeExq,'Wide Interneuron');
% 
% nexttile
% histogram(isolation_distance(idxAll),    'Normalization','probability','BinWidth',1)
% hold on 
% histogram(isolation_distance_exq(idxExq),'Normalization','probability','BinWidth',1)
% hold off
% title('pyramids isolation distance')
% legend('all','exquisite')
% 
% nexttile
% histogram(L_ratio(idxAll),    'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(L_ratio_exq(idxExq),'Normalization','probability','BinWidth',0.1)
% hold off
% title('pyramids L ratio')
% legend('all','exquisite')
% 
% nexttile
% histogram(isolation_distance_exq(idxExq),  'Normalization','probability','BinWidth',1,'FaceColor',"#D95319")
% hold on 
% histogram(isolation_distance_synch,'Normalization','probability','BinWidth',1,'FaceColor',"#7E2F8E")
% hold off
% title('pyramids isolation distance')
% legend('exquisite','exquisite synch')
% 
% nexttile
% histogram(L_ratio_exq(idxExq),  'Normalization','probability','BinWidth',0.1,'FaceColor',"#D95319")
% hold on 
% histogram(L_ratio_synch,'Normalization','probability','BinWidth',0.1,'FaceColor',"#7E2F8E")
% hold off
% title('pyramids L ratio')
% legend('exquisite','exquisite synch')


%%

% nexttile
% histogram(isolation_distance(ratLoopNo == 1),'Normalization','probability','BinWidth',1)
% hold on 
% histogram(isolation_distance_synch(ratLoopNoSynch == 1),'Normalization','probability','BinWidth',1)
% % histogram(isolation_distance(noiseIdxBinary),'Normalization','probability','BinWidth',1)
% hold off
% title('isolation distance')
% 
% nexttile
% histogram(L_ratio(ratLoopNo == 1),'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(L_ratio_synch(ratLoopNoSynch == 1),'Normalization','probability','BinWidth',0.1)
% % histogram(L_ratio(noiseIdxBinary),'Normalization','probability','BinWidth',1)
% hold off
% title('L ratio')
% 
% nexttile
% histogram(isolation_distance(ratLoopNo == 2),'Normalization','probability','BinWidth',1)
% hold on 
% histogram(isolation_distance_synch(ratLoopNoSynch == 2),'Normalization','probability','BinWidth',1)
% % histogram(isolation_distance(noiseIdxBinary),'Normalization','probability','BinWidth',1)
% hold off
% title('isolation distance')
% 
% nexttile
% histogram(L_ratio(ratLoopNo == 2),'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(L_ratio_synch(ratLoopNoSynch == 2),'Normalization','probability','BinWidth',0.1)
% % histogram(L_ratio(noiseIdxBinary),'Normalization','probability','BinWidth',1)
% hold off
% title('L ratio')
% 
% nexttile
% histogram(isolation_distance(ratLoopNo == 3),'Normalization','probability','BinWidth',1)
% hold on 
% histogram(isolation_distance_synch(ratLoopNoSynch == 3),'Normalization','probability','BinWidth',1)
% % histogram(isolation_distance(noiseIdxBinary),'Normalization','probability','BinWidth',1)
% hold off
% title('isolation distance')
% 
% nexttile
% histogram(L_ratio(ratLoopNo == 3),'Normalization','probability','BinWidth',0.1)
% hold on 
% histogram(L_ratio_synch(ratLoopNoSynch == 3),'Normalization','probability','BinWidth',0.1)
% % histogram(L_ratio(noiseIdxBinary),'Normalization','probability','BinWidth',1)
% hold off
% title('L ratio')

%% within cluster analysis - Han

filterFlag     = false; 
noDemeanFlag   = false;
fpass          = 300;
fs             = 3e4;
preLength      = 5*30;
postLength     = 5*30;

figPath = '/media/nasko/WD_BLACK31/ClusterAnaTemp/figures/';

resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
fig_use        = 102;

animalNames = {'N','S','U','AG1','AG2','Roy'};

NspikesPlot = 1000;

for loopRat = 1:6

    if loopRat == 1
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat')  
        neuronsN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/unitsGlobalSnapshots.mat')

        channelSet = {1:16, ...
                      17:32, ... 
                      33:48, ...
                      49:64, ...
                      65:80, ...
                      81:96, ...
                      97:112, ...
                      113:128};

    elseif loopRat == 2
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStatsClusterAna.mat') 
        neuronsS = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/unitsGlobalSnapshots.mat')

        channelSet = {1:10, ...
                      11:21, ...
                      22:32, ... 
                      33:43, ...
                      44:54, ...
                      55:63, ...
                      64:79, ...
                      80:95, ...
                      96:111, ...
                      112:127, ...
                      128:143, ...
                      144:159, ...
                      160:175, ...
                      176:191};

    elseif loopRat == 3
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAna.mat') 
        neuronsU = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
    elseif loopRat == 4
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStatsClusterAna.mat')
        UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
        cell_metrics = load([UtkuPath 'AG_2019-12-23_NSD' '/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);
    elseif loopRat == 5
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStatsClusterAna.mat')
        UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
        cell_metrics = load([UtkuPath 'AG_2019-12-27_NSD' '/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);
    elseif loopRat == 6
        load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')

        spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
        load(spike_data_fullpath, 'spikes')

        bx = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-behavior.mat');

        cellIDs         = 1:size(spikes.RoyMaze1,2);
        neuronTemp      = cell2num({spikes.RoyMaze1.id}');

        cellTypeList = strings(size(spikes.RoyMaze1,2),1);
        cellTypeList([spikes.RoyMaze1.quality] == 8) = "i";
        cellTypeList(neuronTemp(:,2) == 1)           = "noise";
        cellTypeList(strcmp(cellTypeList,""))        = "p"; % double check!!
    end

    for loopPairs = 1:size(pairGroupStatTable,1)

        refShank           = pairGroupStatTable.refShank(loopPairs);
        refChan            = pairGroupStatTable.refChannel(loopPairs);
        refID              = pairGroupStatTable.refNeuronID(loopPairs);
        refSpikeTimes      = pairGroupStatTable.refSpikeTimes{loopPairs};
        refSpikeTimesSynch = pairGroupStatTable.refSynch{loopPairs};
        refCellType        = pairGroupStatTable.refCellExplorerType{loopPairs};
        [~,refSpikeTimesSynchIdx,~] = intersect(refSpikeTimes,refSpikeTimesSynch);

        tarShank           = pairGroupStatTable.tarShank(loopPairs);
        tarChan            = pairGroupStatTable.tarChannel(loopPairs);
        tarID              = pairGroupStatTable.tarNeuronID(loopPairs);
        tarSpikeTimes      = pairGroupStatTable.tarSpikeTimes{loopPairs};
        tarSpikeTimesSynch = pairGroupStatTable.tarSynch{loopPairs};
        tarCellType        = pairGroupStatTable.tarCellExplorerType{loopPairs};
        [~,tarSpikeTimesSynchIdx,~] = intersect(tarSpikeTimes,tarSpikeTimesSynch);

        d                  = pairGroupStatTable.pairDistance(loopPairs);

        tWave = pairGroupStatTable.tWave{loopPairs};

        %% 

        titleStrRef = ['rat ' animalNames{loopRat} ': '...
                        num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')' ' - ' ...
                        num2str(refID) refCellType ' (sh: ' num2str(refShank) ')' ' pair ' num2str(d) '\mum '];

        titleStrTar = ['rat ' animalNames{loopRat} ': '...
                        num2str(refID) refCellType ' (sh: ' num2str(refShank) ')' ' - ' ...
                        num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')' ' pair ' num2str(d) '\mum '];

        saveStr  = ['rat ' animalNames{loopRat} ': '...
                    num2str(refID) refCellType ' (sh ' num2str(refShank) ')' ' - ' ...
                    num2str(tarID) tarCellType ' (sh ' num2str(tarShank) ')' ' pair ' num2str(d) ' microns.jpeg'];

        % [~,tarChan] = max(max(abs(unitsGlobalSnapshots.waveforms{tarID,1}(channelSet{tarShank},:)),[],2));

        if loopRat == 1
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'    num2str(refShank) '.mat']) 
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratN_shank'    num2str(tarShank) '.mat']) 
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find((neuronsN.neuron_ids + 1) == tarID);
            loopIdxRef = find((neuronsN.neuron_ids + 1) == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratN_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 2
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'    num2str(refShank) '.mat']);
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratS_shank'    num2str(tarShank) '.mat']);
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find((neuronsS.neuron_ids + 1) == tarID);
            loopIdxRef = find((neuronsS.neuron_ids + 1) == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratS_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 3
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'    num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratU_shank'    num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find((neuronsU.neuron_ids + 1) == tarID);
            loopIdxRef = find((neuronsU.neuron_ids + 1) == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratU_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 4
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG1_shank' num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find(cell_metrics.cell_metrics.UID == tarID);
            loopIdxRef = find(cell_metrics.cell_metrics.UID == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG1_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 5
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratAG2_shank' num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find(cell_metrics.cell_metrics.UID == tarID);
            loopIdxRef = find(cell_metrics.cell_metrics.UID == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratAG2_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        elseif loopRat == 6
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(refShank) '.mat'])
            clusterValNeuronsRef = clusterValNeuron;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeurons_ratRoy_shank' num2str(tarShank) '.mat'])
            clusterValNeuronsTar = clusterValNeuron;

            loopIdxTar = find(cellIDs == tarID);
            loopIdxRef = find(cellIDs == refID);

            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(refShank) '_loopNeuron' num2str(loopIdxRef) '.mat'])
            clusterValNeuronsPCAvalsRef = clusterValAllChans;
            load(['/media/nasko/WD_BLACK31/ClusterAnaTemp/clusterValNeuronsPvals_ratRoy_shank' num2str(tarShank) '_loopNeuron' num2str(loopIdxTar) '.mat'])
            clusterValNeuronsPCAvalsTar = clusterValAllChans;

            clear clusterValNeuron
            clear clusterValAllChans
        end

        %% gather PC values 

        % PCAchans = clusterValNeuron.templateChannels{clusterValNeuron.cellID == tarID};
        % PCAchansIdx = ~isnan(PCAchans);

        if loopRat == 6
            cluIdxRef  = 1;
            cluIdxTar  = 1;
        else
            cluIdxRef  = clusterValNeuronsRef.clusIdx(find(clusterValNeuronsRef.cellID == refID));
            cluIdxTar  = clusterValNeuronsTar.clusIdx(find(clusterValNeuronsTar.cellID == tarID));
        end

        PCAvalsRef = clusterValNeuronsPCAvalsRef;
        PCAvalsTar = clusterValNeuronsPCAvalsTar;
        % PCAvals = {PCAvals{1,PCAchansIdx}};

        PCAref      = PCAvalsRef{cluIdxRef};
        PCAtar      = PCAvalsTar{cluIdxTar};
        PCAsynchRef = PCAref(refSpikeTimesSynchIdx(refSpikeTimesSynchIdx < size(PCAref,1)),:); % needed cutoff
        PCAsynchTar = PCAtar(tarSpikeTimesSynchIdx(tarSpikeTimesSynchIdx < size(PCAtar,1)),:); % needed cutoff

        %% template indice split 

        spikeTemplatesTar             = clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)};
        spikeTemplatesSynchTar        = clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)}(tarSpikeTimesSynchIdx(tarSpikeTimesSynchIdx < size(PCAtar,1))); % needed cutoff 

        primaryTemplateTar            = mode(spikeTemplatesTar);

        primaryTemplateTarIdx         = find(spikeTemplatesTar == primaryTemplateTar); 
        secondaryTemplateTarIdx       = find(spikeTemplatesTar ~= primaryTemplateTar);

        primaryTemplatesSynchTarIdx   = find(spikeTemplatesSynchTar == primaryTemplateTar); 
        secondaryTemplatesSynchTarIdx = find(spikeTemplatesSynchTar ~= primaryTemplateTar);

        spikeTemplatesRef             = clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)};
        spikeTemplatesSynchRef        = clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)}(refSpikeTimesSynchIdx(refSpikeTimesSynchIdx < size(PCAref,1)));

        primaryTemplateRef            = mode(spikeTemplatesRef);

        primaryTemplateRefIdx         = find(spikeTemplatesRef == primaryTemplateRef); 
        secondaryTemplateRefIdx       = find(spikeTemplatesRef ~= primaryTemplateRef);

        primaryTemplatesSynchRefIdx   = find(spikeTemplatesSynchRef == primaryTemplateRef); 
        secondaryTemplatesSynchRefIdx = find(spikeTemplatesSynchRef ~= primaryTemplateRef);

        %% getting the waveform snippets 

        % timeShift      = 514;
        % chanData       = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarCh) 'highpass300hz.mat']);

        if loopRat == 1
            chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(tarChan) 'highpass300hz.mat']);
            chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(refChan) 'highpass300hz.mat']);
        elseif loopRat == 2
            chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(tarChan) 'highpass300hz.mat']);
            chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(refChan) 'highpass300hz.mat']);
        elseif loopRat == 3
            chanDataTar = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(tarChan) 'highpass300hz.mat']);
            chanDataRef = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(refChan) 'highpass300hz.mat']);
        elseif loopRat == 4
            chanDataRef = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(refShank) '/ch' num2str(refChan) 'highpass300hz.mat']);
            chanDataTar = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(tarShank) '/ch' num2str(tarChan) 'highpass300hz.mat']);
        elseif loopRat == 5
            chanDataRef = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(refShank) '/ch' num2str(refChan) 'highpass300hz.mat']);
            chanDataTar = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(tarShank) '/ch' num2str(tarChan) 'highpass300hz.mat']);
        elseif loopRat == 6
            chanDataRef = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(refChan) 'highpass300hz.mat']);
            chanDataTar = load(['/media/nasko/WD_BLACK31/HiroRawDataPerChannelFiltered/Roy/maze/highpass300hz/ch' num2str(tarChan) 'highpass300hz.mat']);
        end

        if loopRat == 6
            refLat       = round((refSpikeTimes       - chanDataRef.data.onsetTime/1e6)*fs) + 1;
            refSynchLat  = round((refSpikeTimesSynch  - chanDataRef.data.onsetTime/1e6)*fs) + 1;

            tarLat       = round((tarSpikeTimes       - chanDataTar.data.onsetTime/1e6)*fs) + 1;
            tarSynchLat  = round((tarSpikeTimesSynch  - chanDataTar.data.onsetTime/1e6)*fs) + 1;
        else
            refLat       = round(refSpikeTimes*fs)';
            refSynchLat  = round(refSpikeTimesSynch*fs)';

            tarLat       = round(tarSpikeTimes*fs)';
            tarSynchLat  = round(tarSpikeTimesSynch*fs)';
        end

        %%

        refLatPrimaryTemplate        = refLat(primaryTemplateRefIdx);
        refLatSecondaryTemplate      = refLat(secondaryTemplateRefIdx);

        refSynchLatPrimaryTemplate   = refSynchLat(primaryTemplatesSynchRefIdx);
        refSynchLatSecondaryTemplate = refSynchLat(secondaryTemplatesSynchRefIdx);

        tarLatPrimaryTemplate        = tarLat(primaryTemplateTarIdx);
        tarLatSecondaryTemplate      = tarLat(secondaryTemplateTarIdx);

        tarSynchLatPrimaryTemplate   = tarSynchLat(primaryTemplatesSynchTarIdx);
        tarSynchLatSecondaryTemplate = tarSynchLat(secondaryTemplatesSynchTarIdx);

        %%

        if length(refLat) > NspikesPlot
            refLatRandPerm = sort(refLat(randperm(length(refLat),NspikesPlot)));
        else
            refLatRandPerm = refLat;
        end

        if length(refLatPrimaryTemplate) > NspikesPlot
            refLatPrimaryTemplateRandPerm = sort(refLatPrimaryTemplate(randperm(length(refLatPrimaryTemplate),NspikesPlot)));
        else
            refLatPrimaryTemplateRandPerm = refLatPrimaryTemplate;
        end

        if length(refLatSecondaryTemplate) > NspikesPlot
            refLatSecondaryTemplateRandPerm = sort(refLatSecondaryTemplate(randperm(length(refLatSecondaryTemplate),NspikesPlot)));
        else
            refLatSecondaryTemplateRandPerm = refLatSecondaryTemplate;
        end

        if length(refSynchLat) > NspikesPlot
            refSynchLatRandPerm = sort(refSynchLat(randperm(length(refSynchLat),NspikesPlot)));
        else
            refSynchLatRandPerm = refSynchLat;
        end

        if length(refSynchLatPrimaryTemplate) > NspikesPlot
            refSynchLatPrimaryTemplateRandPerm = sort(refSynchLatPrimaryTemplate(randperm(length(refSynchLatPrimaryTemplate),NspikesPlot)));
        else
            refSynchLatPrimaryTemplateRandPerm = refSynchLatPrimaryTemplate;
        end

        if length(refSynchLatSecondaryTemplate) > NspikesPlot
            refSynchLatSecondaryTemplateRandPerm = sort(refSynchLatSecondaryTemplate(randperm(length(refSynchLatSecondaryTemplate),NspikesPlot)));
        else
            refSynchLatSecondaryTemplateRandPerm = refSynchLatSecondaryTemplate;
        end

        if loopRat == 6
            [refSponWaveMean,      refSponWaveforms]                                       = waveformAvg(chanDataRef.data.channel,refLatRandPerm,     preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [refSynchLatWaveMean,  refPeakLatWaveforms]                                    = waveformAvg(chanDataRef.data.channel,refSynchLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
        else
            [refSponWaveMean,      refSponWaveforms]                                       = waveformAvg(chanDataRef.data,        refLatRandPerm,      preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [refSynchLatWaveMean,  refPeakLatWaveforms]                                    = waveformAvg(chanDataRef.data,        refSynchLatRandPerm, preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

            [refSponWaveMeanPrimaryTemplate,      refSponWaveformsPrimaryTemplate]         = waveformAvg(chanDataRef.data,refLatPrimaryTemplateRandPerm,        preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [refSynchLatWaveMeanPrimaryTemplate,  refPeakLatWaveformsPrimaryTemplate]      = waveformAvg(chanDataRef.data,refSynchLatPrimaryTemplateRandPerm,   preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

            [refSponWaveMeanSecondaryTemplate,      refSponWaveformsSecondaryTemplate]     = waveformAvg(chanDataRef.data,refLatSecondaryTemplateRandPerm,      preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [refSynchLatWaveMeanSecondaryTemplate,  refPeakLatWaveformsSecondaryTemplate]  = waveformAvg(chanDataRef.data,refSynchLatSecondaryTemplateRandPerm, preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
        end

        if length(tarLat) > NspikesPlot
            tarLatRandPerm = sort(tarLat(randperm(length(tarLat),NspikesPlot)));
        else
            tarLatRandPerm = tarLat;
        end

        if length(tarLatPrimaryTemplate) > NspikesPlot
            tarLatPrimaryTemplateRandPerm = sort(tarLatPrimaryTemplate(randperm(length(tarLatPrimaryTemplate),NspikesPlot)));
        else
            tarLatPrimaryTemplateRandPerm = tarLatPrimaryTemplate;
        end

        if length(tarLatSecondaryTemplate) > NspikesPlot
            tarLatSecondaryTemplateRandPerm = sort(tarLatSecondaryTemplate(randperm(length(tarLatSecondaryTemplate),NspikesPlot)));
        else
            tarLatSecondaryTemplateRandPerm = tarLatSecondaryTemplate;
        end

        if length(tarSynchLat) > NspikesPlot
            tarSynchLatRandPerm = sort(tarSynchLat(randperm(length(tarSynchLat),NspikesPlot)));
        else
            tarSynchLatRandPerm = tarSynchLat;
        end

        if length(tarSynchLatPrimaryTemplate) > NspikesPlot
            tarSynchLatPrimaryTemplateRandPerm = sort(tarSynchLatPrimaryTemplate(randperm(length(tarSynchLatPrimaryTemplate),NspikesPlot)));
        else
            tarSynchLatPrimaryTemplateRandPerm = tarSynchLatPrimaryTemplate;
        end

        if length(tarSynchLatSecondaryTemplate) > NspikesPlot
            tarSynchLatSecondaryTemplateRandPerm = sort(tarSynchLatSecondaryTemplate(randperm(length(tarSynchLatSecondaryTemplate),NspikesPlot)));
        else
            tarSynchLatSecondaryTemplateRandPerm = tarSynchLatSecondaryTemplate;
        end


        if loopRat == 6
            [tarSponWaveMean,      tarSponWaveforms]     = waveformAvg(chanDataTar.data.channel, tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [tarSynchLatWaveMean,  tarPeakLatWaveforms]  = waveformAvg(chanDataTar.data.channel, tarSynchLatRandPerm,    preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
        else
            [tarSponWaveMean,      tarSponWaveforms]     = waveformAvg(chanDataTar.data,         tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [tarSynchLatWaveMean,  tarPeakLatWaveforms]  = waveformAvg(chanDataTar.data,         tarSynchLatRandPerm,    preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

            [tarSponWaveMeanPrimaryTemplate,      tarSponWaveformsPrimaryTemplate]         = waveformAvg(chanDataTar.data,tarLatPrimaryTemplateRandPerm,        preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [tarSynchLatWaveMeanPrimaryTemplate,  tarPeakLatWaveformsPrimaryTemplate]      = waveformAvg(chanDataTar.data,tarSynchLatPrimaryTemplateRandPerm,   preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

            [tarSponWaveMeanSecondaryTemplate,      tarSponWaveformsSecondaryTemplate]     = waveformAvg(chanDataTar.data,tarLatSecondaryTemplateRandPerm,      preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
            [tarSynchLatWaveMeanSecondaryTemplate,  tarPeakLatWaveformsSecondaryTemplate]  = waveformAvg(chanDataTar.data,tarSynchLatSecondaryTemplateRandPerm, preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
        end


        %%    

        synchPlotColor = 'r';

        figure(102)

        hcomb = figure(102);

        res_type = 'QHD';
        pos = [70 230 (3/5)*1920*0.5 1080*0.6]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

        tiledlayout(6,3,'Padding','none','TileSpacing','compact')

        tR     = pairGroupStatTable.CCGbinLagTimes{loopPairs}*1000;
        ccgR   = pairGroupStatTable.pairRawCCG{loopPairs};
        GSPExc = pairGroupStatTable.GSPExc{loopPairs};

        nexttile(1)
        plot(tR,ccgR,'k','LineWidth',1)
        hold on
        scatter(tR(find(GSPExc)),ccgR(find(GSPExc)),100,'ko','LineWidth',1)
        hold off

        xlim([-1,1])
        ylims = get(gca,'ylim');
        ylabel('Spike Probability')
        xlabel('[ms]')
        title(titleStrTar)
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        nexttile(10)
        plot(tR,flipud(ccgR),'k','LineWidth',1)
        hold on
        scatter(tR(find(flipud(GSPExc))),flipud(ccgR(find(GSPExc))),100,'ko','LineWidth',1)
        hold off

        xlim([-1,1])
        ylims = get(gca,'ylim');
        ylabel('Spike Probability')
        xlabel('[ms]')
        title(titleStrRef)
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off

        %%

        nexttile(2)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         tarSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(tarSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(tarPeakLatWaveformsPrimaryTemplate,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPeakLatWaveformsPrimaryTemplate(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor','r')
        end
        for i = 1:size(tarPeakLatWaveformsSecondaryTemplate,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     tarPeakLatWaveformsSecondaryTemplate(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor','b')
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        title('waveform snippets')
        box off

        nexttile(3)
        plot((-preLength:postLength-1)*(30/1000),tarSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),tarSynchLatWaveMeanPrimaryTemplate,'r','LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),tarSynchLatWaveMeanSecondaryTemplate,'b','LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        title('waveform averages')
        legend('spon spikes','synch largest template','synch all other templates','location','eastoutside')
        box off

        nexttile(11)
        patchline((-preLength:postLength-1)*(30/1000), ... 
                         refSponWaveforms(:,1), ...
                         'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        hold on
        for i = 2:size(refSponWaveforms,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refSponWaveforms(:,i), ...
                     'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
        end
        for i = 1:size(refPeakLatWaveformsPrimaryTemplate,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPeakLatWaveformsPrimaryTemplate(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor','r')
        end
        for i = 1:size(refPeakLatWaveformsSecondaryTemplate,2)
            patchline((-preLength:postLength-1)*(30/1000), ... 
                     refPeakLatWaveformsSecondaryTemplate(:,i), ...
                     'linewidth',1,'edgealpha',0.1,'edgecolor','b')
        end
        hold off
        % ylim([-6 2])
        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        title('waveform snippets')
        box off

        nexttile(12)
        plot((-preLength:postLength-1)*(30/1000),refSponWaveMean,'Color',[.7 .7 .7],'LineWidth',1)
        hold on
        plot((-preLength:postLength-1)*(30/1000),refSynchLatWaveMeanPrimaryTemplate,'r','LineWidth',1)
        plot((-preLength:postLength-1)*(30/1000),refSynchLatWaveMeanSecondaryTemplate,'b','LineWidth',1)
        hold off

        xlim([-1,1])
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        xlabel('[ms]');
        ylabel('[mV]');
        set(gca, 'YDir','reverse')
        legend('spon spikes','synch largest template','synch all other templates','location','eastoutside')
        box off

        %%

        nexttile(4)
        scatter(PCAtar(primaryTemplateTarIdx,1),              PCAtar(primaryTemplateTarIdx,2),              [], [.7 .7 .7])
        hold on
        scatter(PCAtar(secondaryTemplateTarIdx,1),            PCAtar(secondaryTemplateTarIdx,2),            [], [0 0 0])
        scatter(PCAsynchTar(primaryTemplatesSynchTarIdx,1),   PCAsynchTar(primaryTemplatesSynchTarIdx,2),   [], synchPlotColor)
        scatter(PCAsynchTar(secondaryTemplatesSynchTarIdx,1), PCAsynchTar(secondaryTemplatesSynchTarIdx,2), [], [0 0 1])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        % xlim([-2500 2500])
        % ylim([-2500 2500])
        daspect([1 1 1])
        title('PC1 v PC2')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
%             legend('noSynch','synch','location','southoutside')

        nexttile(5)
        scatter(PCAtar(primaryTemplateTarIdx,2),              PCAtar(primaryTemplateTarIdx,3),              [], [.7 .7 .7])
        hold on
        scatter(PCAtar(secondaryTemplateTarIdx,2),            PCAtar(secondaryTemplateTarIdx,3),            [], [0 0 0])
        scatter(PCAsynchTar(primaryTemplatesSynchTarIdx,2),   PCAsynchTar(primaryTemplatesSynchTarIdx,3),   [], synchPlotColor)
        scatter(PCAsynchTar(secondaryTemplatesSynchTarIdx,2), PCAsynchTar(secondaryTemplatesSynchTarIdx,3), [], [0 0 1])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        % xlim([-2500 2500])
        % ylim([-2500 2500])
        daspect([1 1 1])
        title('PC1 v PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
%             legend('noSynch','synch','location','southoutside')
%             title(["ch: " num2str(tarCh)])

        nexttile(6)
        scatter(PCAtar(primaryTemplateTarIdx,1),                    PCAtar(primaryTemplateTarIdx,3),              [], [.7 .7 .7])
        hold on
        scatter(PCAtar(secondaryTemplateTarIdx,1),                  PCAtar(secondaryTemplateTarIdx,3),            [], [0 0 0])
        scatter(PCAsynchTar(primaryTemplatesSynchTarIdx,1),         PCAsynchTar(primaryTemplatesSynchTarIdx,3),   [], synchPlotColor)
        scatter(PCAsynchTar(secondaryTemplatesSynchTarIdx,1),       PCAsynchTar(secondaryTemplatesSynchTarIdx,3), [], [0 0 1])

        hold off
        xlabel('PC2')
        ylabel('PC3')
        % xlim([-2500 2500])
        % ylim([-2500 2500])
        daspect([1 1 1])
        title('PC2 v PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        legend('spon largest template','synch largest template', ...
               'spon all other templates','synch all other templates', ...
               'Location','eastoutside')

        [h1,edges1] = histcounts(spikeTemplatesTar, ... 
            'Normalization','Probability','BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)})]);
        [h2,edges2] = histcounts(spikeTemplatesSynchTar, ... 
            'Normalization','Probability','BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)})]);
        h1counts = histcounts(spikeTemplatesTar, ... 
            'BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)})]);
        h2counts = histcounts(spikeTemplatesSynchTar, ... 
            'BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsTar.spikeTemplates{find(clusterValNeuronsTar.cellID == tarID)})]);

        bfilt = ~((h1 <= 0.01) & (h2 <= 0.01));
        % bfilt = logical(ones(1,length(h1)));

        nexttile([1 3])
        b = bar(categorical(edges1(bfilt)),[h1(bfilt); h2(bfilt)]');

        b(1).FaceColor = 'flat';
        b(2).FaceColor = 'flat';
        for loopColors = 1:sum(bfilt)
            if isscalar(edges1(bfilt))
                b(1).CData(loopColors,:) = [.7 .7 .7];
                b(2).CData(loopColors,:) = [1 0 0];
            elseif loopColors == find(edges1(bfilt) == primaryTemplateTar) 
                b(1).CData(loopColors,:) = [.7 .7 .7];
                b(2).CData(loopColors,:) = [1 0 0];
            else
                b(1).CData(loopColors,:) = [0 0 0];
                b(2).CData(loopColors,:) = [0 0 1];
            end
        end
        ylabel('Probability')
        xlabel('Template #')
        title('cluster template composition (normalized by number of spikes in each template)')
        text(0.3,0.5,{['largest template: ' num2str(primaryTemplateTar)], ...
                        ['largest template non-synch: '   num2str(h1counts(edges1 == primaryTemplateTar)) '/' num2str(sum(h1counts))],...
                        ['largest template synch: '       num2str(h2counts(edges1 == primaryTemplateTar)) '/' num2str(sum(h2counts))],...
                        ['all other template non-synch: ' num2str(sum(h1counts(edges1(1:end-1) ~= primaryTemplateTar))) '/' num2str(sum(h1counts))],...
                        ['all other template synch: '     num2str(sum(h2counts(edges1(1:end-1) ~= primaryTemplateTar))) '/' num2str(sum(h2counts))]},...
                        'FontSize',5,'FontName','Arial')
        % legend('spon spikes','synch spikes','Location','eastoutside')
        box off
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')

        nexttile(13)
        scatter(PCAref(primaryTemplateRefIdx,1),              PCAref(primaryTemplateRefIdx,2),              [], [.7 .7 .7])
        hold on
        scatter(PCAref(secondaryTemplateRefIdx,1),            PCAref(secondaryTemplateRefIdx,2),            [], [0 0 0])
        scatter(PCAsynchRef(primaryTemplatesSynchRefIdx,1),   PCAsynchRef(primaryTemplatesSynchRefIdx,2),   [], synchPlotColor)
        scatter(PCAsynchRef(secondaryTemplatesSynchRefIdx,1), PCAsynchRef(secondaryTemplatesSynchRefIdx,2), [], [0 0 1])
        hold off
        xlabel('PC1')
        ylabel('PC2')
        % xlim([-2500 2500])
        % ylim([-2500 2500])
        daspect([1 1 1])
        title('PC1 v PC2')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
%             legend('noSynch','synch','location','southoutside')

        nexttile(14)
        scatter(PCAref(primaryTemplateRefIdx,2),              PCAref(primaryTemplateRefIdx,3),              [], [.7 .7 .7])
        hold on
        scatter(PCAref(secondaryTemplateRefIdx,2),            PCAref(secondaryTemplateRefIdx,3),            [], [0 0 0])
        scatter(PCAsynchRef(primaryTemplatesSynchRefIdx,2),   PCAsynchRef(primaryTemplatesSynchRefIdx,3),   [], synchPlotColor)
        scatter(PCAsynchRef(secondaryTemplatesSynchRefIdx,2), PCAsynchRef(secondaryTemplatesSynchRefIdx,3), [], [0 0 1])
        hold off
        xlabel('PC1')
        ylabel('PC3')
        % xlim([-2500 2500])
        % ylim([-2500 2500])
        daspect([1 1 1])
        title('PC1 v PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
%             legend('noSynch','synch','location','southoutside')
%             title(["ch: " num2str(tarCh)])

        nexttile(15)
        scatter(PCAref(primaryTemplateRefIdx,1),              PCAref(primaryTemplateRefIdx,3),              [], [.7 .7 .7])
        hold on
        scatter(PCAref(secondaryTemplateRefIdx,1),            PCAref(secondaryTemplateRefIdx,3),            [], [0 0 0])
        scatter(PCAsynchRef(primaryTemplatesSynchRefIdx,1),   PCAsynchRef(primaryTemplatesSynchRefIdx,3),   [], synchPlotColor)
        scatter(PCAsynchRef(secondaryTemplatesSynchRefIdx,1), PCAsynchRef(secondaryTemplatesSynchRefIdx,3), [], [0 0 1])
        hold off
        xlabel('PC2')
        ylabel('PC3')
        % xlim([-2500 2500])
        % ylim([-2500 2500])
        daspect([1 1 1])
        title('PC2 v PC3')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        legend('spon spikes primary','synch spikes primary', ...
               'spon spikes secondary','synch spikes secondary', ...
               'Location','eastoutside')

        nexttile([1 3])
        [h1,edges1] = histcounts(spikeTemplatesRef, ... 
            'Normalization','Probability','BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)})]);
        [h2,edges2] = histcounts(spikeTemplatesSynchRef, ... 
            'Normalization','Probability','BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)})]);
        h1counts = histcounts(spikeTemplatesRef, ... 
            'BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)})]);
        h2counts = histcounts(spikeTemplatesSynchRef, ... 
            'BinWidth',1, ...
            'BinLimits',[0,max(clusterValNeuronsRef.spikeTemplates{find(clusterValNeuronsRef.cellID == refID)})]);

        bfilt = ~((h1 <= 0.01) & (h2 <= 0.01));
        % bfilt = logical(ones(1,length(h1)));

        b = bar(categorical(edges1(bfilt)),[h1(bfilt); h2(bfilt)]');

        b(1).FaceColor = 'flat';
        b(2).FaceColor = 'flat';
        for loopColors = 1:sum(bfilt)
            if isscalar(edges1(bfilt))
                b(1).CData(loopColors,:) = [.7 .7 .7];
                b(2).CData(loopColors,:) = [1 0 0];
            elseif loopColors == find(edges1(bfilt) == primaryTemplateRef) 
                b(1).CData(loopColors,:) = [.7 .7 .7];
                b(2).CData(loopColors,:) = [1 0 0];
            else
                b(1).CData(loopColors,:) = [0 0 0];
                b(2).CData(loopColors,:) = [0 0 1];
            end
        end
        ylabel('Probability')
        xlabel('Template #')
        title('cluster template composition (normalized by number of spikes in each template)')
        text(0.3,0.5,{['largest template: ' num2str(primaryTemplateRef)], ...
                        ['largest template non-synch: '   num2str(h1counts(edges1 == primaryTemplateRef)) '/' num2str(sum(h1counts))],...
                        ['largest template synch: '       num2str(h2counts(edges1 == primaryTemplateRef)) '/' num2str(sum(h2counts))],...
                        ['all other template non-synch: ' num2str(sum(h1counts(edges1(1:end-1) ~= primaryTemplateRef))) '/' num2str(sum(h1counts))],...
                        ['all other template synch: '     num2str(sum(h2counts(edges1(1:end-1) ~= primaryTemplateRef))) '/' num2str(sum(h2counts))]},...
                        'FontSize',5,'FontName','Arial')
        % legend('spon spikes','synch spikes','Location','eastoutside')
        box off
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')


        save_file = fullfile(figPath, saveStr);
        print(fig_use, save_file,'-djpeg',resolution_use);


        %%


        % PCA = [];
        % for loopPCAchans = sum(PCAchansIdx)
        %     PCA = [PCA double(PCAvals{1,loopPCAchans})];
        % end
        % 
        % COM = mean(PCA);
        % 
        % for loopSpikeTimes = 1:size(PCA,1)
        %     D(loopSpikeTimes)   = sqrt(sum((COM - PCA(loopSpikeTimes,:)) .^ 2));
        % end
        % 
        % %%
        % [~,~,idx] = intersect(pairGroupStatTable.tarSynch{loopPairs},pairGroupStatTable.tarSpikeTimes{loopPairs});
        % for loopSpikeTimes = 1:size(idx,1)
        %     Dsynch(loopSpikeTimes)   = sqrt(sum((COM - PCA(idx(loopSpikeTimes),:)) .^ 2));
        % end
        % 
        % histogram(D,     'Normalization','probability','BinWidth',1)
        % hold on
        % histogram(Dsynch,'Normalization','probability','BinWidth',1)
        % hold off

    end

end


% % temp = intersect(cell2num(clusterValNeuron.templateChannels));

%             % find common channels
%             commonChannelsNonNoise = intersect(clusterValNeuron.templateChannels{1},clusterValNeuron.templateChannels{2});
%             for loopNeurons = 3:sum(~noiseIdxBinary)
%                 commonChannelsNonNoise = intersect(commonChannelsNonNoise,clusterValNeuron.templateChannels{loopNeurons});
%             end

%             % check which noise clusters share good unit template channels 
%             removeNoiseUnitIdx = [];
%             for loopNeurons = 1:length(noiseIdx)
%                 if isempty(intersect(commonChannelsNonNoise,clusterValNeuron.templateChannels{noiseIdx(loopNeurons)}))
%                     removeNoiseUnitIdx = [removeNoiseUnitIdx noiseIdx(loopNeurons)];
%                 end
%             end
% 
%             clusterValNeuron(removeNoiseUnitIdx,:) = [];
% 
%             % new noise detection
%             noiseIdxBinary = contains(string(clusterValNeuron.cellID),'999');
%             noiseIdx       = find(noiseIdxBinary);
% 
%             commonChannels = commonChannelsNonNoise;
%             for loopNeurons = 1:length(noiseIdx)
%                 commonChannels = intersect(commonChannels,clusterValNeuron.templateChannels{noiseIdx(loopNeurons)});
%             end
% 
%             if isempty(commonChannels)
%                continue
%             end

% usedShanks = unique(neuronShank);
% for loopShank = 1:length(usedShanks)
%     
%     loopShank 
%     
%     idx = usedShanks(loopShank) == shank_list;
%     
%     featuresTemp      = features(idx,:);
%     cluster_indexTemp = cluster_index(idx,:);
%     
%     isolation_distance = [isolation_distance; IsolationDistance_calc(featuresTemp,cluster_indexTemp)];
%     L_ratio            = [L_ratio;            L_ratio_calc(          featuresTemp,cluster_indexTemp)];
%     
% end 

% %%
% 
% temp = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/similar_templates.npy'));
% 
% temp = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/similar_templates.npy'));
% 
% temp(temp == 0) = NaN;

% temp = nanmean(temp);
% 
% histogram(temp,'Normalization','probability','BinWidth',0.01)

%% global snapshots

% figPath = '/media/nasko/WD_BLACK31/ClusterAnaTemp/figures/';
% 
% resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
% fig_use = 102;
% 
% animalNames = {'N','S','U','AG1','AG2','Roy'};
% 
% for loopRat = 1:6
% 
%     if loopRat == 1
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStatsClusterAna.mat')  
% 
%         channelSet = {1:16, ...
%                       17:32, ... 
%                       33:48, ...
%                       49:64, ...
%                       65:80, ...
%                       81:96, ...
%                       97:112, ...
%                       113:128};
% 
%     elseif loopRat == 2
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStatsClusterAna.mat') 
% 
%         channelSet = {1:10, ...
%                       11:21, ...
%                       22:32, ... 
%                       33:43, ...
%                       44:54, ...
%                       55:63, ...
%                       64:79, ...
%                       80:95, ...
%                       96:111, ...
%                       112:127, ...
%                       128:143, ...
%                       144:159, ...
%                       160:175, ...
%                       176:191};
% 
%     elseif loopRat == 3
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStatsClusterAna.mat') 
% 
%         channelSet = {1:16, ...
%                       17:32, ... 
%                       33:48, ...
%                       49:64, ...
%                       65:80, ...
%                       81:96, ...
%                       97:112, ...
%                       113:128, ...
%                       129:144, ...
%                       145:160, ...
%                       161:176, ...
%                       177:192};
% 
%     elseif loopRat == 4
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStatsClusterAna.mat')
% 
%         channelSet = {1:32, ...
%                       33:64, ... 
%                       65:96, ...
%                       97:128, ...
%                       129:160, ...
%                       161:192};
% 
%     elseif loopRat == 5
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStatsClusterAna.mat')
% 
%         channelSet = {1:32, ...
%                       33:64, ... 
%                       65:96, ...
%                       97:128, ...
%                       129:160, ...
%                       161:192};
% 
%     elseif loopRat == 6
%         load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStatsClusterAna.mat')
% 
%         channelSet = {1:8, ...
%                       9:16, ... 
%                       17:24, ...
%                       25:32, ...
%                       33:40, ...
%                       41:48, ...
%                       48:56, ...
%                       57:64};
%     end
% 
%     tWave = pairGroupStatTable.tWave{1};
% 
%     for loopPairs = 1:size(pairGroupStatTable,1)
% 
%         refShank           = pairGroupStatTable.refShank(loopPairs);
%         refChan            = pairGroupStatTable.refChannel(loopPairs);
%         refID              = pairGroupStatTable.refNeuronID(loopPairs);
%         refSpikeTimes      = pairGroupStatTable.refSpikeTimes{loopPairs};
%         refSpikeTimesSynch = pairGroupStatTable.refSynch{loopPairs};
%         refCellType        = pairGroupStatTable.refCellExplorerType{loopPairs};
%         [~,refSpikeTimesSynchIdx,~] = intersect(refSpikeTimes,refSpikeTimesSynch);
% 
%         tarShank           = pairGroupStatTable.tarShank(loopPairs);
%         tarChan            = pairGroupStatTable.tarChannel(loopPairs);
%         tarID              = pairGroupStatTable.tarNeuronID(loopPairs);
%         tarSpikeTimes      = pairGroupStatTable.tarSpikeTimes{loopPairs};
%         tarSpikeTimesSynch = pairGroupStatTable.tarSynch{loopPairs};
%         tarCellType        = pairGroupStatTable.tarCellExplorerType{loopPairs};
%         [~,tarSpikeTimesSynchIdx,~] = intersect(tarSpikeTimes,tarSpikeTimesSynch);
% 
%         tWave = pairGroupStatTable.tWave{loopPairs};
%         tR    = pairGroupStatTable.CCGbinLagTimes{loopPairs,1};
% 
%         refStr = ['rat ' animalNames{loopRat} ': ' num2str(refID) refCellType ' (sh: ' num2str(refShank) ')'];
%         tarStr = ['rat ' animalNames{loopRat} ': ' num2str(tarID) tarCellType ' (sh: ' num2str(tarShank) ')'];
% 
%         %% exquisite partners 
% 
%         refPartnersIdx      = fliplr(refID == [pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID]);
%         tarPartnersIdx      = fliplr(tarID == [pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID]);
% 
%         refPartnersID       = [pairGroupStatTable.refNeuronID(refPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarNeuronID(refPartnersIdx(:,2))];
%         tarPartnersID       = [pairGroupStatTable.refNeuronID(tarPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarNeuronID(tarPartnersIdx(:,2))];
% 
%         refPartnersChan     = [pairGroupStatTable.refChannel(refPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarChannel(refPartnersIdx(:,2))];
%         tarPartnersChan     = [pairGroupStatTable.refChannel(tarPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarChannel(tarPartnersIdx(:,2))];
% 
%         refPartnersShank    = [pairGroupStatTable.refShank(refPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarShank(refPartnersIdx(:,2))];
%         tarPartnersShank    = [pairGroupStatTable.refShank(tarPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarShank(tarPartnersIdx(:,2))];
% 
%         refPartnersCellType = [pairGroupStatTable.refCellExplorerType(refPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarCellExplorerType(refPartnersIdx(:,2))];
%         tarPartnersCellType = [pairGroupStatTable.refCellExplorerType(tarPartnersIdx(:,1)); ...
%                                pairGroupStatTable.tarCellExplorerType(tarPartnersIdx(:,2))];
% 
%         refPartnersDistance = [pairGroupStatTable.pairDistance(refPartnersIdx(:,1)); ...
%                                pairGroupStatTable.pairDistance(refPartnersIdx(:,2))];
%         tarPartnersDistance = [pairGroupStatTable.pairDistance(tarPartnersIdx(:,1)); ...
%                                pairGroupStatTable.pairDistance(tarPartnersIdx(:,2))];
% 
% 
%         refPartnersCCGs     = [cell2num(pairGroupStatTable.pairRawCCG') ... 
%                                flipud(cell2num(pairGroupStatTable.pairRawCCG'))];
%         tarPartnersCCGs     = [cell2num(pairGroupStatTable.pairRawCCG') ... 
%                                flipud(cell2num(pairGroupStatTable.pairRawCCG'))];
% 
%         refPartnersCCGs(:,~[refPartnersIdx(:,1); refPartnersIdx(:,2)]) = [];
%         tarPartnersCCGs(:,~[tarPartnersIdx(:,1); tarPartnersIdx(:,2)]) = [];
% 
%         % refPartnersCCGs     = {cell2num(pairGroupStatTable.pairRawCCG{refPartnersIdx(:,1),1}) ... 
%         %                        cell2num(pairGroupStatTable.pairRawCCG{refPartnersIdx(:,2),1})};
%         % tarPartnersCCGs     = {cell2num(pairGroupStatTable.pairRawCCG{tarPartnersIdx(:,1),1}) ... 
%         %                        cell2num(pairGroupStatTable.pairRawCCG{tarPartnersIdx(:,2),1})};
% 
%         %%
% 
%         % saveStr  = ['rat ' animalNames{loopRat} ': '...
%         %             num2str(refID) refCellType ' (sh ' num2str(refShank) ')' ' - ' ...
%         %             num2str(tarID) tarCellType ' (sh ' num2str(tarShank) ')' ' pair ' num2str(d) ' microns.jpeg'];
% 
%         if loopRat == 2
%             refWaveforms = pairGroupStatTable.refWaveforms{loopPairs};
%             tarWaveforms = pairGroupStatTable.tarWaveforms{loopPairs};
% 
%             refWaveforms(sum(refWaveforms,2) == 0,:) = [];
%             tarWaveforms(sum(tarWaveforms,2) == 0,:) = [];
%         else
%             refWaveforms = pairGroupStatTable.refWaveforms{loopPairs};
%             tarWaveforms = pairGroupStatTable.tarWaveforms{loopPairs};
%         end
% 
%         %% ref
% 
%         hcomb = figure(102);
% 
%         res_type = 'QHD';
%         pos = [1720 2562 3*560 2*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
%         arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% 
%         tiledlayout(3,size(channelSet,2),'Padding','none','TileSpacing','compact')
% 
%         for loopSets = 1:size(channelSet,2)
% 
%             tempTrough(loopSets) = max(abs(min(refWaveforms(channelSet{loopSets},:)')));
%             tempPeak(loopSets)   = max(abs(max(refWaveforms(channelSet{loopSets},:)')));
% 
%         end
% 
%         tempTrough = -max(tempTrough) - 0.1;
%         tempPeak   =  max(tempPeak)   + 0.1;
% 
%         for loopSets = 1:size(channelSet,2)
% 
%             partnerStrList = {};
%             if sum(loopSets == refPartnersShank) ~= 0
%                 tempIdx = find(loopSets == refPartnersShank);
%                 for loopTemp = 1:length(tempIdx)
% 
%                     partnerStrList = [partnerStrList; [num2str(refPartnersID(tempIdx(loopTemp))) ...
%                                                        refPartnersCellType{tempIdx(loopTemp)} ' (ch: ' ... 
%                                                        num2str(refPartnersChan(tempIdx(loopTemp))) ', dist: ' ...
%                                                        num2str(refPartnersDistance(tempIdx(loopTemp))) ')']];
% 
%                 end
%             end
% 
%             nexttile(loopSets)
%             plot(tWave,refWaveforms(channelSet{loopSets},:)');
%             if sum(loopSets == refPartnersShank) ~= 0
%                 text(-0.9,tempTrough + 0.05,partnerStrList,...
%                      'FontSize',5,'FontName','Arial')
%             end
% 
%             if loopSets == 1
%                 title([refStr ' shank ' num2str(loopSets)])
%             else 
%                 title([' shank ' num2str(loopSets)])
%             end
% 
%             xlim([-1,1])
%             ylim([tempTrough tempPeak])
%             set(gca,'FontSize',5)
%             set(gca,'FontName','Arial')
%             xlabel('[ms]');
%             ylabel('[mV]');
%             set(gca, 'YDir','reverse')
%             box off
% 
%             nexttile(loopSets + size(channelSet,2))
%             plot(tWave,refWaveforms(channelSet{loopSets},:)');
%             xlim([-1,1])
%             set(gca,'FontSize',5)
%             set(gca,'FontName','Arial')
%             xlabel('[ms]');
%             ylabel('[mV]');
%             set(gca, 'YDir','reverse')
%             box off
% 
%             if sum(loopSets == refPartnersShank) ~= 0
%                 nexttile(loopSets + 2*size(channelSet,2))
%                 tempIdx = find(loopSets == refPartnersShank);
% 
%                 stackedplot(tR*1000,refPartnersCCGs(:,tempIdx),'k','DisplayLabels',repmat({"Counts"},1,size(refPartnersCCGs(:,tempIdx),2)))
%                 % stackedplot(tR*1000,refPartnersCCGs,'k','DisplayLabels',repmat({"Counts"},1,size(refPartnersCCGs,2)))
% 
%                 % plot(tR*1000,refPartnersCCGs{1}./length(refSpikeTimes),'k')
% 
%                 % hold on
%                 % for loopTemp = 2:length(tempIdx)
%                 %     stkplt
%                 %     plot(tR*1000,refPartnersCCGs{tempIdx(loopTemp)}./length(refSpikeTimes),'k')
%                 % end
%                 % hold off
% 
%                 xlim([-1,1])
%                 % ylims = get(gca,'ylim');
%                 % if length(tempIdx) > 1
%                 %     set(gca,'ytick',[])
%                 % else 
%                     % ylabel('Spike Probability')
%                 % end 
%                 xlabel('[ms]')
%                 set(gca,'FontSize',5)
%                 set(gca,'FontName','Arial')
%                 % box off
% 
% 
%             end
% 
%         end
% 
%         save_file = fullfile(figPath, [refStr ' waveform global view.jpeg']);
%         print(fig_use, save_file,'-djpeg',resolution_use);
% 
%         close all
% 
%         %% tar
% 
%         hcomb = figure(102);
% 
%         res_type = 'QHD';
%         pos = [1720 2562 3*560 2*420]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
%         arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
% 
%         tiledlayout(3,size(channelSet,2),'Padding','none','TileSpacing','compact')
% 
%         for loopSets = 1:size(channelSet,2)
% 
%             tempTrough(loopSets) = max(abs(min(tarWaveforms(channelSet{loopSets},:)')));
%             tempPeak(loopSets)   = max(abs(max(tarWaveforms(channelSet{loopSets},:)')));
% 
%         end
% 
%         tempTrough = -max(tempTrough) - 0.1;
%         tempPeak   =  max(tempPeak)   + 0.1;
% 
%         for loopSets = 1:size(channelSet,2)
% 
%             partnerStrList = {};
%             if sum(loopSets == tarPartnersShank) ~= 0
%                 tempIdx = find(loopSets == tarPartnersShank);
%                 for loopTemp = 1:length(tempIdx)
% 
%                     partnerStrList = [partnerStrList; [num2str(tarPartnersID(tempIdx(loopTemp))) ...
%                                                        tarPartnersCellType{tempIdx(loopTemp)} ' (ch: ' ... 
%                                                        num2str(tarPartnersChan(tempIdx(loopTemp))) ', dist: ' ...
%                                                        num2str(tarPartnersDistance(tempIdx(loopTemp))) ')']];
% 
%                 end
%             end
% 
%             nexttile
%             plot(tWave,tarWaveforms(channelSet{loopSets},:)');
%             if sum(loopSets == tarPartnersShank) ~= 0
%                 text(-0.9,tempTrough + 0.05,partnerStrList,...
%                      'FontSize',5,'FontName','Arial')
%             end
% 
%             if loopSets == 1
%                 title([tarStr ' shank ' num2str(loopSets)])
%             else 
%                 title([' shank ' num2str(loopSets)])
%             end
% 
%             xlim([-1,1])
%             ylim([tempTrough tempPeak])
%             set(gca,'FontSize',5)
%             set(gca,'FontName','Arial')
%             xlabel('[ms]');
%             ylabel('[mV]');
%             set(gca, 'YDir','reverse')
%             box off
% 
%             nexttile(loopSets + size(channelSet,2))
%             plot(tWave,tarWaveforms(channelSet{loopSets},:)');
%             xlim([-1,1])
%             set(gca,'FontSize',5)
%             set(gca,'FontName','Arial')
%             xlabel('[ms]');
%             ylabel('[mV]');
%             set(gca, 'YDir','reverse')
%             box off
% 
%             if sum(loopSets == tarPartnersShank) ~= 0
%                 nexttile(loopSets + 2*size(channelSet,2))
%                 tempIdx = find(loopSets == tarPartnersShank);
% 
%                 stackedplot(tR*1000,tarPartnersCCGs(:,tempIdx),'k','DisplayLabels',repmat({"Counts"},1,size(tarPartnersCCGs(:,tempIdx),2)))
%                 % stackedplot(tR*1000,tarPartnersCCGs,'k','DisplayLabels',repmat({"Counts"},1,size(tarPartnersCCGs,2)))
%                 % plot(tR*1000,tarPartnersCCGs{1}./length(tarSpikeTimes),'k')
% 
%                 % hold on
%                 % for loopTemp = 2:length(tempIdx)
%                 %     stkplt
%                 %     plot(tR*1000,tarPartnersCCGs{tempIdx(loopTemp)}./length(tarSpikeTimes),'k')
%                 % end
%                 % hold off
% 
%                 xlim([-1,1])
%                 % ylims = get(gca,'ylim');
%                 % if length(tempIdx) > 1
%                 %     set(gca,'ytick',[])
%                 % else 
%                     % ylabel('Spike Probability')
%                 % end 
%                 xlabel('[ms]')
%                 set(gca,'FontSize',5)
%                 set(gca,'FontName','Arial')
%                 % box off
% 
% 
%             end
% 
%         end
% 
%         save_file = fullfile(figPath, [tarStr ' waveform global view.jpeg']);
%         print(fig_use, save_file,'-djpeg',resolution_use);
% 
%         close all
% 
%     end
% 
% end

%% histograms to compare exquisite and millisecond synchrony

% RoySummaryStatsGJ = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ_CA1.mat');
% AGsummaryStatsGJ  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ_CA1.mat');
% NSUsummaryStatsGJ = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ_CA1.mat');
% 
% RoySummaryStatsExq = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStats_CA1.mat');
% AGsummaryStatsExq  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStats_CA1.mat');
% NSUsummaryStatsExq = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStats_CA1.mat');
% 
% BOTsummaryStatsGJ  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ_CA1.mat');
% BOTsummaryStatsExq = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStats_CA1.mat');
% 
% tiledlayout(3,2)
% 
% nexttile(1)
% histogram(100*BOTsummaryStatsGJ.fracGJperInd{1,1}{1,1},'BinWidth',1)
% hold on 
% histogram(100*BOTsummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1},'BinWidth',1)
% hold off
% xlim([0 60])
% legend('ms','exq')
% xlabel('percent')
% ylabel('count')
% title('Allen Int. CA1 ii pairs percentages')
% 
% nexttile(3)
% histogram(100*BOTsummaryStatsGJ.fracGJperInd{2,1}{1,1},'BinWidth',1)
% hold on 
% histogram(100*BOTsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1},'BinWidth',1)
% hold off
% xlim([0 60])
% xlabel('percent')
% ylabel('count')
% title('Allen Int. CA1 pp pairs percentages')
% 
% nexttile(5)
% histogram(100*BOTsummaryStatsGJ.fracGJperInd{2,1}{1,1},'BinWidth',1)
% hold on 
% histogram(100*BOTsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1},'BinWidth',1)
% hold off
% xlim([0 60])
% xlabel('percent')
% ylabel('count')
% title('Allen Int. CA1 ip or pi percentages')
% 
% nexttile(2)
% histogram(100*[     RoySummaryStatsGJ.fracGJperInd{1,1}{1,1}{1,1} ...
%                mean(AGsummaryStatsGJ.fracGJperInd{1,1}{1,1}{1,1}) ...
%                     NSUsummaryStatsGJ.fracGJperInd{1,1}{1,1}{1,1}],'BinWidth',1)
% hold on 
% histogram(100*[     RoySummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1} ...
%                mean(AGsummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1}) ...
%                     NSUsummaryStatsExq.fracExqPerInd{1,1}{1,1}{1,1}],'BinWidth',1)
% hold off
% xlim([0 60])
% xlabel('percent')
% ylabel('count')
% title('rats CA1 ii pairs percentages')
% 
% nexttile(4)
% histogram(100*[     RoySummaryStatsGJ.fracGJperInd{2,1}{1,1}{1,1} ...
%                mean(AGsummaryStatsGJ.fracGJperInd{2,1}{1,1}{1,1}) ...
%                     NSUsummaryStatsGJ.fracGJperInd{2,1}{1,1}{1,1}],'BinWidth',1)
% hold on 
% histogram(100*[     RoySummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1} ...
%                mean(AGsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1}) ...
%                     NSUsummaryStatsExq.fracExqPerInd{2,1}{1,1}{1,1}],'BinWidth',1)
% hold off
% xlim([0 60])
% xlabel('percent')
% ylabel('count')
% title('rats CA1 pp pairs percentages')
% 
% nexttile(6)
% histogram(100*[     RoySummaryStatsGJ.fracGJperInd{3,1}{1,1}{1,1} ...
%                mean(AGsummaryStatsGJ.fracGJperInd{3,1}{1,1}{1,1}) ...
%                     NSUsummaryStatsGJ.fracGJperInd{3,1}{1,1}{1,1}],'BinWidth',1)
% hold on 
% histogram(100*[     RoySummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1} ...
%                mean(AGsummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1}) ...
%                     NSUsummaryStatsExq.fracExqPerInd{3,1}{1,1}{1,1}],'BinWidth',1)
% hold off
% xlim([0 60])
% xlabel('percent')
% ylabel('count')
% title('rats CA1 ip or pi percentages')

%%

function [PC1,PC2,PC3, ...
           spikeTimeCounts] ... 
             = filterPCsforChannel(fetTable, ...
                                   cluTable, ...
                                   spikeTimes, ...
                                   bx, ...
                                   PCsets, ...
                                   chIdx, ...
                                   tarUnitStackSynch, ...
                                   adjustFactor)  
                        
 % filter PCs for channel
    fetTable = fetTable(:,[PCsets(chIdx,:),31]);

    clusList = unique(cluTable); 
    for i = 2:length(clusList)

        % all sessions spike counts
        cluIdx = find(cluTable ==  clusList(i));
        spikeTimeCounts(i) = length(cluIdx); 

        % screen only for spike times in session maze
        mazeBinary = ((spikeTimes > bx.behavior.RoyMaze1.datFrame(2,1)) & (spikeTimes < bx.behavior.RoyMaze1.datFrame(2,2)));

        %

        cluIdxMaze = find((cluTable ==  clusList(i)) & (mazeBinary));

        PC1temp = fetTable(cluIdxMaze,1); PC1{i} = PC1temp{:,:};
        PC2temp = fetTable(cluIdxMaze,2); PC2{i} = PC2temp{:,:};
        PC3temp = fetTable(cluIdxMaze,3); PC3{i} = PC3temp{:,:};

        %

%         if length(cluIdx) == unitSpikeCountAllSess

            spikeTimesUnit        = spikeTimes(cluIdxMaze);
            % spikeTimesUnitSynch   = tarUnitStackSynch{i-1}+bx.behavior.RoyMaze1.datFrame(2,1) - adjustFactor;
            % spikeTimesUnitNoSynch = setdiff(spikeTimesUnit,spikeTimesUnitSynch);
            % 
            % [~,~,cluIdxSynch]     = intersect(spikeTimesUnitSynch,spikeTimes);
            % [~,~,cluIdxNoSynch]   = intersect(spikeTimesUnitNoSynch,spikeTimes);
            % 
            % PC1temp = fetTable(cluIdxSynch,1);   PC1synch{i} = PC1temp{:,:};
            % PC2temp = fetTable(cluIdxSynch,2);   PC2synch{i} = PC2temp{:,:};
            % PC3temp = fetTable(cluIdxSynch,3);   PC3synch{i} = PC3temp{:,:};
            % 
            % PC1temp = fetTable(cluIdxNoSynch,1); PC1noSynch{i} = PC1temp{:,:};
            % PC2temp = fetTable(cluIdxNoSynch,2); PC2noSynch{i} = PC2temp{:,:};
            % PC3temp = fetTable(cluIdxNoSynch,3); PC3noSynch{i} = PC3temp{:,:};

%         end

    end
 
end

