% datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
datapath = '/media/nasko/WD_BLACK31/BOTtemp/';

connType = 'GJ';

animalsList = [715093703
               719161530
               721123822
               732592105
               737581020
               739448407
               742951821
               743475441
               744228101
               746083955
               750332458
               750749662
               751348571
               754312389
               754829445
               755434585
               756029989
               757216464
               757970808
               758798717
               759883607
               760345702
               760693773
               761418226
               762120172
               762602078
               763673393
               766640955
               767871931
               768515987
               771160300
               771990200
               773418906
               774875821
               778240327
               778998620
               779839471
               781842082
               786091066
               787025148
               789848216
               791319847
               793224716
               794812542
               797828357
               798911424
               799864342
               816200189
               819186360
               819701982
               821695405
               829720705
               831882777
               835479236
               839068429
               839557629
               840012044
               847657808];

noJitter              = [];
noJitterII            = [];
noJitterPP            = [];
noJitterIPandPI       = [];
noDoubleJitter        = [];
noDoubleJitterII      = [];
noDoubleJitterPP      = [];
noDoubleJitterIPandPI = [];
           
for loopAnimals = 1:58
    
    display(['animal: '  num2str(animalsList(loopAnimals))])

    %%
    if strcmp(connType,'exq')
        load([datapath 'groupStatsMouse' num2str(animalsList(loopAnimals))])
    elseif strcmp(connType,'GJ')
        load([datapath 'groupStatsMousePutativeGJ' num2str(animalsList(loopAnimals))])
    end
    
    pairGroupStatTable = pairGroupStatsTable;
    clear pairGroupStatsTable

     % remove duplicates
    [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
    pairGroupStatTable = pairGroupStatTable(idx,:);
    
    noJitter(loopAnimals) = size(pairGroupStatTable,1);
    noJitterII(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & ...
                                                      strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow'),1),1);  
    noJitterPP(loopAnimals) = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                      (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
    noJitterPI              = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                       strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow'),1),1);  
    noJitterIP              = size(pairGroupStatTable( strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')  & ...
                                                      (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
    noJitterIPandPI(loopAnimals) = noJitterPI + noJitterIP;
    
    %%
    if strcmp(connType,'exq')
        load([datapath 'groupStatsMouse' num2str(animalsList(loopAnimals)) 'NULL'])
    elseif strcmp(connType,'GJ')
        load([datapath 'groupStatsMousePutativeGJ' num2str(animalsList(loopAnimals)) 'NULL'])
    end
    
    pairGroupStatTable = pairGroupStatsTable;
    clear pairGroupStatsTable
    
    if isempty(pairGroupStatTable)
        noDoubleJitter(loopAnimals)       = 0;
        noDoubleJitterII(loopAnimals)     = 0; 
        noDoubleJitterPP(loopAnimals)     = 0;
        noDoubleJitterIPandPI(loopAnimals) = 0;
    else
         % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
    
        noDoubleJitter(loopAnimals) = size(pairGroupStatTable,1);
        noDoubleJitterII(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & ...
                                                      strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow'),1),1);  
        noDoubleJitterPP(loopAnimals) = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                          (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
        noDoubleJitterPI              = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                           strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow'),1),1);  
        noDoubleJitterIP              = size(pairGroupStatTable( strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')  & ...
                                                          (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
        noDoubleJitterIPandPI(loopAnimals) = noDoubleJitterPI + noDoubleJitterIP;
    end
    
end

figure
histogram(noJitter,'BinWidth',10)
hold on
histogram(noDoubleJitter,'BinWidth',10)
hold off
title('all')

figure
histogram(noJitterII,'BinWidth',10)
hold on
histogram(noDoubleJitterII,'BinWidth',10)
hold off
title('ii')

figure
histogram(noJitterPP,'BinWidth',10)
hold on
histogram(noDoubleJitterPP,'BinWidth',10)
hold off
title('pp')

figure
histogram(noJitterIPandPI,'BinWidth',10)
hold on
histogram(noDoubleJitterIPandPI,'BinWidth',10)
hold off
title('ip and pi')

sum(noJitter)             
sum(noJitterII)          
sum(noJitterPP)           
sum(noJitterIPandPI)      
sum(noDoubleJitter)        
sum(noDoubleJitterII)      
sum(noDoubleJitterPP)     
sum(noDoubleJitterIPandPI) 

% [~,p,~,~] = ttest2(noJitterII,noDoubleJitterII)
% [~,p,~,~] = ttest2(noJitterPP,noDoubleJitterPP)
% [~,p,~,~] = ttest2(noJitterIPandPI,noDoubleJitterIPandPI)

[p,h,stats] = signtest(noJitterII,noDoubleJitterII)
[p,h,stats] = signtest(noJitterPP,noDoubleJitterPP)
[p,h,stats] = signtest(noJitterIPandPI,noDoubleJitterIPandPI)

%%

noJitter              = [];
noJitterII            = [];
noJitterPP            = [];
noJitterIPandPI       = [];
noDoubleJitter        = [];
noDoubleJitterII      = [];
noDoubleJitterPP      = [];
noDoubleJitterIPandPI = [];

% single jitter
if strcmp(connType,'exq')
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoy.mat')
    ratsJitter{1} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/groupStatsRatAGday1.mat')
    ratAG1jitter  = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/groupStatsRatAGday2.mat')
    ratAG2jitter  = pairGroupStatTable;
    ratsJitter{2} = [ratAG1jitter; ratAG2jitter];
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/groupStatsRatN.mat');
    ratsJitter{3} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/groupStatsRatS.mat');
    ratsJitter{4} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/groupStatsRatU.mat');
    ratsJitter{5} = pairGroupStatTable;
elseif strcmp(connType,'GJ')
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoyPutativeGJ.mat')
    ratsJitter{1} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/groupStatsRatAGday1putativeGJ.mat')
    ratAG1jitter  = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/groupStatsRatAGday2putativeGJ.mat')
    ratAG2jitter  = pairGroupStatTable;
    ratsJitter{2} = [ratAG1jitter; ratAG2jitter];
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/groupStatsRatNputativeGJ.mat');
    ratsJitter{3} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/groupStatsRatSputativeGJ.mat');
    ratsJitter{4} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/groupStatsRatUputativeGJ.mat');
    ratsJitter{5} = pairGroupStatTable;
end

% double jitter
if strcmp(connType,'exq')
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoyNULL.mat')
    ratsNull{1} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/nullGroupStatsRatAGday1.mat');
    ratAG1null  = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/nullGroupStatsRatAGday2.mat');
    ratAG2null  = pairGroupStatTable;
    ratsNull{2} = [ratAG1null; ratAG2null];
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/nullGroupStatsRatN.mat');
    ratsNull{3} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/nullGroupStatsRatS.mat');
    ratsNull{4} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/nullGroupStatsRatU.mat');
    ratsNull{5} = pairGroupStatTable;
elseif strcmp(connType,'GJ')
    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/groupStatsRatRoyNULLputativeGJ.mat')
    ratsNull{1} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/nullGroupStatsRatAGday1putativeGJ.mat');
    ratAG1null  = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/nullGroupStatsRatAGday2putativeGJ.mat');
    ratAG2null  = pairGroupStatTable;
    ratsNull{2} = [ratAG1null; ratAG2null];
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/nullGroupStatsRatNputativeGJ.mat');
    ratsNull{3} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/nullGroupStatsRatSputativeGJ.mat');
    ratsNull{4} = pairGroupStatTable;
    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/nullGroupStatsRatUputativeGJ.mat');
    ratsNull{5} = pairGroupStatTable;
end

for loopAnimals = 1:5
    
    display(['animal: '  num2str(loopAnimals)])

    %%
    
    pairGroupStatTable = ratsJitter{loopAnimals};
    
    % remove duplicates
    [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
    pairGroupStatTable = pairGroupStatTable(idx,:);
    
    noJitter(loopAnimals) = size(pairGroupStatTable,1);
    if loopAnimals == 2
        noJitterPP(loopAnimals) = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                          (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
        noJitterII(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & ...
                                                          strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow'),1),1);  
        noJitterPI              = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                           strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow'),1),1);  
        noJitterIP              = size(pairGroupStatTable( strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')  & ...
                                                          (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
        noJitterIPandPI(loopAnimals) = noJitterPI + noJitterIP;
    else 
        noJitterPP(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'p') & ...
                                                          strcmp(pairGroupStatTable.tarCellExplorerType,'p'),1),1);  
        noJitterII(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i') & ...
                                                          strcmp(pairGroupStatTable.tarCellExplorerType,'i'),1),1);  
        noJitterPI              = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'p') & ...
                                                          strcmp(pairGroupStatTable.tarCellExplorerType,'i'),1),1);  
        noJitterIP              = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i') & ...
                                                          strcmp(pairGroupStatTable.tarCellExplorerType,'p'),1),1);  
        noJitterIPandPI(loopAnimals) = noJitterPI + noJitterIP;
    end
    
    %%
    
    pairGroupStatTable = ratsNull{loopAnimals};
    
    if isempty(pairGroupStatTable)
        noDoubleJitter(loopAnimals)        = 0;
        noDoubleJitterII(loopAnimals)      = 0; 
        noDoubleJitterPP(loopAnimals)      = 0;
        noDoubleJitterIPandPI(loopAnimals) = 0;
    else
         % remove duplicates
        [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
        pairGroupStatTable = pairGroupStatTable(idx,:);
    
        noDoubleJitter(loopAnimals) = size(pairGroupStatTable,1);
        if loopAnimals == 2
            noDoubleJitterPP(loopAnimals) = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                                    (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);
            noDoubleJitterII(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & ...
                                                                    strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow'),1),1);  
            noDoubleJitterPI              = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
                                                                     strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow'),1),1);  
            noDoubleJitterIP              = size(pairGroupStatTable( strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')  & ...
                                                                    (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide')),1),1);  
            noDoubleJitterIPandPI(loopAnimals) = noDoubleJitterPI + noDoubleJitterIP;
        else
            noDoubleJitterPP(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'p') & ...
                                                                    strcmp(pairGroupStatTable.tarCellExplorerType,'p'),1),1);   
            noDoubleJitterII(loopAnimals) = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i') & ...
                                                                    strcmp(pairGroupStatTable.tarCellExplorerType,'i'),1),1);   
            noDoubleJitterPI              = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'p') & ...
                                                                    strcmp(pairGroupStatTable.tarCellExplorerType,'i'),1),1);  
            noDoubleJitterIP              = size(pairGroupStatTable(strcmp(pairGroupStatTable.refCellExplorerType,'i') & ...
                                                                    strcmp(pairGroupStatTable.tarCellExplorerType,'p'),1),1); 
            noDoubleJitterIPandPI(loopAnimals) = noDoubleJitterPI + noDoubleJitterIP;
        end
    end
    
end

figure
histogram(noJitter,'BinWidth',10)
hold on
histogram(noDoubleJitter,'BinWidth',10)
hold off
title('all')

figure
histogram(noJitterII,'BinWidth',10)
hold on
histogram(noDoubleJitterII,'BinWidth',10)
hold off
title('ii')

figure
histogram(noJitterPP,'BinWidth',10)
hold on
histogram(noDoubleJitterPP,'BinWidth',10)
hold off
title('pp')

figure
histogram(noJitterIPandPI,'BinWidth',10)
hold on
histogram(noDoubleJitterIPandPI,'BinWidth',10)
hold off
title('ip and pi')

sum(noJitter)             
sum(noJitterII)          
sum(noJitterPP)           
sum(noJitterIPandPI)      
sum(noDoubleJitter)        
sum(noDoubleJitterII)      
sum(noDoubleJitterPP)     
sum(noDoubleJitterIPandPI) 
% [~,p,~,~] = ttest2(noJitterII,noDoubleJitterII)
% [~,p,~,~] = ttest2(noJitterPP,noDoubleJitterPP)
% [~,p,~,~] = ttest2(noJitterIPandPI,noDoubleJitterIPandPI)

[p,h,stats] = signrank(noJitterII,noDoubleJitterII)
[p,h,stats] = signrank(noJitterPP,noDoubleJitterPP)
[p,h,stats] = signrank(noJitterIPandPI,noDoubleJitterIPandPI)

%%

% for loopAnimals = 1:58
%     
%     display(['animal: '  num2str(animalsList(loopAnimals))])
% 
%     %%
%     
%     load([datapath 'AllenInstitute' num2str(animalsList(loopAnimals)) 'processedGroupStats' ])
%    
%      % remove duplicates
%     [~,idx,~] = unique(sort([pairGroupStatTable.refNeuronID pairGroupStatTable.tarNeuronID],2),'rows');
%     pairGroupStatTable = pairGroupStatTable(idx,:);
%     
%     noJitter(loopAnimals) = size(pairGroupStatTable,1);
%     noJitterII(loopAnimals) = size(pairGroupStatTable((strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & ...
%                                                        strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow')) & ...
%                                                        cell2num(pairGroupStatTable.flagExq),1),1);  
%     noJitterPP(loopAnimals) = size(pairGroupStatTable(((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
%                                                        (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide'))) & ...
%                                                         cell2num(pairGroupStatTable.flagExq),1),1);  
%     noJitterPI              = size(pairGroupStatTable(((strcmp(pairGroupStatTable.refCellExplorerType,'i-wide') | strcmp(pairGroupStatTable.tarCellExplorerType,'p')) & ...
%                                                         strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')) & ...
%                                                         cell2num(pairGroupStatTable.flagExq),1),1);   
%     noJitterIP              = size(pairGroupStatTable(( strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')  & ...
%                                                        (strcmp(pairGroupStatTable.refCellExplorerType,'p')      | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide'))) & ...
%                                                         cell2num(pairGroupStatTable.flagExq),1),1);    
%     noJitterIPandPI(loopAnimals) = noJitterPI + noJitterIP;
%   
%     
% end

function pval = randTest(column1,column2)
    for iter = 1:100000
        numNoSwaps(iter) = 58-sum(round(rand(58,1)));
    end
    
    pval = sum(0 == numNoSwaps)/100000;
end