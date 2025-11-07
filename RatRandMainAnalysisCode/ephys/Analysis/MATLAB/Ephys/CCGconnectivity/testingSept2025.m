%% exquisite 

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

datasetID   = 'AllenInstitute';
datapath    = '/media/nasko/WD_BLACK31/BOTtemp/';

pairGroupStatTableExq = table;

for loopAnimals = 1:58

    display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])

    load([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);

    pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
    pairGroupStatTableTemp = pairGroupStatTableTemp(:,[4:5,10:11]);

    pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];
    
end

load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStats.mat')
pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableTemp = pairGroupStatTableTemp(:,[2:3,12:13]);
pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat')
pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableTemp = pairGroupStatTableTemp(:,[5:6,17:18]);
pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStats.mat')
pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableTemp = pairGroupStatTableTemp(:,[5:6,17:18]);
pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStats.mat')
pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableTemp = pairGroupStatTableTemp(:,[2:3,14:15]);
pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')
pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableTemp = pairGroupStatTableTemp(:,[2:3,14:15]);
pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStats.mat')
pairGroupStatTableTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableTemp = pairGroupStatTableTemp(:,[2:3,14:15]);
pairGroupStatTableExq  = [pairGroupStatTableExq; pairGroupStatTableTemp];              

%% ms

datapath    = '/media/nasko/WD_BLACK31/BOTtemp/';

pairGroupStatTableGJ = table;

for loopAnimals = 1:58

    display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])

    load([datapath datasetID num2str(animalsList(loopAnimals)) 'putativeGJprocessedGroupStats.mat']);
    
    pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:); 
    pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[4:5,10:11]);

    pairGroupStatTableGJ     = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];
    
end

load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyputativeGJprocessedGroupStats.mat')
pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[2:3,12:13]);
pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];
    
load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1putativeGJprocessedGroupStats.mat')
pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[5:6,17:18]);
pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];
    
load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2putativeGJprocessedGroupStats.mat')
pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[5:6,17:18]);
pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSputativeGJprocessedGroupStats.mat')
pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[2:3,14:15]);
pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUputativeGJprocessedGroupStats.mat')
pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[2:3,14:15]);
pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNputativeGJprocessedGroupStats.mat')
pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
pairGroupStatTableGJtemp = pairGroupStatTableGJtemp(:,[2:3,14:15]);
pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];

for loopRows = 1:size(pairGroupStatTableGJ,1)
    if strcmp(pairGroupStatTableGJ.refCellExplorerType(loopRows),'i-wide')
        pairGroupStatTableGJ.refCellExplorerType{loopRows} = 'p';
    end
    if strcmp(pairGroupStatTableGJ.tarCellExplorerType(loopRows),'i-wide')
        pairGroupStatTableGJ.tarCellExplorerType{loopRows} = 'p';
    end
    if strcmp(pairGroupStatTableGJ.refCellExplorerType(loopRows),'i-narrow')
        pairGroupStatTableGJ.refCellExplorerType{loopRows} = 'i';
    end
    if strcmp(pairGroupStatTableGJ.tarCellExplorerType(loopRows),'i-narrow')
        pairGroupStatTableGJ.tarCellExplorerType{loopRows} = 'i';
    end
end

ii_GJ_pairs_no    = sum(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i') & strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i'))
pp_GJ_pairs_no    = sum(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p') & strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'))

pi_ip_GJ_pairs_no = sum((strcmp(pairGroupStatTableGJ.refCellExplorerType,'p') & strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i')) | ... 
                            (strcmp(pairGroupStatTableGJ.refCellExplorerType,'i') & strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p')))

%%

(size(intersect(pairGroupStatTableExq,pairGroupStatTableGJ),1)/size(pairGroupStatTableExq,1))*100
(size(intersect(pairGroupStatTableExq,pairGroupStatTableGJ),1)/size(pairGroupStatTableGJ,1))*100