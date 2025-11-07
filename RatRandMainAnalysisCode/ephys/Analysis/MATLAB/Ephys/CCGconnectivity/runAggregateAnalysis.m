% for loopAnimals = 1:58
%     
%     aggregateAnalysisV2('AllenInstitute',loopAnimals,'forward');
%     
% end

clear

anaType = 'exq'; % exq or gj

C1anaOnly = false;

datasetID = 'RatsNSU';

if strcmp(datasetID,'RatRoy')
    datapath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';
    brainRegions   = "CA1";
    
    animalsList    = "RatRoy";   
elseif strcmp(datasetID,'RatAG')
    datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    if C1anaOnly
        brainRegions   = "CA1";
    else
        brainRegions   = ["CA1"
                          "CA3"];
    end
    
    subfolderList  = ["AG_2019-12-23_NSD/"
                      "AG_2019-12-27_NSD/"];
    animalsList    = ["RatAGday1"
                      "RatAGday2"];              

elseif strcmp(datasetID,'RatsNSU')
    datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    brainRegions   = "CA1";

    animalsList = ["RatN"
                   "RatS"
                   "RatU"];

elseif strcmp(datasetID,'Steinmetz')
    datapath = '/home/nasko/CUNY_Work_NPUltraWaveforms/data/';

    brainRegions   = ["ACB"
                      "AId6a"
                      "AIv5"
                      "AIv6a"
                      "AON"
                      "CA1"
                      "CA2"
                      "CA3"
                      "CP"
                      "DG-mo"
                      "DG-po"
                      "DG-sg"
                      "EPd"
                      "EW"
                      "FRP6a"
                      "LD"
                      "LGd-co"
                      "LGd-ip"
                      "LGd-sh"
                      "LP"
                      "MA3"
                      "MB"
                      "MOp6a"
                      "MOp6b"
                      "MOs2/3"
                      "MOs5"
                      "MOs6a"
                      "MOs6b"
                      "MRN"
                      "ND"
                      "NPC"
                      "OLF"
                      "OP"
                      "ORBl1"
                      "ORBl2/3"
                      "ORBl5"
                      "ORBl6a"
                      "ORBl6b"
                      "PAG"
                      "PIR"
                      "PO"
                      "PPT"
                      "ProS"
                      "RL"
                      "RN"
                      "RT"
                      "SCig"
                      "SCiw"
                      "SI"
                      "SNr"
                      "STR"
                      "SUB"
                      "VAL"
                      "VISa4"
                      "VISa5"
                      "VISa6a"
                      "VISa6b"
                      "VISam2/3"
                      "VISam4"
                      "VISam5"
                      "VISp4"
                      "VISp5"
                      "VISp6a"
                      "VISpm2/3"
                      "VISpm4"
                      "VM"
                      "VPL"
                      "VPM"
                      "alv"
                      "bsc"
                      "ccs"
                      "cing"
                      "dhc"
                      "em"
                      "fa"
                      "fp"
                      "or"
                      "root"
                      "scwm"];
    animalsList = "ultra";

elseif strcmp(datasetID,'AllenInstitute')
%     datapath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
    datapath = '/media/nasko/WD_BLACK31/BOTtemp/';

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

    optoTaggedList = [PvalbList 
                      VipList 
                      SstList];
    
    if C1anaOnly
        brainRegions   = "CA1";
    else
        brainRegions   = ["APN"   
                          "CA1"   
                          "CA3"   
                          "DG"    
                          "Eth"   
                          "HPF"   
                          "IntG"  
                          "LGd"   
                          "LGv"   
                          "LP"    
                          "MB"    
                          "MGd"   
                          "MGm"   
                          "MGv"   
                          "MRN"   
                          "NOT"   
                          "PO"    
                          "POL"   
                          "ProS"  
                          "SCig"  
                          "SGN"   
                          "SUB"   
                          "TH"    
                          "VIS"   
                          "VISal" 
                          "VISam" 
                          "VISl"  
                          "VISli" 
                          "VISmma"
                          "VISmmp"
                          "VISp"  
                          "VISpm" 
                          "VISrl" 
                          "VL"    
                          "VPM"   
                          "grey"  ];
    end
 end

screenPairType = ["ii","pp","pi","ip"];

if strcmp(anaType,'exq')
    noPair  = cell(length(screenPairType),1);
    noExq   = cell(length(screenPairType),1);
    fracExq = cell(length(screenPairType),1);
    
    noPairPerInd  = cell(length(screenPairType),1);
    noExqPerInd   = cell(length(screenPairType),1);
    fracExqPerInd = cell(length(screenPairType),1);
elseif strcmp(anaType,'gj')
    noPair  = cell(length(screenPairType),1);
    noGJ    = cell(length(screenPairType),1);
    fracGJ  = cell(length(screenPairType),1);

    noPairPerInd  = cell(length(screenPairType),1);
    noGJperInd   = cell(length(screenPairType),1);
    fracGJperInd = cell(length(screenPairType),1);
end 

%% old code

% for loopPairTypes = 1:length(screenPairType)
%     
%     for loopBrainRegions = 1:length(brainRegions)
%         
%         display(['running brain region: ' brainRegions{loopBrainRegions}])
% 
%         for loopAnimals = 1:58
% 
%             display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])
% 
%             pairGroupStatTable = parLoad([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);
% 
%             % filter brain region
%             brainRegionFilt = strcmp(pairGroupStatTable.brainRegion,brainRegions{loopBrainRegions});
% 
%             pairGroupStatTable = pairGroupStatTable(brainRegionFilt,:);
% 
%             % filter pair type
%             if strcmp(screenPairType{loopPairTypes},'pp')
%                 pairTypeFilt = (strcmp(pairGroupStatTable.refCellExplorerType,'p') | strcmp(pairGroupStatTable.refCellExplorerType,'i-wide')) & ...
%                                (strcmp(pairGroupStatTable.tarCellExplorerType,'p') | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide'));
%             elseif strcmp(screenPairType{loopPairTypes},'ii')
%                 pairTypeFilt = strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow');
%             elseif strcmp(screenPairType{loopPairTypes},'ip')
%                 pairTypeFilt = strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & ...
%                                (strcmp(pairGroupStatTable.tarCellExplorerType,'p') | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide'));
%             elseif strcmp(screenPairType{loopPairTypes},'pi')
%                 pairTypeFilt = (strcmp(pairGroupStatTable.refCellExplorerType,'p') | strcmp(pairGroupStatTable.refCellExplorerType,'i-wide')) & ...
%                                strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow');
%             end
% 
%             pairGroupStatTable = pairGroupStatTable(pairTypeFilt,:);
% 
%             noPair{loopPairTypes}{loopBrainRegions}      = noPair{loopPairTypes}{loopBrainRegions}      + size(pairGroupStatTable,1)/2;
%             noExq{loopPairTypes}{loopBrainRegions}       = noExq{loopPairTypes}{loopBrainRegions}       + sum([pairGroupStatTable.flagExq{:,1}])/2;
%             noExcMonSyp{loopPairTypes}{loopBrainRegions} = noExcMonSyp{loopPairTypes}{loopBrainRegions} + sum([pairGroupStatTable.flagExcMonSyp{:,1}]);
%             noInhMonSyp{loopPairTypes}{loopBrainRegions} = noInhMonSyp{loopPairTypes}{loopBrainRegions} + sum([pairGroupStatTable.flagInhMonSyp{:,1}]);
%             noGJ{loopPairTypes}{loopBrainRegions}        = noGJ{loopPairTypes}{loopBrainRegions}        + sum([pairGroupStatTable.flagGJ{:,1}])/2;
%         end
%     end
% end 

%%
% for loopDataset = 1:1
    for loopPairTypes = 1:length(screenPairType)
        display(['running pair type: ' screenPairType{loopPairTypes}])
        
        % Preallocate temporary variables for each parallel iteration
        temp_noPair  = zeros(length(brainRegions), 1);
        temp_noExq   = zeros(length(brainRegions), 1);
        temp_noGJ    = zeros(length(brainRegions), 1);
        temp_fracExq = zeros(length(brainRegions), 1);
        temp_fracGJ  = zeros(length(brainRegions), 1);

        for loopBrainRegions = 1:length(brainRegions)
            display(['running brain region: ' brainRegions{loopBrainRegions}])

            % Initialize temporary variables for each brain region
            temp_noPair(loopBrainRegions)  = 0;
            temp_fracExq(loopBrainRegions) = 0;
            temp_fracGJ(loopBrainRegions)  = 0;
            temp_fracExq(loopBrainRegions) = 0;
            temp_fracGJ(loopBrainRegions)  = 0;
            
            temp_noPairPerInd  = cell(length(brainRegions),1);
            temp_noExqPerInd   = cell(length(brainRegions),1);
            temp_fracExqPerInd = cell(length(brainRegions),1);
            temp_noGJperInd    = cell(length(brainRegions),1);
            temp_fracGJperInd  = cell(length(brainRegions),1);

            for loopAnimals = 1:length(animalsList)
                
                if strcmp(datasetID,'RatRoy')
                    display(['running dataset: ' datasetID ', animal: '  animalsList{loopAnimals}])
                    if strcmp(anaType,'exq')
                        pairGroupStatTable = parLoad([datapath animalsList{loopAnimals} 'processedGroupStats.mat']);
                    elseif strcmp(anaType,'gj')
                        pairGroupStatTable = parLoad([datapath animalsList{loopAnimals} 'putativeGJprocessedGroupStats.mat']);
                    end 
                elseif strcmp(datasetID,'RatAG')
                    display(['running dataset: ' datasetID ', animal: '  animalsList{loopAnimals}])
                    if strcmp(anaType,'exq')
                        pairGroupStatTable = parLoad([datapath subfolderList{loopAnimals} '/' animalsList{loopAnimals} 'processedGroupStats.mat']);
                    elseif strcmp(anaType,'gj')
                        pairGroupStatTable = parLoad([datapath subfolderList{loopAnimals} '/' animalsList{loopAnimals} 'putativeGJprocessedGroupStats.mat']);
                    end 
                elseif strcmp(datasetID,'RatsNSU')
                    display(['running dataset: ' datasetID ', animal: '  animalsList{loopAnimals}])
                    if strcmp(anaType,'exq')
                        pairGroupStatTable = parLoad([datapath animalsList{loopAnimals} '/' animalsList{loopAnimals} 'processedGroupStats.mat']);
                    elseif strcmp(anaType,'gj')
                        pairGroupStatTable = parLoad([datapath animalsList{loopAnimals} '/' animalsList{loopAnimals} 'putativeGJprocessedGroupStats.mat']);
                    end 
                elseif strcmp(datasetID,'AllenInstitute')
                    display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])
                    if strcmp(anaType,'exq')
                        pairGroupStatTable = parLoad([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);
                    elseif strcmp(anaType,'gj')
                        pairGroupStatTable = parLoad([datapath datasetID num2str(animalsList(loopAnimals)) 'putativeGJprocessedGroupStats.mat']);
                    end 
                elseif strcmp(datasetID,'Steinmetz')
                    display(['running dataset: ' datasetID ', animal: '  animalsList{loopAnimals}])
                    pairGroupStatTable = parLoad([datapath  'SteinmetzprocessedGroupStats.mat']);
                end

                % filter brain region
                try
                    brainRegionFilt = strcmp(pairGroupStatTable.brainRegion,brainRegions{loopBrainRegions});
                catch
                    brainRegionFilt = strcmp(string(pairGroupStatTable.refBrainRegion),brainRegions{loopBrainRegions});
                end
                pairGroupStatTable = pairGroupStatTable(brainRegionFilt,:);

                % filter pair type
                if strcmp(screenPairType{loopPairTypes},'pp')
                    pairTypeFilt = (strcmp(pairGroupStatTable.refCellExplorerType,'p') | strcmp(pairGroupStatTable.refCellExplorerType,'i-wide')) & ...
                                   (strcmp(pairGroupStatTable.tarCellExplorerType,'p') | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide'));
                elseif strcmp(screenPairType{loopPairTypes},'ii')
                    pairTypeFilt = (strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow') & strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow')) | ... 
                                   (strcmp(pairGroupStatTable.refCellExplorerType,'i')        & strcmp(pairGroupStatTable.tarCellExplorerType,'i'));
                elseif strcmp(screenPairType{loopPairTypes},'ip')
                    pairTypeFilt = (strcmp(pairGroupStatTable.refCellExplorerType,'i') | strcmp(pairGroupStatTable.refCellExplorerType,'i-narrow')) & ...
                                   (strcmp(pairGroupStatTable.tarCellExplorerType,'p') | strcmp(pairGroupStatTable.tarCellExplorerType,'i-wide'));
                elseif strcmp(screenPairType{loopPairTypes},'pi')
                    pairTypeFilt = (strcmp(pairGroupStatTable.refCellExplorerType,'p') | strcmp(pairGroupStatTable.refCellExplorerType,'i-wide')) & ...
                                   (strcmp(pairGroupStatTable.tarCellExplorerType,'i') | strcmp(pairGroupStatTable.tarCellExplorerType,'i-narrow'));
                end

                % count opto pairs
                if strcmp(datasetID,'AllenInstitute')
                    if sum(animalsList(loopAnimals) == PvalbList) 
                        pairOptoFilt = (strcmp(pairGroupStatTable.refOptoType5X,'Pvalb') | strcmp(pairGroupStatTable.tarOptoType5X,'Pvalb')); 
                    elseif sum(animalsList(loopAnimals) == VipList) 
                        pairOptoFilt = (strcmp(pairGroupStatTable.refOptoType5X,'Vip')   | strcmp(pairGroupStatTable.tarOptoType5X,'Vip')); 
                    elseif sum(animalsList(loopAnimals) == SstList) 
                        pairOptoFilt = (strcmp(pairGroupStatTable.refOptoType5X,'Sst')   | strcmp(pairGroupStatTable.tarOptoType5X,'Sst')); 
                    end

                    
                end

                pairGroupStatTable = pairGroupStatTable(pairTypeFilt,:);
                
                
                if strcmp(anaType,'exq')
                    
                    if strcmp(screenPairType(loopPairTypes),"ii")
                        normExq = [pairGroupStatTable.nPairCompII([pairGroupStatTable.flagExq{:,1}])];
                    elseif strcmp(screenPairType(loopPairTypes),"pp")
                        normExq = [pairGroupStatTable.nPairCompPP([pairGroupStatTable.flagExq{:,1}])];
                    elseif strcmp(screenPairType(loopPairTypes),"pi") || strcmp(screenPairType(loopPairTypes),"ip")
                        normExq = [pairGroupStatTable.nPairCompIPorPI([pairGroupStatTable.flagExq{:,1}])];
                    end
                    
                    if strcmp(datasetID,'AllenInstitute')
                        uniqueNoProbesNorm = size(unique({pairGroupStatTable.probeID{[pairGroupStatTable.flagExq{:,1}]}}),2);
                    elseif strcmp(datasetID,'RatRoy')
                        uniqueNoProbesNorm = 1;
                    else
                        uniqueNoProbesNorm = size(unique([pairGroupStatTable.refProbe([pairGroupStatTable.flagExq{:,1}])]),1);
                    end

                    % Accumulate results for each brain region
                    temp_noPair(loopBrainRegions)  = temp_noPair(loopBrainRegions)  + size(pairGroupStatTable,1)/2;
                    temp_noExq(loopBrainRegions)   = temp_noExq(loopBrainRegions)   + sum([pairGroupStatTable.flagExq{[pairGroupStatTable.flagExq{:,1}],1}])/2;
                    temp_fracExq(loopBrainRegions) = temp_fracExq(loopBrainRegions) + (sum([pairGroupStatTable.flagExq{[pairGroupStatTable.flagExq{:,1}],1}]./normExq')/uniqueNoProbesNorm)/2;
                    
                    temp_noPairPerInd{loopBrainRegions}  = [temp_noPairPerInd{loopBrainRegions}  size(pairGroupStatTable,1)/2];
                    temp_noExqPerInd{loopBrainRegions}   = [temp_noExqPerInd{loopBrainRegions}   sum([pairGroupStatTable.flagExq{[pairGroupStatTable.flagExq{:,1}],1}])/2];
                    temp_fracExqPerInd{loopBrainRegions} = [temp_fracExqPerInd{loopBrainRegions} (sum([pairGroupStatTable.flagExq{[pairGroupStatTable.flagExq{:,1}],1}]./normExq')/uniqueNoProbesNorm)/2];

                elseif strcmp(anaType,'gj')
                    
                    if strcmp(screenPairType(loopPairTypes),"ii")
                        normGJ  = [pairGroupStatTable.nPairCompII([pairGroupStatTable.flagGJ{:,1}])];
                    elseif strcmp(screenPairType(loopPairTypes),"pp")
                        normGJ  = [pairGroupStatTable.nPairCompPP([pairGroupStatTable.flagGJ{:,1}])];
                    elseif strcmp(screenPairType(loopPairTypes),"pi") || strcmp(screenPairType(loopPairTypes),"ip")
                        normGJ  = [pairGroupStatTable.nPairCompIPorPI([pairGroupStatTable.flagGJ{:,1}])];
                    end

                    if strcmp(datasetID,'AllenInstitute')
                        uniqueNoProbesNorm = size(unique({pairGroupStatTable.probeID{[pairGroupStatTable.flagGJ{:,1}]}}),2);
                    elseif strcmp(datasetID,'RatRoy')
                        uniqueNoProbesNorm = 1;
                    else
                        uniqueNoProbesNorm = size(unique([pairGroupStatTable.refProbe([pairGroupStatTable.flagGJ{:,1}])]),1);
                    end
                    
                    % Accumulate results for each brain region
                    temp_noPair(loopBrainRegions)  = temp_noPair(loopBrainRegions)  + size(pairGroupStatTable,1)/2;
                    temp_noGJ(loopBrainRegions)    = temp_noGJ(loopBrainRegions)    + sum([pairGroupStatTable.flagGJ{ [pairGroupStatTable.flagGJ{:,1}],1}])/2;
                    temp_fracGJ(loopBrainRegions)  = temp_fracGJ(loopBrainRegions)  + (sum([pairGroupStatTable.flagGJ{ [pairGroupStatTable.flagGJ{:,1}],1}]./normGJ')/uniqueNoProbesNorm)/2;
                    
                    temp_noPairPerInd{loopBrainRegions}  = [temp_noPairPerInd{loopBrainRegions}  size(pairGroupStatTable,1)/2];
                    temp_noGJperInd{loopBrainRegions}    = [temp_noGJperInd{loopBrainRegions}    sum([pairGroupStatTable.flagGJ{[pairGroupStatTable.flagGJ{:,1}],1}])/2];
                    temp_fracGJperInd{loopBrainRegions}  = [temp_fracGJperInd{loopBrainRegions}  (sum([pairGroupStatTable.flagGJ{[pairGroupStatTable.flagGJ{:,1}],1}]./normGJ')/uniqueNoProbesNorm)/2];

                end 
                
            end
        end
        
        if strcmp(anaType,'exq')
            % Assign the results back to the cell arrays
            noPair{loopPairTypes}  = temp_noPair;
            noExq{loopPairTypes}   = temp_noExq;
            fracExq{loopPairTypes} = temp_fracExq;

            noPairPerInd{loopPairTypes}{loopBrainRegions}  = temp_noPairPerInd;
            noExqPerInd{loopPairTypes}{loopBrainRegions}   = temp_noExqPerInd;
            fracExqPerInd{loopPairTypes}{loopBrainRegions} = temp_fracExqPerInd;
        elseif strcmp(anaType,'gj')
            % Assign the results back to the cell arrays
            noPair{loopPairTypes}  = temp_noPair;
            noGJ{loopPairTypes}    = temp_noGJ;
            fracGJ{loopPairTypes}  = temp_fracGJ;

            noPanoPairPerIndir{loopPairTypes}{loopBrainRegions}  = temp_noPairPerInd;
            noGJperInd{loopPairTypes}{loopBrainRegions}          = temp_noGJperInd;
            fracGJperInd{loopPairTypes}{loopBrainRegions}        = temp_fracGJperInd;
        end 
        
        
    end
% end
    
for loopPairTypes = 1:length(screenPairType)
    
    if strcmp(anaType,'exq')
        noPairTotal{loopPairTypes}  = sum(noPair{loopPairTypes});
        noExqTotal{loopPairTypes}   = sum(noExq{loopPairTypes});
        fracExqTotal{loopPairTypes} = sum(fracExq{loopPairTypes});
    elseif strcmp(anaType,'gj')
        noPairTotal{loopPairTypes}  = sum(noPair{loopPairTypes});
        noGJTotal{loopPairTypes}    = sum(noGJ{loopPairTypes});
        fracGJTotal{loopPairTypes}  = sum(fracGJ{loopPairTypes});
    end 
    
end

if strcmp(anaType,'exq')
    figure
    bar(categorical(brainRegions),100*[fracExq{1,1}/size(animalsList,1) fracExq{2,1}/size(animalsList,1) fracExq{4,1}/size(animalsList,1)])
    legend(screenPairType{1}, screenPairType{2}, [screenPairType{3} ' or ' screenPairType{4}])
elseif strcmp(anaType,'gj')
    figure
    bar(categorical(brainRegions),100*[fracGJ{1,1}/size(animalsList,1) fracGJ{2,1}/size(animalsList,1) fracGJ{4,1}/size(animalsList,1)])
    legend(screenPairType{1}, screenPairType{2}, [screenPairType{3} ' or ' screenPairType{4}])
end 

if strcmp(anaType,'gj')
    if C1anaOnly
        if strcmp(datasetID,'RatRoy')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ_CA1.mat')
        elseif strcmp(datasetID,'RatAG')            
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ_CA1.mat')
        elseif strcmp(datasetID,'RatsNSU')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ_CA1.mat')
        elseif strcmp(datasetID,'Steinmetz')

        elseif strcmp(datasetID,'AllenInstitute')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ_CA1.mat')
        end
    else
        if strcmp(datasetID,'RatRoy')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ.mat')
        elseif strcmp(datasetID,'RatAG')            
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ.mat')
        elseif strcmp(datasetID,'RatsNSU')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ.mat')
        elseif strcmp(datasetID,'Steinmetz')

        elseif strcmp(datasetID,'AllenInstitute')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ.mat')
        end
    end
end

if strcmp(anaType,'exq')
    if C1anaOnly
        if strcmp(datasetID,'RatRoy')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStats_CA1.mat')
        elseif strcmp(datasetID,'RatAG')            
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStats_CA1.mat')
        elseif strcmp(datasetID,'RatsNSU')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStats_CA1.mat')
        elseif strcmp(datasetID,'Steinmetz')

        elseif strcmp(datasetID,'AllenInstitute')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStats_CA1.mat')
        end
    else
        if strcmp(datasetID,'RatRoy')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStats.mat')
        elseif strcmp(datasetID,'RatAG')            
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStats.mat')
        elseif strcmp(datasetID,'RatsNSU')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStats.mat')
        elseif strcmp(datasetID,'Steinmetz')

        elseif strcmp(datasetID,'AllenInstitute')
            save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStats.mat')
        end
    end
end

% tiledlayout(3,1)
% nexttile
% bar(categorical(brainRegions),noExq{1,1})
% title(screenPairType{1})
% nexttile
% bar(categorical(brainRegions),noExq{2,1})
% title(screenPairType{2})
% nexttile
% bar(categorical(brainRegions),noExq{4,1})
% title([screenPairType{3} ' or ' screenPairType{4}])


function pairGroupStatTable = parLoad(file)

    load(file);
    
end


