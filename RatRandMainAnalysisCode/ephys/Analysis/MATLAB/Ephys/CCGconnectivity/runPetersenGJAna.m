function runPetersenGJAna
    
    filtExqFlag  = false;
    alphaGraph1  = 0.025;
    alphaGraph2  = 0.05;
    
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
%     datapath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
    
    datapath    = '/media/nasko/WD_BLACK31/BOTtemp/';

    
    %% Allen gap junction

    pairGroupStatTableGJ = table;

    for loopAnimals = 1:58

        display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])

        load([datapath datasetID num2str(animalsList(loopAnimals)) 'putativeGJprocessedGroupStats.mat']);
        
        pairGroupStatTableGJtemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
        
        if filtExqFlag
           
            load([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);
            
            pairGroupStatTableExqTemp = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);
            
%             [~,idx,~] = intersect([pairGroupStatTableGJtemp.refNeuronID, pairGroupStatTableGJtemp.tarNeuronID], ... 
%                                   [pairGroupStatTableExqTemp.refNeuronID,pairGroupStatTableExqTemp.tarNeuronID], ...
%                                   'rows');
%                           
%             pairGroupStatTableGJtemp(idx,:) = [];
            
        end
        
        pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTableGJtemp];
        
    end

%     fwhm = findFwhm([pairGroupStatTableGJ.ccgHD{:,1}]);
%     
%     pairGroupStatTableGJ(fwhm < 6,:) = [];
    
    %% rats 

    load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyputativeGJprocessedGroupStats.mat')

    royGJTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
    
    for loopPairs = 1:size(royGJTable,1)
        CCGzscored = zscore(royGJTable.ccgHD{loopPairs});
        if (CCGzscored(106) > 3) && (~(CCGzscored(105) > 1) && ~(CCGzscored(107) > 1))
            royGJTable.sharpZeroLag(loopPairs) = 1;
        else
            royGJTable.sharpZeroLag(loopPairs) = 0;
        end
    end
    
    royGJTableSharpZeroLag   = royGJTable(royGJTable.sharpZeroLag == 1,:);
    royGJTableSharpNoZeroLag = royGJTable(royGJTable.sharpZeroLag == 0,:);
    
    if filtExqFlag
        load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStats.mat')

        royExqTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);

        [~,idx,~] = intersect([royGJTable.refNeuronID,royGJTable.tarNeuronID], ... 
                              [royExqTable.refNeuronID,royExqTable.tarNeuronID], ...
                              'rows');
                          
        royGJTable(idx,:) = [];
    end 
    
%     if ~isempty(royGJTable)
%         fwhm = findFwhm([royGJTable.ccgHD{:,1}]);
% 
%         royGJTable(fwhm < 6,:) = [];
%     end
    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1putativeGJprocessedGroupStats.mat')

    AG1GJTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
    
    for loopPairs = 1:size(AG1GJTable,1)
        CCGzscored = zscore(AG1GJTable.ccgHD{loopPairs});
        if CCGzscored(106) > 3 && (~(CCGzscored(105) > 3) && ~(CCGzscored(107) > 3))
            AG1GJTable.sharpZeroLag(loopPairs) = 1;
        else
            AG1GJTable.sharpZeroLag(loopPairs) = 0;
        end
    end
    
    AG1GJTableSharpZeroLag   = AG1GJTable(AG1GJTable.sharpZeroLag == 1,:);
    AG1GJTableSharpNoZeroLag = AG1GJTable(AG1GJTable.sharpZeroLag == 0,:);
    
    if filtExqFlag
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat')
        
        AG1ExqTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);
       
        [~,idx,~] = intersect([AG1GJTable.refNeuronID,AG1GJTable.tarNeuronID], ... 
                              [AG1ExqTable.refNeuronID,AG1ExqTable.tarNeuronID], ...
                              'rows');
        AG1GJTable(idx,:) = [];
    end
    
%     if ~isempty(AG1GJTable)
%         fwhm = findFwhm([AG1GJTable.ccgHD{:,1}]);
% 
%         AG1GJTable(fwhm < 6,:) = [];
%     end
    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2putativeGJprocessedGroupStats.mat')

    AG2GJTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
    
    for loopPairs = 1:size(AG2GJTable,1)
        CCGzscored = zscore(AG2GJTable.ccgHD{loopPairs});
        if CCGzscored(106) > 3 && (~(CCGzscored(105) > 3) && ~(CCGzscored(107) > 3))
            AG2GJTable.sharpZeroLag(loopPairs) = 1;
        else
            AG2GJTable.sharpZeroLag(loopPairs) = 0;
        end
    end
    
    AG2GJTableSharpZeroLag   = AG2GJTable(AG2GJTable.sharpZeroLag == 1,:);
    AG2GJTableSharpNoZeroLag = AG2GJTable(AG2GJTable.sharpZeroLag == 0,:);
    
    if filtExqFlag
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat')
        
        AG2ExqTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);
       
        [~,idx,~] = intersect([AG2GJTable.refNeuronID,AG2GJTable.tarNeuronID], ... 
                              [AG2ExqTable.refNeuronID,AG2ExqTable.tarNeuronID], ...
                              'rows');
        AG2GJTable(idx,:) = [];
    end
    
%     if ~isempty(AG2GJTable)
%         fwhm = findFwhm([AG2GJTable.ccgHD{:,1}]);
% 
%         AG2GJTable(fwhm < 6,:) = [];
%     end
    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSputativeGJprocessedGroupStats.mat')

    ratSgjTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
    
    for loopPairs = 1:size(ratSgjTable,1)
        CCGzscored = zscore(ratSgjTable.ccgHD{loopPairs});
        if CCGzscored(106) > 3 && (~(CCGzscored(105) > 3) && ~(CCGzscored(107) > 3))
            ratSgjTable.sharpZeroLag(loopPairs) = 1;
        else
            ratSgjTable.sharpZeroLag(loopPairs) = 0;
        end
    end
    
    ratSgjTableSharpZeroLag   = ratSgjTable(ratSgjTable.sharpZeroLag == 1,:);
    ratSgjTableSharpNoZeroLag = ratSgjTable(ratSgjTable.sharpZeroLag == 0,:);
    
    if filtExqFlag
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')
        
        ratSexqTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);
        
        [~,idx,~] = intersect([ratSgjTable.refNeuronID,ratSgjTable.tarNeuronID], ... 
                              [ratSexqTable.refNeuronID,ratSexqTable.tarNeuronID], ...
                              'rows');
        ratSgjTable(idx,:) = [];
    end

%     if ~isempty(ratSgjTable)
%         fwhm = findFwhm([ratSgjTable.ccgHD{:,1}]);
% 
%         ratSgjTable(fwhm < 6,:) = [];
%     end
    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUputativeGJprocessedGroupStats.mat')

    ratUgjTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
    
    for loopPairs = 1:size(ratUgjTable,1)
        CCGzscored = zscore(ratUgjTable.ccgHD{loopPairs});
        if CCGzscored(106) > 3 && (~(CCGzscored(105) > 3) && ~(CCGzscored(107) > 3))
            ratUgjTable.sharpZeroLag(loopPairs) = 1;
        else
            ratUgjTable.sharpZeroLag(loopPairs) = 0;
        end
    end
    
    ratUgjTableSharpZeroLag   = ratUgjTable(ratUgjTable.sharpZeroLag == 1,:);
    ratUgjTableSharpNoZeroLag = ratUgjTable(ratUgjTable.sharpZeroLag == 0,:);
    
    if filtExqFlag
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')
        
        ratUexqTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);
        
        [~,idx,~] = intersect([ratUgjTable.refNeuronID,ratUgjTable.tarNeuronID], ... 
                              [ratUexqTable.refNeuronID,ratUexqTable.tarNeuronID], ...
                              'rows');
        ratUgjTable(idx,:) = [];
    end
    
%     if ~isempty(ratUgjTable)
%         fwhm = findFwhm([ratUgjTable.ccgHD{:,1}]);
% 
%         ratUgjTable(fwhm < 6,:) = [];
%     end
    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNputativeGJprocessedGroupStats.mat')

    ratNgjTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
    
    for loopPairs = 1:size(ratNgjTable,1)
        CCGzscored = zscore(ratNgjTable.ccgHD{loopPairs});
        if CCGzscored(106) > 3 && (~(CCGzscored(105) > 3) && ~(CCGzscored(107) > 3))
            ratNgjTable.sharpZeroLag(loopPairs) = 1;
        else
            ratNgjTable.sharpZeroLag(loopPairs) = 0;
        end
    end
    
    ratNgjTableSharpZeroLag   = ratNgjTable(ratNgjTable.sharpZeroLag == 1,:);
    ratNgjTableSharpNoZeroLag = ratNgjTable(ratNgjTable.sharpZeroLag == 0,:); 
        
    if filtExqFlag
        load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')
        
        ratNexqTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagExq) == 1,:);
        
        [~,idx,~] = intersect([ratNgjTable.refNeuronID,ratNgjTable.tarNeuronID], ... 
                              [ratNexqTable.refNeuronID,ratNexqTable.tarNeuronID], ...
                              'rows');
        ratNgjTable(idx,:) = [];
    end
    
%     if ~isempty(ratNgjTable)
%         fwhm = findFwhm([ratNgjTable.ccgHD{:,1}]);
% 
%         ratNgjTable(fwhm < 6,:) = [];
%     end
    %%

    figure('Renderer', 'painters', 'Position', [1720 2562 2*560*0.125 0.5*420])

    tiledlayout(1,2,'Padding','none','TileSpacing','compact')

%     nexttile
    % plot(pairGroupStatTableGJ.tHD{1,1}*1000,...
    %      mean([pairGroupStatTableGJ.ccgHD{:,1}]/sum([pairGroupStatTableGJ.ccgHD{:,1}]),2),...
    %       'k','LineWidth',4);
    % hold on 
    % for loopPairs = 1:size(pairGroupStatTableGJ,1)
    %     patchline(pairGroupStatTableGJ.tHD{1,1}*1000, ... 
    %               [pairGroupStatTableGJ.ccgHD{loopPairs,1}]/sum([pairGroupStatTableGJ.ccgHD{loopPairs,1}]), ...
    %              'linewidth',1,'edgealpha',0.1)
    % end
    
    tR = royGJTable.tHD{1,1};
    
    mouseCCGs = [[pairGroupStatTableGJ.ccgHD{:,1}]./repmat(cellfun(@length,[pairGroupStatTableGJ.refSpikeTimes]),1,211)'];
    ratCCGs   = [[ratNgjTable.ccgHD{:,1}]./repmat(cellfun(@length,[ratNgjTable.refSpikeTimes]),1,211)', ...
                 [ratSgjTable.ccgHD{:,1}]./repmat(cellfun(@length,[ratSgjTable.refSpikeTimes]),1,211)', ...
                 [ratUgjTable.ccgHD{:,1}]./repmat(cellfun(@length,[ratUgjTable.refSpikeTimes]),1,211)', ...
                 [AG1GJTable.ccgHD{:,1}]./repmat(cellfun(@length, [AG1GJTable.refSpikeTimes]),1,211)', ...
                 [AG2GJTable.ccgHD{:,1}]./repmat(cellfun(@length, [AG2GJTable.refSpikeTimes]),1,211)', ...
                 [royGJTable.ccgHD{:,1}]./repmat(cellfun(@length, [royGJTable.refSpikeTimes]),1,211)'];
             
    % mouseCCGselected = mouseCCGs(:,randperm(size(mouseCCGs,2),100));
    % ratCCGselected   = ratCCGs(:,randperm(size(ratCCGs,2),100));
    
    mouseCCGselectedidx = [1,9,30,59,80,94,111,142,144,191];
    ratCCGselected      = [2,28,51,67,88,104,140,151,162,236];

    % mouseCCGselectedidx = mouseCCGselectedidx([4,2,7,6,5,8,9,10,3,1]);

    mouseCCGselected = mouseCCGs(:,mouseCCGselectedidx);
    ratCCGselected   = ratCCGs(:,ratCCGselected);
    
    % stackedplot(tR,mouseCCGselected,'k','DisplayLabels',repmat({""},1,size(mouseCCGselected,2)))
    % set(gca,'ytick',[])
    % set(gca,'yticklabel',[])
    % stackedplot(tR*1000,refPartnersCCGs(:,tempIdx),'k','DisplayLabels',repmat({"Counts"},1,size(refPartnersCCGs(:,tempIdx),2)))

    nexttile
    plot(tR, mouseCCGselected(:,1),'k','linewidth',1)
    
    hold on
    xline(-1,'--r');
    xline(1, '--r');
    for loopPairs = 2:10 %size(mouseCCGs,2)
        stkplt
        plot(tR,mouseCCGselected(:,loopPairs),'k','linewidth',1)
    end
    stkplt
    plot(tR,zeros(length(tR),1),'Color',[1, 0, 0, 0])
    stkplt
    plot(tR,zeros(length(tR),1),'Color',[1, 0, 0, 0])
    hold off
    
    xlim([-3.5,3.5])
    xlabel('[ms]')
    
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile 
    plot(tR, ratCCGselected(:,1),'k','linewidth',1)
    
    hold on
    xline(-1,'--r');
    xline(1, '--r');
    for loopPairs = 2:10 %size(ratCCGsrand,2)
        stkplt
        plot(tR,ratCCGselected(:,loopPairs),'k','linewidth',1)
    end
    stkplt
    plot(tR,zeros(length(tR),1),'Color',[1, 0, 0, 0])
    stkplt
    plot(tR,zeros(length(tR),1),'Color',[1, 0, 0, 0])
    hold off
    xlim([-3.5,3.5])
    xlabel('[ms]')

    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    %%
    
%     nexttile
%     for loopPairs = 1:20 %size(mouseCCGs,2)
%         patchline(tR, mouseCCGrand(:,loopPairs),'linewidth',1,'edgealpha',0.5)
%     end
%     xline(-1,'--r');
%     xline(1, '--r');
%     
%     nexttile
%     for loopPairs = 1:20 %size(ratCCGs,2)
%         patchline(tR, ratCCGrand(:,loopPairs),'linewidth',1,'edgealpha',0.5)
%     end
%     xline(-1,'--r');
%     xline(1, '--r');
    
    %%
    
%     bar(pairGroupStatTableGJ.tHD{1,1},mean(mouseCCGs,2),'k','LineWidth',1);
%     
%     
%     xline(-1/30,'r')
%     xline(1/30,'r')
%      hold off
%     xlim([-2 2])
% %     ylim([0 0.01])
%     title({'Mouse aggregate','millisecond synchrony'})
%     ylabel('Spike Probability')
%     xlabel('[ms]')
%     set(gca,'FontSize',5)
%     set(gca,'FontName', 'Arial')
%     box off
%     
%     nexttile
%            
%     hold on 
%     bar(royGJTable.tHD{1,1},mean(ratCCGs,2),'k','LineWidth',1);
%            
%     xline(-1/30,'r')
%     xline(1/30,'r')
%     
%     xline(-1/30,'r')
%     xline(1/30,'r')
%      hold off
%     xlim([-2 2])
% %     ylim([0 0.01])
%     title({'Rat aggregate','millisecond synchrony'})
%     ylabel('Spike Probability')
%     xlabel('[ms]')
%     set(gca,'FontSize',5)
%     set(gca,'FontName', 'Arial')
%     box off
    
%     figure
%     plot(mean([[ratNgjTableSharpZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[ratNgjTableSharpZeroLag.refSpikeTimes]),1,211)', ...
%                [ratSgjTableSharpZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[ratSgjTableSharpZeroLag.refSpikeTimes]),1,211)', ...
%                [ratUgjTableSharpZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[ratUgjTableSharpZeroLag.refSpikeTimes]),1,211)', ...
%                [AG1GJTableSharpZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[AG1GJTableSharpZeroLag.refSpikeTimes]),1,211)', ...
%                [AG2GJTableSharpZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[AG2GJTableSharpZeroLag.refSpikeTimes]),1,211)', ...
%                [royGJTableSharpZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[royGJTableSharpZeroLag.refSpikeTimes]),1,211)'],2),...
%                'k','LineWidth',1);    
%            
%     figure 
% %     hold on
%      plot(royGJTable.tHD{1,1}, ...
%          mean([[ratNgjTableSharpNoZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[ratNgjTableSharpNoZeroLag.refSpikeTimes]),1,211)', ...
%                [ratSgjTableSharpNoZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[ratSgjTableSharpNoZeroLag.refSpikeTimes]),1,211)', ...
%                [ratUgjTableSharpNoZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[ratUgjTableSharpNoZeroLag.refSpikeTimes]),1,211)', ...
%                [AG1GJTableSharpNoZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[AG1GJTableSharpNoZeroLag.refSpikeTimes]),1,211)', ...
%                [AG2GJTableSharpNoZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[AG2GJTableSharpNoZeroLag.refSpikeTimes]),1,211)', ...
%                [royGJTableSharpNoZeroLag.ccgHD{:,1}]./repmat(cellfun(@length,[royGJTableSharpNoZeroLag.refSpikeTimes]),1,211)'],2),...
%                'k','LineWidth',1);        
% %     hold off
%            
%     for loopPairs = 1:size(pairGroupStatTableGJ,1)
%         patchline(pairGroupStatTableGJ.tHD{1,1}, ... 
%                  [pairGroupStatTableGJ.ccgHD{loopPairs,1}]/length(pairGroupStatTableGJ.refSpikeTimes{loopPairs}), ...
%                  'linewidth',1,'edgealpha',alphaGraph1)
%     end
%     for loopPairs = 1:size(AG1GJTable,1)
%         patchline(AG1GJTable.tHD{1,1}, ... 
%                  [AG1GJTable.ccgHD{loopPairs,1}]/length(AG1GJTable.refSpikeTimes{loopPairs}), ...
%                  'linewidth',2,'edgealpha',alphaGraph1)
%     end
%     for loopPairs = 1:size(AG2GJTable,1)
%         patchline(AG2GJTable.tHD{1,1}, ... 
%                  [AG2GJTable.ccgHD{loopPairs,1}]/length(AG2GJTable.refSpikeTimes{loopPairs}), ...
%                  'linewidth',2,'edgealpha',alphaGraph1)
%     end
%     for loopPairs = 1:size(ratNgjTable,1)
%         patchline(ratNgjTable.tHD{1,1}, ... 
%                  [ratNgjTable.ccgHD{loopPairs,1}]/length(ratNgjTable.refSpikeTimes{loopPairs}), ...
%                  'linewidth',1,'edgealpha',alphaGraph1)
%     end
%     for loopPairs = 1:size(ratSgjTable,1)
%         patchline(ratSgjTable.tHD{1,1}, ... 
%                  [ratSgjTable.ccgHD{loopPairs,1}]/length(ratSgjTable.refSpikeTimes{loopPairs}), ...
%                  'linewidth',1,'edgealpha',alphaGraph1)
%     end
%     for loopPairs = 1:size(ratUgjTable,1)
%         patchline(ratUgjTable.tHD{1,1}, ... 
%                  [ratUgjTable.ccgHD{loopPairs,1}]/length(ratUgjTable.refSpikeTimes{loopPairs}), ...
%                  'linewidth',1,'edgealpha',alphaGraph1)
%     end
%     for loopPairs = 1:size(royGJTable,1)
%         patchline(royGJTable.tHD{1,1}*1000, ... 
%                  [royGJTable.ccgHD{loopPairs,1}]/length(royGJTable.refSpikeTimes{loopPairs}), ...
%                  'linewidth',2,'edgealpha',0.1)
%     end
%     hold off
%     xlim([-3.5 3.5])
% %     ylim([0 0.01])
%     ylabel('Spike Probability')
%     xlabel('[ms]')
%     set(gca,'FontSize',5)
%     set(gca,'FontName', 'Arial')
%     box off
    
    close all

    hcomb = figure(102);
    
    res_type = 'QHD';
    pos = [70 230 1920*0.125 0.6*1080/3]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    hist3([cell2num([pairGroupStatTableGJ.refTroughToPeakLength;       ...
                     pairGroupStatTableGJ.tarTroughToPeakLength;       ... 
                     ratNgjTable.refTroughToPeak; ...
                     ratNgjTable.tarTroughToPeak; ...
                     ratSgjTable.refTroughToPeak; ...
                     ratSgjTable.tarTroughToPeak; ...
                     ratUgjTable.refTroughToPeak; ...
                     ratUgjTable.tarTroughToPeak; ...
                     AG1GJTable.refTroughToPeak;  ...
                     AG1GJTable.tarTroughToPeak;  ...
                     AG2GJTable.refTroughToPeak;  ...
                     AG2GJTable.tarTroughToPeak]), ...
                    [pairGroupStatTableGJ.refACGtauRise;       ...
                     pairGroupStatTableGJ.tarACGtauRise;       ...
                     ratNgjTable.refACGtauRise; ...
                     ratNgjTable.tarACGtauRise; ...
                     ratSgjTable.refACGtauRise; ...
                     ratSgjTable.tarACGtauRise; ...
                     ratUgjTable.refACGtauRise; ...
                     ratUgjTable.tarACGtauRise; ...
                     AG1GJTable.refACGtauRise;   ...
                     AG1GJTable.tarACGtauRise;   ...
                     AG2GJTable.refACGtauRise;   ...
                     AG2GJTable.tarACGtauRise]],...
            'CDataMode','auto','FaceColor','interp','Edges',{0:1/30:1.5 0:1/5:10},'EdgeColor','none')
        view(2)
        
        % xlim([0 0.8])
        xlim([0.1 0.8])
        ylim([0 10])
        clim([0 150])
        
        colormap turbo
        colorbar('north','color','w')
        ylabel('ACG \tau_r_i_s_e')
        xlabel('Trough-to-peak [ms]')
        set(gca,'FontSize',5)
        set(gca,'FontName','Arial')
        box off
        set(gcf, 'Renderer', 'Painters');
    
%     if filtExqFlag
        hcomb = figure(103);
    
        res_type = 'QHD';
        pos = [70 230 1920*0.125 0.6*1080/3]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
        
        pairGroupStatTableGJ = pairGroupStatTableGJ(randperm(size(pairGroupStatTableGJ,1),500),:);
        
        pyramids = [[cell2num([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p'));       ...
                               pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'));        ...
                               pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide'));   ...
                               pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide'))]); ...
                               ratNgjTable.refTroughToPeak(strcmp(ratNgjTable.refCellExplorerType,'p')); ...
                               ratNgjTable.tarTroughToPeak(strcmp(ratNgjTable.tarCellExplorerType,'p')); ...
                               ratSgjTable.refTroughToPeak(strcmp(ratSgjTable.refCellExplorerType,'p')); ...
                               ratSgjTable.tarTroughToPeak(strcmp(ratSgjTable.tarCellExplorerType,'p')); ...
                               ratUgjTable.refTroughToPeak(strcmp(ratUgjTable.refCellExplorerType,'p')); ...
                               ratUgjTable.tarTroughToPeak(strcmp(ratUgjTable.tarCellExplorerType,'p')); ...
                               AG1GJTable.refTroughToPeak(strcmp(AG1GJTable.refCellExplorerType,'p'));   ...
                               AG1GJTable.tarTroughToPeak(strcmp(AG1GJTable.tarCellExplorerType,'p'));   ...
                               AG1GJTable.refTroughToPeak(strcmp(AG1GJTable.refCellExplorerType,'i-wide'));  ...
                               AG1GJTable.tarTroughToPeak(strcmp(AG1GJTable.tarCellExplorerType,'i-wide'));  ...
                               AG2GJTable.refTroughToPeak(strcmp(AG2GJTable.refCellExplorerType,'p'));   ...
                               AG2GJTable.tarTroughToPeak(strcmp(AG2GJTable.tarCellExplorerType,'p'));   ...
                               AG2GJTable.refTroughToPeak(strcmp(AG2GJTable.refCellExplorerType,'i-wide'));  ...
                               AG2GJTable.tarTroughToPeak(strcmp(AG2GJTable.tarCellExplorerType,'i-wide'))], ... 
                              [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p'));       ...
                               pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'));       ...
                               pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide'));  ...
                               pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide')); ...
                               ratNgjTable.refACGtauRise(strcmp(ratNgjTable.refCellExplorerType,'p')); ...
                               ratNgjTable.tarACGtauRise(strcmp(ratNgjTable.tarCellExplorerType,'p')); ...
                               ratSgjTable.refACGtauRise(strcmp(ratSgjTable.refCellExplorerType,'p')); ...
                               ratSgjTable.tarACGtauRise(strcmp(ratSgjTable.tarCellExplorerType,'p')); ...
                               ratUgjTable.refACGtauRise(strcmp(ratUgjTable.refCellExplorerType,'p')); ...
                               ratUgjTable.tarACGtauRise(strcmp(ratUgjTable.tarCellExplorerType,'p')); ...
                               AG1GJTable.refACGtauRise(strcmp(AG1GJTable.refCellExplorerType,'p'));   ...
                               AG1GJTable.tarACGtauRise(strcmp(AG1GJTable.tarCellExplorerType,'p'));   ...
                               AG1GJTable.refACGtauRise(strcmp(AG1GJTable.refCellExplorerType,'i-wide'));  ...
                               AG1GJTable.tarACGtauRise(strcmp(AG1GJTable.tarCellExplorerType,'i-wide'));  ...
                               AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'p'));   ...
                               AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'p'));   ...
                               AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'i-wide'));  ...
                               AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'i-wide'))]]; 
                           
        interneurons = [[cell2num([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow'));   ...
                                   pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow'))]); ...
                                   ratNgjTable.refTroughToPeak(strcmp(ratNgjTable.refCellExplorerType,'i')); ...
                                   ratNgjTable.tarTroughToPeak(strcmp(ratNgjTable.tarCellExplorerType,'i')); ...
                                   ratSgjTable.refTroughToPeak(strcmp(ratSgjTable.refCellExplorerType,'i')); ...
                                   ratSgjTable.tarTroughToPeak(strcmp(ratSgjTable.tarCellExplorerType,'i')); ...
                                   ratUgjTable.refTroughToPeak(strcmp(ratUgjTable.refCellExplorerType,'i')); ...
                                   ratUgjTable.tarTroughToPeak(strcmp(ratUgjTable.tarCellExplorerType,'i')); ...
                                   AG1GJTable.refTroughToPeak(strcmp(AG1GJTable.refCellExplorerType,'i-narrow'));  ...
                                   AG1GJTable.tarTroughToPeak(strcmp(AG1GJTable.tarCellExplorerType,'i-narrow'));  ...
                                   AG2GJTable.refTroughToPeak(strcmp(AG2GJTable.refCellExplorerType,'i-narrow'));  ...
                                   AG2GJTable.tarTroughToPeak(strcmp(AG2GJTable.tarCellExplorerType,'i-narrow'))], ... 
                                  [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow'));  ...
                                   pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow')); ...
                                   ratNgjTable.refACGtauRise(strcmp(ratNgjTable.refCellExplorerType,'i')); ...
                                   ratNgjTable.tarACGtauRise(strcmp(ratNgjTable.tarCellExplorerType,'i')); ...
                                   ratSgjTable.refACGtauRise(strcmp(ratSgjTable.refCellExplorerType,'i')); ...
                                   ratSgjTable.tarACGtauRise(strcmp(ratSgjTable.tarCellExplorerType,'i')); ...
                                   ratUgjTable.refACGtauRise(strcmp(ratUgjTable.refCellExplorerType,'i')); ...
                                   ratUgjTable.tarACGtauRise(strcmp(ratUgjTable.tarCellExplorerType,'i')); ...
                                   AG1GJTable.refACGtauRise(strcmp(AG1GJTable.refCellExplorerType,'i-narrow'));  ...
                                   AG1GJTable.tarACGtauRise(strcmp(AG1GJTable.tarCellExplorerType,'i-narrow'));  ...
                                   AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'i-narrow'));  ...
                                   AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'i-narrow'))]];
        
        plot(pyramids(:,1),pyramids(:,2),'^','LineWidth',1)
        hold on
        plot(interneurons(:,1),interneurons(:,2),'o','LineWidth',1)
%     else
%         semilogy([cell2num([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p'));       ...
%                             pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'));        ...
%                             pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide'));   ...
%                             pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide'))]); ...
%                             ratNgjTable.refTroughToPeak(strcmp(ratNgjTable.refCellExplorerType,'p')); ...
%                             ratNgjTable.tarTroughToPeak(strcmp(ratNgjTable.tarCellExplorerType,'p')); ...
%                             ratSgjTable.refTroughToPeak(strcmp(ratSgjTable.refCellExplorerType,'p')); ...
%                             ratSgjTable.tarTroughToPeak(strcmp(ratSgjTable.tarCellExplorerType,'p')); ...
%                             ratUgjTable.refTroughToPeak(strcmp(ratUgjTable.refCellExplorerType,'p')); ...
%                             ratUgjTable.tarTroughToPeak(strcmp(ratUgjTable.tarCellExplorerType,'p')); ...
%                             AG1GJTable.refTroughToPeak(strcmp(AG1GJTable.refCellExplorerType,'p'));   ...
%                             AG1GJTable.tarTroughToPeak(strcmp(AG1GJTable.tarCellExplorerType,'p'));   ...
%                             AG1GJTable.refTroughToPeak(strcmp(AG1GJTable.refCellExplorerType,'i-wide'));  ...
%                             AG1GJTable.tarTroughToPeak(strcmp(AG1GJTable.tarCellExplorerType,'i-wide'));  ...
%                             AG2GJTable.refTroughToPeak(strcmp(AG2GJTable.refCellExplorerType,'p'));   ...
%                             AG2GJTable.tarTroughToPeak(strcmp(AG2GJTable.tarCellExplorerType,'p'));   ...
%                             AG2GJTable.refTroughToPeak(strcmp(AG2GJTable.refCellExplorerType,'i-wide'));  ...
%                             AG2GJTable.tarTroughToPeak(strcmp(AG2GJTable.tarCellExplorerType,'i-wide'))], ... 
%                            [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p'));       ...
%                             pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'));       ...
%                             pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide'));  ...
%                             pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide')); ...
%                             ratNgjTable.refACGtauRise(strcmp(ratNgjTable.refCellExplorerType,'p')); ...
%                             ratNgjTable.tarACGtauRise(strcmp(ratNgjTable.tarCellExplorerType,'p')); ...
%                             ratSgjTable.refACGtauRise(strcmp(ratSgjTable.refCellExplorerType,'p')); ...
%                             ratSgjTable.tarACGtauRise(strcmp(ratSgjTable.tarCellExplorerType,'p')); ...
%                             ratUgjTable.refACGtauRise(strcmp(ratUgjTable.refCellExplorerType,'p')); ...
%                             ratUgjTable.tarACGtauRise(strcmp(ratUgjTable.tarCellExplorerType,'p')); ...
%                             AG1GJTable.refACGtauRise(strcmp(AG1GJTable.refCellExplorerType,'p'));   ...
%                             AG1GJTable.tarACGtauRise(strcmp(AG1GJTable.tarCellExplorerType,'p'));   ...
%                             AG1GJTable.refACGtauRise(strcmp(AG1GJTable.refCellExplorerType,'i-wide'));  ...
%                             AG1GJTable.tarACGtauRise(strcmp(AG1GJTable.tarCellExplorerType,'i-wide'));  ...
%                             AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'p'));   ...
%                             AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'p'));   ...
%                             AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'i-wide'));  ...
%                             AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'i-wide'))], ...
%                             '^','LineWidth',1)
%         hold on
%         semilogy([cell2num([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow'));   ...
%                             pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow'))]); ...
%                             ratNgjTable.refTroughToPeak(strcmp(ratNgjTable.refCellExplorerType,'i')); ...
%                             ratNgjTable.tarTroughToPeak(strcmp(ratNgjTable.tarCellExplorerType,'i')); ...
%                             ratSgjTable.refTroughToPeak(strcmp(ratSgjTable.refCellExplorerType,'i')); ...
%                             ratSgjTable.tarTroughToPeak(strcmp(ratSgjTable.tarCellExplorerType,'i')); ...
%                             ratUgjTable.refTroughToPeak(strcmp(ratUgjTable.refCellExplorerType,'i')); ...
%                             ratUgjTable.tarTroughToPeak(strcmp(ratUgjTable.tarCellExplorerType,'i')); ...
%                             AG1GJTable.refTroughToPeak(strcmp(AG1GJTable.refCellExplorerType,'i-narrow'));  ...
%                             AG1GJTable.tarTroughToPeak(strcmp(AG1GJTable.tarCellExplorerType,'i-narrow'));  ...
%                             AG2GJTable.refTroughToPeak(strcmp(AG2GJTable.refCellExplorerType,'i-narrow'));  ...
%                             AG2GJTable.tarTroughToPeak(strcmp(AG2GJTable.tarCellExplorerType,'i-narrow'))], ... 
%                            [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow'));  ...
%                             pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow')); ...
%                             ratNgjTable.refACGtauRise(strcmp(ratNgjTable.refCellExplorerType,'i')); ...
%                             ratNgjTable.tarACGtauRise(strcmp(ratNgjTable.tarCellExplorerType,'i')); ...
%                             ratSgjTable.refACGtauRise(strcmp(ratSgjTable.refCellExplorerType,'i')); ...
%                             ratSgjTable.tarACGtauRise(strcmp(ratSgjTable.tarCellExplorerType,'i')); ...
%                             ratUgjTable.refACGtauRise(strcmp(ratUgjTable.refCellExplorerType,'i')); ...
%                             ratUgjTable.tarACGtauRise(strcmp(ratUgjTable.tarCellExplorerType,'i')); ...
%                             AG1GJTable.refACGtauRise(strcmp(AG1GJTable.refCellExplorerType,'i-narrow'));  ...
%                             AG1GJTable.tarACGtauRise(strcmp(AG1GJTable.tarCellExplorerType,'i-narrow'));  ...
%                             AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'i-narrow'));  ...
%                             AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'i-narrow'))], ...
%                             'o','LineWidth',1)
%     end
    for loopPairs = 1:size(pairGroupStatTableGJ,1)
        patchline(cell2num([pairGroupStatTableGJ.refTroughToPeakLength(loopPairs),pairGroupStatTableGJ.tarTroughToPeakLength(loopPairs)]), ... 
                           [pairGroupStatTableGJ.refACGtauRise(loopPairs),        pairGroupStatTableGJ.tarACGtauRise(loopPairs)], ...
                          'linewidth',1,'edgealpha',alphaGraph2)
    end
    for loopPairs = 1:size(ratNgjTable,1)
        patchline([ratNgjTable.refTroughToPeak(loopPairs),ratNgjTable.tarTroughToPeak(loopPairs)], ... 
                  [ratNgjTable.refACGtauRise(loopPairs),  ratNgjTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',alphaGraph2)
    end
    for loopPairs = 1:size(ratSgjTable,1)
        patchline([ratSgjTable.refTroughToPeak(loopPairs),ratSgjTable.tarTroughToPeak(loopPairs)], ... 
                  [ratSgjTable.refACGtauRise(loopPairs),  ratSgjTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',alphaGraph2)
    end
    for loopPairs = 1:size(ratUgjTable,1)
        patchline([ratUgjTable.refTroughToPeak(loopPairs),ratUgjTable.tarTroughToPeak(loopPairs)], ... 
                  [ratUgjTable.refACGtauRise(loopPairs),  ratUgjTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',alphaGraph2)
    end
    if ~filtExqFlag
        for loopPairs = 1:size(AG1GJTable,1)
            patchline([AG1GJTable.refTroughToPeak(loopPairs),AG1GJTable.tarTroughToPeak(loopPairs)], ... 
                      [AG1GJTable.refACGtauRise(loopPairs),AG1GJTable.tarACGtauRise(loopPairs)], ...
                      'linewidth',2,'edgealpha',alphaGraph2)
        end
    end
    for loopPairs = 1:size(AG2GJTable,1)
        patchline([AG2GJTable.refTroughToPeak(loopPairs),AG2GJTable.tarTroughToPeak(loopPairs)], ... 
                  [AG2GJTable.refACGtauRise(loopPairs),AG2GJTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',2,'edgealpha',alphaGraph2)
    end
    hold off
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
    legend('p','i','Location','best')
%     ylim([0.05 30])
%     xlim([0.1 2])
    ylim([0 10])
    % xlim([0 1])
    xlim([0.1 0.8])
    ylabel('ACG \tau_r_i_s_e')
    xlabel('Trough-to-peak [ms]')
    set(gca,'FontSize',5)
    set(gca, 'FontName', 'Arial')
    set(gcf, 'Renderer', 'Painters');
    box off
    
    %% count pairs
    
    % GJlabels = table;
    % 
    % GJlabels.refCellExplorerType = [pairGroupStatTableGJ.refCellExplorerType
    %                                 ratNgjTable.refCellExplorerType
    %                                 ratSgjTable.refCellExplorerType
    %                                 ratUgjTable.refCellExplorerType
    %                                 royGJTable.refCellExplorerType
    %                                 AG1GJTable.refCellExplorerType
    %                                 AG2GJTable.refCellExplorerType
    %                                 royGJTable.refCellExplorerType];
    % 
    % GJlabels.tarCellExplorerType = [pairGroupStatTableGJ.tarCellExplorerType
    %                                 ratNgjTable.tarCellExplorerType
    %                                 ratSgjTable.tarCellExplorerType
    %                                 ratUgjTable.tarCellExplorerType
    %                                 royGJTable.tarCellExplorerType
    %                                 AG1GJTable.tarCellExplorerType
    %                                 AG2GJTable.tarCellExplorerType
    %                                 royGJTable.tarCellExplorerType];
    % 
    % for loopRows = 1:size(GJlabels,1)
    %     if strcmp(GJlabels.refCellExplorerType(loopRows),'i-wide')
    %         GJlabels.refCellExplorerType{loopRows} = 'p';
    %     end
    %     if strcmp(GJlabels.tarCellExplorerType(loopRows),'i-wide')
    %         GJlabels.tarCellExplorerType{loopRows} = 'p';
    %     end
    %     if strcmp(GJlabels.refCellExplorerType(loopRows),'i-narrow')
    %         GJlabels.refCellExplorerType{loopRows} = 'i';
    %     end
    %     if strcmp(GJlabels.tarCellExplorerType(loopRows),'i-narrow')
    %         GJlabels.tarCellExplorerType{loopRows} = 'i';
    %     end
    % end
    % 
    % ii_GJ_pairs_no    = sum(strcmp(GJlabels.refCellExplorerType,'i') & strcmp(GJlabels.tarCellExplorerType,'i'))
    % pp_GJ_pairs_no    = sum(strcmp(GJlabels.refCellExplorerType,'p') & strcmp(GJlabels.tarCellExplorerType,'p'))
    % 
    % pi_ip_GJ_pairs_no = sum((strcmp(GJlabels.refCellExplorerType,'p') & strcmp(GJlabels.tarCellExplorerType,'i')) | ... 
    %                         (strcmp(GJlabels.refCellExplorerType,'i') & strcmp(GJlabels.tarCellExplorerType,'p')))
    
end

function fwhm = findFwhm(ccgs)

    for loopCCGs = 1:size(ccgs,2)
        
        % Find the half max value.
        halfMax = (min(ccgs(:,loopCCGs)) + max(ccgs(:,loopCCGs))) / 2;
        % Find where the data first drops below half the max.
        index1 = find(ccgs(:,loopCCGs) >= halfMax, 1, 'first');
        % Find where the data last rises above half the max.
        index2 = find(ccgs(:,loopCCGs) >= halfMax, 1, 'last');
        
        fwhm(loopCCGs) = index2-index1 + 1; % FWHM in indexes.
    end
end




