function runPetersenAna

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
    datapath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';

%     pairGroupStatTableGJ = table;

%     for loopAnimals = 1:58
%     
%         display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])
%     
%         load([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);
%         
%         pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:)];
%         
%     end
%     
%     tiledlayout(1,2)
%     
%     nexttile
%     plot(pairGroupStatTableGJ.CCGbinLagTimes{1,1}*1000,...
%          [            pairGroupStatTableGJ.pairRawCCG{cell2num(pairGroupStatTableGJ.flagGJ(1:end)) == 1,1}],...
%           'k','LineWidth',1);
%     xlim([-3.5 3.5])
%     ylabel('counts')
%     xlabel('Time [msec]')
%     
%     nexttile
%     semilogy([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p')); pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'))], ... 
%              [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p')); pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'))], ...
%              'o','LineWidth',2)
%     
%     hold on
%     semilogy([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow')); pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow'))], ... 
%              [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow')); pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow'))], ...
%              'o','LineWidth',2)
%     semilogy([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide')); pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide'))], ... 
%              [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide')); pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide'))], ...
%              'o','LineWidth',2)
%     for loopPairs = 1:size(pairGroupStatTableGJ,1)
%         semilogy([pairGroupStatTableGJ.refTroughToPeakLength(loopPairs),pairGroupStatTableGJ.tarTroughToPeakLength(loopPairs)], ... 
%                  [pairGroupStatTableGJ.refACGtauRise(loopPairs),pairGroupStatTableGJ.tarACGtauRise(loopPairs)], ...
%                  'k-','LineWidth',2)
%     end
%     hold off
%     
%     legend('p','i-narrow','i-wide')
%     xlim([0 0.6])
%     ylabel('ACG \tau_r_i_s_e')
%     xlabel('trough-to-peak [ms]')

    %% Allen gap junction

    pairGroupStatTableGJ = table;

    for loopAnimals = 1:58

        display(['running dataset: ' datasetID ', animal: '  num2str(animalsList(loopAnimals))])

        load([datapath datasetID num2str(animalsList(loopAnimals)) 'processedGroupStats.mat']);

        pairGroupStatTableGJ = [pairGroupStatTableGJ; pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:)];

    end

    pairGroupStatTableGJ(cell2num(pairGroupStatTableGJ.flagExq == 1),:) = [];

    fwhm = findFwhm([pairGroupStatTableGJ.pairRawCCG{:,1}]);
    
    pairGroupStatTableGJ(fwhm < 6,:) = [];
    
    %% rats 

%     load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStats.mat')
% 
%     royGJTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);
% 
%     royGJTable(cell2num(royGJTable.flagExq == 1),:) = [];
%     
%     fwhm = findFwhm([royGJTable.pairRawCCG{:,1}]);
%     
%     royGJTable(fwhm < 6,:) = [];

    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat')

    AG1GJTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);

    AG1GJTable(cell2num(AG1GJTable.flagExq == 1),:) = [];
    
    fwhm = findFwhm([AG1GJTable.pairRawCCG{:,1}]);
    
    AG1GJTable(fwhm < 6,:) = [];

    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStats.mat')

    AG2GJTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);

    AG2GJTable(cell2num(AG2GJTable.flagExq == 1),:) = [];
    
    fwhm = findFwhm([AG2GJTable.pairRawCCG{:,1}]);
    
    AG2GJTable(fwhm < 6,:) = [];

    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')

    ratSgjTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);

    ratSgjTable(cell2num(ratSgjTable.flagExq == 1),:) = [];
    
    fwhm = findFwhm([ratSgjTable.pairRawCCG{:,1}]);
    
    ratSgjTable(fwhm < 6,:) = [];
    
    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStats.mat')

    ratUgjTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);

    ratUgjTable(cell2num(ratUgjTable.flagExq == 1),:) = [];
    
    fwhm = findFwhm([ratUgjTable.pairRawCCG{:,1}]);
    
    ratUgjTable(fwhm < 6,:) = [];

    %%

    load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStats.mat')

    ratNgjTable = pairGroupStatTable(cell2num(pairGroupStatTable.flagGJ(1:size(pairGroupStatTable,1)/2)) == 1,:);

    ratNgjTable(cell2num(ratNgjTable.flagExq == 1),:) = [];

    fwhm = findFwhm([ratNgjTable.pairRawCCG{:,1}]);
    
    ratNgjTable(fwhm < 6,:) = [];
    
    %%

    figure('Renderer', 'painters', 'Position', [1720 2562 2*560*0.5 420*0.5])

    tiledlayout(1,2,'Padding','none','TileSpacing','compact')

    nexttile
    % plot(pairGroupStatTableGJ.CCGbinLagTimes{1,1}*1000,...
    %      mean([pairGroupStatTableGJ.pairRawCCG{:,1}]/sum([pairGroupStatTableGJ.pairRawCCG{:,1}]),2),...
    %       'k','LineWidth',4);
    % hold on 
    % for loopPairs = 1:size(pairGroupStatTableGJ,1)
    %     patchline(pairGroupStatTableGJ.CCGbinLagTimes{1,1}*1000, ... 
    %               [pairGroupStatTableGJ.pairRawCCG{loopPairs,1}]/sum([pairGroupStatTableGJ.pairRawCCG{loopPairs,1}]), ...
    %              'linewidth',1,'edgealpha',0.1)
    % end
    plot(pairGroupStatTableGJ.CCGbinLagTimes{1,1}*1000,...
         mean([[pairGroupStatTableGJ.pairRawCCG{:,1}]./repmat(cellfun(@length,[pairGroupStatTableGJ.refSpikeTimes]),1,211)', ...
               [ratNgjTable.pairRawCCG{:,1}]./repmat(cellfun(@length,[ratNgjTable.refSpikeTimes]),1,211)', ...
               [ratSgjTable.pairRawCCG{:,1}]./repmat(cellfun(@length,[ratSgjTable.refSpikeTimes]),1,211)', ...
               [ratUgjTable.pairRawCCG{:,1}]./repmat(cellfun(@length,[ratUgjTable.refSpikeTimes]),1,211)' ...
               ],2),...
               'k','LineWidth',1);
%                [AG2GJTable.pairRawCCG{:,1}]./repmat(cellfun(@length,[AG2GJTable.refSpikeTimes]),1,211)', ...
%                [AG1GJTable.pairRawCCG{:,1}], ...
%                [royGJTable.pairRawCCG{:,1}]...
    hold on 
    for loopPairs = 1:size(pairGroupStatTableGJ,1)
        patchline(pairGroupStatTableGJ.CCGbinLagTimes{1,1}*1000, ... 
                 [pairGroupStatTableGJ.pairRawCCG{loopPairs,1}]/length(pairGroupStatTableGJ.refSpikeTimes{loopPairs}), ...
                 'linewidth',1,'edgealpha',0.1)
    end
%     for loopPairs = 1:size(AG1GJTable,1)
%         patchline(AG1GJTable.CCGbinLagTimes{1,1}*1000, ... 
%                  [AG1GJTable.pairRawCCG{loopPairs,1}], ...
%                  'linewidth',2,'edgealpha',0.1)
%     end
%     for loopPairs = 1:size(AG2GJTable,1)
%         patchline(AG2GJTable.CCGbinLagTimes{1,1}*1000, ... 
%                  [AG2GJTable.pairRawCCG{loopPairs,1}], ...
%                  'linewidth',2,'edgealpha',0.1)
%     end
    for loopPairs = 1:size(ratNgjTable,1)
        patchline(ratNgjTable.CCGbinLagTimes{1,1}*1000, ... 
                 [ratNgjTable.pairRawCCG{loopPairs,1}]/length(ratNgjTable.refSpikeTimes{loopPairs}), ...
                 'linewidth',1,'edgealpha',0.1)
    end
    for loopPairs = 1:size(ratSgjTable,1)
        patchline(ratSgjTable.CCGbinLagTimes{1,1}*1000, ... 
                 [ratSgjTable.pairRawCCG{loopPairs,1}]/length(ratSgjTable.refSpikeTimes{loopPairs}), ...
                 'linewidth',1,'edgealpha',0.1)
    end
    for loopPairs = 1:size(ratUgjTable,1)
        patchline(ratUgjTable.CCGbinLagTimes{1,1}*1000, ... 
                 [ratUgjTable.pairRawCCG{loopPairs,1}]/length(ratUgjTable.refSpikeTimes{loopPairs}), ...
                 'linewidth',1,'edgealpha',0.1)
    end
%     for loopPairs = 1:size(royGJTable,1)
%         patchline(royGJTable.CCGbinLagTimes{1,1}*1000, ... 
%                  [royGJTable.pairRawCCG{loopPairs,1}], ...
%                  'linewidth',2,'edgealpha',0.1)
%     end
    hold off
    xlim([-3.5 3.5])
    ylabel('Probability')
    xlabel('[ms]')
    set(gca,'FontSize',5)
    set(gca,'FontName', 'Arial')
    box off

    nexttile
    semilogy([cell2num([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p'));       ...
                        pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'));        ...
                        pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide'));   ...
                        pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide'))]); ...
                        ratNgjTable.refTroughToPeak(strcmp(ratNgjTable.refCellExplorerType,'p')); ...
                        ratNgjTable.tarTroughToPeak(strcmp(ratNgjTable.tarCellExplorerType,'p')); ...
                        ratSgjTable.refTroughToPeak(strcmp(ratSgjTable.refCellExplorerType,'p')); ...
                        ratSgjTable.tarTroughToPeak(strcmp(ratSgjTable.tarCellExplorerType,'p')); ...
                        ratUgjTable.refTroughToPeak(strcmp(ratUgjTable.refCellExplorerType,'p')); ...
                        ratUgjTable.tarTroughToPeak(strcmp(ratUgjTable.tarCellExplorerType,'p'))], ... 
                       [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'p'));       ...
                        pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'p'));       ...
                        pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-wide'));  ...
                        pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-wide')); ...
                        ratNgjTable.refACGtauRise(strcmp(ratNgjTable.refCellExplorerType,'p')); ...
                        ratNgjTable.tarACGtauRise(strcmp(ratNgjTable.tarCellExplorerType,'p')); ...
                        ratSgjTable.refACGtauRise(strcmp(ratSgjTable.refCellExplorerType,'p')); ...
                        ratSgjTable.tarACGtauRise(strcmp(ratSgjTable.tarCellExplorerType,'p')); ...
                        ratUgjTable.refACGtauRise(strcmp(ratUgjTable.refCellExplorerType,'p')); ...
                        ratUgjTable.tarACGtauRise(strcmp(ratUgjTable.tarCellExplorerType,'p'))], ...
                        '^','LineWidth',1)
    hold on
%     semilogy(cell2num([AG2GJTable.refTroughToPeakLength(strcmp(AG2GJTable.refCellExplorerType,'p'));       ...
%                        AG2GJTable.tarTroughToPeakLength(strcmp(AG2GJTable.tarCellExplorerType,'p'));        ...
%                        AG2GJTable.refTroughToPeakLength(strcmp(AG2GJTable.refCellExplorerType,'i-wide'));  ...
%                        AG2GJTable.tarTroughToPeakLength(strcmp(AG2GJTable.tarCellExplorerType,'i-wide'))]), ... 
%                       [AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'p'));       ...
%                        AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'p'));        ...
%                        AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'i-wide'));  ...
%                        AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'i-wide'))], ...
%                      '^','LineWidth',4)
    semilogy([cell2num([pairGroupStatTableGJ.refTroughToPeakLength(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow'));   ...
                        pairGroupStatTableGJ.tarTroughToPeakLength(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow'))]); ...
                        ratNgjTable.refTroughToPeak(strcmp(ratNgjTable.refCellExplorerType,'i')); ...
                        ratNgjTable.tarTroughToPeak(strcmp(ratNgjTable.tarCellExplorerType,'i')); ...
                        ratSgjTable.refTroughToPeak(strcmp(ratSgjTable.refCellExplorerType,'i')); ...
                        ratSgjTable.tarTroughToPeak(strcmp(ratSgjTable.tarCellExplorerType,'i')); ...
                        ratUgjTable.refTroughToPeak(strcmp(ratUgjTable.refCellExplorerType,'i')); ...
                        ratUgjTable.tarTroughToPeak(strcmp(ratUgjTable.tarCellExplorerType,'i'))], ... 
                       [pairGroupStatTableGJ.refACGtauRise(strcmp(pairGroupStatTableGJ.refCellExplorerType,'i-narrow'));  ...
                        pairGroupStatTableGJ.tarACGtauRise(strcmp(pairGroupStatTableGJ.tarCellExplorerType,'i-narrow')); ...
                        ratNgjTable.refACGtauRise(strcmp(ratNgjTable.refCellExplorerType,'i')); ...
                        ratNgjTable.tarACGtauRise(strcmp(ratNgjTable.tarCellExplorerType,'i')); ...
                        ratSgjTable.refACGtauRise(strcmp(ratSgjTable.refCellExplorerType,'i')); ...
                        ratSgjTable.tarACGtauRise(strcmp(ratSgjTable.tarCellExplorerType,'i')); ...
                        ratUgjTable.refACGtauRise(strcmp(ratUgjTable.refCellExplorerType,'i')); ...
                        ratUgjTable.tarACGtauRise(strcmp(ratUgjTable.tarCellExplorerType,'i'))], ...
                        'o','LineWidth',1)
%     semilogy(cell2num([AG2GJTable.refTroughToPeakLength(strcmp(AG2GJTable.refCellExplorerType,'i-narrow'));  ...
%                        AG2GJTable.tarTroughToPeakLength(strcmp(AG2GJTable.tarCellExplorerType,'i-narrow'))]), ... 
%                       [AG2GJTable.refACGtauRise(strcmp(AG2GJTable.refCellExplorerType,'i-narrow'));  ...
%                        AG2GJTable.tarACGtauRise(strcmp(AG2GJTable.tarCellExplorerType,'i-narrow'))], ...
%                       'o','LineWidth',4)
    for loopPairs = 1:size(pairGroupStatTableGJ,1)
        patchline(cell2num([pairGroupStatTableGJ.refTroughToPeakLength(loopPairs),pairGroupStatTableGJ.tarTroughToPeakLength(loopPairs)]), ... 
                           [pairGroupStatTableGJ.refACGtauRise(loopPairs),        pairGroupStatTableGJ.tarACGtauRise(loopPairs)], ...
                          'linewidth',1,'edgealpha',0.5)
    end
    for loopPairs = 1:size(ratNgjTable,1)
        patchline([ratNgjTable.refTroughToPeak(loopPairs),ratNgjTable.tarTroughToPeak(loopPairs)], ... 
                  [ratNgjTable.refACGtauRise(loopPairs),  ratNgjTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',0.5)
    end
    for loopPairs = 1:size(ratSgjTable,1)
        patchline([ratSgjTable.refTroughToPeak(loopPairs),ratSgjTable.tarTroughToPeak(loopPairs)], ... 
                  [ratSgjTable.refACGtauRise(loopPairs),  ratSgjTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',0.5)
    end
    for loopPairs = 1:size(ratUgjTable,1)
        patchline([ratUgjTable.refTroughToPeak(loopPairs),ratUgjTable.tarTroughToPeak(loopPairs)], ... 
                  [ratUgjTable.refACGtauRise(loopPairs),  ratUgjTable.tarACGtauRise(loopPairs)], ...
                  'linewidth',1,'edgealpha',0.5)
    end
%     for loopPairs = 1:size(AG2GJTable,1)
%         patchline(cell2num([AG2GJTable.refTroughToPeakLength(loopPairs),AG2GJTable.tarTroughToPeakLength(loopPairs)]), ... 
%                            [AG2GJTable.refACGtauRise(loopPairs),AG2GJTable.tarACGtauRise(loopPairs)], ...
%                           'linewidth',2,'edgealpha',0.5)
%     end
    hold off
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    legend('p','i','Location','NorthWest')
    ylim([0.05 30])
    xlim([0.1 2])
    ylabel('ACG \tau_r_i_s_e')
    xlabel('Trough-to-peak [ms]')
    set(gca,'FontSize',5)
    set(gca, 'FontName', 'Arial')
    set(gcf, 'Renderer', 'Painters');
    box off
    

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




