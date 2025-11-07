function plotNeuropixel
    
    subsetFlag = 1;

    animalID = 750332458;
%     probeID = 757904554;

%     animalID = 715093703;
%     probeID = 810755803;

    dataPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';

%     load([dataPath ['neuronsMouse' num2str(animalID) 'probe' num2str(probeID) '.mat']])
    load([dataPath ['neuronsMouse' num2str(animalID) '.mat']])


    dataJscale1 = load([dataPath 'BOT_mouse' num2str(animalID) '_jscale1_alpha5_pairs.mat']);
    dataJscale5 = load([dataPath 'BOT_mouse' num2str(animalID) '_jscale5_alpha5_pairs.mat']);

%     ExcPairs = [dataJscale1.pairs.ExcPairs(:,1:2); dataJscale5.pairs.ExcPairs(:,1:2)];
%     InhPairs = [dataJscale1.pairs.InhPairs(:,1:2); dataJscale5.pairs.InhPairs(:,1:2)];
%     GapPairs = [dataJscale1.pairs.GapPairs(:,1:2); dataJscale5.pairs.GapPairs(:,1:2)];
% 
%     [~,idxUnique,~] = unique(ExcPairs,'rows');
%     ExcPairs(idxUnique,:) = [];
% 
%     [~,idxUnique,~] = unique(InhPairs,'rows');
%     InhPairs(idxUnique,:) = [];
% 
%     [~,idxUnique,~] = unique(GapPairs,'rows');
%     GapPairs(idxUnique,:) = [];
% 
%     allPairs = [ExcPairs; InhPairs; GapPairs];
%     allPairs = [ExcPairs; InhPairs];
%     [~,idxUnique,~] = unique(allPairs,'rows');
%     allPairs = allPairs(idxUnique,:);
    
    %% specific subset
    
    if subsetFlag
        
        allPairs = [951809078 951809195];
        
%         allPairs = [950954941 950955053; 
%                     950950270 950950116;
%                     950928518 950934268;
%                     950928179 950934268;
%                     950928179 950927775];
                
%         allPairs = [950948881 950948853;
%                     950948631 950948026;
%                     950955965 950955844];
        
    end
        
    %%

    jscale         = 1;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    tWave          = (-20:61)/30;
    
    figPath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/figures';

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    for i = 1:size(allPairs,1)
        
        display([num2str(100*(i/size(allPairs,1))) '% done'])
        
        idx1    = find(neurons.unitID == allPairs(i,1));
        idx2    = find(neurons.unitID == allPairs(i,2));
        
        res1    = neurons.spikeTimes{1,idx1};
        res2    = neurons.spikeTimes{1,idx2};

        cellID1 = allPairs(i,1);
        cellID2 = allPairs(i,2);

        ch1     = double(neurons.peakChan(idx1));
        ch2     = double(neurons.peakChan(idx2));

        probe1  = char(neurons.probeName{idx1});
        probe2  = char(neurons.probeName{idx2});
                
        region1 = char(neurons.brainRegion{idx1});
        region2 = char(neurons.brainRegion{idx2});
        
        refAvgWaves = neurons.avgWaveforms{idx1};
        tarAvgWaves = neurons.avgWaveforms{idx2};
        
        tiledlayout(96,4)
        neuropixelChIdx = fliplr(1:size(refAvgWaves,1));
        
        for j = 1:size(refAvgWaves,1)
            
            j
            
            nexttile(neuropixelChIdx(j));
            
            plot(tWave,refAvgWaves(j,:)/1000,'LineWidth',2)
            hold on
            plot(tWave,tarAvgWaves(j,:)/1000,'LineWidth',2)
            hold off
%             ylabel('revere order voltage [mV]')
            xlim([-1 1])
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
        end
            
            
%             nexttile
%             plot(tWave,refAvgWaveMaxChan/1000,'LineWidth',2)
%             hold on
%             plot(tWave,refAvgWaveOpsChan/1000,'LineWidth',2)
%             hold off
%             set(gca, 'YDir','reverse')
%             xlabel('time [ms]')
%             ylabel('revere order voltage [mV]')
%             title(['Reference e-spike: ' num2str(cellID1) ' (ch: ' num2str(ch1) ', probe: ' probe1 ', region: ' region1 ')'])
%             legend('ref. activity at ref. max chan','ref. activity at tar. max chan')
%             xlim([-3.5 3.5])
% %             xlim([-1 1])
% 
%             nexttile
%             plot(tWave,tarAvgWaveMaxChan/1000,'LineWidth',2)
%             hold on
%             plot(tWave,refAvgWaveOpsChan/1000,'LineWidth',2)
%             hold off
%             set(gca, 'YDir','reverse')
%             set(gca, 'XDir','reverse')
%             xlabel('reverse order time [ms]')
%             ylabel('revere order voltage [mV]')
%             title(['Target e-spike: ' num2str(cellID2) ' (ch: ' num2str(ch2) ', probe: ' probe2 ', region: ' region2 ')'])
%             legend('tar. activity at tar. max chan','tar. activity at ref. max chan')
%             xlim([-3.5 3.5])
% %             xlim([-1 1])
%             
%             save_file = fullfile(figPath, pairStr);
%             print(fig_use, save_file,'-djpeg',resolution_use);
% 
%             close all
% 
%             hcomb = figure(102);
%             arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
%         end
        
        
    end
    
end