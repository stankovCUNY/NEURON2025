function UtkuRunCCGsWaveMap(dayNo)

    UtkuData      = loadUtkuV2(dayNo); % NOTE: spike time data is in minutes!!!
    pairType      = 'exquisite';
    filteredPairs = UtkuScreenJitterV2(dayNo,pairType);
    
    %%

    jscale         = 5;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    fpass          = 300;
%     binSize        = 1/5000;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = false;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    Nwaveforms     = 100;
    filterFlag     = false;
    
    chSimpleLayout = ... 
        {[15  25  8   16  26  7   17  27  6   18  28  5   19  29  4   20  30  3   21  31  2   22  32  1   14  24  9   13  23  10  12  11],...
         [47  57  40  48  58  39  49  59  38  50  60  37  51  61  36  52  62  35  53  63  34  54  64  33  46  56  41  45  55  42  44  43],...
         [79  89  72  80  90  71  81  91  70  82  92  69  83  93  68  84  94  67  85  95  66  86  96  65  78  88  73  77  87  74  76  75],...
         [111 121 104 112 122 103 113 123 102 114 124 101 115 125 100 116 126 99  117 127 98  118 128 97  110 120 105 109 119 106 108 107],...
         [143 153 136 144 154 135 145 155 134 146 156 133 147 157 132 148 158 131 149 159 130 150 160 129 142 152 137 141 151 138 140 139],...
         [175 185 168 176 186 167 177 187 166 178 188 165 179 189 164 180 190 163 181 191 162 182 192 161 174 184 169 173 183 170 172 171],...
         [207 217 200 208 218 199 209 219 198 210 220 197 211 221 196 212 222 195 213 223 194 214 224 193 206 216 201 205 215 202 204 203],...
         [239 249 232 240 250 231 241 251 230 242 252 229 243 253 228 244 254 227 245 255 226 246 256 225 238 248 233 237 247 234 236 235],...
         };
    
%     chTiledCoors = {[1 2],3,4,[5 6],7,8,[9 10],11,12,[13 14],15,16,[17 18],19,20,[21 22],23,24,[25 26],27,28,[29 30],31,32,[33 34],35,36,[37 38],39,40,[41,42],[43,44]};
    chTiledCoors = {[1 2],3,4,[1 2],7,8,[1 2],11,12,[1 2],15,16,[1 2],19,20,[1 2],23,24,[1 2],27,28,[1 2],31,32,[1 2],35,36,[1 2],39,40,[1,2],[1,2]};
  
    preLength  = 36;
    postLength = 36;
    
    % organizing and plotting variabls
%     waveTimeShort   = UtkuData.cell_metrics.waveforms.time{1,1};     % 48 time points
    waveTimeLong    = UtkuData.cell_metrics.waveforms.time_all{1,1}; % 72 time points
    electrodeGroups = UtkuData.cell_metrics.general.electrodeGroups; 

    figpath     = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    UtkuPath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    day = str2double(dayNo);
    if day == 1
        datapath = [UtkuPath 'AG_2019-12-23_NSD' '/'];
        RawDataPath = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/' 'AG_2019-12_23_NSD' '/'];
    elseif day == 2
        datapath = [UtkuPath 'AG_2019-12-27_NSD' '/'];
        RawDataPath = ['/media/nasko/WD_BLACK3/UtkuFilteredDataPerChannel/AG/' 'AG_2019-12_27_NSD' '/'];
    end 

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    for i = 46:size(filteredPairs,1)
        tic
        
        res1    = UtkuData.cell_metrics.spikes.times{1,filteredPairs(i,1)};
        res2    = UtkuData.cell_metrics.spikes.times{1,filteredPairs(i,2)};

        cellID1 = UtkuData.cell_metrics.UID(filteredPairs(i,1));
        cellID2 = UtkuData.cell_metrics.UID(filteredPairs(i,2));

        ch1     = UtkuData.cell_metrics.maxWaveformCh1(filteredPairs(i,1));
        ch2     = UtkuData.cell_metrics.maxWaveformCh1(filteredPairs(i,2));

        clu1    = UtkuData.cell_metrics.cluID(filteredPairs(i,1));
        clu2    = UtkuData.cell_metrics.cluID(filteredPairs(i,2));
        
        shank1  = UtkuData.cell_metrics.shankID(filteredPairs(i,1));
        shank2  = UtkuData.cell_metrics.shankID(filteredPairs(i,2));
        
        % channel index set
        channelIndexSet = chSimpleLayout{shank1};
        
        cell1type = UtkuData.cell_metrics.putativeCellType{1,filteredPairs(i,1)};
        cell2type = UtkuData.cell_metrics.putativeCellType{1,filteredPairs(i,2)};
        
        region1 = UtkuData.cell_metrics.brainRegion{1,filteredPairs(i,1)}; 
        region2 = UtkuData.cell_metrics.brainRegion{1,filteredPairs(i,2)};
        
        cell1waveform = UtkuData.cell_metrics.waveforms.filt_all{1,filteredPairs(i,1)};
        cell2waveform = UtkuData.cell_metrics.waveforms.filt_all{1,filteredPairs(i,2)};
        
        waveform1Idx = find(electrodeGroups{1,shank1} == ch1);
        waveform2Idx = find(electrodeGroups{1,shank2} == ch2);
        
        % some error handling: use on-shank pairs only and different
        % channels
        if shank1 ~= shank2
            break
        elseif ch1 == ch2
            break
        end
        
%         spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
%         spikeTimeIndxCell2 = round((res2-onsetTime/1e6)*fs) + 1;
        
        if strcmp(cell1type,'Pyramidal Cell')
            cell1type = 'p';
        elseif strcmp(cell1type,'Narrow Interneuron')
            cell1type = 'i-narrow';
        elseif strcmp(cell1type,'Wide Interneuron')
            cell1type = 'i-wide';
        end
        
        if strcmp(cell2type,'Pyramidal Cell')
            cell2type = 'p';
        elseif strcmp(cell2type,'Narrow Interneuron')
            cell2type = 'i-narrow';
        elseif strcmp(cell2type,'Wide Interneuron')
            cell2type = 'i-wide';
        end
        
        tiledlayout(22,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
        
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);
        
        % removing the synchronous spikes                
        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(res1, res2, ones(211,1),duration);
                        
        res1sync = [];
        res2sync = [];

        for k = 1:211

            if (GSPExc(k) == 1) || (GSPInh(k) == 1)
                res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
                res2sync = [res2sync; SyncSpBinAll{k}(:,2)];
            end

        end

        % remove sync spikes for waveforms
        res1nosync = setdiff(res1,res1sync);
        res2nosync = setdiff(res2,res2sync);
        
%         if Nwaveforms < length(res1nosync)
%             res1nosync = res1nosync(randperm(length(res1nosync),Nwaveforms));
%         end
%         if Nwaveforms < length(res2nosync)
%             res2nosync = res2nosync(randperm(length(res2nosync),Nwaveforms));
%         end
        
        spikeTimeIndxCell1 = round((res1nosync)*fs) + 1; 
        spikeTimeIndxCell2 = round((res2nosync)*fs) + 1;                
                        
        % some data prep for plotting
        pairCCGscreen.pair(i,:) = [cellID1,cellID2];
        pairCCGscreen.GSPExc{i} = GSPExc;
        pairCCGscreen.GSPInh{i} = GSPInh;

%             pairStr = ['AG' dayNo ' - ' num2str(cellID1) cell1type ' (ch: ' num2str(ch1) ', sh: ' num2str(shank1) ', clu: ' num2str(clu1) ', region: ' region1 ')' ' v ' ...
%                                         num2str(cellID2) cell2type ' (ch: ' num2str(ch2) ', sh: ' num2str(shank2) ', clu: ' num2str(clu2) ', region: ' region2 ')'];

%             pairStr = [num2str(cellID1) cell1type ' v ' ...
%                        num2str(cellID2) cell2type];
% 
%             title(pairStr)
% 
%             ylims = get(gca,'ylim');
% 
%             if any(GSPExc)
%                 hold on;
%                 plot(tR(GSPExc == 1)*1000, 0.95*ylims(2), 'k^');
%             end
%             if any(GSPInh)
%                 hold on;
%                 plot(tR(GSPInh == 1)*1000, 0.95*ylims(2),'kv');
%             end
%             xlim([-(duration/2)*1000 (duration/2)*1000])
        
        tic    
        for j = 1:32
                        
            nexttile(chTiledCoors{j})

            chanDataAtCell1 = load([RawDataPath '/shank' num2str(shank1) '/ch' num2str(channelIndexSet(j)) 'highpass300hz.mat' ]);
            [chTimeSeriesCell1, chWaveformsCell1] = waveformAvg(double(chanDataAtCell1.data),spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);
            clear chanDataAtCell1
            
            chanDataAtCell2 = load([RawDataPath '/shank' num2str(shank2) '/ch' num2str(channelIndexSet(j)) 'highpass300hz.mat' ]);
            [chTimeSeriesCell2, chWaveformsCell2] = waveformAvg(double(chanDataAtCell2.data),spikeTimeIndxCell2,preLength,postLength,fpass,fs,filterFlag);
            clear chanDataAtCell2
            
            plot(waveTimeLong,chTimeSeriesCell1,'LineWidth',3,'Color','#0072BD')
            hold on
            plot(waveTimeLong,chTimeSeriesCell2,'LineWidth',3,'Color','#D95319')
            plot(waveTimeLong,chTimeSeriesCell1+chTimeSeriesCell2,'LineWidth',3,'Color','k')
            hold off
            
            set(gca, 'YDir','reverse')
            xlabel('Time [ms]')
            
            
            if length(chTiledCoors{j}) == 1
                xlim([-1.2 1.2])
            elseif length(chTiledCoors{j}) == 2
                xlim([-2.4 2.4])
                ylabel('revere order voltage [mV]')
                legend(['unit ' num2str(cellID1) cell1type], ...
                       ['unit ' num2str(cellID2) cell2type], ...
                       ['sum'])
            end
            
            if channelIndexSet(j) == ch1
                title(['sh: ' num2str(shank1) ' ch: ' num2str(channelIndexSet(j))],'Color','#0072BD')
            elseif channelIndexSet(j) == ch2
                title(['sh: ' num2str(shank1) ' ch: ' num2str(channelIndexSet(j))],'Color','#D95319')
            else
                title(['sh: ' num2str(shank1) ' ch: ' num2str(channelIndexSet(j))])
            end
        end
        toc
        
%             save_file = fullfile(figpath, ['waveform ana - ' pairStr]);
%             print(fig_use, save_file,'-djpeg',resolution_use);

        close all

        hcomb = figure(102);
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

        toc
        
    end
    
end