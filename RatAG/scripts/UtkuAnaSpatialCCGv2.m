clear; close all;

dayNo = '1';

UtkuData      = loadUtkuV2(dayNo); % NOTE: spike time data is in minutes!!!
filteredPairs = UtkuScreenJitterV2(dayNo,'all');

% constants
fs               = 30000;
BinSize          = 1/fs;
Duration         = 0.03;

fig_use        = 102;
resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.

UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

day = str2double(dayNo);
if day == 1
    datapath = [UtkuPath 'AG_2019-12-23_NSD' '/'];
elseif day == 2
    datapath = [UtkuPath 'AG_2019-12-27_NSD' '/'];
end 

% 
chSimpleLayout = { ...
     [139; 171; 203; 235:  11;  43;  75; 107]  ...                                          %1
     [140; 172; 204; 236;  12;  44;  76; 108], ...                                          %2
     [138; 151; 170; 183; 202; 215; 234; 247;  10;  23;  42;  55;  74;  87; 106; 119], ...  %3
     [141; 173; 205; 237;  13;  45;  77; 109], ...                                          %4
     [137; 152; 169; 184; 201; 216; 233; 248;   9;  24;  41;  56;  73;  88; 105; 120], ...  %5 
     [142; 174; 206; 238;  14;  46;  78; 110], ...                                          %6
     [129; 160; 161; 192; 193; 224; 225; 256;   1;  32;  33;  64;  65;  96;  97; 128], ...  %7
     [150; 182; 214; 246;  22;  54;  86; 118], ...                                          %8
     [130; 159; 162; 191; 194; 223; 226; 255;   2;  31;  34;  63;  66;  95;  98; 127], ...  %9 
     [149; 181; 213; 245;  21;  53;  85; 117], ...                                          %10
     [131; 158; 163; 190; 195; 222; 227; 254;   3;  30;  35;  62;  67;  94;  99; 126], ...  %11
     [148; 180; 212; 244;  20;  52;  84; 116], ...                                          %12
     [132; 157; 164; 189; 196; 221; 228; 253;   4;  29;  36;  61;  68;  93; 100; 125], ...  %13
     [147; 179; 211; 243;  19;  51;  83; 115], ...                                          %14    
     [133; 156; 165; 188; 197; 220; 229; 252;   5;  28;  37;  60;  69;  92; 101; 124], ...  %15
     [146; 178; 210; 242;  18;  50;  82; 114], ...                                          %16
     [134; 155; 166; 187; 198; 219; 230; 251;   6;  27;  38;  59;  70;  91; 102; 123], ...  %17
     [145; 177; 209; 241;  17;  49;  81; 113], ...                                          %18
     [135; 154; 167; 186; 199; 218; 231; 250;   7;  26;  39;  58;  71;  90; 103; 122], ...  %19
     [144; 176; 208; 240;  16;  48;  80; 112], ...                                          %20
     [136; 153; 168; 185; 200; 217; 232; 249;   8;  25;  40;  57;  72;  89; 104; 121], ...  %21
     [143; 175; 207; 239;  15;  47;  79; 111], ...                                          %22
     };

chIdx = { ...
    [  1;   3;   5;   7;   9;  11;  13;  15], ... %1
    [ 17;  19;  21;  23;  25;  27;  29;  31], ... %2
    [ 33:48]', ...                                %3
    [ 49;  51;  53;  55;  57;  59;  61;  63], ... %4
    [ 65:80]', ...                                %5
    [ 81;  83;  85;  87;  89;  91;  93;  95], ... %6
    [ 97:112]', ...                               %7
    [113; 115; 117; 119; 121; 123; 125; 127], ... %8
    [129:144]', ...                               %9
    [145; 147; 149; 151; 153; 155; 157; 159], ... %10
    [161:176]', ...                               %11
    [177; 179; 181; 183; 185; 187; 189; 191], ... %12
    [193:208]', ...                               %13
    [209; 211; 213; 215; 217; 219; 221; 223], ... %14
    [225:240]', ...                               %15
    [241; 243; 245; 247; 249; 251; 253; 255], ... %16
    [257:272]', ...                               %17
    [273; 275; 277; 279; 281; 283; 285; 287], ... %18
    [289:304], ...                                %19
    [305; 307; 309; 311; 313; 315; 317; 319], ... %20
    [321:336], ...                                %21
    [337; 339; 341; 343; 345; 347; 349; 351]  ... %22
    };
 
% Monitor specific plot settings.
screensize = get(0,'screensize');
% initiate figure
hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 1440 2560]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

tiledlayout(22,16, 'Padding', 'none', 'TileSpacing', 'compact'); 

templateTime = (-30:30)/30;

% unique neurons
refNeuron = unique(filteredPairs(:));

for i = 46:length(refNeuron)
    
    %
    ref  = UtkuData.cell_metrics.spikes.times{1,refNeuron(i)};
    res1 = ref; % convert to seconds
    
    cluID       = UtkuData.cell_metrics.UID(refNeuron(i)); % matlab to python indexing?
    shankID     = UtkuData.cell_metrics.shankID(refNeuron(i));
    chIDref     = UtkuData.cell_metrics.maxWaveformCh(refNeuron(i));

    %%
    
    colOne  = (refNeuron(i) == filteredPairs(:,1));
    colTwo  = (refNeuron(i) == filteredPairs(:,2));
    colBoth = (colOne | colTwo);
    
    colIdxFlipped = find(colTwo);
    
    copyNetwork = filteredPairs;
    
    for j = 1:length(colIdxFlipped)
        copyNetwork(colIdxFlipped(j),:) = fliplr(filteredPairs(colIdxFlipped(j),:));
    end
    
    refNetwork = copyNetwork(find(colBoth),:);
%     refNetwork = [refNetwork; refNeuron(i) refNeuron(i)]; % add ACG
    
    %%

    if shankID == 1
        templateChIdx = chIDref;
    elseif shankID == 2
        templateChIdx = chIDref - 32;   
    elseif shankID == 3
        templateChIdx = chIDref - 64;
    elseif shankID == 4
        templateChIdx = chIDref - 96;
    end
    
    clear dataNetPlot

    for j = 1:size(refNetwork,1)
        
        chIDtar    = UtkuData.cell_metrics.maxWaveformCh(refNetwork(j,2));
        for k = 1:size(chSimpleLayout,2)
            temp = chSimpleLayout{k}; % loop shank rows

            if sum(temp == chIDtar)
               dataNetPlot(j).rowIdx1 = k; 
               dataNetPlot(j).rowIdx2 = find(temp == chIDtar);
               
               if sum(k == [1,2,4,6,8,10,12,14,16,18,20,22]); dataNetPlot(j).plotSize = [1,2]; end
               if sum(k == [3,5,7,9,11,13,15,17,19,21]);      dataNetPlot(j).plotSize = [1,1]; end
               
               break
            end
        end
        
        %%
        
        res2 = UtkuData.cell_metrics.spikes.times{refNetwork(j,2)};

        dataNetPlot(j).cell1ID      = UtkuData.cell_metrics.UID(refNeuron(i));
        dataNetPlot(j).cell2ID      = UtkuData.cell_metrics.UID(refNetwork(j,2));

        dataNetPlot(j).cell1shank   = UtkuData.cell_metrics.shankID(refNeuron(i));
        dataNetPlot(j).cell2shank   = UtkuData.cell_metrics.shankID(refNetwork(j,2));

        dataNetPlot(j).cell1channel = UtkuData.cell_metrics.maxWaveformCh1(refNeuron(i)); 
        dataNetPlot(j).cell2channel = UtkuData.cell_metrics.maxWaveformCh1(refNetwork(j,2)); 
        
        cell1type                   = UtkuData.cell_metrics.putativeCellType{1,refNeuron(i)};
        cell2type                   = UtkuData.cell_metrics.putativeCellType{1,refNetwork(j,2)};
        
        if strcmp(cell1type,'Pyramidal Cell')
            cell1type = 'p';
        elseif strcmp(cell1type,'Narrow Interneuron') || strcmp(cell1type,'Wide Interneuron')
            cell1type = 'i';
        end
        
        dataNetPlot(j).cell1type    = cell1type;
        dataNetPlot(j).cell2type    = cell2type;
        
        dataNetPlot(j).cell1spikeTimes = res1;
        dataNetPlot(j).cell2spikeTimes = res2;
        
        dataNetPlot(j).cell1waveform = UtkuData.cell_metrics.waveforms.raw{1,refNeuron(i)};
        dataNetPlot(j).cell1waveTime = UtkuData.cell_metrics.waveforms.time{1,refNeuron(i)};
        
        [ccgR, tR] = CCG([res1;res2],[ones(size(res1));2*ones(size(res2))], ...
            'binSize', BinSize, 'duration', Duration, 'Fs', 1/fs,...
            'norm', 'counts');

        ccg = ccgR(:,1,2)/sum(ccgR(:,1,2));
%         ccg = ccgR(:,1,2);
        if length(res1) == length(res2)
            ccg(round(length(ccg)/2)) = 0;
        end
        dataNetPlot(j).CCG = ccg;
    end
    
    targetCh =[dataNetPlot.cell2channel];
    uniqueCh = unique(targetCh);
    
    for j = 1:length(uniqueCh)
                
        clear legendStr
        
        occurances = find(targetCh == uniqueCh(j));
        nexttile(chIdx{dataNetPlot(occurances(1)).rowIdx1}(dataNetPlot(occurances(1)).rowIdx2), ...
                 dataNetPlot(occurances(1)).plotSize)
        plot(tR*1000,dataNetPlot(occurances(1)).CCG)
        if length(occurances) > 1
            hold on
            for k = 2:length(occurances)
                 plot(tR*1000,dataNetPlot(occurances(k)).CCG)
            end
            hold off
            
            legendStr{1} = ['tar'  num2str(dataNetPlot(occurances(1)).cell2ID)];
            for k = 2:length(occurances)
                legendStr{k} = ['tar'  num2str(dataNetPlot(occurances(k)).cell2ID)];
            end
        else 
            legendStr = ['tar'  num2str(dataNetPlot(occurances(1)).cell2ID)];
        end
%         
%         yyaxis right
%         plot(dataNetPlot(j).cell1waveTime,dataNetPlot(j).cell1waveform)
%         set(gca, 'YDir','reverse')
        
        if dataNetPlot(occurances(1)).plotSize(2) == 1
            xlim([-1.5 1.5])
        elseif dataNetPlot(occurances(1)).plotSize(2) == 2
            xlim([-3 3])
        end
%         legend(legendStr,'Location','southoutside','FontSize',8)
        
%         ylabel('prob.')
%         xlabel('[ms]')
        
        if chIDref == chSimpleLayout{dataNetPlot(occurances(1)).rowIdx1}(dataNetPlot(occurances(1)).rowIdx2)
            title(['\color{red}' sprintf('sh:%d ch:%d',dataNetPlot(occurances(1)).cell2shank, ...
                chSimpleLayout{dataNetPlot(occurances(1)).rowIdx1}(dataNetPlot(occurances(1)).rowIdx2))])
        else
            title(sprintf('sh:%d ch:%d',dataNetPlot(occurances(1)).cell2shank, ...
                chSimpleLayout{dataNetPlot(occurances(1)).rowIdx1}(dataNetPlot(occurances(1)).rowIdx2)))
        end
    end
    
%     save_file = fullfile(datapath, ['sig_spatialCCGs_Utku_neuron_' num2str(refNeuron(i))]);
%     print(fig_use, save_file,'-djpeg',resolution_use);
    
    close all
    
    hcomb = figure(102);
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    
    tiledlayout(22,16, 'Padding', 'none', 'TileSpacing', 'compact'); 
end