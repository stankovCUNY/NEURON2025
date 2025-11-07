function BapunFixMergeSpikeOffset

    fs = 3e4;
    t  = (-35:36)*(1000/fs);

    data_dir = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    animalIDs  = {'N','S','U'};
    
    for loopAnimal = 1:length(animalIDs)
        
        RatID = animalIDs{loopAnimal};
        disp(['looping rat: ' RatID])
        
        if strcmp(RatID,'N')

            chanID_xml = [13  2   14  1   15  0   12  3   11  4   10  5   9   6   8   7   ...
                          29  18  30  17  31  16  28  19  27  20  26  21  25  22  24  23  ...
                          45  34  46  33  47  32  44  35  43  36  42  37  41  38  40  39  ...
                          61  50  62  49  63  48  60  51  59  52  58  53  57  54  56  55  ...
                          77  66  78  65  79  64  76  67  75  68  74  69  73  70  72  71  ...
                          93  82  94  81  95  80  92  83  91  84  90  85  89  86  88  87  ...
                          109 98  110 97  111 96  108 99  107 100 106 101 105 102 104 103 ...
                          125 114 126 113 127 112 124 115 123 116 122 117 121 118 120 119
                          ];  

            neurons     = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
            snippetsDir = '/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/RatN/snippetsPeakChan/';
            saveName    = 'RatN_adjustedSpikeTimes.mat';

            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/unitsGlobalSnapshots.mat');

        elseif strcmp(RatID,'S')

            chanID_xml = [19  3   26  20  5   8   14  0   4   25     ...
                          2   17  16  6   24  21  10  30  22  9   1  ...
                          27  18  31  12  23  28  15  11  29  7   13 ...
                          36  48  58  56  45  32  53  47  33  40  50 ...
                          34  41  52  54  43  37  57  60  49  44  63 ...
                          46  61  59  39  35  62  42  38  55         ...
                          77  66  78  65  79  64  76  67  75  68  74  69  73  70  72  71 ...
                          93  82  94  81  95  80  92  83  91  84  90  85  89  86  88  87 ...
                          109 98  110 97  111 96  108 99  107 100 106 101 105 102 104 103 ...
                          125 114 126 113 127 112 124 115 123 116 122 117 121 118 120 119 ...
                          141 130 142 129 143 128 140 131 139 132 138 133 137 134 136 135 ...
                          157 146 158 145 159 144 156 147 155 148 154 149 153 150 152 151 ...
                          173 162 174 161 175 160 172 163 171 164 170 165 169 166 168 167 ...
                          189 178 190 177 191 176 188 179 187 180 186 181 185 182 184 183];

            neurons     = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
            snippetsDir =      '/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/RatS/snippetsPeakChan/';
            saveName    = 'RatS_adjustedSpikeTimes.mat';

            load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/unitsGlobalSnapshots.mat');

        elseif strcmp(RatID,'U')

            chanID_xml = [1:16, ...
                          17:32, ... 
                          33:48, ...
                          49:64, ...
                          65:80, ...
                          81:96, ...
                          97:112, ...
                          113:128, ...
                          129:144, ...
                          145:160, ...
                          161:176, ...
                          177:192];

            neurons     = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
            snippetsDir =      '/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/RatU/snippetsPeakChan/';
            saveName    = 'RatU_adjustedSpikeTimes.mat';

        end

        adjustedSpikeTimes = table;

        for loopUnits = 1:length(neurons.neuron_ids) 

            cellID     = neurons.neuron_ids(loopUnits) + 1;
    %         peakChan   = neurons.peak_channels(loopUnits) + 1;
            if strcmp(RatID,'N') || strcmp(RatID,'S')
                [~,peakChan] = max(max(abs(unitsGlobalSnapshots.waveforms{loopUnits,1}),[],2));
            else
                peakChan = find(double(neurons.peak_channels(loopUnits) + 1) == chanID_xml);
            end
            spikeTimes = neurons.spiketrains{loopUnits};

            if exist([snippetsDir 'unit' num2str(cellID) 'peakCh' num2str(peakChan) 'snippets.mat'],'file') > 0

                load([snippetsDir 'unit' num2str(cellID) 'peakCh' num2str(peakChan) 'snippets.mat']);

                if strcmp(neurons.neuron_type(loopUnits,:),'mua  ')
                    continue
                end

                %% determine positive or negative spike
                [~,idxFeature] = max(abs(snippetsData.waveformAvg));
                if snippetsData.waveformAvg(idxFeature) < 0
                    [~,idx] = min(snippetsData.waveforms);
                elseif snippetsData.waveformAvg(idxFeature) > 0
                    [~,idx] = max(snippetsData.waveforms);
                end

                modeWavePeakIdx = mode(idx);

                offset = idx - modeWavePeakIdx;
                
                loopUnits
                
                adjustedSpikeTimes.cellID(loopUnits)     = cellID;
                
                % quick fix sometimes there are more spike times than waveforms (need to debug why)
                if length(spikeTimes) > length(offset)
                    spikeTimes = spikeTimes(1:length(offset));
                end
                
                adjustedSpikeTimes.spikeTimes{loopUnits} = spikeTimes - offset*(1e3/fs);
                                
            end

        end

        save([data_dir 'Rat' RatID '/' saveName],'adjustedSpikeTimes')
    end
end