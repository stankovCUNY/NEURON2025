function chanData = UtkuLoadUnfiltChanData(dayNo,shank,channel)

    day = str2double(dayNo);
    if day == 1
        filepath = ['/media/nasko/WD_BLACK/UtkuRawMultiArray/AG/AG_2019-12_23_NSD/shank' num2str(shank) '/'];
    elseif day == 2
        filepath = ['/media/nasko/WD_BLACK/UtkuRawMultiArray/AG/AG_2019-12_27_NSD/shank' num2str(shank) '/'];
    end 
    
    load([filepath 'rawDataCh' num2str(channel) '.mat'],'data');
    chanData = single(data);
end