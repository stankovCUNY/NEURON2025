function chanData = UtkuLoad300hz(dayNo,shank,channel)

    day = str2double(dayNo);
    if day == 1
        filepath = ['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(shank) '/'];
    elseif day == 2
        filepath = ['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(shank) '/'];
    end 
    
    load([filepath 'ch' num2str(channel) 'highpass300hz.mat'],'data');
    chanData = single(data);
end