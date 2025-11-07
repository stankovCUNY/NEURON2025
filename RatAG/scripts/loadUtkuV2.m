function cell_metrics = loadUtkuV2(day)
    
    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    day = str2double(day);
    if day == 1
        cell_metrics = load([UtkuPath 'AG_2019-12-23_NSD' '/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);
    elseif day == 2
        cell_metrics = load([UtkuPath 'AG_2019-12-27_NSD' '/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);
    end    

end