function noPairCompTable = countPairsHiro

    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";

    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')

    spikesFiltered = spikes.RoyMaze1(:,find([spikes.RoyMaze1.StablePrePost]'));

    for loopUnits = 1:size(spikesFiltered,2)

        temp1 = spikesFiltered(loopUnits).id;
        temp2 = spikesFiltered(loopUnits).quality;

        unitShankIDvector(loopUnits)      = temp1(1);
        unitShankQualityVector(loopUnits) = temp2;
    end

    for loopShanks = 1:8

       unitsOnShank(loopShanks)        = sum(unitShankIDvector == loopShanks);
       unitsOnOTherShanks(loopShanks)  = sum(~(unitShankIDvector == loopShanks));
       
       pUnitsOnShank(loopShanks)       = sum((unitShankIDvector == loopShanks)  & ~(unitShankQualityVector == 8));
       pUnitsOnOTherShanks(loopShanks) = sum(~(unitShankIDvector == loopShanks) & ~(unitShankQualityVector == 8));
       
       iUnitsOnShank(loopShanks)       = sum((unitShankIDvector == loopShanks)  & (unitShankQualityVector == 8));
       iUnitsOnOTherShanks(loopShanks) = sum(~(unitShankIDvector == loopShanks) & (unitShankQualityVector == 8));
       
    end

    noPairCompTable = table;

    noPairCompTable.animalID            = 'Roy';
    noPairCompTable.probeName           = '1';
    noPairCompTable.brainRegion         = 'CA1';
    noPairCompTable.nPairComparisons    = unitsOnShank*unitsOnOTherShanks';
    noPairCompTable.nPairCompPP         = pUnitsOnShank*pUnitsOnOTherShanks';
    noPairCompTable.nPairCompII         = iUnitsOnShank*iUnitsOnOTherShanks';
    noPairCompTable.nPairCompIPorPI     = unitsOnShank*unitsOnOTherShanks' - pUnitsOnShank*pUnitsOnOTherShanks' - iUnitsOnShank*iUnitsOnOTherShanks';
    
end
    
