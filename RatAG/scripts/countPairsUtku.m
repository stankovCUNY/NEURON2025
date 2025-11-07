function noPairCompTable = countPairsUtku

    noPairCompTable = table;

    %%

    noPairCompTable.animalID(1:4,:)   = "AG1";

    UtkuData = loadUtkuV2('1');

    noPairCompTable.probe(1,:)            = 1;
    noPairCompTable.brainRegion(1,:)      = "CA1";  
    n = sum(strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID <= 4));
    noPairCompTable.nPairComparisons(1,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        (strcmp(UtkuData.cell_metrics.putativeCellType,"Pyramidal Cell") | strcmp(UtkuData.cell_metrics.putativeCellType,"Wide Interneuron")));
    noPairCompTable.nPairCompPP(1,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        strcmp(UtkuData.cell_metrics.putativeCellType,"Narrow Interneuron"));
    noPairCompTable.nPairCompII(1,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(1,:) = noPairCompTable.nPairComparisons(1,:) - ...
                                           noPairCompTable.nPairCompPP(1,:) - ...
                                           noPairCompTable.nPairCompII(1,:);
    %%
    
    noPairCompTable.probe(3,:)            = 1;
    noPairCompTable.brainRegion(3,:)      = "CA3"; 
    n = sum(strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID <= 4));
    noPairCompTable.nPairComparisons(3,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        (strcmp(UtkuData.cell_metrics.putativeCellType,"Pyramidal Cell") | strcmp(UtkuData.cell_metrics.putativeCellType,"Wide Interneuron")));
    noPairCompTable.nPairCompPP(3,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        strcmp(UtkuData.cell_metrics.putativeCellType,"Narrow Interneuron"));
    noPairCompTable.nPairCompII(3,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(3,:) = noPairCompTable.nPairComparisons(3,:) - ...
                                           noPairCompTable.nPairCompPP(3,:) - ...
                                           noPairCompTable.nPairCompII(3,:);
    %%
    
    noPairCompTable.probe(2,:)            = 2;
    noPairCompTable.brainRegion(2,:)      = "CA1";
    n = sum(strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID > 4));
    noPairCompTable.nPairComparisons(2,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID > 4)) & ...
        (strcmp(UtkuData.cell_metrics.putativeCellType,"Pyramidal Cell") | strcmp(UtkuData.cell_metrics.putativeCellType,"Wide Interneuron")));
    noPairCompTable.nPairCompPP(2,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID > 4)) & ...
        strcmp(UtkuData.cell_metrics.putativeCellType,"Narrow Interneuron"));
    noPairCompTable.nPairCompII(2,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(2,:) = noPairCompTable.nPairComparisons(2,:) - ...
                                           noPairCompTable.nPairCompPP(2,:) - ...
                                           noPairCompTable.nPairCompII(2,:);
    
    %%
    
    noPairCompTable.probe(4,:)            = 2;
    noPairCompTable.brainRegion(4,:)      = "CA3"; 
    n = sum(strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID > 4));
    noPairCompTable.nPairComparisons(4,:) = (n*(n-1))/2;

    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID > 4)) & ...
        (strcmp(UtkuData.cell_metrics.putativeCellType,"Pyramidal Cell") | strcmp(UtkuData.cell_metrics.putativeCellType,"Wide Interneuron")));
    noPairCompTable.nPairCompPP(4,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID > 4)) & ...
        strcmp(UtkuData.cell_metrics.putativeCellType,"Narrow Interneuron"));
    noPairCompTable.nPairCompII(4,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(4,:) = noPairCompTable.nPairComparisons(4,:) - ...
                                           noPairCompTable.nPairCompPP(4,:) - ...
                                           noPairCompTable.nPairCompII(4,:);
    
    %%

    noPairCompTable.animalID([5,6],:) = "AG2";

    UtkuData = loadUtkuV2('2');

    noPairCompTable.probe(5,:)            = 1;
    noPairCompTable.brainRegion(5,:)      = "CA1";  
    n = sum(strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID <= 4));
    noPairCompTable.nPairComparisons(5,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        (strcmp(UtkuData.cell_metrics.putativeCellType,"Pyramidal Cell") | strcmp(UtkuData.cell_metrics.putativeCellType,"Wide Interneuron")));
    noPairCompTable.nPairCompPP(5,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA1") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        strcmp(UtkuData.cell_metrics.putativeCellType,"Narrow Interneuron"));
    noPairCompTable.nPairCompII(5,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(5,:) = noPairCompTable.nPairComparisons(5,:) - ...
                                           noPairCompTable.nPairCompPP(5,:) - ...
                                           noPairCompTable.nPairCompII(5,:);

    noPairCompTable.probe(6,:)            = 1;
    noPairCompTable.brainRegion(6,:)      = "CA3";
    n = sum(strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID <= 4));
    noPairCompTable.nPairComparisons(6,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        (strcmp(UtkuData.cell_metrics.putativeCellType,"Pyramidal Cell") | strcmp(UtkuData.cell_metrics.putativeCellType,"Wide Interneuron")));
    noPairCompTable.nPairCompPP(6,:) = (n*(n-1))/2;
    
    n = sum((strcmp(UtkuData.cell_metrics.brainRegion,"CA3") & (UtkuData.cell_metrics.shankID <= 4)) & ...
        strcmp(UtkuData.cell_metrics.putativeCellType,"Narrow Interneuron"));
    noPairCompTable.nPairCompII(6,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(6,:) = noPairCompTable.nPairComparisons(6,:) - ...
                                           noPairCompTable.nPairCompPP(6,:) - ...
                                           noPairCompTable.nPairCompII(6,:);

end
