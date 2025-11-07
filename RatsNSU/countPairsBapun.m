function noPairCompTable = countPairsBapun

    noPairCompTable = table;

    %% Rat N

    noPairCompTable.animalID(1,:) = "N";

    neuronsN    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');

    noPairCompTable.probe(1,:)            = 1;
    noPairCompTable.brainRegion(1,:)      = "CA1";
    n = sum(~(sum(neuronsN.neuron_type == 'mua  ',2) == 5));
    noPairCompTable.nPairComparisons(1,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsN.neuron_type == 'pyr  ',2) == 5));
    noPairCompTable.nPairCompPP(1,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsN.neuron_type == 'inter',2) == 5));
    noPairCompTable.nPairCompII(1,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(1,:) = noPairCompTable.nPairComparisons(1,:) - ...
                                           noPairCompTable.nPairCompPP(1,:) - ...
                                           noPairCompTable.nPairCompII(1,:);
    
    %% Rat S

    noPairCompTable.animalID([2,3],:) = "S";

    neuronsS    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');

    noPairCompTable.probe(2,:)            = 1;
    noPairCompTable.brainRegion(2,:)      = "CA1";
    n = sum(~(sum(neuronsS.neuron_type == 'mua  ',2) == 5) & (neuronsS.shank_ids <= 5)');
    noPairCompTable.nPairComparisons(2,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsS.neuron_type == 'pyr  ',2) == 5) & (neuronsS.shank_ids <= 5)');
    noPairCompTable.nPairCompPP(2,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsS.neuron_type == 'inter',2) == 5) & (neuronsS.shank_ids <= 5)');
    noPairCompTable.nPairCompII(2,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(2,:) = noPairCompTable.nPairComparisons(2,:) - ...
                                           noPairCompTable.nPairCompPP(2,:) - ...
                                           noPairCompTable.nPairCompII(2,:);
    
    noPairCompTable.probe(3,:)            = 2;
    noPairCompTable.brainRegion(3,:)      = "CA1";
    n = sum(~(sum(neuronsS.neuron_type == 'mua  ',2) == 5) & (neuronsS.shank_ids > 5)');
    noPairCompTable.nPairComparisons(3,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsS.neuron_type == 'pyr  ',2) == 5) & (neuronsS.shank_ids > 5)');
    noPairCompTable.nPairCompPP(3,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsS.neuron_type == 'inter',2) == 5) & (neuronsS.shank_ids > 5)');
    noPairCompTable.nPairCompII(3,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(3,:) = noPairCompTable.nPairComparisons(3,:) - ...
                                           noPairCompTable.nPairCompPP(3,:) - ...
                                           noPairCompTable.nPairCompII(3,:);
    
    %% Rat U

    noPairCompTable.animalID([4,5],:) = "U";

    neuronsU    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');

    noPairCompTable.probe(4,:)            = 1;
    noPairCompTable.brainRegion(4,:)      = "CA1";
    n = sum(~(sum(neuronsU.neuron_type == 'mua  ',2) == 5) & (neuronsU.shank_ids <= 7)');
    noPairCompTable.nPairComparisons(4,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsU.neuron_type == 'pyr  ',2) == 5) & (neuronsU.shank_ids <= 7)');
    noPairCompTable.nPairCompPP(4,:) = (n*(n-1))/2;

    n = sum((sum(neuronsU.neuron_type == 'inter',2) == 5) & (neuronsU.shank_ids <= 7)');
    noPairCompTable.nPairCompII(4,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(4,:) = noPairCompTable.nPairComparisons(4,:) - ...
                                           noPairCompTable.nPairCompPP(4,:) - ...
                                           noPairCompTable.nPairCompII(4,:);

    noPairCompTable.probe(5,:)            = 2;
    noPairCompTable.brainRegion(5,:)      = "CA1";
    n = sum(~(sum(neuronsU.neuron_type == 'mua  ',2) == 5) & (neuronsU.shank_ids > 7)');
    noPairCompTable.nPairComparisons(5,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsU.neuron_type == 'pyr  ',2) == 5) & (neuronsU.shank_ids > 7)');
    noPairCompTable.nPairCompPP(5,:) = (n*(n-1))/2;
    
    n = sum((sum(neuronsU.neuron_type == 'inter',2) == 5) & (neuronsU.shank_ids > 7)');
    noPairCompTable.nPairCompII(5,:) = (n*(n-1))/2;
    
    noPairCompTable.nPairCompIPorPI(5,:) = noPairCompTable.nPairComparisons(5,:) - ...
                                           noPairCompTable.nPairCompPP(5,:) - ...
                                           noPairCompTable.nPairCompII(5,:);

end
