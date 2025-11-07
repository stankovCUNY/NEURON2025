function noPairCompTable = countPairsBOT

    % copied code from https://www.mathworks.com/matlabcentral/fileexchange/90900-brain-observatory-toolbox
    
    savePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
    
    % obtain all sessions in the ephys dataset
    sessions = bot.listSessions('ephys');
    
    desiredSessions = sessions;
    animalsList = desiredSessions.id;
    
%     for loopAnimal = 1:length(animalsList)
%         
%         display([num2str((loopAnimal/length(animalsList))*100) '% done'])
%         
%         animalID = double(animalsList(loopAnimal));
%         session = bot.getSessions(desiredSessions(desiredSessions.id == animalID, :));
% 
%     end

    %% 

    noPairCompTable = table;

    for loopAnimal = 1:length(animalsList)
        
        display([num2str((loopAnimal/length(animalsList))*100) '% done'])
        
        animalID = double(animalsList(loopAnimal));
        
%         session = bot.getSessions(desiredSessions(desiredSessions.id == animalID, :));
%         units = session.units;
        
        load([savePath ['neuronsMouse' num2str(animalID) '.mat']])

        probeList = unique(string(neurons.probeName));

        for loopProbe = 1:length(probeList)
            neuronsProbe = neurons(string(neurons.probeName) == probeList(loopProbe),:);

            brainRegionList = unique(string(neuronsProbe.brainRegion));
            for loopBrainRegion = 1:length(brainRegionList)
                
                tempNoPairCompTable = table;
                
                % number of unique pairs
                n = size(neuronsProbe(string(neuronsProbe.brainRegion) == brainRegionList(loopBrainRegion),:),1);
                nPairComparisons = (n*(n-1))/2; 
                
                % number of pp pairs
                n = size(neuronsProbe(string(neuronsProbe.brainRegion) == brainRegionList(loopBrainRegion) & ...
                                     (strcmp(neuronsProbe.putativeCellType,'p') | strcmp(neuronsProbe.putativeCellType,'i-wide')),:),1);
                nPairCompPP = (n*(n-1))/2;
                
                % number of ii pairs
                n = size(neuronsProbe(string(neuronsProbe.brainRegion) == brainRegionList(loopBrainRegion) & ...
                                     strcmp(neuronsProbe.putativeCellType,'i-narrow'),:),1);
                nPairCompII = (n*(n-1))/2;
                
                % number of ip or pi pairs
                nPairCompIPorPI = nPairComparisons - ...
                                  nPairCompPP - ...
                                  nPairCompII;

                tempNoPairCompTable.animalID         = animalID;
                tempNoPairCompTable.probeName        = probeList(loopProbe);
                tempNoPairCompTable.brainRegion      = brainRegionList(loopBrainRegion);
                tempNoPairCompTable.nPairComparisons = nPairComparisons;
                tempNoPairCompTable.nPairCompPP      = nPairCompPP;
                tempNoPairCompTable.nPairCompII      = nPairCompII;
                tempNoPairCompTable.nPairCompIPorPI  = nPairCompIPorPI;

                noPairCompTable = [noPairCompTable; tempNoPairCompTable];
            end
        end

    end
    
    save([savePath 'noPairCompTable.mat'],'noPairCompTable')
end