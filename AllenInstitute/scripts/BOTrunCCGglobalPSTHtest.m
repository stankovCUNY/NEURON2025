function BOTrunCCGglobalPSTHtest
     
     % jitter constants
     jscale         = 1;
     alpha_name     = 5;
     duration       = 0.007;
     fs             = 30000;
     binSize        = 1/fs;
     fig_use        = 102;
     njitter        = 500;
     alpha          = 0.05;
     for_grant      = false;
     plotFlag       = true;
     plotOptoCCG    = false;
     plotUnionCCG   = false;
     resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
     tWave          = (-19:62)/30;
     
     % paths
     dataPath      = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
     figPath       = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/figures';
     unionClusPath = '/media/nasko/WD_BLACK4/unionClusters/';     

     % obtain all sessions in the ephys dataset
     sessions = bot.listSessions('ephys');
    
     % obtain all sessions with the desired session_type and then select the first one
%      desiredSessions = sessions(sessions.session_type == 'brain_observatory_1.1', :); % use Brain Observatory for now
%      desiredSessions = sessions(sessions.session_type == 'functional_connectivity', :); % use Brain Observatory for now
     desiredSessions = sessions;

     animalsList = desiredSessions.id;
    
     for loopAnimal = 1:length(animalsList)
        
        animalID = animalsList(loopAnimal);
         
        display(['animal: ' num2str(animalID)])
        
        if ~exist([figPath '/' num2str(animalID)],'dir')
            mkdir([figPath '/' num2str(animalID)])
        end
        
        load([dataPath ['neuronsMouse' num2str(animalsList(loopAnimal)) '.mat']])
        
        allSpikeTimes = sort(cell2num([neurons.spikeTimes]));

        [N,edges] = histcounts(allSpikeTimes,'BinWidth',0.001,'BinLimits',[100,200]);
     end
     

end