function BOTscreenConv(nullFlag)
    
%     animalID = 750332458;
%     probeID = 757904554;
    
%     animalID = 715093703;
%     probeID = 810755803;
    
    % obtain all sessions in the ephys dataset
    sessions = bot.listSessions('ephys');
    
    % obtain all sessions with the desired session_type and then select the first one
%     desiredSessions = sessions(sessions.session_type == 'brain_observatory_1.1', :); % use Brain Observatory for now
%     desiredSessions = sessions(sessions.session_type == 'functional_connectivity', :); % use Brain Observatory for now

    animalsList = sessions.id;

    for loopAnimal = 1:length(animalsList)
        
        tic
        
        display(['animal: ' num2str(animalsList(loopAnimal))])
        
        dataPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
        
        animalID = double(animalsList(loopAnimal));
        
        % probe is in animal 831882777 and 839557629 have error with bot.getChannels
%         if animalID == 831882777 
%             continue
%         elseif animalID == 839557629 
%             continue
%         end
        
        load([dataPath ['neuronsMouse' num2str(animalID) '.mat']])
        
        % probe is in animal 763673393 and probe E
%         if (animalID == 763673393) 
%             neurons(423:548,:) = [];
%         end
        
        probeString = string(neurons.probeName);
        probeList   = unique(probeString);
        
        for loopProbe = 1:length(probeList)
            
            display(['probe: ' num2str(probeList(loopProbe))])
            
            probeFilter  = strcmp(probeList(loopProbe),probeString);    
            neuronsProbe = neurons(probeFilter,:);

            brainRegionString = string(neuronsProbe.brainRegion);
            brainRegionList   = unique(brainRegionString);
        
            for loopBrainRegions = 1:length(brainRegionList)
                
                display(['region: ' brainRegionList{loopBrainRegions}])
                
                brainRegionFilter  = strcmp(brainRegionList(loopBrainRegions),brainRegionString); 
                neuronsBrainRegion = neuronsProbe(brainRegionFilter,:); 
                
                spikeTimes = []; % spiket for Eran
                spikeInd   = []; % spikeind for Eran

                for i = 1:length(neuronsBrainRegion.unitID) 
                    spikeTimes = [neuronsBrainRegion.spikeTimes{i}; spikeTimes];
                    spikeInd   = [ones(length(neuronsBrainRegion.spikeTimes{i}),1)*(double(neuronsBrainRegion.unitID(i))); spikeInd];
                end
                
                if nullFlag
                    one_ms = 0.001;
                    delta = 5*one_ms;
                    spikeTimes = floor(spikeTimes/delta)*delta + rand(size(spikeTimes))*delta;
                end
                
                cellID     = double(neuronsBrainRegion.unitID)';             % Cells for Eran
                probeID    = double(cell2num(neuronsBrainRegion.probeID))';	 % shank for Eran, set as 1 bc of neuropixel
                SampleRate = 30000;                              % SampleRate for Eran
                alpha      = 0.05;                               % alpah for Eran
                wintype    = 'gauss';                            % wintype for Eran

                %% jscale 1
                jscale = 1; % jscale for Eran

                [ExcPairs, ...
                 InhPairs, ...
                 GapPairs, ...
                 Rzero] = EranConv_group_forUtku(spikeTimes', ...
                                                 spikeInd', ...
                                                 cellID, ...
                                                 SampleRate, ...
                                                 jscale, ...
                                                 alpha,...
                                                 probeID, ...
                                                 wintype);

                pairs.ExcPairs   = ExcPairs;
                pairs.InhPairs   = InhPairs;
                pairs.GapPairs   = GapPairs;
                pairs.RZero      = Rzero;
                pairs.jscale     = jscale;
                pairs.nullFlag   = nullFlag;
                pairs.spikeTimes = spikeTimes;
                pairs.spikeInd   = spikeInd;
                
                if nullFlag
                    save(fullfile(dataPath,['BOTmouse' num2str(animalID) 'NULL' ...
                                num2str(probeList(loopProbe)) ... 
                                'region' brainRegionList{loopBrainRegions} ...
                                '_jscale' num2str(jscale) '_alpha' ...
                                num2str(round(alpha*100)) '_pairs']), ...
                                'pairs', 'jscale', 'alpha')
                else
                    save(fullfile(dataPath,['BOTmouse' num2str(animalID) ...
                                num2str(probeList(loopProbe)) ... 
                                'region' brainRegionList{loopBrainRegions} ...
                                '_jscale' num2str(jscale) '_alpha' ...
                                num2str(round(alpha*100)) '_pairs']), ...
                                'pairs', 'jscale', 'alpha')
                end
                
                %% jscale 5
                jscale = 5;  % jscale for Eran

                [ExcPairs, ...
                 InhPairs, ...
                 GapPairs, ...
                 Rzero] = EranConv_group_forUtku(spikeTimes', ...
                                                 spikeInd', ...
                                                 cellID, ...
                                                 SampleRate, ...
                                                 jscale, ...
                                                 alpha,...
                                                 probeID, ...
                                                 wintype);

                pairs.ExcPairs   = ExcPairs;
                pairs.InhPairs   = InhPairs;
                pairs.GapPairs   = GapPairs;
                pairs.RZero      = Rzero;
                pairs.jscale     = jscale;
                pairs.nullFlag   = nullFlag;
                pairs.spikeTimes = spikeTimes;
                pairs.spikeInd   = spikeInd;
                
                if nullFlag
                    save(fullfile(dataPath,['BOTmouse' num2str(animalID) 'NULL' ...
                                num2str(probeList(loopProbe)) ... 
                                'region' brainRegionList{loopBrainRegions} ...
                                '_jscale' num2str(jscale) '_alpha' ...
                                num2str(round(alpha*100)) '_pairs']), ...
                                'pairs', 'jscale', 'alpha')
                else
                    save(fullfile(dataPath,['BOTmouse' num2str(animalID) ...
                                num2str(probeList(loopProbe)) ... 
                                'region' brainRegionList{loopBrainRegions} ...
                                '_jscale' num2str(jscale) '_alpha' ...
                                num2str(round(alpha*100)) '_pairs']), ...
                                'pairs', 'jscale', 'alpha')
                end
                        
            end
        end
        
        toc
    end
    
 end