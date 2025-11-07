function prepBOTdata

    % copied code from https://www.mathworks.com/matlabcentral/fileexchange/90900-brain-observatory-toolbox
    
    savePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
    
    % obtain all sessions in the ephys dataset
    sessions = bot.listSessions('ephys');
    
    % obtain all sessions with the desired session_type and then select the first one
%     desiredSessions = sessions(sessions.session_type == 'brain_observatory_1.1', :); % use Brain Observatory for now
%     desiredSessions = sessions(sessions.session_type == 'functional_connectivity', :); 
    desiredSessions = sessions;

    animalsList = desiredSessions.id;
    
    for loopAnimal = 1:length(animalsList)
        
        display([num2str((loopAnimal/length(animalsList))*100) '% done'])
        
        animalID = double(animalsList(loopAnimal));
        session = bot.getSessions(desiredSessions(desiredSessions.id == animalID, :));

    end
    
%     animalID = 715093703; % pick animal
%     animalID = 750332458; % pick animal
%     sess = bot.getSessions(desiredSessions(desiredSessions.id == animalID, :)); 
    
    % obtain all probes in the session
%     probes = sess.probes;
    
    % obtain the brain region that each probe examines
    %disp(probes.ephys_structure_acronyms)

    %% single probe
%     % pick probe
%     probeID = 810755803;
%     probe = bot.getProbes(probes(probes.id == probeID, :));
%     
%     % select all units of the probe
%     probeUnits = probe.units;
%     
%     % select all units in the CA1/CA3 using logical indexing
%     probeUnits = probeUnits(probeUnits.ephys_structure_acronym == 'CA1' | ...
%                             probeUnits.ephys_structure_acronym == 'CA3', :);
%                         
%     % obtain `unit` objects
%     for i = 1:size(probeUnits,1)
%         
%         i
%         tic
%         unit = bot.getUnits(probeUnits(i, :));
%         
%         neurons.unitID(i)      = unit.id;
%         neurons.peakChan(i)    = probeUnits.peak_channel(double(probeUnits.id) == double(unit.id));
%         neurons.probeName{i}   = probeUnits.probe_name(double(probeUnits.id) == double(unit.id));
%         neurons.probeID{i}     = probeUnits.ephys_probe_id(double(probeUnits.id) == double(unit.id));
%         neurons.animalID{i}    = probeUnits.ephys_session_id(double(probeUnits.id) == double(unit.id));
%         neurons.brainRegion{i} = probeUnits.ephys_structure_acronym(double(probeUnits.id) == double(unit.id));
%         neurons.spikeTimes{i}  = sess.spike_times.spike_times{double(sess.spike_times{:,1}) == double(unit.id),1};
%         toc
%     end
%     
%     saveName = ['neuronsMouse' num2str(animalID) 'probe' num2str(probeID) '.mat'];
%     
%     save([savePath saveName],'neurons')

    %% all probes
    
    for loopAnimal = 24:length(animalsList)
        
        display([num2str((loopAnimal/length(animalsList))*100) '% done'])
        
        animalID = double(animalsList(loopAnimal));
        
        % probe is in animal 839557629 have error with bot.getChannels
%         if animalID == 839557629 
%             continue
%         end
        
        session = bot.getSessions(desiredSessions(desiredSessions.id == animalID, :));
        units = session.units;
        
        invalid_times = session.invalid_times;
        invalid_times_resorted = table;
        for loopInvTimes = 1:size(invalid_times,1)
            invalid_times_resorted.start_time(loopInvTimes) = invalid_times.start_time(loopInvTimes);
            invalid_times_resorted.stop_time(loopInvTimes)  = invalid_times.stop_time(loopInvTimes);
            invalid_times_resorted.probeID{loopInvTimes}    = invalid_times.tags{loopInvTimes,1}(3);
        end
        
    %     allUnits = allUnits(allUnits.ephys_structure_acronym == 'CA1' | ...
    %                         allUnits.ephys_structure_acronym == 'CA3', :);

    %     allUnits = allUnits((allUnits.ephys_structure_acronym == 'MB' | allUnits.ephys_structure_acronym == 'LGd') | ... 
    %                         (allUnits.ephys_structure_acronym == 'VISl' | allUnits.ephys_structure_acronym == 'VISrl'), :);

        channelsTable = session.channels;
%         channelUnderStudy = bot.getChannels(channelsTable(1, :));
        
        % optotag
        if strcmp(sessions.full_genotype{sessions.id == animalID},"Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt")
            load([savePath 'Sst_mouse_' num2str(animalID) '_cre_pos_units.mat' ]);
            creLabel= 'Sst';
        elseif strcmp(sessions.full_genotype{sessions.id == animalID},"Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt")
            load([savePath 'Pvalb_mouse_' num2str(animalID) '_cre_pos_units.mat' ]);
            creLabel = 'Pvalb';
        elseif strcmp(sessions.full_genotype{sessions.id == animalID},"Vip-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt")
            load([savePath 'Vip_mouse_' num2str(animalID) '_cre_pos_units.mat' ]);
            creLabel = 'Vip';
        elseif strcmp(sessions.full_genotype{sessions.id == animalID},"wt/wt")
            cre_pos_units = [];
        end
        
        % initiate neurons table
        neurons = table;
        
        % obtain `unit` objects
        for loopUnit = 1:size(units,1)

            loopUnit
             
            % probe is in animal 763673393 and probe E
%             if (animalID == 763673393) && ((loopUnit > 423) && (loopUnit < 548))
%                 continue
%             elseif (animalID == 831882777) && ((loopUnit > 513) && (loopUnit < 658)) % probe is in animal 831882777 and probe F
%                 continue
%             end
            
            tic
                  
%             unit = bot.getUnits(units(loopUnit, :));
            unitID = units.id(loopUnit);
            fs     = units.probe_sampling_rate(loopUnit);
            
            % hacky error handling when listed fs is 0Hz. All probes so far
            % are sampled at 30kHz
            if fs == 0
               fs = 3e4;
            end
            
            % assign optotag
            if ~isempty(find(double(cre_pos_units) == double(unitID)))
                neurons.putativeOptotagType{loopUnit} = creLabel;
            else
                neurons.putativeOptotagType{loopUnit} = 'No opto label';
            end
                
            spikeTimes         = session.spike_times.spike_times{double(session.spike_times{:,1}) == double(unitID),1};
            
%             avgWaveforms = channelUnderStudy.session.mean_waveforms{double(channelUnderStudy.session.mean_waveforms.unit_id) == double(unitID),2}{1,1};\
            avgWaveforms = session.mean_waveforms.waveform_mean{find(double(double(session.mean_waveforms.unit_id) == double(unitID))),1};
            probeName    = units.probe_name(double(units.id) == double(unitID));
            peakChan     = units.peak_channel(double(units.id) == double(unitID)) + 1;
%             [~,troughTime] = min(avgWaveforms(peakChan,:));
%             [~,troughToPeakLength] = max(avgWaveforms(peakChan,troughTime:end));
%             troughToPeakLength = (troughToPeakLength/fs)*1000; % convert from samples to ms
            troughToPeakLength = units.waveform_duration(double(units.id) == double(unitID)); % better results
            
            % remove spike during invalid times 
            if ~isempty(invalid_times_resorted)
                if ~isempty(find(strcmp(string(invalid_times_resorted.probeID),string(probeName)))) 
                    indTimesIdx = find(strcmp(string(invalid_times_resorted.probeID),string(probeName)));

                    for loopInvalidTimes = 1:length(indTimesIdx)
                        spikeTimes(...
                                  (spikeTimes > invalid_times_resorted.start_time(indTimesIdx(loopInvalidTimes))) &...
                                  (spikeTimes < invalid_times_resorted.stop_time(indTimesIdx(loopInvalidTimes)))...
                                  ) = [];
                    end
                end
            end
                
            neurons.putativeCellType{loopUnit}  = cellClass(spikeTimes,troughToPeakLength,fs);
            neurons.spikeTimes{loopUnit}        = spikeTimes;
            neurons.unitID(loopUnit)            = unitID;
            neurons.peakChan(loopUnit)          = peakChan;
            neurons.probeName{loopUnit}         = probeName;
            neurons.probeID{loopUnit}           = units.ephys_probe_id(double(units.id) == double(unitID));
            neurons.animalID{loopUnit}          = units.ephys_session_id(double(units.id) == double(unitID));
            neurons.brainRegion{loopUnit}       = units.ephys_structure_acronym(double(units.id) == double(unitID));
            neurons.probeVertPos(loopUnit)      = units.probe_vertical_position(double(units.id) == double(unitID));
            neurons.probeHortPos(loopUnit)      = units.probe_horizontal_position(double(units.id) == double(unitID));
            neurons.avgWaveforms{loopUnit}      = avgWaveforms;
            
            neurons.Lratio(loopUnit)            = units.L_ratio(double(units.id) == double(unitID));
            neurons.Dprime(loopUnit)            = units.d_prime(double(units.id) == double(unitID));
            neurons.firingRate(loopUnit)        = units.firing_rate(double(units.id) == double(unitID));
            neurons.ISIviolations(loopUnit)     = units.isi_violations(double(units.id) == double(unitID));
            neurons.isolationDistance(loopUnit) = units.isolation_distance(double(units.id) == double(unitID));
            neurons.snr(loopUnit)               = units.snr(double(units.id) == double(unitID));
            neurons.amplitudeCutoff(loopUnit)   = units.amplitude_cutoff(double(units.id) == double(unitID));
            neurons.presenceRatio(loopUnit)     = units.presence_ratio(double(units.id) == double(unitID));
            
            toc

        end

        saveName = ['neuronsMouse' num2str(animalID) 'qualityMetrics.mat'];

        save([savePath saveName],'neurons')
    end
end