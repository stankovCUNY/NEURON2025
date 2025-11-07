function prepBOToptotaggingData

    geneStr = 'Vip';

    % copied code from https://www.mathworks.com/matlabcentral/fileexchange/90900-brain-observatory-toolbox
    
    figPath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/figures/';
    dataPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_BOT/data/';
    
    % obtain all sessions in the ephys dataset
    sessions = bot.listSessions('ephys');
    
    if strcmp(geneStr,"Sst")
        % Sst
        sessions = sessions(strcmp(sessions.full_genotype,"Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt"),:);
    elseif strcmp(geneStr,"Vip")
        % Vip
        sessions = sessions(strcmp(sessions.full_genotype,"Vip-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt"),:);
    elseif strcmp(geneStr,"Pvalb")
        % Pvalb
        sessions = sessions(strcmp(sessions.full_genotype,"Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt"),:);
    end
        
    desiredSessions = sessions;
    animalsList = desiredSessions.id;
    
    for loopAnimal = 1:length(animalsList)
        
        display([num2str((loopAnimal/length(animalsList))*100) '% done'])
        
        animalID = double(animalsList(loopAnimal));
        
         % probe is in animal 831882777 and 839557629 have error with bot.getChannels
%         if animalID == 831882777 
%             continue
%         elseif animalID == 839557629 
%             continue
%         end
        
        session = bot.getSessions(desiredSessions(desiredSessions.id == animalID, :));
        units = session.units;
        
        %% optotagging
        load([dataPath 'optoStimData_' geneStr '_' num2str(animalID) '.mat'])
        
        idx = find((optoStimData.duration > 0.009) & (optoStimData.duration < 0.02));
                
        trials.start_time    = optoStimData.start_time(idx);
        trials.stop_time     = optoStimData.stop_time(idx);
        trials.condition     = optoStimData.condition(idx,:);
        trials.level         = optoStimData.level(idx);
        trials.duration      = optoStimData.duration(idx);
        trials.stimulus_name = optoStimData.stimulus_name(idx,:);
        
        time_resolution = 0.0005; % 0.5 ms bins;
        
        bin_edges = -0.01:time_resolution:0.025;
        
        dataStruct = optotagging_spike_counts(bin_edges, trials, units, session);
        
        plot_optotagging_response(dataStruct,bin_edges,geneStr,animalID,figPath)
        
        %%
        
        % between -0.01 and -0.002 ms
        baseline = dataStruct.spikeCounts(:,(dataStruct.timeRelativeToStimulusOnset >= -0.01) & (dataStruct.timeRelativeToStimulusOnset <= -0.002),:);
        baselineAvgFRperTrial  = squeeze(mean(baseline,2)/time_resolution);
        baselineAvgFR          = mean(baselineAvgFRperTrial);
        
        % between 0.001 and 0.009 ms
        evoked   = dataStruct.spikeCounts(:,(dataStruct.timeRelativeToStimulusOnset >= 0.001) & (dataStruct.timeRelativeToStimulusOnset <= 0.009),:);
        evokedAvgFRperTrial    = squeeze(mean(evoked,2)/time_resolution);
        evokedAvgFR            = mean(evokedAvgFRperTrial);
        
        %%
        
        hcomb = figure(102);
        fig_use = 102;
        pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
        resolution_use = '-r300';
        
        titleStr = [geneStr ' mouse ' num2str(animalID) ' evoked rate scatter plot'];
        saveStr  = [geneStr '_mouse_' num2str(animalID) '_evoked_rate_scatter_plot'];
        
        axis_limit = 250;
        scatter(baselineAvgFR,evokedAvgFR)
        hold on
        line([0 axis_limit],[0 axis_limit],  'LineStyle','--','LineWidth',1,'Color','k') % unitary line
        line([0 axis_limit],[0 axis_limit*2],'LineStyle','--','LineWidth',1,'Color','r') % cutoff line
        hold off
        
        xlim([0 axis_limit])
        ylim([0 axis_limit])
        xlabel('Baseline rate (Hz)')
        ylabel('Evoked rate (Hz)')
        title(titleStr)
        
        save_file = fullfile([figPath num2str(animalID) '/'], saveStr);
        print(fig_use,save_file,'-djpeg',resolution_use)
        
        close all
        
        %%
        
        cre_pos_units = dataStruct.unitID((evokedAvgFR./(baselineAvgFR + 1)) > 2);
        
        saveStr  = [geneStr '_mouse_' num2str(animalID) '_cre_pos_units'];
        save([figPath saveStr '.mat'],'cre_pos_units');
        
    end
    
end

function dataStruct = optotagging_spike_counts(bin_edges, trials, units, session)

    time_resolution = mean(diff(bin_edges));
    
    spike_matrix = zeros(length(trials.start_time), length(bin_edges), size(units,1));
    
    for loopUnit = 1:size(units,1)
        
        loopUnit
        
         % probe is in animal 763673393 and probe E
%         if (session.id == 763673393) && ((loopUnit > 423) && (loopUnit < 548))
%             continue
%         elseif (session.id == 831882777) && ((loopUnit > 513) && (loopUnit < 658)) % probe is in animal 831882777 and probe F
%             continue
%         elseif (session.id == 839557629) && ((loopUnit > 0) && (loopUnit < 451)) % probe is in animal 839557629 and all probes 
%             continue    
%         end
        
        try
            unit = bot.getUnits(units(loopUnit, :));
            spike_times = session.spike_times.spike_times{double(session.spike_times{:,1}) == double(unit.id),1};
        catch
            spike_times = session.spike_times.spike_times{session.spike_times.unit_id == double(units.id(loopUnit, :))}; 
        end
        
        for loopTrials = 1:length(trials.start_time)
                        
            in_range = (spike_times > (trials.start_time(loopTrials) + bin_edges(1))) & ...
                       (spike_times < (trials.start_time(loopTrials) + bin_edges(end)));
            
            binned_times = ceil((spike_times(in_range) - (trials.start_time(loopTrials) + bin_edges(1)))/time_resolution);
            spike_matrix(loopTrials,binned_times,loopUnit) = 1;       
        end
    end
    
    dataStruct.spikeCounts                 = spike_matrix;
    dataStruct.trialID                     = 1:length(trials.start_time); 
    dataStruct.unitID                      = units.id;
    dataStruct.timeRelativeToStimulusOnset = bin_edges;
    
end 

function plot_optotagging_response(dataStruct,bin_edges,geneStr,animalID,figPath)
    
    hcomb = figure(102);
    fig_use = 102;
    pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));
    resolution_use = '-r300';
    
    titleStr = [geneStr ' mouse ' num2str(animalID) ' optotagging response'];
    saveStr  = [geneStr '_mouse_' num2str(animalID) '_optotagging_response'];
    
    time_resolution = mean(diff(dataStruct.timeRelativeToStimulusOnset));

    imagesc(squeeze(mean(dataStruct.spikeCounts)/time_resolution)');
    
    hold on
    line([22 22],[0 length(dataStruct.unitID)],'LineStyle','--','LineWidth',1,'Color','w') % at bound 0.0005 ms
    line([40 40],[0 length(dataStruct.unitID)],'LineStyle','--','LineWidth',1,'Color','w') % at bound 0.0095 ms
    hold off
    
    set(gca, 'XTick', (1:7)*10+1, 'XTickLabel', bin_edges((1:7)*10+1))
    xlabel('Time (s)')
    ylabel('Unit #')
    
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0, 300]);
    title(titleStr)
    
    save_file = fullfile([figPath num2str(animalID) '/'], saveStr);
    print(fig_use,save_file,'-djpeg',resolution_use)
    
    close all
    
end