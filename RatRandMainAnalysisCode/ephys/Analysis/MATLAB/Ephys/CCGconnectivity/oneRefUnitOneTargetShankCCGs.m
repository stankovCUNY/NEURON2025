function [tarUnitStackSynch] = oneRefUnitOneTargetShankCCGs(resRef,tarShankID,binOfInterestBackup) 
    
    jscale         = 1;
    duration       = 0.007;
    fig_use        = 102;
    njitter        = 1000;
    alpha          = 0.25;
    for_grant      = false;
    fs             = 30000;
    binSize        = 1/fs;

    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";

    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-basics.mat'),   'basics');
    
    % reference shank
%     for i = 1:size(spikes.RoyMaze1,2)
%         if size(spikes.RoyMaze1(i).time,2) == refTotalNoOfSpikes
%             
%         end
%     end
    
    tarAllUnits = [];

    % target shank
    for i = 1:size(spikes.RoyMaze1,2)
        if spikes.RoyMaze1(i).id(1) == tarShankID
            temp = spikes.RoyMaze1(i).time;
            temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
            tarUnitStack{i} = temp;
            tarAllUnits  = [tarAllUnits temp];
        else
            tarUnitStack{i} = [];
        end
        tarUnitStack = tarUnitStack(~cellfun('isempty',tarUnitStack));
    end
    
    for i = 1:size(tarUnitStack,2)
        
        % figure(i)
        
        resTar = tarUnitStack{i}'/1e6;
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resRef,resTar,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', i, ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',true);
        % title(num2str(i))
        
        %%
        
        [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(resRef, resTar, ones(length(LSPExc),1),duration);
        
        resRefSync = [];
        resTarSync = [];
        
        % error handling if statements
        if (sum(GSPExc) == 0) && ~isempty(binOfInterestBackup)
            if isempty(SyncSpBinAll{binOfInterestBackup})
                resTarSync = [];
            else
                resTarSync = SyncSpBinAll{binOfInterestBackup}(:,2);
            end
        else 
            for k = 1:length(GSPExc)
    
                if (GSPExc(k) == 1) %|| (GSPInh(k) == 1)
    %                 res1sync = [res1sync; SyncSpBinAll{k}(:,1)];
                    resTarSync = [resTarSync; SyncSpBinAll{k}(:,2)];
                end
    
            end
        end
        
        tarUnitStackSynch{i} = resTarSync;
        
    end
    
end